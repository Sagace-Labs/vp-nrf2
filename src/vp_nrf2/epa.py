"""Rebuild the training-only EPA AEID 1110 table from its official export.

The v4.3 activity calls are paired with EPA's archived chemical identity SDF;
unresolved identities and Inconclusive calls do not become labels.
"""

from __future__ import annotations

import io
import re
import tempfile
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

import pandas as pd
import requests
from rdkit import Chem, RDLogger

from vp_core import dataset
from vp_core.standardise import standardise_many
from vp_nrf2 import data
from vp_nrf2.are_ensemble import external_only

ACTIVITY_URL = "https://clowder.edap-cluster.com/api/files/6a85a260e4b0731a36bd04cc"
IDENTITY_URL = "https://clowder.edap-cluster.com/api/files/6114f600e4b0856fdc65865c"


def _download(url: str, path: Path) -> None:
    with requests.get(url, stream=True, timeout=120) as response:
        response.raise_for_status()
        with path.open("wb") as output:
            for chunk in response.iter_content(1024 * 1024):
                output.write(chunk)


def _worksheet(archive: zipfile.ZipFile, sheet_number: int) -> list[list[str]]:
    """Read sparse XLSX cells without shifting absent columns."""
    ns = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    strings = []
    if "xl/sharedStrings.xml" in archive.namelist():
        root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
        strings = ["".join(item.itertext()) for item in root.findall("m:si", ns)]
    sheet = ET.fromstring(archive.read(f"xl/worksheets/sheet{sheet_number}.xml"))
    rows = []
    for row in sheet.findall(".//m:row", ns):
        values = {}
        for cell in row.findall("m:c", ns):
            column = 0
            for char in re.match(r"[A-Z]+", cell.attrib["r"]).group():
                column = 26 * column + ord(char) - ord("A") + 1
            value = cell.find("m:v", ns)
            content = "" if value is None else value.text
            if cell.get("t") == "s":
                content = strings[int(content)]
            elif cell.get("t") == "inlineStr":
                content = "".join(cell.find("m:is", ns).itertext())
            values[column - 1] = content
        rows.append([values.get(i, "") for i in range(max(values, default=-1) + 1)])
    return rows


def _identity_map(path: Path) -> dict[str, str]:
    RDLogger.DisableLog("rdApp.*")
    with zipfile.ZipFile(path) as archive:
        sdf = next(name for name in archive.namelist() if name.endswith(".sdf"))
        result = {}
        for mol in Chem.ForwardSDMolSupplier(io.BytesIO(archive.read(sdf))):
            if mol is not None and mol.HasProp("DSSTox_Substance_Id"):
                result[mol.GetProp("DSSTox_Substance_Id")] = Chem.MolToSmiles(mol)
    return result


def rebuild(activity: Path, identity: Path) -> pd.DataFrame:
    """Parse AEID 1110 calls, standardise structures and remove primary identities."""
    with zipfile.ZipFile(activity) as archive:
        name = next(name for name in archive.namelist() if "aeid1110_" in name.lower())
        with zipfile.ZipFile(io.BytesIO(archive.read(name))) as workbook:
            rows = _worksheet(workbook, 2)
    expected = ["PUBCHEM_RESULT_TAG", "PUBCHEM_EXT_DATASOURCE_REGID", "PUBCHEM_ACTIVITY_OUTCOME"]
    if rows[0][:3] != expected:
        raise ValueError("Unexpected EPA PubChem export schema")
    records = pd.DataFrame(rows[6:], columns=rows[0])
    records["dtxsid"] = records.PUBCHEM_ACTIVITY_URL.str.split("/").str[-1]
    records["label"] = records.PUBCHEM_ACTIVITY_OUTCOME.map({"1": 0, "2": 1})
    records["smiles_raw"] = records.dtxsid.map(_identity_map(identity))
    called = records.dropna(subset=["smiles_raw", "label"]).copy()
    called["smiles"], called["inchikey"] = standardise_many(called.smiles_raw.tolist())
    rows = []
    for key, group in called.dropna(subset=["smiles", "inchikey"]).groupby("inchikey"):
        share = float(group.label.mean())
        if share != 0.5:
            rows.append({"inchikey": key, "smiles": group.smiles.iloc[0],
                         "label": int(share > 0.5), "n_calls": len(group)})
    full = pd.DataFrame(rows).sort_values("inchikey").reset_index(drop=True)
    return external_only(full, data.load())


def fetch() -> pd.DataFrame:
    """Explicitly download and rebuild the training-only standardized table."""
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        activity, identity = root / "activity.zip", root / "identity.zip"
        _download(ACTIVITY_URL, activity)
        _download(IDENTITY_URL, identity)
        table = rebuild(activity, identity)
    dataset.write_table(table, data.EPA_PATH, labels=("label",))
    return table
