import json
import ssl
import urllib.request
import urllib.error
from typing import Optional

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


def _fetch_json(url: str, timeout: int = 10) -> Optional[dict]:
    try:
        ctx = ssl._create_unverified_context()
        req = urllib.request.Request(url, headers={"User-Agent": "ChemStructGen/1.0"})
        with urllib.request.urlopen(req, timeout=timeout, context=ctx) as resp:
            return json.loads(resp.read().decode())
    except Exception:
        return None


def name_to_smiles(name: str) -> Optional[str]:
    from urllib.parse import quote
    url = f"{PUBCHEM_BASE}/compound/name/{quote(name)}/property/CanonicalSMILES,ConnectivitySMILES/JSON"
    data = _fetch_json(url)
    if data and "PropertyTable" in data:
        props = data["PropertyTable"].get("Properties", [])
        if props:
            return props[0].get("CanonicalSMILES") or props[0].get("ConnectivitySMILES")
    return None


def name_to_inchi(name: str) -> Optional[str]:
    url = f"{PUBCHEM_BASE}/compound/name/{urllib.parse.quote(name)}/property/InChI/JSON"
    data = _fetch_json(url)
    if data and "PropertyTable" in data:
        props = data["PropertyTable"].get("Properties", [])
        if props:
            return props[0].get("InChI")
    return None


def name_to_info(name: str) -> Optional[dict]:
    url = f"{PUBCHEM_BASE}/compound/name/{urllib.parse.quote(name)}/property/MolecularFormula,MolecularWeight,CanonicalSMILES,InChI/JSON"
    data = _fetch_json(url)
    if data and "PropertyTable" in data:
        props = data["PropertyTable"].get("Properties", [])
        if props:
            return props[0]
    return None
