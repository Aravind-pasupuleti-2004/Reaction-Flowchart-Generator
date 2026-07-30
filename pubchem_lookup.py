import json
import re
import ssl
import urllib.request
import urllib.error
from urllib.parse import quote
from typing import Optional

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
CACTUS_BASE = "https://cactus.nci.nih.gov/chemical/structure"

_STEREO_PREFIX = re.compile(r"^\([\d,]*[RSEZ][\d,]*\)[- ]?")


def _strip_stereo(name: str) -> str:
    return _STEREO_PREFIX.sub("", name).strip()


def _fetch_json(url: str, timeout: int = 10) -> Optional[dict]:
    try:
        ctx = ssl._create_unverified_context()
        req = urllib.request.Request(url, headers={"User-Agent": "ChemStructGen/1.0"})
        with urllib.request.urlopen(req, timeout=timeout, context=ctx) as resp:
            return json.loads(resp.read().decode())
    except Exception:
        return None


def _fetch_text(url: str, timeout: int = 15) -> Optional[str]:
    try:
        ctx = ssl._create_unverified_context()
        req = urllib.request.Request(url, headers={"User-Agent": "ChemStructGen/1.0"})
        with urllib.request.urlopen(req, timeout=timeout, context=ctx) as resp:
            return resp.read().decode().strip()
    except Exception:
        return None


def _try_pubchem(name: str, properties: str) -> Optional[str]:
    url = f"{PUBCHEM_BASE}/compound/name/{quote(name)}/property/{properties}/JSON"
    data = _fetch_json(url)
    if data and "PropertyTable" in data:
        props = data["PropertyTable"].get("Properties", [])
        if props:
            return props[0].get("CanonicalSMILES") or props[0].get("ConnectivitySMILES")
    return None


def _try_cactus(name: str) -> Optional[str]:
    url = f"{CACTUS_BASE}/{quote(name)}/smiles"
    return _fetch_text(url)


def name_to_smiles(name: str) -> Optional[str]:
    result = _try_pubchem(name, "CanonicalSMILES,ConnectivitySMILES")
    if result:
        return result
    stripped = _strip_stereo(name)
    if stripped != name:
        result = _try_pubchem(stripped, "CanonicalSMILES,ConnectivitySMILES")
        if result:
            return result
    result = _try_cactus(name)
    if result:
        return result
    if stripped != name:
        return _try_cactus(stripped)
    return None


def name_to_inchi(name: str) -> Optional[str]:
    url = f"{PUBCHEM_BASE}/compound/name/{quote(name)}/property/InChI/JSON"
    data = _fetch_json(url)
    if data and "PropertyTable" in data:
        props = data["PropertyTable"].get("Properties", [])
        if props:
            return props[0].get("InChI")
    return None


def name_to_info(name: str) -> Optional[dict]:
    url = f"{PUBCHEM_BASE}/compound/name/{quote(name)}/property/MolecularFormula,MolecularWeight,CanonicalSMILES,InChI/JSON"
    data = _fetch_json(url)
    if data and "PropertyTable" in data:
        props = data["PropertyTable"].get("Properties", [])
        if props:
            return props[0]
    return None
