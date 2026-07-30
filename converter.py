import traceback
from typing import Optional, Tuple, Dict
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from utils import detect_input_type, logger
from pubchem_lookup import name_to_smiles
from validation import compute_diagnostics, compare_molecules, generate_report


_OPSIN_AVAILABLE = None

def _try_opsin(name: str) -> Optional[str]:
    global _OPSIN_AVAILABLE
    if _OPSIN_AVAILABLE is False:
        return None
    try:
        import pyopsin
        _OPSIN_AVAILABLE = True
        opsin = pyopsin.PyOpsin()
        result = opsin.to_smiles_single(name)
        if result and result != "None":
            return result
    except ImportError:
        try:
            import opsin
            _OPSIN_AVAILABLE = True
            smiles = opsin.name_to_smiles(name)
            if smiles:
                return smiles
        except ImportError:
            _OPSIN_AVAILABLE = False
            logger.info("OPSIN not available (requires Java). Falling back to PubChem.")
        except Exception as exc:
            logger.debug(f"OPSIN (legacy) failed for '{name}': {exc}")
    except Exception as exc:
        msg = str(exc)
        if "JVMNotFoundException" in msg or "JAVA_HOME" in msg:
            logger.info("Java required for pyopsin. Install Java or use SMILES input directly.")
        else:
            logger.debug(f"pyopsin failed for '{name}': {exc}")
    return None


def _try_pubchem(name: str) -> Optional[str]:
    try:
        smiles = name_to_smiles(name)
        if smiles:
            return smiles
    except Exception as exc:
        logger.debug(f"PubChem lookup failed for '{name}': {exc}")
    return None


def convert_to_mol(
    input_text: str,
) -> Tuple[Optional[Chem.Mol], str, str]:
    input_text = input_text.strip()
    if not input_text:
        return None, "unknown", "Empty input"

    detected_type = detect_input_type(input_text)
    mol = None
    status = "pending"

    if detected_type == "smiles":
        mol = Chem.MolFromSmiles(input_text)
        if mol is None:
            status = "failed"
            return None, "smiles", "Invalid SMILES string"
        status = "success"

    elif detected_type == "inchi":
        mol = Chem.MolFromInchi(input_text)
        if mol is None:
            status = "failed"
            return None, "inchi", "Invalid InChI string"
        status = "success"

    elif detected_type == "inchikey":
        try:
            mol = Chem.MolFromInchiKey(input_text)
            if mol is None:
                raise ValueError("Invalid InChIKey")
            status = "success"
        except Exception:
            status = "failed"
            return None, "inchikey", "Invalid InChIKey"

    elif detected_type == "iupac_or_common":
        smiles = _try_opsin(input_text)
        if smiles:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                detected_type = "iupac"
                status = "success"
            else:
                mol = None

        if mol is None:
            logger.debug(f"OPSIN failed for '{input_text}', trying PubChem...")
            smiles = _try_pubchem(input_text)
            if smiles:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    detected_type = "common_name"
                    status = "success"
                else:
                    mol = None

        if mol is None:
            status = "failed"
            return None, "iupac_or_common", f"Cannot resolve name: '{input_text}'"

    else:
        status = "failed"
        return None, "unknown", f"Unrecognized input type: '{input_text}'"

    if mol is not None:
        try:
            Chem.SanitizeMol(mol)
        except Exception as exc:
            logger.warning(f"Sanitization warning for '{input_text}': {exc}")

    return mol, detected_type, status


def get_molecule_info(mol: Chem.Mol, input_text: str, detected_type: str) -> dict:
    if mol is None:
        return {}

    info = {
        "original_input": input_text,
        "detected_type": detected_type,
    }

    try:
        info["canonical_smiles"] = Chem.MolToSmiles(mol, canonical=True)
    except Exception:
        info["canonical_smiles"] = "N/A"

    try:
        info["molecular_formula"] = rdMolDescriptors.CalcMolFormula(mol)
    except Exception:
        info["molecular_formula"] = "N/A"

    try:
        info["molecular_weight"] = round(Descriptors.MolWt(mol), 4)
    except Exception:
        info["molecular_weight"] = "N/A"

    try:
        info["exact_mass"] = round(Descriptors.ExactMolWt(mol), 4)
    except Exception:
        info["exact_mass"] = "N/A"

    try:
        info["inchi"] = Chem.MolToInchi(mol)
    except Exception:
        info["inchi"] = "N/A"

    try:
        info["inchikey"] = Chem.MolToInchiKey(mol)
    except Exception:
        info["inchikey"] = "N/A"

    try:
        info["num_atoms"] = mol.GetNumAtoms()
    except Exception:
        info["num_atoms"] = "N/A"

    try:
        info["num_heavy_atoms"] = mol.GetNumHeavyAtoms()
    except Exception:
        info["num_heavy_atoms"] = "N/A"

    try:
        info["num_rings"] = rdMolDescriptors.CalcNumRings(mol)
    except Exception:
        info["num_rings"] = "N/A"

    try:
        info["num_aromatic_rings"] = rdMolDescriptors.CalcNumAromaticRings(mol)
    except Exception:
        info["num_aromatic_rings"] = "N/A"

    try:
        info["num_stereocenters"] = Chem.rdMolDescriptors.CalcNumAtomStereoCenters(mol)
    except Exception:
        info["num_stereocenters"] = "N/A"

    try:
        info["logp"] = round(Descriptors.MolLogP(mol), 3)
    except Exception:
        info["logp"] = "N/A"

    try:
        diag = compute_diagnostics(mol)
        info["diagnostics"] = diag
        info["report"] = generate_report(mol, diag)
    except Exception:
        info["diagnostics"] = {}
        info["report"] = "Diagnostics unavailable."

    return info
