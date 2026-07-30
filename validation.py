from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, rdMolChemicalFeatures
from rdkit.Chem.rdchem import BondType
from typing import Optional, Dict, List, Tuple


def compute_diagnostics(mol: Chem.Mol) -> Dict:
    if mol is None:
        return {}

    try:
        Chem.SanitizeMol(mol)
    except Exception:
        pass

    diagnostics = {}

    diagnostics["formula"] = rdMolDescriptors.CalcMolFormula(mol)
    diagnostics["mol_weight"] = round(Descriptors.MolWt(mol), 4)
    diagnostics["exact_mass"] = round(Descriptors.ExactMolWt(mol), 4)
    diagnostics["num_atoms"] = mol.GetNumAtoms()
    diagnostics["num_heavy_atoms"] = mol.GetNumHeavyAtoms()
    diagnostics["num_rings"] = rdMolDescriptors.CalcNumRings(mol)
    diagnostics["num_aromatic_rings"] = rdMolDescriptors.CalcNumAromaticRings(mol)
    diagnostics["num_aliphatic_rings"] = rdMolDescriptors.CalcNumAliphaticRings(mol)
    diagnostics["num_rotatable_bonds"] = rdMolDescriptors.CalcNumRotatableBonds(mol)
    diagnostics["num_hba"] = rdMolDescriptors.CalcNumHBA(mol)
    diagnostics["num_hbd"] = rdMolDescriptors.CalcNumHBD(mol)
    diagnostics["num_stereocenters"] = Chem.rdMolDescriptors.CalcNumAtomStereoCenters(mol)

    try:
        diagnostics["inchi"] = Chem.MolToInchi(mol)
    except Exception:
        diagnostics["inchi"] = "N/A"
    try:
        diagnostics["inchikey"] = Chem.MolToInchiKey(mol)
    except Exception:
        diagnostics["inchikey"] = "N/A"

    try:
        diagnostics["canonical_smiles"] = Chem.MolToSmiles(mol, canonical=True)
    except Exception:
        diagnostics["canonical_smiles"] = "N/A"

    elem_counts = {}
    for atom in mol.GetAtoms():
        sym = atom.GetSymbol()
        elem_counts[sym] = elem_counts.get(sym, 0) + 1
    diagnostics["elemental_composition"] = elem_counts

    functional_groups = _detect_functional_groups(mol)
    diagnostics["functional_groups"] = functional_groups

    diagnostics["logp"] = round(Descriptors.MolLogP(mol), 3)

    return diagnostics


def _detect_functional_groups(mol: Chem.Mol) -> Dict[str, bool]:
    groups = {}
    smarts_patterns = {
        "carboxylic_acid": "[CX3](=O)[OX2H1]",
        "ester": "[CX3](=O)[OX2][CX4]",
        "ketone": "[CX3](=O)[#6]",
        "aldehyde": "[CX3H1](=O)[#6]",
        "alcohol": "[OX2H]",
        "ether": "[OD2]([#6])[#6]",
        "primary_amine": "[NH2]",
        "secondary_amine": "[NH1]",
        "tertiary_amine": "[N]([#6])([#6])[#6]",
        "amide": "[CX3](=O)[NX3]",
        "nitro": "[N+](=O)[O-]",
        "sulfoxide": "[SX3](=O)",
        "sulfone": "[SX4](=O)(=O)",
        "halogen": "[F,Cl,Br,I]",
        "aromatic_ring": "a",
        "oxo": "[CX3]=[OX1]",
        "methoxy": "[OX2]C",
        "lactam": "[NR]1[CX3](=O)[C,c][C,c][C,c][C,c]1",
    }
    for name, smarts in smarts_patterns.items():
        try:
            pat = Chem.MolFromSmarts(smarts)
            if pat:
                groups[name] = mol.HasSubstructMatch(pat)
            else:
                groups[name] = False
        except Exception:
            groups[name] = False
    return groups


def compare_molecules(
    generated_mol: Chem.Mol,
    reference_smiles: Optional[str] = None,
    reference_formula: Optional[str] = None,
    reference_inchikey: Optional[str] = None,
) -> Dict:
    results = {
        "smiles_match": None,
        "formula_match": None,
        "inchikey_match": None,
        "weight_match": None,
        "issues": [],
        "confidence": 0.0,
    }

    gen_smiles = Chem.MolToSmiles(generated_mol, canonical=True)
    gen_formula = rdMolDescriptors.CalcMolFormula(generated_mol)
    gen_inchikey = Chem.MolToInchiKey(generated_mol)
    gen_weight = Descriptors.MolWt(generated_mol)

    checks_passed = 0
    checks_total = 0

    if reference_smiles:
        checks_total += 1
        ref_mol = Chem.MolFromSmiles(reference_smiles)
        if ref_mol:
            ref_smiles = Chem.MolToSmiles(ref_mol, canonical=True)
            if gen_smiles == ref_smiles:
                results["smiles_match"] = True
                checks_passed += 1
            else:
                results["smiles_match"] = False
                results["issues"].append(
                    f"Canonical SMILES mismatch: generated={gen_smiles}, expected={ref_smiles}"
                )
        else:
            results["smiles_match"] = "invalid_reference"

    if reference_formula:
        checks_total += 1
        if gen_formula == reference_formula:
            results["formula_match"] = True
            checks_passed += 1
        else:
            results["formula_match"] = False
            results["issues"].append(
                f"Formula mismatch: generated={gen_formula}, expected={reference_formula}"
            )

    if reference_inchikey:
        checks_total += 1
        if gen_inchikey == reference_inchikey:
            results["inchikey_match"] = True
            checks_passed += 1
        else:
            results["inchikey_match"] = False
            results["issues"].append(
                f"InChIKey mismatch: generated={gen_inchikey}, expected={reference_inchikey}"
            )

    if checks_total > 0:
        results["confidence"] = round(checks_passed / checks_total * 100, 1)
    else:
        results["confidence"] = None

    return results


def generate_report(mol: Chem.Mol, diagnostics: Dict) -> str:
    if mol is None:
        return "No molecule generated."

    lines = []
    lines.append("=== STRUCTURAL VALIDATION REPORT ===")
    lines.append("")
    lines.append(f"Canonical SMILES: {diagnostics.get('canonical_smiles', 'N/A')}")
    lines.append(f"Molecular Formula: {diagnostics.get('formula', 'N/A')}")
    lines.append(f"Molecular Weight: {diagnostics.get('mol_weight', 'N/A')}")
    lines.append(f"Exact Mass: {diagnostics.get('exact_mass', 'N/A')}")
    lines.append(f"InChI: {diagnostics.get('inchi', 'N/A')}")
    lines.append(f"InChIKey: {diagnostics.get('inchikey', 'N/A')}")
    lines.append("")

    lines.append("--- Composition ---")
    elem = diagnostics.get("elemental_composition", {})
    for sym in sorted(elem.keys()):
        lines.append(f"  {sym}: {elem[sym]}")
    lines.append("")

    lines.append("--- Topology ---")
    lines.append(f"  Atoms: {diagnostics.get('num_atoms', 'N/A')} (heavy: {diagnostics.get('num_heavy_atoms', 'N/A')})")
    lines.append(f"  Rings: {diagnostics.get('num_rings', 'N/A')} (aromatic: {diagnostics.get('num_aromatic_rings', 'N/A')}, aliphatic: {diagnostics.get('num_aliphatic_rings', 'N/A')})")
    lines.append(f"  Rotatable bonds: {diagnostics.get('num_rotatable_bonds', 'N/A')}")
    lines.append(f"  Stereocenters: {diagnostics.get('num_stereocenters', 'N/A')}")
    lines.append("")

    lines.append("--- Physicochemical ---")
    lines.append(f"  LogP: {diagnostics.get('logp', 'N/A')}")
    lines.append(f"  H-Bond acceptors: {diagnostics.get('num_hba', 'N/A')}")
    lines.append(f"  H-Bond donors: {diagnostics.get('num_hbd', 'N/A')}")
    lines.append("")

    lines.append("--- Functional Groups ---")
    fg = diagnostics.get("functional_groups", {})
    detected_fgs = [name for name, present in fg.items() if present]
    if detected_fgs:
        for name in detected_fgs:
            lines.append(f"  [x] {name}")
    else:
        lines.append("  (none detected)")
    lines.append("")

    return "\n".join(lines)
