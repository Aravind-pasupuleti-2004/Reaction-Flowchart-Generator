import numpy as np
from typing import List, Optional, Tuple
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from reaction_form import ReactionData, ReactantInput, ProductInput


MORGAN_RADIUS = 2
MORGAN_BITS = 2048

CONDITION_FEATURE_KEYS = [
    "temperature", "time", "pressure", "ph", "concentration"
]


def _morgan_fingerprint(mol: Chem.Mol, radius: int = MORGAN_RADIUS, nbits: int = MORGAN_BITS) -> np.ndarray:
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)
    arr = np.zeros((nbits,), dtype=np.float32)
    for i in fp.GetOnBits():
        arr[i] = 1.0
    return arr


def _compute_molecular_descriptors(mol: Chem.Mol) -> np.ndarray:
    desc = [
        Descriptors.MolWt(mol),
        Descriptors.MolLogP(mol),
        rdMolDescriptors.CalcNumHBA(mol),
        rdMolDescriptors.CalcNumHBD(mol),
        rdMolDescriptors.CalcNumRotatableBonds(mol),
        rdMolDescriptors.CalcNumRings(mol),
        rdMolDescriptors.CalcNumAromaticRings(mol),
        mol.GetNumAtoms(),
        mol.GetNumHeavyAtoms(),
        Descriptors.TPSA(mol),
    ]
    return np.array(desc, dtype=np.float32)


def prepare_feature_vector(reaction: ReactionData) -> Optional[np.ndarray]:
    valid_reactants = [r for r in reaction.reactants if r.is_valid()]
    if not valid_reactants or not reaction.product or not reaction.product.is_valid():
        return None

    reactant_fp = np.zeros((MORGAN_BITS,), dtype=np.float32)
    reactant_descs = []
    for r in valid_reactants:
        fp = _morgan_fingerprint(r.mol)
        reactant_fp = np.maximum(reactant_fp, fp)
        reactant_descs.append(_compute_molecular_descriptors(r.mol))

    if reactant_descs:
        avg_reactant_desc = np.mean(reactant_descs, axis=0)
    else:
        avg_reactant_desc = np.zeros(10, dtype=np.float32)

    product_fp = _morgan_fingerprint(reaction.product.mol)
    product_desc = _compute_molecular_descriptors(reaction.product.mol)

    cond = reaction.conditions
    condition_vec = np.array([
        cond.temperature if cond.temperature_unit == "C" else cond.temperature * 1.8 + 32,
        cond.time if cond.time_unit == "h" else cond.time / 60,
        cond.pressure or 1.0,
        cond.ph or 7.0,
        cond.concentration or 0.1,
    ], dtype=np.float32)

    feature_vector = np.concatenate([
        reactant_fp,
        product_fp,
        avg_reactant_desc,
        product_desc,
        condition_vec,
    ])

    return feature_vector


def get_feature_dimension() -> int:
    return MORGAN_BITS * 2 + 10 * 2 + len(CONDITION_FEATURE_KEYS)


def get_feature_names() -> List[str]:
    names = []
    for i in range(MORGAN_BITS):
        names.append(f"reactant_fp_{i}")
    for i in range(MORGAN_BITS):
        names.append(f"product_fp_{i}")
    desc_names = ["MolWt", "MolLogP", "HBA", "HBD", "RotBonds", "Rings", "AroRings", "Atoms", "HeavyAtoms", "TPSA"]
    for n in desc_names:
        names.append(f"avg_reactant_{n}")
    for n in desc_names:
        names.append(f"product_{n}")
    names.extend(CONDITION_FEATURE_KEYS)
    return names
