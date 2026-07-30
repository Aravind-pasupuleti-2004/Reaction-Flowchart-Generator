import re
from dataclasses import dataclass, field
from typing import List, Optional, Tuple, Dict
from converter import convert_to_mol, get_molecule_info
from utils import logger


@dataclass
class ReactantInput:
    name: str
    mol: Optional[object] = None
    smiles: str = ""
    formula: str = ""
    mw: float = 0.0
    detected_type: str = ""
    conversion_status: str = "pending"
    error: str = ""
    equivalents: float = 1.0

    def is_valid(self) -> bool:
        return self.mol is not None and self.conversion_status == "success"


@dataclass
class ProductInput:
    name: str
    mol: Optional[object] = None
    smiles: str = ""
    formula: str = ""
    mw: float = 0.0
    detected_type: str = ""
    conversion_status: str = "pending"
    error: str = ""

    def is_valid(self) -> bool:
        return self.mol is not None and self.conversion_status == "success"


@dataclass
class ReactionConditions:
    solvent: str = ""
    catalyst: str = ""
    reagent: str = ""
    temperature: float = 25.0
    temperature_unit: str = "C"
    time: float = 1.0
    time_unit: str = "h"
    pressure: Optional[float] = None
    pressure_unit: str = "bar"
    ph: Optional[float] = None
    concentration: Optional[float] = None
    concentration_unit: str = "M"


@dataclass
class StageData:
    """A single synthetic stage within a multi-stage reaction sequence."""
    stage_number: int = 1
    reactants: List[ReactantInput] = field(default_factory=list)
    product: Optional[ProductInput] = None
    conditions: ReactionConditions = field(default_factory=ReactionConditions)
    validation_errors: List[str] = field(default_factory=list)

    def validate(self) -> bool:
        self.validation_errors.clear()
        if not self.reactants:
            self.validation_errors.append(f"Stage {self.stage_number}: At least one reactant is required.")
        else:
            for i, r in enumerate(self.reactants):
                if not r.name or not r.name.strip():
                    self.validation_errors.append(f"Stage {self.stage_number}, Reactant {i+1}: name is empty.")
                elif r.conversion_status == "failed":
                    self.validation_errors.append(f"Stage {self.stage_number}, Reactant {i+1} ('{r.name}'): {r.error}")
        if self.product is None or not self.product.name.strip():
            self.validation_errors.append(f"Stage {self.stage_number}: Product name is required.")
        elif self.product.conversion_status == "failed":
            self.validation_errors.append(f"Stage {self.stage_number}, Product ('{self.product.name}'): {self.product.error}")
        cond = self.conditions
        if cond.temperature < -273.15:
            self.validation_errors.append(f"Stage {self.stage_number}: Invalid temperature (below absolute zero).")
        if cond.time <= 0:
            self.validation_errors.append(f"Stage {self.stage_number}: Invalid reaction time (must be > 0).")
        if cond.ph is not None and (cond.ph < 0 or cond.ph > 14):
            self.validation_errors.append(f"Stage {self.stage_number}: Invalid pH (must be 0-14).")
        return len(self.validation_errors) == 0

    def to_export_dict(self) -> Dict:
        return {
            "stage_number": self.stage_number,
            "reactants": [
                {
                    "name": r.name,
                    "smiles": r.smiles or "",
                    "formula": r.formula or "",
                    "mw": r.mw,
                    "mol": r.mol,
                }
                for r in self.reactants if r.is_valid()
            ],
            "product": {
                "name": self.product.name,
                "smiles": self.product.smiles or "",
                "formula": self.product.formula or "",
                "mw": self.product.mw,
                "mol": self.product.mol,
            } if self.product and self.product.is_valid() else None,
            "conditions": {
                "solvent": self.conditions.solvent or "",
                "catalyst": self.conditions.catalyst or "",
                "reagent": self.conditions.reagent or "",
                "temperature": self.conditions.temperature,
                "time": self.conditions.time,
                "pressure": self.conditions.pressure,
                "ph": self.conditions.ph,
                "concentration": self.conditions.concentration,
            },
        }


@dataclass
class ReactionData:
    reactants: List[ReactantInput] = field(default_factory=list)
    product: Optional[ProductInput] = None
    conditions: ReactionConditions = field(default_factory=ReactionConditions)
    validation_errors: List[str] = field(default_factory=list)
    is_valid: bool = False

    def validate(self) -> bool:
        self.validation_errors.clear()

        if not self.reactants:
            self.validation_errors.append("At least one reactant is required.")
        else:
            for i, r in enumerate(self.reactants):
                if not r.name or not r.name.strip():
                    self.validation_errors.append(f"Reactant {i+1} name is empty.")
                elif r.conversion_status == "failed":
                    self.validation_errors.append(f"Reactant {i+1} ('{r.name}'): {r.error}")

        if self.product is None or not self.product.name.strip():
            self.validation_errors.append("Product name is required.")
        elif self.product.conversion_status == "failed":
            self.validation_errors.append(f"Product ('{self.product.name}'): {self.product.error}")

        cond = self.conditions
        if cond.temperature < -273.15:
            self.validation_errors.append(f"Invalid temperature: {cond.temperature} (below absolute zero).")
        if cond.time <= 0:
            self.validation_errors.append(f"Invalid reaction time: {cond.time} (must be > 0).")
        if cond.ph is not None and (cond.ph < 0 or cond.ph > 14):
            self.validation_errors.append(f"Invalid pH: {cond.ph} (must be 0-14).")

        self.is_valid = len(self.validation_errors) == 0
        return self.is_valid


def process_reactant(name: str) -> ReactantInput:
    entry = ReactantInput(name=name)
    if not name or not name.strip():
        entry.conversion_status = "failed"
        entry.error = "Empty name"
        return entry

    try:
        mol, dtype, status = convert_to_mol(name)
        if status == "success" and mol is not None:
            info = get_molecule_info(mol, name, dtype)
            entry.mol = mol
            entry.smiles = info.get("canonical_smiles", "")
            entry.formula = info.get("molecular_formula", "")
            entry.mw = info.get("molecular_weight", 0)
            entry.detected_type = dtype
            entry.conversion_status = "success"
        else:
            entry.conversion_status = "failed"
            entry.error = f"Cannot convert: '{name}'"
    except Exception as exc:
        entry.conversion_status = "failed"
        entry.error = str(exc)

    return entry


def process_product(name: str) -> ProductInput:
    entry = ProductInput(name=name)
    if not name or not name.strip():
        entry.conversion_status = "failed"
        entry.error = "Empty name"
        return entry

    try:
        mol, dtype, status = convert_to_mol(name)
        if status == "success" and mol is not None:
            info = get_molecule_info(mol, name, dtype)
            entry.mol = mol
            entry.smiles = info.get("canonical_smiles", "")
            entry.formula = info.get("molecular_formula", "")
            entry.mw = info.get("molecular_weight", 0)
            entry.detected_type = dtype
            entry.conversion_status = "success"
        else:
            entry.conversion_status = "failed"
            entry.error = f"Cannot convert: '{name}'"
    except Exception as exc:
        entry.conversion_status = "failed"
        entry.error = str(exc)

    return entry
