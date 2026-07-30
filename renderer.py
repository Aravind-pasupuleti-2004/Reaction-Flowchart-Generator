import io
import base64
from typing import Optional
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D


def _prepare_mol(mol: Chem.Mol) -> Chem.Mol:
    try:
        from rdkit.Chem import rdDepictor
        rdDepictor.SetPreferCoordGen(True)
    except Exception:
        pass
    return mol


def mol_to_png_base64(mol: Chem.Mol, size: tuple = (600, 400)) -> Optional[str]:
    if mol is None:
        return None
    try:
        mol = _prepare_mol(mol)
        img = Draw.MolToImage(mol, size=size)
        buf = io.BytesIO()
        img.save(buf, format="PNG")
        buf.seek(0)
        return base64.b64encode(buf.read()).decode("utf-8")
    except Exception:
        return None


def mol_to_png_bytes(mol: Chem.Mol, size: tuple = (600, 400)) -> Optional[bytes]:
    if mol is None:
        return None
    try:
        mol = _prepare_mol(mol)
        img = Draw.MolToImage(mol, size=size)
        buf = io.BytesIO()
        img.save(buf, format="PNG")
        buf.seek(0)
        return buf.read()
    except Exception:
        return None


def mol_to_svg(mol: Chem.Mol, size: tuple = (600, 400)) -> Optional[str]:
    if mol is None:
        return None
    try:
        mol = _prepare_mol(mol)
        drawer = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
        opts = drawer.drawOptions()
        opts.legendFontSize = 14
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()
    except Exception:
        return None


def mol_to_molfile(mol: Chem.Mol) -> Optional[str]:
    if mol is None:
        return None
    try:
        return Chem.MolToMolBlock(mol)
    except Exception:
        return None


def save_png(mol: Chem.Mol, filepath: str, size: tuple = (600, 400)) -> bool:
    try:
        data = mol_to_png_bytes(mol, size)
        if data:
            with open(filepath, "wb") as f:
                f.write(data)
            return True
        return False
    except Exception:
        return False
