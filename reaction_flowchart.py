from typing import List, Dict, Optional
import base64
import io


def _strip_indent(s: str) -> str:
    return "\n".join(
        line.lstrip() if line.strip() else ""
        for line in s.split("\n")
    )


def _mol_to_b64(mol) -> Optional[str]:
    from rdkit.Chem import Draw
    try:
        img = Draw.MolToImage(mol, size=(280, 200))
        buf = io.BytesIO()
        img.save(buf, format="PNG")
        return base64.b64encode(buf.getvalue()).decode("utf-8")
    except Exception:
        return None


def build_reaction_scheme_html(
    reactants: List[Dict],
    product: Optional[Dict],
    conditions: Dict,
) -> str:
    r"""Generate HTML for a chemistry journal-style reaction flowchart.

    Parameters
    ----------
    reactants : list of dict
        Each with keys: name, smiles, formula, mol (RDKit Mol object optional)
    product : dict or None
        Keys: name, smiles, formula, mol (optional)
    conditions : dict
        Keys: solvent, catalyst, reagent, temperature, time, pressure, ph, concentration
    yield_percent : float, optional
    confidence : float, optional
    """
    cond_lines = []
    if conditions.get("solvent"):
        cond_lines.append(conditions["solvent"])
    if conditions.get("catalyst"):
        cond_lines.append(conditions["catalyst"])
    if conditions.get("reagent"):
        cond_lines.append(conditions["reagent"])
    temp = conditions.get("temperature")
    time_val = conditions.get("time")
    cond_parts = []
    if temp is not None:
        cond_parts.append(f"{temp}°C")
    if time_val is not None:
        cond_parts.append(f"{time_val} h")
    if cond_parts:
        cond_lines.append(", ".join(cond_parts))
    extra = []
    if conditions.get("pressure"):
        extra.append(f"P = {conditions['pressure']} bar")
    if conditions.get("ph"):
        extra.append(f"pH = {conditions['ph']}")
    if conditions.get("concentration"):
        extra.append(f"[{conditions['concentration']} M]")
    cond_text = "<br>".join(cond_lines) if cond_lines else ""
    extra_text = "<br>".join(extra) if extra else ""

    def _render_molecule(mol_data: Dict, label: str, width: int = 280, height: int = 200) -> str:
        mol = mol_data.get("mol")
        b64 = _mol_to_b64(mol) if mol else None
        name = mol_data.get("name", label)
        smi = mol_data.get("smiles", "")
        formula = mol_data.get("formula", "")
        if b64:
            return f"""\
            <div class="rxn-mol" style="width:{width}px;">
              <div class="rxn-mol-label">{label}</div>
              <img src="data:image/png;base64,{b64}" style="width:{width}px;height:{height}px;object-fit:contain;">
              <div class="rxn-smiles">{name}</div>
              <div class="rxn-smiles" style="color:#6b7280;font-size:11px;">{formula}</div>
            </div>"""
        return f"""\
        <div class="rxn-mol" style="width:{width}px;">
          <div class="rxn-mol-label">{label}</div>
          <div class="rxn-noimg">No structure</div>
          <div class="rxn-smiles">{name}</div>
          <div class="rxn-smiles" style="color:#6b7280;font-size:11px;">{formula}</div>
        </div>"""

    rxn_html = '<div class="rxn-row">'
    for i, r in enumerate(reactants):
        rxn_html += _render_molecule(r, f"Reactant {i+1}")
        if i < len(reactants) - 1:
            rxn_html += '<div class="rxn-plus">+</div>'
    rxn_html += "</div>"

    arrow_html = f"""\
    <div class="rxn-arrow-section">
      <div class="rxn-conditions">
    {cond_text}
    {("<br>" + extra_text) if extra_text and cond_text else extra_text}
      </div>
      <div class="rxn-arrow-body">
        <div class="rxn-arrow-line"></div>
        <div class="rxn-arrowhead">▶</div>
      </div>
    </div>"""

    prod_html = '<div class="rxn-row">'
    if product:
        prod_html += _render_molecule(product, "Product")
    prod_html += "</div>"

    html = f"""\
<div class="rxn-container">
  <div class="rxn-title">Reaction Process Flowchart</div>
  <div class="rxn-scheme">
{rxn_html}
{arrow_html}
{prod_html}
  </div>
</div>

<style>
  .rxn-container {{
    background: #111827;
    border: 1px solid #1f2937;
    border-radius: 12px;
    padding: 24px;
    margin: 16px 0;
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
  }}
  .rxn-title {{
    font-size: 18px;
    font-weight: 600;
    color: #f3f4f6;
    margin-bottom: 20px;
    padding-bottom: 12px;
    border-bottom: 1px solid #1f2937;
  }}
  .rxn-scheme {{
    display: flex;
    flex-direction: column;
    align-items: center;
    gap: 6px;
  }}
  .rxn-row {{
    display: flex;
    align-items: center;
    justify-content: center;
    gap: 16px;
    flex-wrap: wrap;
  }}
  .rxn-mol {{
    background: #1a1a2e;
    border: 1px solid #2d2d5e;
    border-radius: 10px;
    padding: 12px;
    text-align: center;
    transition: all 0.2s ease;
  }}
  .rxn-mol:hover {{
    border-color: #3b82f6;
    box-shadow: 0 0 12px rgba(59,130,246,0.2);
  }}
  .rxn-mol-label {{
    font-size: 12px;
    font-weight: 600;
    color: #9ca3af;
    text-transform: uppercase;
    letter-spacing: 0.5px;
    margin-bottom: 6px;
  }}
  .rxn-noimg {{
    width: 280px;
    height: 200px;
    display: flex;
    align-items: center;
    justify-content: center;
    color: #4b5563;
    font-size: 13px;
    border: 1px dashed #374151;
    border-radius: 6px;
  }}
  .rxn-smiles {{
    font-size: 12px;
    color: #d1d5db;
    margin-top: 6px;
    word-break: break-all;
  }}
  .rxn-plus {{
    font-size: 28px;
    font-weight: 300;
    color: #6b7280;
    padding: 0 8px;
  }}
  .rxn-arrow-section {{
    display: flex;
    flex-direction: column;
    align-items: center;
    padding: 8px 0;
    width: 100%;
    max-width: 600px;
  }}
  .rxn-conditions {{
    background: #1f2937;
    border: 1px solid #374151;
    border-radius: 8px;
    padding: 10px 20px;
    font-size: 13px;
    color: #e5e7eb;
    text-align: center;
    line-height: 1.6;
    margin-bottom: 6px;
    min-width: 200px;
  }}
  .rxn-arrow-body {{
    display: flex;
    align-items: center;
    width: 80%;
    max-width: 400px;
  }}
  .rxn-arrow-line {{
    flex: 1;
    height: 2px;
    background: #4b5563;
  }}
  .rxn-arrowhead {{
    font-size: 18px;
    color: #6b7280;
    margin-left: -2px;
  }}
</style>
"""

    return _strip_indent(html)
