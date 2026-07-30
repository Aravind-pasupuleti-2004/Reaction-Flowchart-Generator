from typing import List, Dict, Optional


def _strip_indent(s: str) -> str:
    return "\n".join(
        line.lstrip() if line.strip() else ""
        for line in s.split("\n")
    )


STATUS_COLORS = {
    "pending": {"bg": "#1a1a2e", "border": "#2d2d5e", "text": "#6b7280", "glow": "none"},
    "running": {"bg": "#0f2b3d", "border": "#2563eb", "text": "#60a5fa", "glow": "0 0 12px rgba(37,99,235,0.4)"},
    "completed": {"bg": "#0a2e1a", "border": "#16a34a", "text": "#86efac", "glow": "0 0 12px rgba(22,163,74,0.3)"},
    "failed": {"bg": "#2e0a0a", "border": "#dc2626", "text": "#fca5a5", "glow": "0 0 12px rgba(220,38,38,0.3)"},
}

STATUS_ICONS = {
    "pending": "○",
    "running": "◎",
    "completed": "●",
    "failed": "●",
}

NODE_ICONS = {
    "inputs": "📝",
    "validation": "✅",
    "opsin": "🧪",
    "rdkit": "⚛️",
    "features": "⚙️",
    "results": "📊",
}

DEFAULT_STEPS = [
    {"id": "inputs", "label": "User Inputs", "detail": "Reactant(s), Product, Solvent, Catalyst, Temperature, Time"},
    {"id": "validation", "label": "Input Validation", "detail": "Required fields, IUPAC conversion, numeric bounds"},
    {"id": "opsin", "label": "IUPAC → Molecular Structure (OPSIN)", "detail": ""},
    {"id": "rdkit", "label": "RDKit Processing", "detail": "Canonical SMILES, InChI, Formula, MW"},
    {"id": "features", "label": "Feature Engineering", "detail": "Morgan fingerprints, molecular descriptors, condition vector"},
    {"id": "results", "label": "Prediction Results", "detail": "Molecular Structures, Properties"},
]


def build_workflow_html(
    steps: List[Dict],
    num_structures: int = 0,
    status: str = "pending",
    expanded_node: Optional[str] = None,
) -> str:
    nodes_html = ""
    for i, step in enumerate(steps):
        sid = step["id"]
        label = step["label"]
        detail = step.get("detail", "")
        step_status = step.get("status", "pending")
        colors = STATUS_COLORS.get(step_status, STATUS_COLORS["pending"])
        icon = NODE_ICONS.get(sid, "•")
        status_icon = STATUS_ICONS.get(step_status, "○")
        is_expanded = expanded_node == sid

        detail_html = ""
        if is_expanded and detail:
            detail_html = f"""\
            <div class="node-detail">
              <div class="detail-arrow"></div>
              <div class="detail-content">{detail}</div>
            </div>"""

        nodes_html += f"""\
        <div class="workflow-node" id="node-{sid}">
          <div class="node-box" style="background:{colors['bg']};border-color:{colors['border']};box-shadow:{colors['glow']}">
            <div class="node-status" style="color:{colors['text']}">{status_icon}</div>
            <div class="node-icon">{icon}</div>
            <div class="node-label" style="color:{colors['text']}">{label}</div>
          </div>
          {detail_html}
        </div>"""

        if i < len(steps) - 1:
            nodes_html += f"""\
            <div class="workflow-arrow">
              <div class="arrow-shaft"></div>
              <div class="arrow-head">▼</div>
            </div>"""

    html = f"""\
<div class="workflow-container">
  <div class="workflow-title">Processing Pipeline</div>
  <div class="workflow-pipeline">
{nodes_html}
  </div>
</div>

<style>
  .workflow-container {{
    background: #111827;
    border: 1px solid #1f2937;
    border-radius: 12px;
    padding: 24px;
    margin: 16px 0;
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
  }}
  .workflow-title {{
    font-size: 18px;
    font-weight: 600;
    color: #f3f4f6;
    margin-bottom: 20px;
    padding-bottom: 12px;
    border-bottom: 1px solid #1f2937;
  }}
  .workflow-pipeline {{
    display: flex;
    flex-direction: column;
    align-items: center;
    gap: 0;
  }}
  .workflow-node {{
    width: 100%;
    max-width: 480px;
  }}
  .node-box {{
    display: flex;
    align-items: center;
    gap: 12px;
    padding: 14px 18px;
    border-radius: 10px;
    border: 1.5px solid;
    transition: all 0.3s ease;
    cursor: default;
  }}
  .node-box:hover {{
    filter: brightness(1.2);
    transform: translateX(4px);
  }}
  .node-status {{
    font-size: 12px;
    width: 16px;
    text-align: center;
    flex-shrink: 0;
  }}
  .node-icon {{
    font-size: 20px;
    flex-shrink: 0;
  }}
  .node-label {{
    font-size: 14px;
    font-weight: 500;
    line-height: 1.3;
  }}
  .workflow-arrow {{
    display: flex;
    flex-direction: column;
    align-items: center;
    padding: 2px 0;
    color: #374151;
  }}
  .arrow-shaft {{
    width: 2px;
    height: 16px;
    background: #374151;
  }}
  .arrow-head {{
    font-size: 10px;
    line-height: 1;
    color: #4b5563;
  }}
  .node-detail {{
    margin-top: 6px;
    margin-left: 24px;
    animation: fadeIn 0.2s ease;
  }}
  .detail-arrow {{
    width: 0;
    height: 0;
    border-left: 6px solid transparent;
    border-right: 6px solid transparent;
    border-bottom: 6px solid #1f2937;
    margin-left: 12px;
  }}
  .detail-content {{
    background: #1f2937;
    border: 1px solid #374151;
    border-radius: 8px;
    padding: 10px 14px;
    font-size: 12px;
    color: #d1d5db;
    line-height: 1.5;
    font-family: 'SF Mono', 'Fira Code', monospace;
  }}
  @keyframes fadeIn {{
    from {{ opacity: 0; transform: translateY(-4px); }}
    to {{ opacity: 1; transform: translateY(0); }}
  }}
</style>
"""

    return _strip_indent(html)


def make_steps_from_status(
    inputs_ok: bool = False,
    validation_ok: bool = False,
    opsin_ok: bool = False,
    rdkit_ok: bool = False,
    features_ok: bool = False,
    results_ok: bool = False,
    details: Optional[Dict[str, str]] = None,
) -> List[Dict]:
    if details is None:
        details = {}

    status_order = []
    for ok in [inputs_ok, validation_ok, opsin_ok, rdkit_ok, features_ok, results_ok]:
        status_order.append(ok)

    def _status(idx: int) -> str:
        if idx >= len(status_order):
            return "pending"
        return "completed" if status_order[idx] else "pending"

    steps = DEFAULT_STEPS.copy()
    for i, step in enumerate(steps):
        step["status"] = _status(i)
        sid = step["id"]
        if sid in details:
            step["detail"] = details[sid]

    return steps
