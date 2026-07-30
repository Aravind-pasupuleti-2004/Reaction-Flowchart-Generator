import sys
import os
import streamlit as st
from rdkit import Chem

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from renderer import mol_to_png_base64, mol_to_svg, mol_to_molfile
from converter import get_molecule_info
from reaction_form import (
    ReactionData, ReactionConditions, ReactantInput, ProductInput, StageData,
    process_reactant, process_product,
)
from validation import compute_diagnostics
from workflow_chart import build_workflow_html, make_steps_from_status, DEFAULT_STEPS
from reaction_flowchart import build_reaction_scheme_html
from export_flowchart import (
    export_svg_bytes, export_pdf_bytes, export_png_bytes,
    export_docx_bytes, export_structure_svg_bytes, export_structure_pdf_bytes,
    export_structure_png_bytes, export_structure_docx_bytes,
    build_reaction_scheme_svg,
    export_multi_stage_svg_bytes, export_multi_stage_pdf_bytes,
    export_multi_stage_png_bytes, export_multi_stage_docx_bytes,
)

st.set_page_config(page_title="Reaction Flowchart Generator", page_icon="⚗️", layout="wide")
st.title("Reaction Flowchart Generator")
st.markdown("Enter reaction details to generate professional flowcharts and Word documents.")


def render_molecule_card(mol, label: str, info: dict, cols):
    with cols[0]:
        png = mol_to_png_base64(mol, size=(400, 250))
        if png:
            st.image(f"data:image/png;base64,{png}", use_container_width=True)
    with cols[1]:
        st.markdown(f"**{label}**")
        st.code(info.get("canonical_smiles", "N/A"), language="text")
        props = [
            ("Formula", info.get("molecular_formula", "N/A")),
            ("MW", f"{info.get('molecular_weight', 'N/A')} g/mol"),
            ("InChIKey", info.get("inchikey", "N/A")),
        ]
        for k, v in props:
            st.caption(f"{k}: {v}")


# ============================================================
# SESSION STATE
# ============================================================
if "num_stages" not in st.session_state:
    st.session_state.num_stages = 1
if "show_results" not in st.session_state:
    st.session_state.show_results = False
if "stage_data" not in st.session_state:
    st.session_state.stage_data = []


# ============================================================
# INPUT
# ============================================================
if not st.session_state.show_results:
    st.number_input(
        "Number of Stages",
        min_value=1, max_value=10, value=st.session_state.num_stages,
        key="num_stages_input",
        on_change=lambda: (
            setattr(st.session_state, "num_stages", st.session_state.num_stages_input),
            None,
        ),
    )

for i in range(st.session_state.num_stages):
    # ensure stage_data list length matches
    while len(st.session_state.stage_data) <= i:
        st.session_state.stage_data.append(None)

    already_predicted = st.session_state.stage_data[i] is not None
    expanded = not already_predicted

    with st.expander(f"Stage {i + 1}", expanded=expanded):
        if already_predicted:
            sd: StageData = st.session_state.stage_data[i]
            rnames = " + ".join(r.name for r in sd.reactants if r.is_valid())
            st.success(f"**Processed** — {rnames} → {sd.product.name if sd.product else '?'}")
            if st.button(f"Edit Stage {i + 1}", key=f"edit_stage_{i}"):
                st.session_state.stage_data[i] = None
                st.rerun()
            continue

        with st.form(f"stage_form_{i}"):
            st.subheader(f"Reactants — Stage {i + 1}")

            reactant_names = []
            reactant_count_key = f"reactant_count_{i}"
            if reactant_count_key not in st.session_state:
                st.session_state[reactant_count_key] = 1

            for j in range(st.session_state[reactant_count_key]):
                col_r1, col_r2 = st.columns([4, 1])
                with col_r1:
                    name = st.text_input(
                        f"Reactant {j + 1}",
                        key=f"stage_{i}_reactant_{j}",
                        placeholder="IUPAC name or SMILES",
                    )
                    reactant_names.append(name)
                with col_r2:
                    eq = st.number_input(
                        "Eq", min_value=0.0, value=1.0, step=0.1,
                        key=f"stage_{i}_eq_{j}",
                    )

            col_add, col_rem = st.columns([1, 5])
            with col_add:
                if st.form_submit_button("+ Add Reactant", type="secondary"):
                    st.session_state[reactant_count_key] += 1
                    st.rerun()
            with col_rem:
                if st.session_state[reactant_count_key] > 1 and st.form_submit_button("− Remove", type="secondary"):
                    st.session_state[reactant_count_key] -= 1
                    st.rerun()

            st.divider()
            st.subheader("Product")
            product_name = st.text_input(
                "Product", placeholder="IUPAC name or SMILES",
                key=f"stage_{i}_product",
            )

            st.divider()
            st.subheader("Reaction Conditions")
            col_c1, col_c2 = st.columns(2)
            with col_c1:
                solvent = st.text_input("Solvent(s)", placeholder="e.g. ethanol, water", key=f"stage_{i}_solvent")
                catalyst = st.text_input("Catalyst(s)", placeholder="e.g. Pd/C, H2SO4", key=f"stage_{i}_catalyst")
                reagent = st.text_input("Reagent(s)", placeholder="e.g. NaBH4", key=f"stage_{i}_reagent")
            with col_c2:
                temp = st.number_input("Temperature (°C)", value=25.0, step=5.0, key=f"stage_{i}_temp")
                time_val = st.number_input("Time (hours)", value=1.0, min_value=0.1, step=0.5, key=f"stage_{i}_time")
                col_p1, col_p2, col_p3 = st.columns(3)
                with col_p1:
                    pressure = st.text_input("Pressure (bar) (opt)", value="", key=f"stage_{i}_pressure")
                with col_p2:
                    ph = st.text_input("pH (opt)", value="", key=f"stage_{i}_ph")
                with col_p3:
                    conc = st.text_input("Conc. (M) (opt)", value="", key=f"stage_{i}_conc")

            st.divider()
            submitted = st.form_submit_button(
                f"Process Stage {i + 1}", type="primary", use_container_width=True
            )

            if submitted:
                errors = []
                if not any(r.strip() for r in reactant_names if r.strip()):
                    errors.append(f"Stage {i + 1}: At least one reactant is required.")
                if not product_name.strip():
                    errors.append(f"Stage {i + 1}: Product is required.")
                if errors:
                    for e in errors:
                        st.error(e)
                    st.stop()

                with st.spinner(f"Converting IUPAC names for Stage {i + 1}..."):
                    reactants = []
                    for j, rname in enumerate(reactant_names):
                        eq = st.session_state.get(f"stage_{i}_eq_{j}", 1.0)
                        for part in rname.split("+"):
                            part = part.strip()
                            if part:
                                r = process_reactant(part)
                                r.equivalents = eq
                                reactants.append(r)
                    product = process_product(product_name.strip())

                conversion_ok = True
                for r in reactants:
                    if not r.is_valid():
                        st.error(f"Stage {i + 1}, Reactant '{r.name}': {r.error}")
                        conversion_ok = False
                if not product.is_valid():
                    st.error(f"Stage {i + 1}, Product '{product.name}': {product.error}")
                    conversion_ok = False
                if not conversion_ok:
                    st.stop()

                conditions = ReactionConditions(
                    solvent=solvent, catalyst=catalyst, reagent=reagent,
                    temperature=temp, time=time_val,
                    pressure=float(pressure) if pressure.strip() else None,
                    ph=float(ph) if ph.strip() else None,
                    concentration=float(conc) if conc.strip() else None,
                )

                stage_obj = StageData(
                    stage_number=i + 1,
                    reactants=reactants,
                    product=product,
                    conditions=conditions,
                )

                if not stage_obj.validate():
                    for e in stage_obj.validation_errors:
                        st.error(e)
                    st.stop()

                st.session_state.stage_data[i] = stage_obj
                st.session_state.show_results = True
                st.rerun()


# ============================================================
# RESULTS
# ============================================================
completed_stages = [sd for sd in st.session_state.stage_data if sd is not None]
if completed_stages:
    st.divider()
    st.header("Results")

    # ── per-stage results ──
    for idx, sd in enumerate(completed_stages):
        with st.container():
            st.subheader(f"Stage {sd.stage_number} of {len(completed_stages)}")

            rnames = " + ".join(r.name for r in sd.reactants if r.is_valid())
            st.markdown(f"**Reactants:** {rnames}")
            st.markdown(f"**Product:** {sd.product.name if sd.product else '?'}")
            c = sd.conditions
            if c.solvent:
                st.markdown(f"**Solvent:** {c.solvent}")
            if c.catalyst:
                st.markdown(f"**Catalyst:** {c.catalyst}")
            if c.reagent:
                st.markdown(f"**Reagent:** {c.reagent}")
            st.markdown(f"**Temperature:** {c.temperature} °C")
            st.markdown(f"**Time:** {c.time} h")

            # molecular structures for this stage
            st.markdown("**Reactants**")
            for r in sd.reactants:
                if r.is_valid():
                    info = get_molecule_info(r.mol, r.name, r.detected_type)
                    with st.container():
                        cols = st.columns([1, 1])
                        render_molecule_card(r.mol, f"Reactant: {r.name}", info, cols)
                        st.caption(f"Eq: {r.equivalents}")
                        dcol = st.columns([1, 1, 4])
                        with dcol[0]:
                            st.download_button(
                                label="SVG",
                                data=export_structure_svg_bytes(r.mol, r.name, info),
                                file_name=f"stage_{sd.stage_number}_{r.name}.svg",
                                mime="image/svg+xml",
                                use_container_width=True,
                            )
                        with dcol[1]:
                            st.download_button(
                                label="Word",
                                data=export_structure_docx_bytes(r.mol, r.name, info),
                                file_name=f"stage_{sd.stage_number}_{r.name}.docx",
                                mime="application/vnd.openxmlformats-officedocument.wordprocessingml.document",
                                use_container_width=True,
                            )

            st.markdown("**Product**")
            if sd.product and sd.product.is_valid():
                info = get_molecule_info(sd.product.mol, sd.product.name, sd.product.detected_type)
                with st.container():
                    cols = st.columns([1, 1])
                    render_molecule_card(sd.product.mol, f"Product: {sd.product.name}", info, cols)
                    dcol = st.columns([1, 1, 4])
                    with dcol[0]:
                        st.download_button(
                            label="SVG",
                            data=export_structure_svg_bytes(sd.product.mol, sd.product.name, info),
                            file_name=f"stage_{sd.stage_number}_product.svg",
                            mime="image/svg+xml",
                            use_container_width=True,
                        )
                    with dcol[1]:
                        st.download_button(
                            label="Word",
                            data=export_structure_docx_bytes(sd.product.mol, sd.product.name, info),
                            file_name=f"stage_{sd.stage_number}_product.docx",
                            mime="application/vnd.openxmlformats-officedocument.wordprocessingml.document",
                            use_container_width=True,
                        )

            if idx < len(completed_stages) - 1:
                st.divider()

    # ── per-stage flowchart + yield ──
    st.divider()
    st.subheader("Reaction Process Flowchart")

    for idx, sd in enumerate(completed_stages):
        with st.container():
            st.markdown(f"**Stage {sd.stage_number}**")
            rxn_reactants = [
                {"name": r.name, "smiles": r.smiles or "", "formula": r.formula or "", "mw": r.mw, "mol": r.mol}
                for r in sd.reactants if r.is_valid()
            ]
            rxn_product = {
                "name": sd.product.name, "smiles": sd.product.smiles or "",
                "formula": sd.product.formula or "", "mw": sd.product.mw, "mol": sd.product.mol,
            } if sd.product and sd.product.is_valid() else None
            rxn_conditions = {
                "solvent": sd.conditions.solvent or "", "catalyst": sd.conditions.catalyst or "",
                "reagent": sd.conditions.reagent or "", "temperature": sd.conditions.temperature,
                "time": sd.conditions.time, "pressure": sd.conditions.pressure,
                "ph": sd.conditions.ph, "concentration": sd.conditions.concentration,
            }

            scheme_html = build_reaction_scheme_html(
                reactants=rxn_reactants, product=rxn_product, conditions=rxn_conditions,
            )
            st.markdown(scheme_html, unsafe_allow_html=True)

    # ── combined export (all stages in one document) ──
    st.divider()
    st.subheader("Export All Stages (Combined Document)")

    export_data = [sd.to_export_dict() for sd in completed_stages]
    is_multi = len(export_data) > 1

    def _safe_export(fn, *args, **kwargs):
        try:
            return fn(*args, **kwargs)
        except Exception as e:
            st.error(f"Export failed: {e}")
            return b""

    if not is_multi:
        sd = completed_stages[0]
        rxn_reactants = [
            {"name": r.name, "smiles": r.smiles or "", "formula": r.formula or "", "mw": r.mw, "mol": r.mol}
            for r in sd.reactants if r.is_valid()
        ]
        rxn_product = {
            "name": sd.product.name, "smiles": sd.product.smiles or "",
            "formula": sd.product.formula or "", "mw": sd.product.mw, "mol": sd.product.mol,
        } if sd.product and sd.product.is_valid() else None
        rxn_conditions = {
            "solvent": sd.conditions.solvent or "", "catalyst": sd.conditions.catalyst or "",
            "reagent": sd.conditions.reagent or "", "temperature": sd.conditions.temperature,
            "time": sd.conditions.time, "pressure": sd.conditions.pressure,
            "ph": sd.conditions.ph, "concentration": sd.conditions.concentration,
        }
        export_cols = st.columns(4)
        with export_cols[0]:
            st.download_button(
                label="Download PDF",
                data=_safe_export(export_pdf_bytes, rxn_reactants, rxn_product, rxn_conditions),
                file_name="reaction_flowchart.pdf", mime="application/pdf",
                use_container_width=True,
            )
        with export_cols[1]:
            st.download_button(
                label="Download PNG",
                data=_safe_export(export_png_bytes, rxn_reactants, rxn_product, rxn_conditions, dpi=300),
                file_name="reaction_flowchart.png", mime="image/png",
                use_container_width=True,
            )
        with export_cols[2]:
            st.download_button(
                label="Download SVG",
                data=_safe_export(export_svg_bytes, rxn_reactants, rxn_product, rxn_conditions),
                file_name="reaction_flowchart.svg", mime="image/svg+xml",
                use_container_width=True,
            )
        with export_cols[3]:
            st.download_button(
                label="Download Word",
                data=_safe_export(export_docx_bytes, rxn_reactants, rxn_product, rxn_conditions),
                file_name="reaction_flowchart.docx",
                mime="application/vnd.openxmlformats-officedocument.wordprocessingml.document",
                use_container_width=True,
            )
    else:
        export_cols = st.columns(4)
        with export_cols[0]:
            st.download_button(
                label="Download PDF",
                data=_safe_export(export_multi_stage_pdf_bytes, export_data),
                file_name="reaction_flowchart.pdf",
                mime="application/pdf",
                use_container_width=True,
            )
        with export_cols[1]:
            st.download_button(
                label="Download PNG",
                data=_safe_export(export_multi_stage_png_bytes, export_data, dpi=300),
                file_name="reaction_flowchart.png",
                mime="image/png",
                use_container_width=True,
            )
        with export_cols[2]:
            st.download_button(
                label="Download SVG",
                data=_safe_export(export_multi_stage_svg_bytes, export_data),
                file_name="reaction_flowchart.svg",
                mime="image/svg+xml",
                use_container_width=True,
            )
        with export_cols[3]:
            st.download_button(
                label="Download Word",
                data=_safe_export(export_multi_stage_docx_bytes, export_data),
                file_name="reaction_flowchart.docx",
                mime="application/vnd.openxmlformats-officedocument.wordprocessingml.document",
                use_container_width=True,
            )

    # ── Software Processing Pipeline (last stage) ──
    with st.expander("Software Processing Pipeline"):
        last = completed_stages[-1]
        all_compounds = [(r.name, r.mol, f"Reactant {j+1}") for j, r in enumerate(last.reactants) if r.is_valid()]
        if last.product and last.product.is_valid():
            all_compounds.append((last.product.name, last.product.mol, "Product"))

        pipe_details = {
            "inputs": f"Reactants: {' + '.join(r.name for r in last.reactants if r.is_valid())}<br>Product: {last.product.name if last.product else ''}<br>Solvent: {last.conditions.solvent or 'N/A'}<br>Catalyst: {last.conditions.catalyst or 'N/A'}<br>Temperature: {last.conditions.temperature} °C<br>Time: {last.conditions.time} h",
            "validation": "All required fields present.<br>All IUPAC names resolved to valid molecules.",
            "opsin": "<br>".join(f"{n}: {Chem.MolToSmiles(m) if m else '?'}" for n, m, _ in all_compounds if m),
            "rdkit": "<br>".join(f"{n}: {Chem.rdMolDescriptors.CalcMolFormula(m) if m else '?'}, {round(Chem.Descriptors.MolWt(m), 2) if m else 0} g/mol" for n, m, _ in all_compounds if m),
            "features": f"Feature vector dimension: 2048×2 + 10×2 + 5 = {2048*2+10*2+5}",
            "results": f"Stages: {len(completed_stages)}",
        }
        pipe_steps = make_steps_from_status(
            inputs_ok=True, validation_ok=True, opsin_ok=True,
            rdkit_ok=True, features_ok=True, results_ok=True,
            details=pipe_details,
        )
        pipe_html = build_workflow_html(pipe_steps)
        st.markdown(pipe_html, unsafe_allow_html=True)

    if st.button("Start New Process", type="secondary"):
        st.session_state.show_results = False
        st.session_state.stage_data = []
        st.session_state.num_stages = 1
        st.rerun()

st.markdown("---")
st.caption("Powered by RDKit, OPSIN, and scikit-learn. Fully offline, no cloud services required.")
