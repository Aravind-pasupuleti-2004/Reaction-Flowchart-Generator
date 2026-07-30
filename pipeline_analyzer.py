import re
import json
import ssl
import urllib.request
import urllib.parse
import urllib.error
from typing import Optional, Dict, List, Tuple, Any
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem
from converter import convert_to_mol, get_molecule_info
from renderer import mol_to_png_base64, mol_to_svg
from utils import logger

_OPSIN_ENGINE = None


def _get_opsin_engine():
    global _OPSIN_ENGINE
    if _OPSIN_ENGINE is not None:
        return _OPSIN_ENGINE
    try:
        import pyopsin
        _OPSIN_ENGINE = pyopsin.PyOpsin()
        return _OPSIN_ENGINE
    except Exception:
        return None


# ─────────────────────────────────────────────
# STAGE 1: Preprocessing / Normalization
# ─────────────────────────────────────────────
def stage1_preprocess(name: str) -> Dict:
    result = {
        "stage": "1 - Preprocessing / Normalization",
        "input": name,
        "output": "",
        "warnings": [],
        "errors": [],
        "passed": False,
    }

    if not name or not name.strip():
        result["errors"].append("Empty input")
        return result

    cleaned = name.strip()

    cleaned = cleaned.replace("\u2019", "'").replace("\u2018", "'")
    cleaned = cleaned.replace("\u201c", '"').replace("\u201d", '"')

    cleaned = re.sub(r"\s+", " ", cleaned)

    result["output"] = cleaned
    result["passed"] = True
    if cleaned != name.strip():
        result["warnings"].append(f"Normalized whitespace/unicode: '{name.strip()}' -> '{cleaned}'")

    return result


# ─────────────────────────────────────────────
# STAGE 2: IUPAC Name Scanning
# ─────────────────────────────────────────────
_DIFFICULT_PATTERNS = [
    ("stereochemistry", r"\b[(\[]\d*[RSEZ][,\s)\]]|\bcis\b|\btrans\b", False),
    ("fused_ring", r"\[[\d',]+:[\d',]+:[\d',]+\]", False),
    ("bridgehead_locant", r"\b\d+[a-g]\b", False),
    ("spiro", r"\bspiro\b", False),
    ("bicyclo", r"\bbicyclo\b", False),
    ("tricyclo", r"\btricyclo\b", False),
    ("heterocycle_prefix", r"\b(imidazo|pyrido|pyrazino|pyrimidino|pyrrolo|furo|thieno|oxazolo|thiazolo|piperidino|piperazino|morpholino)\b", False),
    ("multicomponent", r"\b(?:hydrochloride|hydrate|monohydrate|dihydrate|hemihydrate)\b", False, "warning"),
    ("composition", r"\bcompd?\b|with\b.*\b(?:and|or)\b", False, "warning"),
]

_DIFFICULT_STEREO_ADJACENT = r"[(\[]\d*[RSEZ][,\s)\]]"

_KNOWN_FUSED_RINGS = [
    r"pyrido\[[\d',]+:[\d',]+:[\d',]+\]pyrazino",
    r"pyrido\[[\d',]+:[\d',]+:[\d',]+\]",
    r"pyrazino\[[\d',]+:[\d',]+:[\d',]+\]",
    r"imidazo\[[\d',]+:[\d',]+:[\d',]+\]",
    r"pyrrolo\[[\d',]+:[\d',]+:[\d',]+\]",
    r"furo\[[\d',]+:[\d',]+:[\d',]+\]",
    r"thieno\[[\d',]+:[\d',]+:[\d',]+\]",
    r"oxazolo\[[\d',]+:[\d',]+:[\d',]+\]",
    r"thiazolo\[[\d',]+:[\d',]+:[\d',]+\]",
    r"\[[\d',]+:[\d',]+:[\d',]+\]",
]


def stage2_scan_iupac(name: str) -> Dict:
    result = {
        "stage": "2 - IUPAC Name Difficulty Scan",
        "input": name,
        "output": name,
        "warnings": [],
        "errors": [],
        "passed": True,
        "difficulty_score": 0,
        "difficulty_level": "simple",
        "detected_constructs": [],
        "stereochemistry_found": False,
        "fused_rings_found": False,
        "bridgehead_found": False,
        "skipped_tokens": [],
    }

    if re.match(r"^InChI=|[BJKLOPQRSTUVY]\d{6,}", name):
        result["warnings"].append("Input is InChI/InChIKey or SMILES, not an IUPAC name. Skipping IUPAC scanning.")
        result["output"] = name
        return result

    if re.match(r"^[A-Za-z0-9@+\-\[\]()\\\/%#$.=,:;]+$", name) and len(name) < 200:
        mol_check = Chem.MolFromSmiles(name)
        if mol_check is not None:
            result["warnings"].append("Input appears to be SMILES, not an IUPAC name. Skipping IUPAC scanning.")
            result["output"] = name
            return result

    for entry in _DIFFICULT_PATTERNS:
        construct_name = entry[0]
        pattern = entry[1]
        case_insensitive = entry[2] if len(entry) > 2 else True
        severity = entry[3] if len(entry) > 3 else "error"
        flags = re.IGNORECASE if case_insensitive else 0
        matches = re.findall(pattern, name, flags)
        if matches:
            result["detected_constructs"].append({"construct": construct_name, "matches": matches, "severity": severity})
            result["difficulty_score"] += len(matches)

            if construct_name == "stereochemistry":
                result["stereochemistry_found"] = True
            if construct_name == "fused_ring":
                result["fused_rings_found"] = True
            if construct_name == "bridgehead_locant":
                result["bridgehead_found"] = True

            if severity == "error":
                result["warnings"].append(f"Complex construct detected: {construct_name} ({', '.join(matches)})")

    for pattern in _KNOWN_FUSED_RINGS:
        if re.search(pattern, name, re.IGNORECASE):
            result["detected_constructs"].append({"construct": "fused_heterocycle_system", "matches": [name], "severity": "error"})
            result["fused_rings_found"] = True
            result["difficulty_score"] += 2

    if result["difficulty_score"] == 0:
        result["difficulty_level"] = "simple"
    elif result["difficulty_score"] <= 2:
        result["difficulty_level"] = "moderate"
    elif result["difficulty_score"] <= 4:
        result["difficulty_level"] = "complex"
    else:
        result["difficulty_level"] = "very complex"

    return result


# ─────────────────────────────────────────────
# STAGE 3: OPSIN Parsing
# ─────────────────────────────────────────────
def stage3_opsin_parse(name: str, scan_result: Dict) -> Dict:
    result = {
        "stage": "3 - OPSIN Parsing",
        "input": name,
        "output": None,
        "output_type": None,
        "warnings": [],
        "errors": [],
        "passed": False,
        "stereochemistry_preserved": None,
        "parsing_exception": None,
        "raw_output": None,
        "tokens_parsed": None,
    }

    engine = _get_opsin_engine()
    if engine is None:
        result["errors"].append("OPSIN engine not available (pyopsin not installed or Java not found)")
        result["passed"] = False
        return result

    try:
        raw = engine.to_smiles_single(name)
        result["raw_output"] = raw

        if raw and raw != "None":
            result["output"] = raw
            result["output_type"] = "smiles"
            result["passed"] = True
            result["tokens_parsed"] = "full"
        else:
            result["errors"].append("OPSIN returned empty result")
            result["passed"] = False
            return result

    except Exception as exc:
        result["parsing_exception"] = str(exc)
        result["errors"].append(f"OPSIN parsing exception: {exc}")
        result["passed"] = False
        return result

    if raw and scan_result.get("stereochemistry_found"):
        mol_check = Chem.MolFromSmiles(raw)
        if mol_check:
            n_stereo = Chem.rdMolDescriptors.CalcNumAtomStereoCenters(mol_check)
            if n_stereo > 0:
                result["stereochemistry_preserved"] = True
            elif n_stereo == 0:
                result["stereochemistry_preserved"] = False
                result["warnings"].append("Stereochemistry specified in name but OPSIN output has 0 stereocenters")

    return result


# ─────────────────────────────────────────────
# STAGE 4: SMILES Validation
# ─────────────────────────────────────────────
def stage4_smiles_validate(smiles: str) -> Dict:
    result = {
        "stage": "4 - SMILES Validation",
        "input": smiles,
        "output": None,
        "warnings": [],
        "errors": [],
        "passed": False,
        "canonical_smiles": None,
        "rdkit_valid": False,
    }

    if not smiles:
        result["errors"].append("No SMILES to validate")
        return result

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        result["errors"].append("RDKit could not parse SMILES")
        return result

    result["rdkit_valid"] = True

    try:
        canonical = Chem.MolToSmiles(mol, canonical=True)
        result["canonical_smiles"] = canonical
        if canonical != smiles:
            result["warnings"].append(f"SMILES was not canonical. Input: {smiles}, Canonical: {canonical}")
    except Exception as exc:
        result["errors"].append(f"Canonicalization failed: {exc}")
        return result

    result["output"] = canonical
    result["passed"] = True
    return result


# ─────────────────────────────────────────────
# STAGE 5: External Reference Comparison
# ─────────────────────────────────────────────
_PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
_CACTUS_BASE = "https://cactus.nci.nih.gov/chemical/structure"


def _pubchem_lookup(name: str) -> Optional[Dict]:
    try:
        ctx = ssl._create_unverified_context()
        url = f"{_PUBCHEM_BASE}/compound/name/{urllib.parse.quote(name)}/property/CanonicalSMILES,MolecularFormula,InChI,InChIKey,MolecularWeight/JSON"
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=15, context=ctx) as resp:
            data = json.loads(resp.read().decode())
            if "PropertyTable" in data:
                props = data["PropertyTable"]["Properties"][0]
                return {
                    "smiles": props.get("CanonicalSMILES"),
                    "formula": props.get("MolecularFormula"),
                    "inchi": props.get("InChI"),
                    "inchikey": props.get("InChIKey"),
                    "mw": props.get("MolecularWeight"),
                    "source": "PubChem",
                }
    except Exception:
        return None


def _cactus_lookup(name: str) -> Optional[Dict]:
    try:
        ctx = ssl._create_unverified_context()
        url = f"{_CACTUS_BASE}/{urllib.parse.quote(name)}/smiles"
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=15, context=ctx) as resp:
            smiles = resp.read().decode().strip()
            if smiles:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    return {
                        "smiles": Chem.MolToSmiles(mol, canonical=True),
                        "formula": rdMolDescriptors.CalcMolFormula(mol),
                        "inchikey": Chem.MolToInchiKey(mol),
                        "source": "NCI CACTUS",
                    }
    except Exception:
        return None


def stage5_external_reference(name: str, generated_smiles: str) -> Dict:
    result = {
        "stage": "5 - External Reference Comparison",
        "input": name,
        "output": None,
        "warnings": [],
        "errors": [],
        "passed": False,
        "reference_found": False,
        "reference_data": None,
        "smiles_match": None,
        "formula_match": None,
        "inchikey_match": None,
        "differences": [],
    }

    ref = _pubchem_lookup(name)
    if ref is None:
        ref = _cactus_lookup(name)

    if ref is None:
        result["warnings"].append("No external reference found (offline or compound not in database)")
        result["passed"] = True
        return result

    result["reference_found"] = True
    result["reference_data"] = ref
    result["passed"] = True

    gen_mol = Chem.MolFromSmiles(generated_smiles)
    ref_mol = Chem.MolFromSmiles(ref["smiles"]) if ref.get("smiles") else None

    if gen_mol and ref_mol:
        gen_canon = Chem.MolToSmiles(gen_mol, canonical=True)
        ref_canon = Chem.MolToSmiles(ref_mol, canonical=True)
        result["smiles_match"] = (gen_canon == ref_canon)
        if not result["smiles_match"]:
            result["differences"].append({
                "field": "canonical_smiles",
                "generated": gen_canon,
                "reference": ref_canon,
            })
    elif ref.get("smiles"):
        result["smiles_match"] = (generated_smiles == ref["smiles"])

    gen_formula = rdMolDescriptors.CalcMolFormula(gen_mol) if gen_mol else None
    if gen_formula and ref.get("formula"):
        result["formula_match"] = (gen_formula == ref["formula"])
        if not result["formula_match"]:
            result["differences"].append({
                "field": "formula",
                "generated": gen_formula,
                "reference": ref["formula"],
            })

    gen_ik = Chem.MolToInchiKey(gen_mol) if gen_mol else None
    if gen_ik and ref.get("inchikey"):
        result["inchikey_match"] = (gen_ik == ref["inchikey"])
        if not result["inchikey_match"]:
            result["differences"].append({
                "field": "inchikey",
                "generated": gen_ik,
                "reference": ref["inchikey"],
            })

    return result


# ─────────────────────────────────────────────
# STAGE 6: RDKit Molecular Graph Analysis
# ─────────────────────────────────────────────
def stage6_molecular_graph(smiles: str) -> Dict:
    result = {
        "stage": "6 - RDKit Molecular Graph",
        "input": smiles,
        "output": None,
        "warnings": [],
        "errors": [],
        "passed": False,
        "formula": None,
        "mw": None,
        "num_atoms": None,
        "num_rings": None,
        "num_stereocenters": None,
        "inchi": None,
        "inchikey": None,
        "has_2d_coords": False,
    }

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        result["errors"].append("RDKit could not build molecule from SMILES")
        return result

    try:
        Chem.SanitizeMol(mol)
    except Exception as exc:
        result["warnings"].append(f"Sanitization warning: {exc}")

    result["formula"] = rdMolDescriptors.CalcMolFormula(mol)
    result["mw"] = round(Descriptors.MolWt(mol), 4)
    result["num_atoms"] = mol.GetNumAtoms()
    result["num_heavy_atoms"] = mol.GetNumHeavyAtoms()
    result["num_rings"] = rdMolDescriptors.CalcNumRings(mol)
    result["num_aromatic_rings"] = rdMolDescriptors.CalcNumAromaticRings(mol)
    result["num_stereocenters"] = Chem.rdMolDescriptors.CalcNumAtomStereoCenters(mol)

    try:
        result["inchi"] = Chem.MolToInchi(mol)
    except Exception:
        result["inchi"] = "N/A"
    try:
        result["inchikey"] = Chem.MolToInchiKey(mol)
    except Exception:
        result["inchikey"] = "N/A"

    try:
        conf = mol.GetConformer()
        if conf is not None:
            result["has_2d_coords"] = True
    except Exception:
        result["has_2d_coords"] = False

    result["output"] = result["inchikey"]
    result["passed"] = True
    return result


# ─────────────────────────────────────────────
# STAGE 7: 2D Coordinate Generation
# ─────────────────────────────────────────────
def stage7_2d_coords(smiles: str) -> Dict:
    result = {
        "stage": "7 - 2D Coordinate Generation",
        "input": smiles,
        "output": None,
        "warnings": [],
        "errors": [],
        "passed": False,
        "coords_generated": False,
        "coordgen_used": False,
        "depiction_issues": [],
    }

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        result["errors"].append("No molecule to generate coordinates for")
        return result

    try:
        from rdkit.Chem import rdDepictor
        rdDepictor.SetPreferCoordGen(True)
        result["coordgen_used"] = True
    except Exception:
        pass

    try:
        AllChem.Compute2DCoords(mol)
        result["coords_generated"] = True
    except Exception as exc:
        result["errors"].append(f"Coordinate generation failed: {exc}")
        return result

    try:
        conf = mol.GetConformer()
        positions = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
        result["output"] = mol
    except Exception as exc:
        result["warnings"].append(f"Coordinate inspection failed: {exc}")

    result["passed"] = True
    return result


# ─────────────────────────────────────────────
# STAGE 8: Image Rendering Analysis
# ─────────────────────────────────────────────
def stage8_rendering(smiles: str) -> Dict:
    result = {
        "stage": "8 - Image Rendering",
        "input": smiles,
        "output": None,
        "warnings": [],
        "errors": [],
        "passed": False,
        "png_generated": False,
        "svg_generated": False,
        "png_size_bytes": None,
        "svg_size_chars": None,
    }

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        result["errors"].append("No molecule to render")
        return result

    try:
        AllChem.Compute2DCoords(mol)
    except Exception:
        pass

    try:
        svg = mol_to_svg(mol)
        if svg:
            result["svg_generated"] = True
            result["svg_size_chars"] = len(svg)
    except Exception as exc:
        result["errors"].append(f"SVG rendering failed: {exc}")

    try:
        png = mol_to_png_base64(mol)
        if png:
            result["png_generated"] = True
            result["png_size_bytes"] = len(png)
    except Exception as exc:
        result["errors"].append(f"PNG rendering failed: {exc}")

    result["passed"] = result["png_generated"] or result["svg_generated"]
    return result


# ─────────────────────────────────────────────
# FULL PIPELINE DIAGNOSTIC
# ─────────────────────────────────────────────
def run_pipeline_diagnostics(name: str) -> Dict:
    report = {
        "input": name,
        "stages": [],
        "overall_status": "unknown",
        "failure_stage": None,
        "failure_reason": None,
        "recommendation": None,
    }

    s1 = stage1_preprocess(name)
    report["stages"].append(s1)

    normalized = s1["output"]

    from utils import detect_input_type
    input_type = detect_input_type(normalized)

    s2 = stage2_scan_iupac(normalized)
    report["stages"].append(s2)

    if input_type in ("smiles", "inchi"):
        from rdkit import Chem
        if input_type == "inchi":
            inchi_mol = Chem.MolFromInchi(normalized)
            inchi_smiles = Chem.MolToSmiles(inchi_mol, canonical=True) if inchi_mol else None
            stage_output = inchi_smiles or normalized
        else:
            stage_output = normalized

        s3 = {
            "stage": "3 - OPSIN Parsing",
            "input": normalized,
            "output": stage_output,
            "output_type": input_type,
            "warnings": [f"Input detected as {input_type.upper()}, not IUPAC name. OPSIN skipped."],
            "errors": [],
            "passed": True,
            "stereochemistry_preserved": None,
            "parsing_exception": None,
            "raw_output": None,
            "tokens_parsed": None,
        }
        opsin_smiles = stage_output
    else:
        s3 = stage3_opsin_parse(normalized, s2)
        opsin_smiles = s3.get("output")

    report["stages"].append(s3)
    if opsin_smiles:
        s4 = stage4_smiles_validate(opsin_smiles)
    else:
        s4 = {
            "stage": "4 - SMILES Validation",
            "input": None,
            "output": None,
            "warnings": [],
            "errors": ["Skipped - no SMILES from OPSIN"],
            "passed": False,
        }
    report["stages"].append(s4)

    canonical_smiles = s4.get("output") or s4.get("canonical_smiles")
    if canonical_smiles:
        s5 = stage5_external_reference(normalized, canonical_smiles)
    else:
        s5 = {
            "stage": "5 - External Reference Comparison",
            "input": normalized,
            "output": None,
            "warnings": ["Skipped - no valid SMILES to compare"],
            "errors": [],
            "passed": True,
            "reference_found": False,
            "reference_data": None,
        }
    report["stages"].append(s5)

    if canonical_smiles:
        s6 = stage6_molecular_graph(canonical_smiles)
    else:
        s6 = {
            "stage": "6 - RDKit Molecular Graph",
            "input": None,
            "output": None,
            "warnings": [],
            "errors": ["Skipped - no valid SMILES"],
            "passed": False,
        }
    report["stages"].append(s6)

    if canonical_smiles:
        s7 = stage7_2d_coords(canonical_smiles)
    else:
        s7 = {
            "stage": "7 - 2D Coordinate Generation",
            "input": None,
            "output": None,
            "warnings": [],
            "errors": ["Skipped - no valid SMILES"],
            "passed": False,
        }
    report["stages"].append(s7)

    if canonical_smiles:
        s8 = stage8_rendering(canonical_smiles)
    else:
        s8 = {
            "stage": "8 - Image Rendering",
            "input": None,
            "output": None,
            "warnings": [],
            "errors": ["Skipped - no valid SMILES"],
            "passed": False,
        }
    report["stages"].append(s8)

    report = _compute_overall_status(report)
    report = _generate_recommendation(report)

    return report


def _compute_overall_status(report: Dict) -> Dict:
    stages = report["stages"]
    failed_stages = [s for s in stages if s.get("passed") is False and "Skipped" not in str(s.get("errors", []))]
    warning_stages = [s for s in stages if s.get("warnings") and s.get("passed") is not False]

    failed_stages_only = [s for s in failed_stages if not any("Skipped" in e for e in s.get("errors", []))]

    if not failed_stages_only:
        if warning_stages:
            report["overall_status"] = "passed_with_warnings"
        else:
            report["overall_status"] = "passed"
    else:
        report["overall_status"] = "failed"
        first_fail = failed_stages_only[0]
        report["failure_stage"] = first_fail["stage"]
        report["failure_reason"] = "; ".join(first_fail.get("errors", []))

    return report


def _generate_recommendation(report: Dict) -> Dict:
    s2 = next((s for s in report["stages"] if s["stage"].startswith("2")), None)
    s3 = next((s for s in report["stages"] if s["stage"].startswith("3")), None)
    s5 = next((s for s in report["stages"] if s["stage"].startswith("5")), None)

    recommendations = []

    if report["overall_status"] == "passed":
        if s5 and s5.get("reference_found") and s5.get("smiles_match") is True:
            recommendations.append("Pipeline successful. Generated structure matches external reference.")
        elif s5 and s5.get("reference_found") and s5.get("smiles_match") is False:
            recommendations.append(
                f"OPSIN parsed the name but the result differs from the external reference. "
                f"Generated SMILES: {s5.get('differences', [{}])[0].get('generated', '?')}, "
                f"Reference SMILES: {s5.get('differences', [{}])[0].get('reference', '?')}. "
                f"Verify the IUPAC name or use the reference SMILES directly."
            )
        elif s5 and not s5.get("reference_found"):
            if s2 and s2.get("difficulty_level") in ("complex", "very complex"):
                recommendations.append(
                    f"Pipeline completed but no external reference found for cross-checking. "
                    f"The IUPAC name contains complex constructs ({s2.get('difficulty_level')}). "
                    f"Review the generated structure manually."
                )
            else:
                recommendations.append("Pipeline completed successfully. No reference available for verification.")
        else:
            recommendations.append("Pipeline completed successfully.")

    elif report["overall_status"] == "passed_with_warnings":
        warnings = []
        for s in report["stages"]:
            if s.get("warnings"):
                warnings.extend(s["warnings"])
        rec = f"Pipeline completed with warnings: {'; '.join(warnings[:2])}. "
        if s5 and s5.get("reference_found") and s5.get("smiles_match") is False:
            rec += "Reference comparison reveals structural differences."
        recommendations.append(rec)

    elif report["overall_status"] == "failed":
        stage_name = report.get("failure_stage", "unknown")
        reason = report.get("failure_reason", "unknown")
        recommendations.append(f"Pipeline failed at {stage_name}: {reason}")

        if "OPSIN" in stage_name:
            if s2 and s2.get("difficulty_level") in ("complex", "very complex"):
                recommendations.append(
                    f"The IUPAC name contains {s2.get('difficulty_level')} constructs "
                    f"(stereochemistry: {s2.get('stereochemistry_found')}, "
                    f"fused rings: {s2.get('fused_rings_found')}). "
                    f"OPSIN may not fully support these features. "
                    f"Try using the SMILES directly."
                )
            recommendations.append("Install a Java runtime (JRE) and ensure JAVA_HOME is set.")

    report["recommendation"] = " ".join(recommendations)
    return report


def summarize_pipeline_report(report: Dict) -> str:
    lines = []
    lines.append("=" * 60)
    lines.append("PIPELINE DIAGNOSTIC REPORT")
    lines.append("=" * 60)
    lines.append(f"Input: {report['input']}")
    lines.append(f"Status: {report['overall_status']}")
    if report.get("failure_stage"):
        lines.append(f"Failure at: {report['failure_stage']}")
        lines.append(f"Reason: {report['failure_reason']}")
    lines.append("")

    for stage in report["stages"]:
        status = "PASS" if stage.get("passed") else "FAIL"
        status = "SKIP" if any("Skipped" in e for e in stage.get("errors", [])) else status
        status = "WARN" if (stage.get("passed") and stage.get("warnings")) else status
        lines.append(f"  [{status}] {stage['stage']}")
        for w in stage.get("warnings", []):
            lines.append(f"         WARNING: {w}")
        for e in stage.get("errors", []):
            if "Skipped" not in e:
                lines.append(f"         ERROR: {e}")

    lines.append("")
    if report.get("recommendation"):
        lines.append(f"Recommendation: {report['recommendation']}")
    lines.append("=" * 60)

    return "\n".join(lines)
