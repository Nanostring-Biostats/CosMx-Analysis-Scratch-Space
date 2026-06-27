from __future__ import annotations

from copy import deepcopy
from typing import Dict, Iterable, List, Optional

import pandas as pd

from ..datamodels import MarkerProgram
from ..inspect import describe_hierarchy
from .schema import _mp
from .matching import match_builtin_hierarchy, suggest_builtin_hierarchies
from .factories import (
    _immune_core,
    _solid_tissue_core,
    _tissue_origin_screen,
    _brain_core,
    _tme_core,
    _colon_tme,
    _kidney_tme,
    _tonsil_tme,
    _breast_tme,
    _pancreas_tme,
    _liver_tme,
    _lung_tme,
    _skin_tme,
)

BUILTIN_HIERARCHIES: Dict[str, Dict[str, object]] = {
    "immune_core": {"factory": _immune_core, "description": "Pan-tissue immune hierarchy.", "category": "immune", "tissue_scope": "pan-tissue"},
    "solid_tissue_core": {"factory": _solid_tissue_core, "description": "Major non-immune solid tissue lineages.", "category": "solid_tissue", "tissue_scope": "pan-solid-tissue"},
    "tissue_origin_screen": {"factory": _tissue_origin_screen, "description": "Shallow tissue-origin screening hierarchy for detect_tissue_type; not intended as a final cell-type hierarchy.", "category": "tissue_detection", "tissue_scope": "pan-tissue"},
    "brain_core": {"factory": _brain_core, "description": "Core neuro hierarchy for brain-like samples.", "category": "brain", "tissue_scope": "brain"},
    "tme_core": {"factory": _tme_core, "description": "Combined solid tissue plus immune hierarchy for mixed tissue microenvironments.", "category": "tme", "tissue_scope": "pan-solid-tissue"},
    "colon_tme": {"factory": _colon_tme, "description": "Colon-flavored tissue microenvironment hierarchy.", "category": "tme", "tissue_scope": "colon"},
    "kidney_tme": {"factory": _kidney_tme, "description": "Kidney-flavored tissue microenvironment hierarchy.", "category": "tme", "tissue_scope": "kidney"},
    "tonsil_tme": {"factory": _tonsil_tme, "description": "Tonsil-flavored tissue microenvironment hierarchy.", "category": "tme", "tissue_scope": "tonsil"},
    "breast_tme": {"factory": _breast_tme, "description": "Breast-flavored tissue microenvironment hierarchy.", "category": "tme", "tissue_scope": "breast"},
    "pancreas_tme": {"factory": _pancreas_tme, "description": "Pancreas-flavored tissue microenvironment hierarchy.", "category": "tme", "tissue_scope": "pancreas"},
    "liver_tme": {"factory": _liver_tme, "description": "Liver-flavored tissue microenvironment hierarchy.", "category": "tme", "tissue_scope": "liver"},
    "lung_tme": {"factory": _lung_tme, "description": "Lung-flavored tissue microenvironment hierarchy.", "category": "tme", "tissue_scope": "lung"},
    "skin_tme": {"factory": _skin_tme, "description": "Skin-flavored tissue microenvironment hierarchy.", "category": "tme", "tissue_scope": "skin"},
}




def _clone_program(node: MarkerProgram) -> MarkerProgram:
    return _mp(
        name=node.name,
        pos=list(node.positive_markers),
        neg=list(node.negative_markers),
        children=[_clone_program(c) for c in node.children],
        description=node.description,
        aliases=list(node.aliases),
        metadata=dict(node.metadata),
    )


def _program_panel_coverage(node: MarkerProgram, gene_set: set[str]) -> dict:
    pos = [g for g in node.positive_markers]
    neg = [g for g in node.negative_markers]
    pos_present = [g for g in pos if g in gene_set]
    neg_present = [g for g in neg if g in gene_set]
    pos_fraction = (len(pos_present) / len(pos)) if len(pos) else 1.0
    neg_fraction = (len(neg_present) / len(neg)) if len(neg) else 1.0
    return {
        "positive_total": len(pos),
        "positive_present": len(pos_present),
        "positive_fraction": pos_fraction,
        "negative_total": len(neg),
        "negative_present": len(neg_present),
        "negative_fraction": neg_fraction,
    }


def _panel_tier_from_gene_count(n_genes: int) -> str:
    if n_genes <= 1500:
        return "small"
    if n_genes <= 9000:
        return "medium"
    return "large"


def _compile_program_for_panel(
    node: MarkerProgram,
    gene_set: set[str],
    panel_tier: str,
    preserve_major_immune_subtypes: bool,
    report_rows: list[dict],
    level: int = 1,
) -> Optional[MarkerProgram]:
    cov = _program_panel_coverage(node, gene_set)
    meta = dict(node.metadata)
    role = str(meta.get("role", "major")).lower()
    requires_ancestor_support = bool(meta.get("requires_ancestor_support", True))
    fallback_to_parent_allowed = bool(meta.get("fallback_to_parent", True))
    is_immune = str(meta.get("category", "")).lower() == "immune" or node.name in {"Immune", "T cell", "B cell", "NK cell", "Myeloid", "Macrophage", "Monocyte", "Dendritic cell", "Mast cell", "Neutrophil", "Plasma cell", "CD4 T cell", "CD8 T cell", "Treg"}
    major_immune_names = {
        "Immune", "T cell", "CD4 T cell", "CD8 T cell", "NK cell", "B cell", "Plasma cell",
        "Myeloid", "Macrophage", "Monocyte", "Dendritic cell", "Mast cell", "Neutrophil"
    }

    compiled = _clone_program(node)
    compiled.positive_markers = [g for g in compiled.positive_markers if g in gene_set]
    compiled.negative_markers = [g for g in compiled.negative_markers if g in gene_set]

    disable = False
    fallback_to_parent = False

    if preserve_major_immune_subtypes and is_immune and node.name in major_immune_names:
        disable = len(compiled.positive_markers) < 2
    else:
        if panel_tier == "small":
            min_frac_major, fb_frac_major = 0.40, 0.60
            min_frac_optional, fb_frac_optional = 0.55, 0.75
        elif panel_tier == "medium":
            min_frac_major, fb_frac_major = 0.30, 0.50
            min_frac_optional, fb_frac_optional = 0.40, 0.60
        else:
            min_frac_major, fb_frac_major = 0.20, 0.35
            min_frac_optional, fb_frac_optional = 0.30, 0.45

        if role == "optional":
            disable = len(compiled.positive_markers) < 2 or cov["positive_fraction"] < min_frac_optional
            fallback_to_parent = cov["positive_fraction"] < fb_frac_optional
        else:
            disable = len(compiled.positive_markers) < 2 or cov["positive_fraction"] < min_frac_major
            fallback_to_parent = cov["positive_fraction"] < fb_frac_major

    compiled_children = []
    for child in node.children:
        child_compiled = _compile_program_for_panel(
            child,
            gene_set=gene_set,
            panel_tier=panel_tier,
            preserve_major_immune_subtypes=preserve_major_immune_subtypes,
            report_rows=report_rows,
            level=level + 1,
        )
        if child_compiled is not None:
            compiled_children.append(child_compiled)
    compiled.children = compiled_children

    # Keep structurally important parents if they retain any children.
    if disable and compiled.children and requires_ancestor_support:
        disable = False
        fallback_to_parent = True

    if fallback_to_parent and not fallback_to_parent_allowed:
        fallback_to_parent = False

    status = "disabled" if disable else ("fallback_to_parent" if fallback_to_parent else "full")
    report_rows.append({
        "name": node.name,
        "level": level,
        "panel_tier": panel_tier,
        "category": meta.get("category"),
        "role": role,
        "lineage_module": meta.get("lineage_module"),
        "requires_ancestor_support": requires_ancestor_support,
        "fallback_to_parent_allowed": fallback_to_parent_allowed,
        "positive_markers_total": cov["positive_total"],
        "positive_markers_present": cov["positive_present"],
        "positive_markers_present_fraction": cov["positive_fraction"],
        "negative_markers_total": cov["negative_total"],
        "negative_markers_present": cov["negative_present"],
        "negative_markers_present_fraction": cov["negative_fraction"],
        "compiled_positive_markers": list(compiled.positive_markers),
        "compiled_negative_markers": list(compiled.negative_markers),
        "status": status,
    })

    if disable:
        return None
    if fallback_to_parent and not compiled.children:
        # leaf survives as a shallow leaf; downstream can treat this as parent fallback semantics
        compiled.metadata = dict(compiled.metadata)
        compiled.metadata["compiled_status"] = "fallback_to_parent"
    return compiled


def compile_builtin_hierarchy_for_panel(
    gene_list,
    builtin_name: str,
    panel_name: Optional[str] = None,
    preserve_major_immune_subtypes: bool = True,
):
    """
    Compile a built-in hierarchy against an input RNA panel gene list.

    The compiler automatically prunes or simplifies weakly covered nodes while
    preserving major immune subtyping across small, medium, and large panels.
    Returns a tuple: (compiled_roots, compilation_report).
    """
    roots = get_builtin_hierarchy(builtin_name, copy=True)
    gene_set = {str(g).upper() for g in list(gene_list) if pd.notna(g) and str(g).strip() != ""}
    report_rows = []
    panel_tier = _panel_tier_from_gene_count(len(gene_set))
    compiled = []
    for root in roots:
        out = _compile_program_for_panel(
            root,
            gene_set=gene_set,
            panel_tier=panel_tier,
            preserve_major_immune_subtypes=preserve_major_immune_subtypes,
            report_rows=report_rows,
            level=1,
        )
        if out is not None:
            compiled.append(out)
    report = pd.DataFrame(report_rows)
    if not report.empty:
        report["builtin_name"] = builtin_name
        report["panel_name"] = panel_name
        report["panel_gene_count"] = len(gene_set)
    return compiled, report

def _flatten(programs: Iterable[MarkerProgram]) -> List[MarkerProgram]:
    out: List[MarkerProgram] = []
    def visit(node: MarkerProgram):
        out.append(node)
        for child in node.children:
            visit(child)
    for root in programs:
        visit(root)
    return out


def list_builtin_hierarchies() -> pd.DataFrame:
    rows = []
    for name, spec in BUILTIN_HIERARCHIES.items():
        roots = spec["factory"]()
        rows.append({
            "name": name,
            "category": spec["category"],
            "tissue_scope": spec["tissue_scope"],
            "description": spec["description"],
            "n_root_programs": len(roots),
            "n_total_programs": len(_flatten(roots)),
        })
    return pd.DataFrame(rows).sort_values(["category", "name"]).reset_index(drop=True)


def get_builtin_hierarchy(name: str, copy: bool = True) -> List[MarkerProgram]:
    if name not in BUILTIN_HIERARCHIES:
        raise KeyError(f"Unknown built-in hierarchy: {name}")
    roots = BUILTIN_HIERARCHIES[name]["factory"]()
    return deepcopy(roots) if copy else roots


def describe_builtin_hierarchy(name: str) -> pd.DataFrame:
    roots = get_builtin_hierarchy(name, copy=True)
    df = describe_hierarchy(roots)
    if not df.empty:
        df["pack"] = name
    return df


def list_builtin_programs(category: Optional[str] = None, tissue_scope: Optional[str] = None) -> pd.DataFrame:
    rows = []
    for pack_name, spec in BUILTIN_HIERARCHIES.items():
        for node in _flatten(spec["factory"]()):
            meta = dict(node.metadata)
            rows.append({
                "pack": pack_name,
                "name": node.name,
                "description": node.description,
                "category": meta.get("category", spec["category"]),
                "tissue_scope": meta.get("tissue_scope", spec["tissue_scope"]),
                "role": meta.get("role"),
                "lineage_module": meta.get("lineage_module"),
                "requires_ancestor_support": meta.get("requires_ancestor_support"),
                "fallback_to_parent": meta.get("fallback_to_parent"),
                "n_positive_markers": len(node.positive_markers),
                "n_negative_markers": len(node.negative_markers),
                "n_children": len(node.children),
                "aliases": list(node.aliases),
            })
    df = pd.DataFrame(rows)
    if category is not None:
        df = df[df["category"] == category]
    if tissue_scope is not None:
        df = df[df["tissue_scope"] == tissue_scope]
    return df.sort_values(["pack", "name"]).reset_index(drop=True)


def get_builtin_programs(pack_name: str, names: Optional[List[str]] = None) -> List[MarkerProgram]:
    roots = get_builtin_hierarchy(pack_name, copy=True)
    flat = _flatten(roots)
    if names is None:
        return flat
    wanted = set(names)
    return [node.clone() for node in flat if node.name in wanted]


def compose_hierarchy(*pack_names: str, deduplicate_roots: bool = True) -> List[MarkerProgram]:
    roots: List[MarkerProgram] = []
    seen = set()
    for pack_name in pack_names:
        for root in get_builtin_hierarchy(pack_name, copy=True):
            if deduplicate_roots and root.name in seen:
                continue
            seen.add(root.name)
            roots.append(root)
    return roots
