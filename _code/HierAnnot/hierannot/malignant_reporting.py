from __future__ import annotations

from typing import Iterable, Optional

import pandas as pd

from .builtin.matching import _normalize_lineage_text, is_immune_like_lineage
from .constants import (
    ANNOT_PATH_SEPARATOR,
    TUMOR_REPORTABLE_GUARD_LINEAGES,
    PROGRAM_REPORT_BLOCK_OFF_TOKENS,
    PROGRAM_REPORT_BLOCK_PRESETS,
    TUMOR_REPORTABLE_BUILTIN_LINEAGES,
)


def _validate_program_report_block_preset(value):
    """Validate and normalize the program-report blocking control.

    Semantics are intentionally reporting-oriented:
    - None, "none", or "off": do not block program reporting.
    - "immune_like": block immune-like final normal calls.
    - "tumor_reportable": report only on curated tumor-reportable built-in hierarchy branches.
    - "lineage_aware": use the selected program's metadata["report_on_lineages"].
    - a sequence of strings: block rows whose final normal label/path contains
      any exact node in the sequence.
    """
    if value is None:
        return None
    if isinstance(value, str):
        text = value.strip()
        if text == "" or text.lower() in PROGRAM_REPORT_BLOCK_OFF_TOKENS:
            return None
        low = text.lower()
        if low in PROGRAM_REPORT_BLOCK_PRESETS:
            return low
        return text
    if isinstance(value, (list, tuple, set)):
        if not all(isinstance(x, str) for x in value):
            raise TypeError(
                "program_report_block_preset must be None, a string preset, or a list/tuple/set of strings"
            )
        cleaned = []
        for x in value:
            nx = str(x).strip()
            if nx:
                cleaned.append(nx)
        return cleaned
    raise TypeError(
        "program_report_block_preset must be None, a string preset, or a list/tuple/set of strings"
    )


def _normalize_node_set(values: Optional[Iterable[str]]) -> set[str]:
    if values is None:
        return set()
    out = set()
    for v in values:
        if v is None:
            continue
        nv = _normalize_lineage_text(str(v))
        if nv:
            out.add(nv)
    return out


def _normalize_report_on_lineages(value) -> list[str]:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return []
    if isinstance(value, str):
        parts = [x.strip() for x in value.split(";") if x.strip()]
        if len(parts) == 1 and "," in parts[0]:
            parts = [x.strip() for x in parts[0].split(",") if x.strip()]
        return parts
    if isinstance(value, (list, tuple, set)):
        return [str(x).strip() for x in value if str(x).strip()]
    return []


def _add_scalar_node_candidate(node_set: set[str], value) -> None:
    if value is None or pd.isna(value):
        return
    norm = _normalize_lineage_text(str(value))
    if norm:
        node_set.add(norm)


def _add_path_node_candidates(node_set: set[str], value) -> None:
    if value is None or pd.isna(value):
        return
    path = str(value).strip()
    if not path:
        return
    for node in [x.strip() for x in path.split(ANNOT_PATH_SEPARATOR) if x.strip()]:
        norm = _normalize_lineage_text(node)
        if norm:
            node_set.add(norm)


def _collect_row_node_candidates(row) -> set[str]:
    """Collect final normal hierarchy nodes for program-report block checks.

    Program reporting should be blocked by the final normal hierarchy call, not
    by diagnostic compartment hints or downstream labels.  This intentionally
    ignores fields such as ``annot_primary_compartment``,
    ``annot_integrated_label``, and ``annot_export_label`` because those can
    disagree with the final routed normal label or make blocking circular.
    """
    node_set: set[str] = set()

    for c in ["annot_final_top_branch", "annot_label"]:
        if c in row.index:
            _add_scalar_node_candidate(node_set, row[c])

    for c in ["annot_path"]:
        if c in row.index:
            _add_path_node_candidates(node_set, row[c])

    if not node_set:
        for c in ["annot_decision_label"]:
            if c in row.index:
                _add_scalar_node_candidate(node_set, row[c])
        for c in ["annot_decision_path"]:
            if c in row.index:
                _add_path_node_candidates(node_set, row[c])

    return node_set


def _node_set_matches_any(node_set: set[str], targets: Iterable[str]) -> bool:
    target_set = _normalize_node_set(targets)
    if not target_set:
        return False
    for node in node_set:
        for target in target_set:
            if node == target or f" {target} " in f" {node} " or f" {node} " in f" {target} ":
                return True
    return False


def _is_tumor_reportable_builtin_lineage(node_set: set[str]) -> bool:
    """Return True when final normal nodes match tumor-reportable lineages.

    The ``tumor_reportable`` preset primarily uses a curated built-in hierarchy
    allowlist plus a few high-confidence generic epithelial/parenchymal names.
    Immune-like labels are guarded by ``is_immune_like_lineage()`` and a short
    additional guard list blocks common stromal, vascular, and mural labels, so
    labels such as ``Tumor-associated macrophage`` or
    ``Cancer-associated fibroblast`` are not accidentally allowed by a broad
    epithelial/parenchymal phrase elsewhere in the final path. Users working
    with custom hierarchies can pass an explicit block list or use
    ``lineage_aware`` with program-level ``report_on_lineages`` metadata.
    """
    if any(is_immune_like_lineage(x) for x in node_set if x):
        return False
    if _node_set_matches_any(node_set, TUMOR_REPORTABLE_GUARD_LINEAGES):
        return False
    return _node_set_matches_any(node_set, TUMOR_REPORTABLE_BUILTIN_LINEAGES)



def _immune_like_from_builtin_text(summary: pd.DataFrame) -> pd.Series:
    def blocked(row) -> bool:
        node_set = _collect_row_node_candidates(row)
        return any(is_immune_like_lineage(x) for x in node_set if x)

    return summary.apply(blocked, axis=1).astype(bool)


def _tumor_reportable_from_builtin_lineages(summary: pd.DataFrame) -> pd.Series:
    def blocked(row) -> bool:
        node_set = _collect_row_node_candidates(row)
        if not node_set:
            return False
        return not _is_tumor_reportable_builtin_lineage(node_set)

    return summary.apply(blocked, axis=1).astype(bool)


def _blocked_by_explicit_nodes(summary: pd.DataFrame, blocked_nodes: Optional[Iterable[str]]) -> pd.Series:
    blocked = _normalize_node_set(blocked_nodes)
    if not blocked:
        return pd.Series(False, index=summary.index)

    def hit(row) -> bool:
        node_set = _collect_row_node_candidates(row)
        return any(node in blocked for node in node_set)

    return summary.apply(hit, axis=1).astype(bool)


def _lineage_aware_reportability(summary: pd.DataFrame) -> tuple[pd.Series, pd.Series, pd.Series]:
    """Return blocked, matched, and report-on-lineages strings for selected programs.

    A program with no metadata["report_on_lineages"] is considered broadly
    reportable. If the selected program defines report_on_lineages, the final
    normal label/path must contain at least one of those lineages.
    """
    blocked = pd.Series(False, index=summary.index, dtype=bool)
    matched = pd.Series(True, index=summary.index, dtype=bool)
    report_on = pd.Series("", index=summary.index, dtype=object)

    def eval_row(row):
        vals = _normalize_report_on_lineages(row.get("annot_malignant_program_report_on_lineages", None))
        if not vals:
            return False, True, ""
        node_set = _collect_row_node_candidates(row)
        ok = _node_set_matches_any(node_set, vals)
        return (not ok), bool(ok), ";".join(vals)

    if len(summary) > 0:
        evaluated = summary.apply(eval_row, axis=1)
        blocked = evaluated.map(lambda x: bool(x[0])).astype(bool)
        matched = evaluated.map(lambda x: bool(x[1])).astype(bool)
        report_on = evaluated.map(lambda x: x[2])
    return blocked, matched, report_on


def resolve_program_report_block_preset(summary: pd.DataFrame, program_report_block_preset=None) -> pd.DataFrame:
    """Resolve report blocking for the flat tumor/auxiliary program track.

    Returns a DataFrame with row-aligned diagnostic columns:
    ``blocked``, ``reason``, ``global_blocked``, ``program_lineage_blocked``,
    ``program_lineage_matched``, and ``program_report_on_lineages``.
    """
    preset = _validate_program_report_block_preset(program_report_block_preset)
    out = pd.DataFrame(index=summary.index)
    global_blocked = pd.Series(False, index=summary.index, dtype=bool)
    lineage_blocked = pd.Series(False, index=summary.index, dtype=bool)
    lineage_matched = pd.Series(True, index=summary.index, dtype=bool)
    report_on = pd.Series("", index=summary.index, dtype=object)
    reasons = pd.Series("", index=summary.index, dtype=object)

    if preset is None:
        pass
    elif isinstance(preset, str):
        if preset == "immune_like":
            global_blocked = _immune_like_from_builtin_text(summary)
            reasons.loc[global_blocked] = "blocked_by_immune_like_preset"
        elif preset == "tumor_reportable":
            global_blocked = _tumor_reportable_from_builtin_lineages(summary)
            reasons.loc[global_blocked] = "blocked_by_tumor_reportable_preset"
        elif preset == "lineage_aware":
            lineage_blocked, lineage_matched, report_on = _lineage_aware_reportability(summary)
            reasons.loc[lineage_blocked] = "blocked_by_program_lineage"
        else:
            global_blocked = _blocked_by_explicit_nodes(summary, [preset])
            reasons.loc[global_blocked] = "blocked_by_node_list"
    elif isinstance(preset, (list, tuple, set)):
        normalized = _validate_program_report_block_preset(preset)
        if normalized is None:
            pass
        else:
            normalized_lower = [str(x).strip().lower() for x in normalized]
            explicit = [x for x in normalized if str(x).strip().lower() not in {"immune_like", "tumor_reportable", "lineage_aware"}]
            if "immune_like" in normalized_lower:
                mask = _immune_like_from_builtin_text(summary)
                global_blocked = global_blocked | mask
                reasons.loc[mask] = "blocked_by_immune_like_preset"
            if "tumor_reportable" in normalized_lower:
                mask = _tumor_reportable_from_builtin_lineages(summary)
                global_blocked = global_blocked | mask
                reasons.loc[mask] = "blocked_by_tumor_reportable_preset"
            if "lineage_aware" in normalized_lower:
                lb, lm, ro = _lineage_aware_reportability(summary)
                lineage_blocked = lineage_blocked | lb
                lineage_matched = lm
                report_on = ro
                reasons.loc[lb] = "blocked_by_program_lineage"
            if explicit:
                mask = _blocked_by_explicit_nodes(summary, explicit)
                global_blocked = global_blocked | mask
                reasons.loc[mask] = "blocked_by_node_list"

    blocked = (global_blocked | lineage_blocked).astype(bool)
    reasons.loc[~blocked] = ""
    out["blocked"] = blocked
    out["reason"] = reasons
    out["global_blocked"] = global_blocked.astype(bool)
    out["program_lineage_blocked"] = lineage_blocked.astype(bool)
    out["program_lineage_matched"] = lineage_matched.astype(bool)
    out["program_report_on_lineages"] = report_on
    return out

