from __future__ import annotations

from .constants import ANNOT_PATH_SEPARATOR
from .utils import _sanitize_label_text
from .hierarchy import decide_cluster_annotation_from_scores
from .malignant_reporting import _validate_program_report_block_preset

from typing import Iterable, Optional, Sequence, List, Dict, Set, Union

import numpy as np
import re
import pandas as pd
import warnings

def _select_annotation_summary_source(result, source: str = "normal") -> pd.DataFrame:
    source = str(source)
    if source == "integrated":
        summary = getattr(result, "integrated_annotations", None)
        if summary is not None and len(summary) > 0:
            return summary.copy()
    return result.cluster_annotations.copy()


def _build_per_level_summary(result) -> pd.DataFrame:
    """Build a wide per-cluster summary from result.level_scores."""
    level_scores = getattr(result, "level_scores", None)
    if level_scores is None or len(level_scores) == 0:
        return pd.DataFrame(columns=["cluster_id"])

    df = level_scores.copy()
    df["cluster"] = df["cluster"].astype(str)
    score_col = "decision_score" if "decision_score" in df.columns else "score"
    best_idx = df.groupby(["cluster", "level"], sort=False)[score_col].idxmax().dropna()
    best_idx = best_idx.astype(int) if len(best_idx) else best_idx
    best = df.loc[best_idx].copy() if len(best_idx) else pd.DataFrame(columns=df.columns)

    wide = pd.DataFrame({"cluster_id": sorted(df["cluster"].unique(), key=lambda x: x)})
    level_values = sorted(best["level"].dropna().unique())
    for level in level_values:
        sub = best[best["level"] == level].copy()
        rename_map = {
            "label": f"annot_l{int(level)}_label",
            score_col: f"annot_l{int(level)}_score",
            "parent_label": f"annot_l{int(level)}_parent_label",
            "positive_score": f"annot_l{int(level)}_positive_score",
            "negative_score": f"annot_l{int(level)}_negative_score",
            "marker_source": f"annot_l{int(level)}_marker_source",
            "score_status": f"annot_l{int(level)}_score_status",
            "branch_supported_raw_score": f"annot_l{int(level)}_branch_supported_raw_score",
            "markers_present": f"annot_l{int(level)}_markers_present",
            "markers_total": f"annot_l{int(level)}_markers_total",
            "markers_present_fraction": f"annot_l{int(level)}_markers_present_fraction",
            "markers_detection_support_fraction": f"annot_l{int(level)}_markers_detection_support_fraction",
            "negative_markers_present": f"annot_l{int(level)}_negative_markers_present",
            "negative_markers_total": f"annot_l{int(level)}_negative_markers_total",
            "negative_markers_present_fraction": f"annot_l{int(level)}_negative_markers_present_fraction",
            "used_fallback_marker_set": f"annot_l{int(level)}_used_fallback_marker_set",
            "control_fallback_used": f"annot_l{int(level)}_control_fallback_used",
            "median_controls_per_positive_marker": f"annot_l{int(level)}_median_controls_per_positive_marker",
            "median_controls_per_negative_marker": f"annot_l{int(level)}_median_controls_per_negative_marker",
            "sibling_specificity_score": f"annot_l{int(level)}_sibling_specificity_score",
            "program_robust_zscore": f"annot_l{int(level)}_program_robust_zscore",
        }
        if score_col != "score" and "score" in sub.columns:
            rename_map["score"] = f"annot_l{int(level)}_raw_score"
        keep_cols = ["cluster", *[c for c in rename_map if c in sub.columns]]
        sub = sub[keep_cols].rename(columns={k: v for k, v in rename_map.items() if k in keep_cols})
        sub = sub.rename(columns={"cluster": "cluster_id"})
        wide = wide.merge(sub, on="cluster_id", how="left")
    return wide


def _normalize_confidence_levels(levels: Optional[Iterable[str]]) -> set[str]:
    if levels is None:
        return {"low", "none"}
    return {str(x).strip().lower() for x in levels}


def _rerun_malignant_summary_from_scores(
    result,
    malignant_status_score_threshold: float = 0.35,
    malignant_raw_score_threshold: float = 0.15,
    malignant_margin_threshold: float = 0.15,
):
    malignant_scores = getattr(result, "malignant_scores", None)
    if malignant_scores is None or len(malignant_scores) == 0:
        return None
    from .malignant_scoring import _format_malignant_annotations
    return _format_malignant_annotations(
        malignant_scores=malignant_scores,
        score_threshold=float(malignant_status_score_threshold),
        raw_score_threshold=float(malignant_raw_score_threshold),
        margin_threshold=float(malignant_margin_threshold),
    )

def _rerun_decision_summary_from_scores(
    result,
    score_threshold: float = 0.05,
    margin_threshold: float = 0.02,
    branch_support_parent_weight: float = 0.65,
    branch_support_descendant_weight: float = 0.35,
    level1_branch_support_parent_weight: float = 0.60,
    level1_branch_support_descendant_weight: float = 0.40,
    branch_routing_raw_weight: float = 0.60,
    weak_raw_score_threshold: float = 0.10,
    min_branch_supported_raw_score: float = 0.0,
    min_leaf_raw_score: float = -0.2,
):
    compiled = getattr(result, "compiled_programs", None)
    all_scores = getattr(result, "all_scores", None)
    if compiled is None or all_scores is None or len(compiled) == 0 or len(all_scores) == 0:
        return None
    rows = []
    clusters = pd.Index(all_scores["cluster"].astype(str).unique())
    for cluster in clusters:
        ann = decide_cluster_annotation_from_scores(
            cluster=str(cluster),
            root_programs=compiled,
            all_scores=all_scores,
            score_threshold=score_threshold,
            margin_threshold=margin_threshold,
            branch_support_parent_weight=branch_support_parent_weight,
            branch_support_descendant_weight=branch_support_descendant_weight,
            level1_branch_support_parent_weight=level1_branch_support_parent_weight,
            level1_branch_support_descendant_weight=level1_branch_support_descendant_weight,
            branch_routing_raw_weight=branch_routing_raw_weight,
            weak_raw_score_threshold=weak_raw_score_threshold,
            min_branch_supported_raw_score=float(min_branch_supported_raw_score),
            min_leaf_raw_score=min_leaf_raw_score,
        )
        rows.append(ann.__dict__)
    if not rows:
        return None
    raw = pd.DataFrame(rows)
    expanded = pd.DataFrame(raw)
    # expand dict columns same style as pipeline
    for col in ["level_labels", "level_scores", "level_raw_scores", "level_branch_supported_raw_scores", "level_margins", "level_confidence"]:
        if col in expanded.columns:
            mapped = expanded[col].apply(lambda x: x if isinstance(x, dict) else {})
            wide = pd.DataFrame(list(mapped))
            rename_prefix = {
                "level_labels": "_label",
                "level_scores": "_score",
                "level_raw_scores": "_raw_score",
                "level_branch_supported_raw_scores": "_branch_supported_raw_score",
                "level_margins": "_margin",
                "level_confidence": "_confidence",
            }[col]
            wide.columns = [str(c) + rename_prefix for c in wide.columns]
            expanded = pd.concat([expanded.drop(columns=[col]), wide], axis=1)
    out = pd.DataFrame({
        "cluster_id": expanded["cluster"].astype(str),
        "annot_label": expanded["final_label"],
        "annot_label_with_cluster": expanded["final_label"].astype(str).map(_sanitize_label_text) + "_" + expanded["cluster"].astype(str),
        "annot_path": expanded["final_path"],
        "annot_level": expanded["final_level"],
        "annot_status": expanded["status"],
        "annot_stop_reason": expanded["stop_reason"],
        "annot_confidence": expanded["confidence"],
        "annot_score": expanded["final_score"],
        "annot_raw_score": expanded["final_raw_score"],
        "annot_branch_supported_raw_score": expanded["final_branch_supported_raw_score"],
        "annot_margin": expanded["final_margin"],
        "annot_best_any_level_label": expanded["best_label_any_level"],
        "annot_best_any_level_score": expanded["best_score_any_level"],
        "annot_low_evidence": expanded["confidence"].astype(str).isin(["low", "none"]) | expanded["final_score"].isna() | expanded["final_level"].fillna(0).astype(int).eq(0),
        "annot_low_absolute_support": pd.to_numeric(expanded["final_branch_supported_raw_score"], errors="coerce").lt(float(min_branch_supported_raw_score)).fillna(True),
        "annot_deeper_level_unresolved": expanded["stop_reason"].astype(str).isin(["below_score_threshold", "ambiguous_sibling_margin", "below_absolute_support"]) & expanded["final_level"].fillna(0).astype(int).gt(0),
        "annot_branch_conflict": False,
        "annot_decision_label": expanded["decision_label"],
        "annot_decision_path": expanded["decision_path"],
        "annot_decision_level": expanded["decision_level"],
        "annot_decision_score": expanded["decision_score"],
        "annot_decision_raw_score": expanded["decision_raw_score"],
        "annot_decision_branch_supported_raw_score": expanded["decision_branch_supported_raw_score"],
        "annot_backoff_applied": expanded["backoff_applied"],
        "annot_backoff_steps": expanded["backoff_steps"],
        "annot_backoff_reason": expanded["backoff_reason"],
    })
    keep_cols = [c for c in expanded.columns if c.startswith("level_")]
    if keep_cols:
        rename = {}
        for c in keep_cols:
            m = re.match(r"level_(\d+)_(.+)", c)
            if m:
                lvl, suffix = m.groups()
                rename[c] = f"annot_l{int(lvl)}_{suffix}"
        out = out.join(expanded[keep_cols].rename(columns=rename))
    # derive top-level competition from all_scores
    l1 = all_scores[all_scores["level"] == 1].copy()
    if not l1.empty:
        score_col = "decision_score" if "decision_score" in l1.columns else "score"
        purity_rows = []
        for cluster, sub in l1.groupby("cluster", sort=False):
            sub = sub.dropna(subset=[score_col]).sort_values(score_col, ascending=False)
            top1 = sub.iloc[0] if len(sub) >= 1 else None
            top2 = sub.iloc[1] if len(sub) >= 2 else None
            gap = float(top1[score_col] - top2[score_col]) if top1 is not None and top2 is not None else np.nan
            purity_rows.append({
                "cluster_id": str(cluster),
                "annot_primary_compartment": str(top1["label"]) if top1 is not None else np.nan,
                "annot_secondary_compartment": str(top2["label"]) if top2 is not None else np.nan,
                "annot_primary_compartment_score": float(top1[score_col]) if top1 is not None else np.nan,
                "annot_secondary_compartment_score": float(top2[score_col]) if top2 is not None else np.nan,
                "annot_top_level_gap": gap,
            })
        purity_df = pd.DataFrame(purity_rows)
        out = out.merge(purity_df, on="cluster_id", how="left")
    return out


def _flatten_compiled_programs_for_export(programs) -> List:
    out = []
    def visit(node):
        out.append(node)
        for child in getattr(node, "children", []) or []:
            visit(child)
    for root in programs or []:
        visit(root)
    return out


def _relationship_maps_for_export(programs) -> tuple[Dict[str, Set[str]], Dict[str, Set[str]]]:
    flat = _flatten_compiled_programs_for_export(programs)
    by_name = {p.name: p for p in flat}
    children_lookup = {p.name: [c.name for c in getattr(p, "children", []) or []] for p in flat}

    descendants: Dict[str, Set[str]] = {}
    def gather_desc(name: str) -> Set[str]:
        out: Set[str] = set()
        for child_name in children_lookup.get(name, []):
            out.add(child_name)
            out |= gather_desc(child_name)
        return out

    for name in by_name:
        descendants[name] = gather_desc(name)

    ancestors: Dict[str, Set[str]] = {name: set() for name in by_name}
    for parent_name, child_names in children_lookup.items():
        for child_name in child_names:
            ancestors.setdefault(child_name, set()).add(parent_name)
    changed = True
    while changed:
        changed = False
        for name in list(ancestors):
            expanded = set(ancestors[name])
            for anc in list(ancestors[name]):
                expanded |= ancestors.get(anc, set())
            if expanded != ancestors[name]:
                ancestors[name] = expanded
                changed = True
    return ancestors, descendants



def _attach_final_call_competition(
    summary: pd.DataFrame,
    all_scores: pd.DataFrame,
    compiled_programs,
    *,
    cluster_id_column: str = "cluster_id",
    final_label_column: str = "annot_label",
    final_score_column: str = "annot_score",
    path_margin_column: str = "annot_margin",
    overwrite: bool = True,
) -> pd.DataFrame:
    if summary is None or summary.empty:
        return summary
    out = summary.copy()
    if all_scores is None or len(all_scores) == 0:
        if "annot_path_margin" not in out.columns:
            out["annot_path_margin"] = out.get(path_margin_column, np.nan)
        if "annot_second_best_call" not in out.columns:
            out["annot_second_best_call"] = np.nan
        if "annot_second_best_call_score" not in out.columns:
            out["annot_second_best_call_score"] = np.nan
        if "annot_final_call_margin" not in out.columns:
            out["annot_final_call_margin"] = np.nan
        if "annot_final_call_margin" in out.columns:
            out["annot_margin"] = out["annot_final_call_margin"]
        return out

    ancestors, descendants = _relationship_maps_for_export(compiled_programs)
    score_col = "branch_support_score" if "branch_support_score" in all_scores.columns else ("decision_score" if "decision_score" in all_scores.columns else "score")

    records = []
    for _, row in out.iterrows():
        cluster = str(row[cluster_id_column]) if cluster_id_column in out.columns else str(row.get("cluster"))
        final_label = str(row.get(final_label_column, row.get("final_label", "Unresolved")))
        final_score = pd.to_numeric(pd.Series([row.get(final_score_column, np.nan)]), errors="coerce").iloc[0]

        sub = all_scores[all_scores["cluster"].astype(str) == cluster].copy()
        if sub.empty:
            records.append({
                cluster_id_column: cluster,
                "annot_second_best_call": np.nan,
                "annot_second_best_call_score": np.nan,
                "annot_final_call_margin": np.nan,
            })
            continue

        sub[score_col] = pd.to_numeric(sub[score_col], errors="coerce")
        sub = sub[np.isfinite(sub[score_col])].copy()
        if sub.empty:
            records.append({
                cluster_id_column: cluster,
                "annot_second_best_call": np.nan,
                "annot_second_best_call_score": np.nan,
                "annot_final_call_margin": np.nan,
            })
            continue

        excluded = {final_label} | ancestors.get(final_label, set()) | descendants.get(final_label, set())
        alt = sub[~sub["label"].astype(str).isin(excluded)].sort_values(score_col, ascending=False)
        if alt.empty:
            second_label = np.nan
            second_score = np.nan
            final_margin = np.nan
        else:
            second = alt.iloc[0]
            second_label = str(second["label"])
            second_score = float(second[score_col])
            final_margin = float(final_score - second_score) if np.isfinite(final_score) and np.isfinite(second_score) else np.nan

        records.append({
            cluster_id_column: cluster,
            "annot_second_best_call": second_label,
            "annot_second_best_call_score": second_score,
            "annot_final_call_margin": final_margin,
        })

    comp = pd.DataFrame(records)
    refresh_cols = ["annot_second_best_call", "annot_second_best_call_score", "annot_final_call_margin"]
    out = out.copy()
    if cluster_id_column in out.columns:
        out[cluster_id_column] = out[cluster_id_column].astype(str)
    if cluster_id_column in comp.columns:
        comp[cluster_id_column] = comp[cluster_id_column].astype(str)
    if overwrite:
        existing_refresh = [c for c in refresh_cols if c in out.columns]
        if existing_refresh:
            out = out.drop(columns=existing_refresh)
    out = out.merge(comp, on=cluster_id_column, how="left")
    if "annot_path_margin" not in out.columns:
        out["annot_path_margin"] = out.get(path_margin_column, np.nan)
    if "annot_final_call_margin" in out.columns:
        out["annot_margin"] = out["annot_final_call_margin"]
    else:
        out["annot_final_call_margin"] = np.nan
        out["annot_margin"] = np.nan
    return out


def _compute_final_call_competition(summary: pd.DataFrame, result) -> pd.DataFrame:
    all_scores = getattr(result, "all_scores", None)
    compiled = getattr(result, "compiled_programs", None)
    return _attach_final_call_competition(
        summary=summary,
        all_scores=all_scores,
        compiled_programs=compiled,
        cluster_id_column="cluster_id",
        final_label_column="annot_label",
        final_score_column="annot_score",
        path_margin_column="annot_margin",
        overwrite=True,
    )


def _attach_malignant_columns(summary: pd.DataFrame, malignant_summary: Optional[pd.DataFrame]) -> pd.DataFrame:
    """Attach malignant annotation columns onto ``summary`` without touching integrated columns."""
    out = summary.copy()
    if malignant_summary is None or len(malignant_summary) == 0:
        return out
    extra = malignant_summary.copy()
    key = "cluster_id"
    if key not in out.columns or key not in extra.columns:
        return out
    out[key] = out[key].astype(str)
    extra[key] = extra[key].astype(str)
    wanted = [c for c in extra.columns if c != key and c.startswith("annot_malignant_")]
    if not wanted:
        return out
    overlap = [c for c in wanted if c in out.columns]
    if overlap:
        out = out.drop(columns=overlap)
    out = out.merge(extra[[key] + wanted], on=key, how="left")
    return out


def _resolve_normal_summary_for_export(
    result,
    *,
    rerun_decision: bool,
    score_threshold: float,
    margin_threshold: float,
    branch_child_rescue_weight: float,
    level1_branch_child_rescue_weight: float,
    branch_routing_raw_weight: float,
    weak_raw_score_threshold: float,
    min_branch_supported_raw_score: float,
    min_leaf_raw_score: float,
) -> Optional[pd.DataFrame]:
    """Return the normal-summary source for export, optionally rerouting from stored scores."""
    if rerun_decision:
        return _rerun_decision_summary_from_scores(
            result,
            score_threshold=score_threshold,
            margin_threshold=margin_threshold,
            branch_support_parent_weight=(1.0 - float(branch_child_rescue_weight)),
            branch_support_descendant_weight=float(branch_child_rescue_weight),
            level1_branch_support_parent_weight=(1.0 - float(level1_branch_child_rescue_weight)),
            level1_branch_support_descendant_weight=float(level1_branch_child_rescue_weight),
            branch_routing_raw_weight=branch_routing_raw_weight,
            weak_raw_score_threshold=weak_raw_score_threshold,
            min_branch_supported_raw_score=float(min_branch_supported_raw_score),
            min_leaf_raw_score=float(min_leaf_raw_score),
        )
    cluster_summary = _select_annotation_summary_source(result, source="normal")
    per_level_summary = _build_per_level_summary(result)
    if per_level_summary is not None and len(per_level_summary) > 0:
        if "cluster_id" in cluster_summary.columns and "cluster_id" in per_level_summary.columns:
            per_level_summary = per_level_summary.drop(columns=["cluster_id"])
        return cluster_summary.join(per_level_summary, how="left")
    return cluster_summary.copy()


def _cached_dataframe(obj, attr: str) -> Optional[pd.DataFrame]:
    """Return a defensive copy of a cached result table when it exists."""
    df = getattr(obj, attr, None)
    if isinstance(df, pd.DataFrame) and len(df) > 0:
        return df.copy()
    return None


def _merge_cached_summary_with_normal(cached_summary: pd.DataFrame, normal_summary: Optional[pd.DataFrame]) -> pd.DataFrame:
    """Keep cached integrated columns authoritative while adding normal-only extras."""
    out = cached_summary.copy()
    if normal_summary is None or len(normal_summary) == 0:
        return out
    if "cluster_id" not in out.columns or "cluster_id" not in normal_summary.columns:
        return out
    out["cluster_id"] = out["cluster_id"].astype(str)
    normal_extra = normal_summary.copy()
    normal_extra["cluster_id"] = normal_extra["cluster_id"].astype(str)
    extra_cols = [c for c in normal_extra.columns if c != "cluster_id" and c not in out.columns]
    if extra_cols:
        out = out.merge(normal_extra[["cluster_id"] + extra_cols], on="cluster_id", how="left")
    return out


def _resolve_malignant_summary_for_export(
    result,
    *,
    rerun_decision: bool,
    malignant_integration_mode: str,
    malignant_status_score_threshold: float,
    malignant_raw_score_threshold: float,
    malignant_margin_threshold: float = 0.15,
) -> Optional[pd.DataFrame]:
    """Return malignant annotations for export.

    By default this uses the malignant annotations cached by the fitted pipeline.
    Threshold-aware malignant reformatting is only performed when
    ``rerun_decision=True``.
    """
    if str(malignant_integration_mode).lower() == "off":
        return None
    if not rerun_decision:
        cached = _cached_dataframe(result, "malignant_annotations")
        if cached is not None:
            return cached
        integrated = _cached_dataframe(result, "integrated_annotations")
        if integrated is not None and "cluster_id" in integrated.columns:
            malignant_cols = [c for c in integrated.columns if c == "cluster_id" or c.startswith("annot_malignant_")]
            if len(malignant_cols) > 1:
                return integrated[malignant_cols].copy()
        return None
    return _rerun_malignant_summary_from_scores(
        result,
        malignant_status_score_threshold=float(malignant_status_score_threshold),
        malignant_raw_score_threshold=float(malignant_raw_score_threshold),
        malignant_margin_threshold=float(malignant_margin_threshold),
    )


def _resolve_base_summary_for_export(
    result,
    *,
    rerun_decision: bool,
    score_threshold: float,
    margin_threshold: float,
    branch_child_rescue_weight: float,
    level1_branch_child_rescue_weight: float,
    branch_routing_raw_weight: float,
    weak_raw_score_threshold: float,
    min_branch_supported_raw_score: float,
    min_leaf_raw_score: float,
    normal_strong_score_threshold: float,
    normal_strong_raw_threshold: float,
    malignant_integration_mode: str,
    malignant_status_score_threshold: float,
    malignant_raw_score_threshold: float,
    malignant_normal_raw_delta_threshold: Optional[float],
    program_report_block_preset=None,
):
    """Resolve the base summary used for export.

    With the default ``rerun_decision=False``, export uses the normal,
    malignant, and integrated annotation tables cached in ``result``. Passing
    threshold arguments to the export helper does not recompute malignant or
    integrated decisions in that mode. Set ``rerun_decision=True`` to rebuild
    normal routing, malignant status formatting, and normal/malignant
    integration from the stored score tables.
    """
    normal_summary = _resolve_normal_summary_for_export(
        result,
        rerun_decision=rerun_decision,
        score_threshold=score_threshold,
        margin_threshold=margin_threshold,
        branch_child_rescue_weight=branch_child_rescue_weight,
        level1_branch_child_rescue_weight=level1_branch_child_rescue_weight,
        branch_routing_raw_weight=branch_routing_raw_weight,
        weak_raw_score_threshold=weak_raw_score_threshold,
        min_branch_supported_raw_score=min_branch_supported_raw_score,
        min_leaf_raw_score=min_leaf_raw_score,
    )
    malignant_summary = _resolve_malignant_summary_for_export(
        result,
        rerun_decision=rerun_decision,
        malignant_integration_mode=malignant_integration_mode,
        malignant_status_score_threshold=malignant_status_score_threshold,
        malignant_raw_score_threshold=malignant_raw_score_threshold,
        malignant_margin_threshold=margin_threshold,
    )

    if normal_summary is None:
        return None, None, None

    if str(malignant_integration_mode).lower() == "off":
        return normal_summary.copy(), normal_summary, malignant_summary

    if not rerun_decision:
        cached_integrated = _cached_dataframe(result, "integrated_annotations")
        if cached_integrated is not None:
            summary = _merge_cached_summary_with_normal(cached_integrated, normal_summary)
            if malignant_summary is not None:
                summary = _attach_malignant_columns(summary, malignant_summary)
            return summary, normal_summary, malignant_summary

        summary = normal_summary.copy()
        if malignant_summary is not None:
            summary = _attach_malignant_columns(summary, malignant_summary)
        return summary, normal_summary, malignant_summary

    from .integration.integrate_tracks import _integrate_normal_and_malignant_annotations
    summary = _integrate_normal_and_malignant_annotations(
        normal_summary.copy(),
        malignant_summary,
        normal_strong_score_threshold=float(normal_strong_score_threshold),
        normal_strong_raw_threshold=float(normal_strong_raw_threshold),
        malignant_status_score_threshold=float(malignant_status_score_threshold),
        malignant_raw_score_threshold=float(malignant_raw_score_threshold),
        malignant_normal_raw_delta_threshold=(None if malignant_normal_raw_delta_threshold is None else float(malignant_normal_raw_delta_threshold)),
        malignant_integration_mode=str(malignant_integration_mode),
        program_report_block_preset=program_report_block_preset,
    )
    summary = _attach_malignant_columns(summary, malignant_summary)
    return summary, normal_summary, malignant_summary


def _attach_export_diagnostics(summary: pd.DataFrame, result) -> pd.DataFrame:
    """Attach export-time diagnostics that depend on cached score tables but do not mutate ``result``."""
    out = _compute_final_call_competition(summary, result)
    out = _attach_best_any_level_raw_score(out, result)
    return out


def _reorder_export_summary_columns(summary: pd.DataFrame, malignant_integration_mode: str) -> pd.DataFrame:
    front_cols = [
        "cluster_id",
        "annot_export_label",
        "annot_export_label_with_cluster",
        "annot_export_status",
        "annot_export_reason",
        "annot_export_malignant_flag",
        "annot_export_tumor_like_label",
    ]
    if str(malignant_integration_mode).lower() != "off":
        front_cols.extend([
            "annot_malignant_name",
            "annot_malignant_label",
            "annot_malignant_label_concise",
            "annot_malignant_status",
            "annot_malignant_reason",
            "annot_integrated_label",
            "annot_integrated_status",
            "annot_integrated_reason",
            "annot_integrated_source",
        ])
    ordered_front = [c for c in front_cols if c in summary.columns]
    remaining = [c for c in summary.columns if c not in ordered_front]
    return summary[ordered_front + remaining]


def _apply_export_view(summary: pd.DataFrame, malignant_integration_mode: str, export_view: str) -> pd.DataFrame:
    """Return either a compact export table or the full diagnostic table."""
    ordered = _reorder_export_summary_columns(summary, malignant_integration_mode)
    view = str(export_view or "compact").strip().lower()
    if view in {"diagnostic", "diagnostics", "full", "all"}:
        return ordered
    if view == "minimal":
        cols = [
            "cluster_id",
            "annot_export_label",
            "annot_export_label_with_cluster",
            "annot_export_status",
            "annot_export_reason",
            "annot_export_malignant_flag",
        ]
        return ordered[[c for c in cols if c in ordered.columns]]
    if view != "compact":
        raise ValueError("export_view must be one of 'compact', 'minimal', or 'diagnostic'")
    compact_cols = [
        "cluster_id",
        "annot_export_label",
        "annot_export_label_with_cluster",
        "annot_export_status",
        "annot_export_reason",
        "annot_export_malignant_flag",
        "annot_export_tumor_like_label",
        "annot_label",
        "annot_label_with_cluster",
        "annot_status",
        "annot_stop_reason",
        "annot_path",
        "annot_level",
        "annot_confidence",
        "annot_score",
        "annot_raw_score",
        "annot_branch_supported_raw_score",
        "annot_margin",
        "annot_final_call_margin",
        "annot_second_best_call",
        "annot_second_best_call_score",
        "annot_unresolved_candidate_label",
        "annot_unresolved_candidate_score",
        "annot_unresolved_candidate_raw_score",
    ]
    if str(malignant_integration_mode).lower() != "off":
        compact_cols.extend([
            "annot_integrated_label",
            "annot_integrated_status",
            "annot_integrated_reason",
            "annot_integrated_source",
            "annot_integrated_normal_state",
            "annot_integrated_malignant_state",
            "annot_malignant_name",
            "annot_malignant_label",
            "annot_malignant_label_concise",
            "annot_malignant_status",
            "annot_malignant_reason",
            "annot_malignant_status_score",
            "annot_malignant_raw_score",
            "annot_malignant_decision_score",
            "annot_malignant_specificity_score",
            "annot_malignant_margin",
            "annot_malignant_normal_raw_delta",
            "annot_malignant_normal_raw_delta_pass",
            "annot_malignant_reporting_blocked",
        ])
    cols = [c for c in compact_cols if c in ordered.columns]
    # Keep per-level summary columns in compact output because they are commonly
    # joined back to spatial/single-cell metadata and used for downstream checks.
    per_level_cols = [c for c in ordered.columns if re.match(r"annot_l\d+_", c) and c not in cols]
    # Keep any additional export-prefixed columns added in future releases.
    extra_export_cols = [c for c in ordered.columns if c.startswith("annot_export_") and c not in cols]
    return ordered[cols + per_level_cols + extra_export_cols]


def _build_tumor_like_export_labels(
    summary: pd.DataFrame,
    normal_label: pd.Series,
    *,
    tumor_like_prefix: str = "tumor_like",
) -> pd.Series:
    """Build concise tumor-like labels while preserving the requested normal-label suffix."""
    prefix = str(tumor_like_prefix).strip() or "tumor_like"
    idx = summary.index
    state = pd.Series("", index=idx, dtype=object)
    for col in ["annot_malignant_label_concise", "annot_malignant_label", "annot_malignant_name"]:
        if col in summary.columns:
            vals = summary[col].astype(str).replace({"nan": "", "None": "", "none": ""})
            state = state.where(state.astype(str).str.len() > 0, vals)
    normal = normal_label.astype(str).replace({"": "unknown", "nan": "unknown", "None": "unknown", "none": "unknown"})
    generic = state.astype(str).str.strip().str.lower().isin({"", "nan", "none", "malignant", "unspecified", "tumor_like"})
    out = prefix + "_" + state.astype(str).str.strip() + "." + normal
    out.loc[generic] = prefix + "." + normal.loc[generic]
    return out




def _resolve_export_malignant_integration_mode(result, malignant_integration_mode, *, rerun_decision: bool = False) -> str:
    """Resolve export-time malignant integration mode.

    When ``rerun_decision=False``, the fitted pipeline configuration stored in
    ``result.resolved_config`` is authoritative whenever available. This keeps
    export summaries from silently recomputing malignant integration under new
    arguments. When ``rerun_decision=True``, explicit arguments are used, with
    the fitted configuration as the fallback for ``None``/``"auto"``.
    """
    cfg = dict(getattr(result, "resolved_config", {}) or {})
    cfg_mode = cfg.get("malignant_integration_mode")
    raw_mode = malignant_integration_mode
    raw_mode_is_auto = raw_mode is None or str(raw_mode).strip().lower() in {"auto", "inherit", "from_result"}

    if not rerun_decision and cfg_mode is not None:
        if not raw_mode_is_auto and str(raw_mode).strip().lower() != str(cfg_mode).strip().lower():
            warnings.warn(
                "malignant_integration_mode was supplied to export but rerun_decision=False; "
                "using result.resolved_config['malignant_integration_mode'] instead. "
                "Set rerun_decision=True to recompute integration under new controls.",
                UserWarning,
            )
        mode = cfg_mode
    elif raw_mode_is_auto:
        mode = cfg_mode
        if mode is None:
            integrated = getattr(result, "integrated_annotations", None)
            if isinstance(integrated, pd.DataFrame) and not integrated.empty:
                if "annot_integrated_status" in integrated.columns:
                    statuses = integrated["annot_integrated_status"].astype(str).str.lower()
                    if statuses.eq("tumor_like").any():
                        mode = "integrate"
                if mode is None and "annot_integrated_label" in integrated.columns:
                    labels = integrated["annot_integrated_label"].astype(str).str.lower()
                    if labels.str.startswith("tumor_like").any():
                        mode = "integrate"
        if mode is None:
            malignant_scores = getattr(result, "malignant_scores", None)
            malignant_annotations = getattr(result, "malignant_annotations", None)
            has_malignant = (malignant_scores is not None and len(malignant_scores) > 0) or (isinstance(malignant_annotations, pd.DataFrame) and len(malignant_annotations) > 0)
            mode = "flag_only" if has_malignant else "off"
    else:
        mode = raw_mode

    mode = str(mode).strip().lower()
    if mode not in {"off", "flag_only", "integrate"}:
        raise ValueError(f"Unsupported malignant_integration_mode: {mode}")
    return mode


def _resolve_export_label_source(export_label_source: str, malignant_integration_mode: str) -> str:
    """Resolve how tumor-like integrated calls should affect annot_export_label.

    ``auto`` is intentionally kept as its own mode: in integrate mode it uses a
    hybrid export policy where resolved normal labels become tumor-like integrated
    labels, mixed/candidate normal safeguards stay primary, and unknown normal
    labels fall back to ``tumor_like_*.unknown``.
    """
    source = str(export_label_source or "auto").strip().lower()
    aliases = {
        "normal": "normal_priority",
        "normal-priority": "normal_priority",
        "normal_priority": "normal_priority",
        "normal_primary": "normal_priority",
        "integrated": "integrated",
        "integrated_label": "integrated",
        "integration": "integrated",
        "auto": "auto",
    }
    if source not in aliases:
        raise ValueError("export_label_source must be one of 'auto', 'normal_priority', or 'integrated'")
    return aliases[source]

def rerun_cluster_annotation_result_from_scores(
    result,
    *,
    score_threshold: float = 0.05,
    margin_threshold: float = 0.02,
    branch_child_rescue_weight: float = 0.35,
    level1_branch_child_rescue_weight: float = 0.40,
    branch_routing_raw_weight: float = 0.60,
    weak_raw_score_threshold: float = 0.10,
    min_branch_supported_raw_score: float = 0.0,
    min_leaf_raw_score: float = -0.2,
    normal_strong_score_threshold: float = 0.35,
    normal_strong_raw_threshold: float = 0.10,
    malignant_integration_mode: str = "flag_only",
    malignant_status_score_threshold: float = 0.35,
    malignant_raw_score_threshold: float = 0.15,
    malignant_normal_raw_delta_threshold: Optional[float] = None,
    program_report_block_preset: Optional[Union[str, Sequence[str]]] = "tumor_reportable",
):
    """Create a new result object by rerouting existing score tables under new cutoffs.

    This function is non-mutating: the input ``result`` is not modified. A new
    ``HierAnnotResult`` is returned so callers can inspect or save a rerouted
    version of the analysis.
    """
    malignant_integration_mode = str(malignant_integration_mode).lower()
    if malignant_integration_mode not in {"off", "flag_only", "integrate"}:
        raise ValueError(f"Unsupported malignant_integration_mode: {malignant_integration_mode}")
    if malignant_integration_mode != "off":
        malignant_scores = getattr(result, "malignant_scores", None)
        if malignant_scores is None or len(malignant_scores) == 0:
            raise ValueError("Malignant integration requested but no malignant_scores found in result")
        program_report_block_preset = _validate_program_report_block_preset(program_report_block_preset)

    from .datamodels import HierAnnotResult

    summary, normal_summary, malignant_summary = _resolve_base_summary_for_export(
        result,
        rerun_decision=True,
        score_threshold=score_threshold,
        margin_threshold=margin_threshold,
        branch_child_rescue_weight=branch_child_rescue_weight,
        level1_branch_child_rescue_weight=level1_branch_child_rescue_weight,
        branch_routing_raw_weight=branch_routing_raw_weight,
        weak_raw_score_threshold=weak_raw_score_threshold,
        min_branch_supported_raw_score=min_branch_supported_raw_score,
        min_leaf_raw_score=min_leaf_raw_score,
        normal_strong_score_threshold=normal_strong_score_threshold,
        normal_strong_raw_threshold=normal_strong_raw_threshold,
        malignant_integration_mode=malignant_integration_mode,
        malignant_status_score_threshold=malignant_status_score_threshold,
        malignant_raw_score_threshold=malignant_raw_score_threshold,
        malignant_normal_raw_delta_threshold=malignant_normal_raw_delta_threshold,
        program_report_block_preset=program_report_block_preset,
    )
    if normal_summary is None:
        raise ValueError("Unable to rerun annotations from the provided result object")

    resolved_config = dict(getattr(result, "resolved_config", {}) or {})
    resolved_config.update({
        "score_threshold": float(score_threshold),
        "margin_threshold": float(margin_threshold),
        "branch_child_rescue_weight": float(branch_child_rescue_weight),
        "level1_branch_child_rescue_weight": float(level1_branch_child_rescue_weight),
        "branch_routing_raw_weight": float(branch_routing_raw_weight),
        "weak_raw_score_threshold": float(weak_raw_score_threshold),
        "min_branch_supported_raw_score": float(min_branch_supported_raw_score),
        "min_leaf_raw_score": float(min_leaf_raw_score),
        "normal_strong_score_threshold": float(normal_strong_score_threshold),
        "normal_strong_raw_threshold": float(normal_strong_raw_threshold),
        "malignant_integration_mode": str(malignant_integration_mode),
        "malignant_status_score_threshold": float(malignant_status_score_threshold),
        "malignant_raw_score_threshold": float(malignant_raw_score_threshold),
        "malignant_normal_raw_delta_threshold": (None if malignant_normal_raw_delta_threshold is None else float(malignant_normal_raw_delta_threshold)),
        "program_report_block_preset": program_report_block_preset,
    })
    metadata = dict(getattr(result, "metadata", {}) or {})
    metadata["derived_from_rerun"] = True

    return HierAnnotResult(
        cluster_annotations=normal_summary.copy(),
        level_scores=getattr(result, "level_scores", None).copy() if getattr(result, "level_scores", None) is not None else None,
        all_scores=getattr(result, "all_scores", None).copy() if getattr(result, "all_scores", None) is not None else None,
        malignant_scores=getattr(result, "malignant_scores", None).copy() if getattr(result, "malignant_scores", None) is not None else None,
        malignant_annotations=malignant_summary.copy(),
        integrated_annotations=summary.copy(),
        normalized_matrix=getattr(result, "normalized_matrix", None),
        control_gene_map=dict(getattr(result, "control_gene_map", {}) or {}),
        compiled_programs=list(getattr(result, "compiled_programs", []) or []),
        compilation_report=getattr(result, "compilation_report", None),
        diagnostics_summary=getattr(result, "diagnostics_summary", None),
        resolved_config=resolved_config,
        metadata=metadata,
        hierarchy=getattr(result, "hierarchy", None),
        malignant_programs=list(getattr(result, "malignant_programs", []) or []) or None,
    )

def _attach_best_any_level_raw_score(summary: pd.DataFrame, result) -> pd.DataFrame:
    """Attach annot_best_any_level_raw_score from result.all_scores when possible."""
    if summary is None or summary.empty:
        return summary
    if "annot_best_any_level_raw_score" in summary.columns:
        return summary
    required = {"cluster_id", "annot_best_any_level_label"}
    if not required.issubset(set(summary.columns)):
        return summary

    all_scores = getattr(result, "all_scores", None)
    if all_scores is None or len(all_scores) == 0:
        out = summary.copy()
        out["annot_best_any_level_raw_score"] = np.nan
        return out
    if not {"cluster", "label", "score"}.issubset(set(all_scores.columns)):
        out = summary.copy()
        out["annot_best_any_level_raw_score"] = np.nan
        return out

    out = summary.copy()
    lookup = all_scores[["cluster", "label", "score"]].copy()
    lookup["cluster"] = lookup["cluster"].astype(str)
    lookup["label"] = lookup["label"].astype(str)
    lookup = lookup.rename(columns={
        "cluster": "cluster_id",
        "label": "annot_best_any_level_label",
        "score": "annot_best_any_level_raw_score",
    })
    lookup = lookup.drop_duplicates(subset=["cluster_id", "annot_best_any_level_label"], keep="first")

    out["cluster_id"] = out["cluster_id"].astype(str)
    out["annot_best_any_level_label"] = out["annot_best_any_level_label"].astype(str)

    overlap = [c for c in ["annot_best_any_level_raw_score"] if c in out.columns]
    if overlap:
        out = out.drop(columns=overlap)
    out = out.merge(lookup, on=["cluster_id", "annot_best_any_level_label"], how="left")
    return out



def _build_export_labels(
    summary: pd.DataFrame,
    label_source_column: str,
    label_with_cluster_source_column: str,
    label_sep: str = "_",
    sanitize_label: bool = True,
    unknown_on_low_confidence: bool = True,
    unknown_confidence_levels: Optional[Iterable[str]] = None,
    unknown_label: str = "unknown",
    unknown_min_score: Optional[float] = None,
    unknown_max_margin: Optional[float] = None,
    unknown_score_column: str = "annot_score",
    unknown_margin_column: str = "annot_final_call_margin",
) -> pd.DataFrame:
    out = summary.copy()
    if label_source_column not in out.columns:
        raise KeyError(f"label_source_column '{label_source_column}' not found in cluster summary")
    levels = _normalize_confidence_levels(unknown_confidence_levels)
    export_label = out[label_source_column].astype(str).copy()
    if unknown_on_low_confidence and "annot_confidence" in out.columns:
        mask = out["annot_confidence"].astype(str).str.lower().isin(levels)
        export_label.loc[mask] = unknown_label
    if unknown_min_score is not None and unknown_score_column in out.columns:
        score_mask = pd.to_numeric(out[unknown_score_column], errors="coerce") < float(unknown_min_score)
        export_label.loc[score_mask.fillna(False)] = unknown_label
    if unknown_max_margin is not None and unknown_margin_column in out.columns:
        margin_mask = pd.to_numeric(out[unknown_margin_column], errors="coerce") < float(unknown_max_margin)
        export_label.loc[margin_mask.fillna(False)] = unknown_label
    out["annot_export_label"] = export_label
    if label_with_cluster_source_column in out.columns:
        cluster_ids = out["cluster_id"].astype(str)
        base = export_label.copy()
        if sanitize_label:
            base = base.map(_sanitize_label_text)
        out["annot_export_label_with_cluster"] = base + label_sep + cluster_ids
    return out


def make_cluster_annotation_export_summary(
    result,
    label_with_cluster: bool = True,
    export_view: str = "compact",
    unknown_on_low_confidence: bool = True,
    unknown_confidence_levels=("low", "none"),
    unknown_on_branch_conflict: bool = False,
    unknown_on_low_absolute_support: bool = False,
    unknown_branch_raw_threshold: Optional[float] = None,
    mixed_on_parent_mixing: bool = True,
    mixed_label_prefix: str = "mixed",
    mixed_branch_separator: str = ".",
    mixed_min_score: float = 0.35,
    unresolved_candidate_min_score: float = 1.0,
    unresolved_candidate_min_raw_score: float = 0.2,
    unresolved_candidate_label_prefix: str = "candidate",
    unknown_label: str = "unknown",
    unknown_min_score=None,
    unknown_max_margin=None,
    unknown_score_column: str = "annot_score",
    unknown_margin_column: str = "annot_final_call_margin",
    rescue_unknown_with_blocked_candidates: bool = True,
    normal_strong_score_threshold: float = 0.35,
    normal_strong_raw_threshold: float = 0.10,
    malignant_integration_mode: Optional[str] = None,
    export_label_source: str = "auto",
    malignant_status_score_threshold: float = 0.35,
    malignant_raw_score_threshold: float = 0.15,
    malignant_normal_raw_delta_threshold: Optional[float] = None,
    program_report_block_preset: Optional[Union[str, Sequence[str]]] = "tumor_reportable",
    rerun_decision: bool = False,
    score_threshold: float = 0.05,
    margin_threshold: float = 0.02,
    branch_child_rescue_weight: float = 0.35,
    level1_branch_child_rescue_weight: float = 0.40,
    branch_routing_raw_weight: float = 0.60,
    weak_raw_score_threshold: float = 0.10,
    min_branch_supported_raw_score: float = 0.0,
    min_leaf_raw_score: float = -0.2,
):
    """
    Build an export-ready cluster annotation summary from a HierAnnot result.

    This is the main public workflow helper for turning `result.cluster_annotations`
    into a flat cluster-level table with export labels, export status, optional
    rerun-based decision refresh, and extra diagnostics.
    
    This function is non-mutating. It returns a derived export table but does
    not update the input ``result``. With the default ``rerun_decision=False``,
    it uses the normal, malignant, and integrated annotation decisions cached in
    the fitted result and applies export-only masking/label formatting on top.
    When ``rerun_decision=True``, the function reroutes the normal hierarchy
    from cached score tables, reformats malignant annotations, recomputes
    normal/malignant integration under the supplied controls, and then applies
    export logic on top of that derived summary.

    Parameters controlling export masking
    ------------------------------------
    label_with_cluster
        If True, also create `annot_export_label_with_cluster` by appending the
        cluster id to the export label.

    export_view
        Controls output width. `"compact"` returns join-ready labels plus key
        normal, malignant, and integrated diagnostics; `"diagnostic"` returns
        the full derived table; `"minimal"` returns only the core export columns.

    unknown_on_low_confidence
        If True, mark clusters as unknown when `annot_low_evidence` is True.

    unknown_confidence_levels
        Reserved compatibility argument. The current export helper relies on the
        precomputed low-evidence flag in `cluster_annotations`.

    unknown_on_branch_conflict
        If True, mark clusters as unknown when `annot_branch_conflict` is True.

    unknown_on_low_absolute_support
        If True, mark clusters as unknown when the precomputed
        `annot_low_absolute_support` flag is True.

    unknown_branch_raw_threshold
        Optional export-time threshold on `annot_branch_supported_raw_score`.
        When provided, clusters below this threshold are exported as unknown
        regardless of `unknown_on_low_absolute_support`.

    unknown_label
        Export label used for unknown clusters.

    unknown_min_score
        Optional export-time lower bound on `unknown_score_column`. Clusters
        below this threshold are exported as unknown.

    unknown_max_margin
        Optional export-time upper bound on `unknown_margin_column`. Clusters
        below this ambiguity margin are exported as unknown.

    unknown_score_column
        Column used with `unknown_min_score`. Default is `annot_score`.

    unknown_margin_column
        Column used with `unknown_max_margin`. Default is
        `annot_final_call_margin`.

    rescue_unknown_with_blocked_candidates
        If True, allow blocked strong candidates to rescue unknown clusters.

    Parameters controlling mixed export labels
    -----------------------------------------
    mixed_on_parent_mixing
        If True, enable mixed-label export based on final-call competition.

    mixed_label_prefix, mixed_branch_separator
        Formatting controls for mixed export labels.

    mixed_min_score
        Minimum score required for both the final call and second-best call
        before a mixed label is emitted.

    Parameters controlling blocked-candidate hints
    ---------------------------------------------
    unresolved_candidate_min_score
        Minimum best-any-level score required to surface a blocked strong
        candidate label. This hint is only applied to unresolved-at-root cases
        and only when `annot_best_any_level_label`,
        `annot_best_any_level_score`, and `annot_best_any_level_raw_score`
        are all present.

    unresolved_candidate_min_raw_score
        Minimum raw score required for the blocked strong candidate.

    unresolved_candidate_label_prefix
        Prefix used for short blocked-candidate labels such as
        `candidate_plasma_cell`.

    Parameters controlling integration confidence
    ---------------------------------------------
    normal_strong_score_threshold, normal_strong_raw_threshold
        Integration/export confidence thresholds for the selected normal-track
        hierarchy label. They do not affect normal hierarchy routing.

    malignant_status_score_threshold, malignant_raw_score_threshold
        Integration/export confidence thresholds for flat program support. For
        tumor integration, only strong `reporting_role="status"` programs
        establish tumor-like identity; strong state/modifier programs remain
        visible in `result.malignant_annotations` but do not integrate by
        themselves.

    Parameters controlling malignant integration
    -------------------------------------------
    malignant_integration_mode
        Controls how malignant results are combined with the normal track,
        for example `"flag_only"`, `"off"` or `"integrate"`. The default
        ``None`` inherits the mode recorded in ``result.resolved_config`` when
        available, so export summaries follow the mode used by
        ``HierAnnotPipeline.fit_score`` unless explicitly overridden.

    export_label_source
        Controls whether tumor-like integrated calls overwrite the primary
        `annot_export_label`. The default `"auto"` is hybrid: in integrate mode,
        resolved normal labels become tumor-like integrated labels, mixed and
        rescued-candidate normal labels stay primary, and unknown normal labels
        fall back to `tumor_like_*.unknown`. Use `"normal_priority"` to keep
        normal/mixed/candidate labels primary whenever available while carrying
        tumor-like labels in `annot_export_tumor_like_label`; use
        `"integrated"` to force reportable tumor-like calls to become the
        primary export label even for mixed/candidate rows.

    malignant_normal_raw_delta_threshold
        Optional integration-time contrast threshold. When provided and
        `rerun_decision=True`, tumor-like integration additionally requires
        malignant raw enrichment to exceed normal-track raw evidence by at least
        this amount. With the default `rerun_decision=False`, the cached
        pipeline-time raw-delta decision is respected instead of recomputed.

    program_report_block_preset
        Controls which normal-track nodes or program-lineage combinations are
        blocked from tumor-like reporting when `rerun_decision=True`. With the
        default `rerun_decision=False`, the cached pipeline-time blocking
        decision is respected. Defaults to `"tumor_reportable"` when decisions are rerun. 
        Use `None` or `"off"` for no blocklist, `"immune_like"`, `"tumor_reportable"`, 
        `"lineage_aware"`, or a list of exact node/preset names. The `"tumor_reportable"` 
        preset is curated for built-in hierarchies; use an explicit node list for 
        custom hierarchy names.
    
    Parameters controlling rerun
    ---------------------------
    rerun_decision
        If False, use the existing normal, malignant, and integrated decisions
        stored in the result object. Thresholds and blocklist controls passed to
        the export helper are not used to recompute malignant integration in this
        mode. If True, rerun hierarchical decision logic and recompute malignant
        integration using the thresholds and routing weights passed to this
        function.

    score_threshold, margin_threshold
        Decision thresholds used only when `rerun_decision=True`.

    branch_child_rescue_weight, level1_branch_child_rescue_weight
        Descendant-rescue weights used only when `rerun_decision=True`.

    branch_routing_raw_weight, weak_raw_score_threshold
        Local evidence and gating controls used only when `rerun_decision=True`.

    min_branch_supported_raw_score
        Minimum branch-supported raw score required by the rerun decision logic.

    min_leaf_raw_score
        Minimum local raw score required for the selected node when
        `rerun_decision=True`. This is a node-level gate, distinct from
        branch-supported rescue and `min_branch_supported_raw_score`. The
        default `-0.2` allows routing to continue through mildly weak parent
        nodes when child-supported branch evidence is strong.

    Notes
    -----
    Export masking is applied on top of the stored or rerun integration table.
    With the default `export_label_source="auto"`,
    `malignant_integration_mode="integrate"` promotes tumor-like labels for
    resolved normal calls, preserves mixed and rescued-candidate labels as the
    primary export label, and uses `tumor_like_*.unknown` only when the normal
    side has no usable label. Tumor-like labels are always retained in
    `annot_export_tumor_like_label` when malignant evidence is reportable.
    """
    malignant_integration_mode = _resolve_export_malignant_integration_mode(result, malignant_integration_mode, rerun_decision=rerun_decision)
    export_label_source = _resolve_export_label_source(export_label_source, malignant_integration_mode)
    if malignant_integration_mode != "off":
        malignant_scores = getattr(result, "malignant_scores", None)
        malignant_annotations = getattr(result, "malignant_annotations", None)
        integrated_annotations = getattr(result, "integrated_annotations", None)
        has_malignant_tables = (
            (malignant_scores is not None and len(malignant_scores) > 0)
            or (isinstance(malignant_annotations, pd.DataFrame) and len(malignant_annotations) > 0)
            or (isinstance(integrated_annotations, pd.DataFrame) and len(integrated_annotations) > 0 and any(c.startswith("annot_malignant_") for c in integrated_annotations.columns))
        )
        if not has_malignant_tables:
            warnings.warn("malignant_integration_mode is not 'off' but no cached malignant annotations or scores are present in the result")
        if rerun_decision:
            program_report_block_preset = _validate_program_report_block_preset(program_report_block_preset)

    summary, normal_summary, malignant_summary = _resolve_base_summary_for_export(
        result,
        rerun_decision=rerun_decision,
        score_threshold=score_threshold,
        margin_threshold=margin_threshold,
        branch_child_rescue_weight=branch_child_rescue_weight,
        level1_branch_child_rescue_weight=level1_branch_child_rescue_weight,
        branch_routing_raw_weight=branch_routing_raw_weight,
        weak_raw_score_threshold=weak_raw_score_threshold,
        min_branch_supported_raw_score=min_branch_supported_raw_score,
        min_leaf_raw_score=min_leaf_raw_score,
        normal_strong_score_threshold=normal_strong_score_threshold,
        normal_strong_raw_threshold=normal_strong_raw_threshold,
        malignant_integration_mode=malignant_integration_mode,
        malignant_status_score_threshold=malignant_status_score_threshold,
        malignant_raw_score_threshold=malignant_raw_score_threshold,
        malignant_normal_raw_delta_threshold=malignant_normal_raw_delta_threshold,
        program_report_block_preset=program_report_block_preset,
    )
    if summary is None:
        raise ValueError("Unable to resolve a cluster annotation summary from the provided result")

    summary = _attach_export_diagnostics(summary, result)
    use_integrated_source = malignant_integration_mode != "off"

    label_source_column = "annot_label" if "annot_label" in summary.columns else "annot_integrated_label"
    export_labels = _build_export_labels(
        summary,
        label_source_column=label_source_column,
        label_with_cluster_source_column="annot_label_with_cluster",
        unknown_on_low_confidence=unknown_on_low_confidence,
        unknown_confidence_levels=unknown_confidence_levels,
        unknown_label=unknown_label,
        unknown_min_score=unknown_min_score,
        unknown_max_margin=unknown_max_margin,
        unknown_score_column=unknown_score_column,
        unknown_margin_column=unknown_margin_column,
    )

    export_status = pd.Series("resolved", index=summary.index, dtype=object)
    export_reason = pd.Series("", index=summary.index, dtype=object)

    malignant_state = (
        summary["annot_integrated_malignant_state"].astype(str).str.lower()
        if "annot_integrated_malignant_state" in summary.columns
        else (summary["annot_malignant_status"].astype(str).str.lower() if "annot_malignant_status" in summary.columns else pd.Series("none", index=summary.index))
    )
    blocked = summary["annot_malignant_reporting_blocked"].fillna(False).astype(bool) if "annot_malignant_reporting_blocked" in summary.columns else pd.Series(False, index=summary.index)
    raw_delta_pass = summary["annot_malignant_normal_raw_delta_pass"].fillna(False).astype(bool) if "annot_malignant_normal_raw_delta_pass" in summary.columns else pd.Series(True, index=summary.index)
    integrated_status = summary["annot_integrated_status"].astype(str).str.lower() if "annot_integrated_status" in summary.columns else pd.Series("", index=summary.index)
    integrated_label = summary["annot_integrated_label"].astype(str) if "annot_integrated_label" in summary.columns else pd.Series("", index=summary.index)

    cached_tumor_like_mask = integrated_status.eq("tumor_like") | integrated_label.str.lower().str.startswith("tumor_like")
    cached_reportable_mask = malignant_state.eq("strong") & ~blocked & raw_delta_pass
    malignant_reportable_mask = use_integrated_source & (cached_tumor_like_mask | cached_reportable_mask)
    tumor_like_integrated_mask = use_integrated_source & (str(malignant_integration_mode).lower() == "integrate") & (cached_tumor_like_mask | cached_reportable_mask)
    export_labels["annot_export_malignant_flag"] = malignant_reportable_mask.astype(bool)
    export_labels["annot_export_tumor_like_label"] = pd.Series(np.nan, index=summary.index, dtype=object)

    unknown_mask = pd.Series(False, index=summary.index)
    if unknown_on_low_confidence and "annot_low_evidence" in summary.columns:
        unknown_mask = unknown_mask | summary["annot_low_evidence"].fillna(False).astype(bool)
    if unknown_on_branch_conflict and "annot_branch_conflict" in summary.columns:
        unknown_mask = unknown_mask | summary["annot_branch_conflict"].fillna(False).astype(bool)
    if unknown_on_low_absolute_support and "annot_low_absolute_support" in summary.columns:
        unknown_mask = unknown_mask | summary["annot_low_absolute_support"].fillna(False).astype(bool)
    if unknown_branch_raw_threshold is not None and "annot_branch_supported_raw_score" in summary.columns:
        unknown_mask = unknown_mask | pd.to_numeric(summary["annot_branch_supported_raw_score"], errors="coerce").lt(float(unknown_branch_raw_threshold)).fillna(False)
    if unknown_min_score is not None and unknown_score_column in summary.columns:
        unknown_mask = unknown_mask | pd.to_numeric(summary[unknown_score_column], errors="coerce").lt(float(unknown_min_score)).fillna(False)
    if unknown_max_margin is not None and unknown_margin_column in summary.columns:
        unknown_mask = unknown_mask | pd.to_numeric(summary[unknown_margin_column], errors="coerce").lt(float(unknown_max_margin)).fillna(False)

    mixed_mask = pd.Series(False, index=summary.index)
    required_mixed_cols = {"annot_final_call_margin", "annot_score", "annot_second_best_call_score", "annot_second_best_call"}
    if mixed_on_parent_mixing and required_mixed_cols.issubset(set(summary.columns)):
        final_margin = pd.to_numeric(summary["annot_final_call_margin"], errors="coerce")
        final_score = pd.to_numeric(summary["annot_score"], errors="coerce")
        second_score = pd.to_numeric(summary["annot_second_best_call_score"], errors="coerce")
        final_top = summary["annot_final_top_branch"].astype(str) if "annot_final_top_branch" in summary.columns else pd.Series("", index=summary.index)
        second_top = summary["annot_second_best_call"].astype(str).str.split(ANNOT_PATH_SEPARATOR).str[0]
        effective_margin_threshold = float(margin_threshold)
        mixed_mask = (
            final_margin.notna()
            & final_margin.le(effective_margin_threshold)
            & final_score.notna()
            & second_score.notna()
            & final_score.ge(float(mixed_min_score))
            & second_score.ge(float(mixed_min_score))
            & final_top.ne(second_top)
        )
        if "annot_parent_mixing" in summary.columns:
            mixed_mask = mixed_mask | summary["annot_parent_mixing"].fillna(False).astype(bool)

    # resolved by default; mixed can apply, but unknown explicitly overrides mixed at the end
    if mixed_mask.any():
        primary = summary["annot_label"].astype(str).map(_sanitize_label_text).str.lower()
        secondary = summary["annot_second_best_call"].astype(str).map(_sanitize_label_text).str.lower()
        mixed_label = mixed_label_prefix + "_" + primary + mixed_branch_separator + secondary
        export_status.loc[mixed_mask] = "mixed"
        export_reason.loc[mixed_mask] = "final_call_competition"
        export_labels.loc[mixed_mask, "annot_export_label"] = mixed_label[mixed_mask]
        if "annot_export_label_with_cluster" in export_labels.columns:
            cid = summary["cluster_id"].astype(str) if "cluster_id" in summary.columns else summary.index.astype(str)
            export_labels.loc[mixed_mask, "annot_export_label_with_cluster"] = mixed_label[mixed_mask] + "_" + cid[mixed_mask]

    export_status.loc[unknown_mask] = "unknown"
    export_reason.loc[unknown_mask] = "unknown_rule"
    export_labels.loc[unknown_mask, "annot_export_label"] = unknown_label
    if "annot_export_label_with_cluster" in export_labels.columns:
        cid = summary["cluster_id"].astype(str) if "cluster_id" in summary.columns else summary.index.astype(str)
        export_labels.loc[unknown_mask, "annot_export_label_with_cluster"] = unknown_label + "_" + cid[unknown_mask]

    # Surface blocked strong candidates only for true unresolved-at-root cases.
    required_candidate_cols = {"annot_best_any_level_label", "annot_best_any_level_score", "annot_best_any_level_raw_score"}
    if required_candidate_cols.issubset(set(summary.columns)):
        best_any_label = summary["annot_best_any_level_label"].astype(str)
        best_any_score = pd.to_numeric(summary["annot_best_any_level_score"], errors="coerce")
        best_any_raw = pd.to_numeric(summary["annot_best_any_level_raw_score"], errors="coerce")
        annot_status = summary["annot_status"].astype(str).str.lower() if "annot_status" in summary.columns else pd.Series("", index=summary.index)
        annot_label = summary["annot_label"].astype(str).str.lower() if "annot_label" in summary.columns else pd.Series("", index=summary.index)
        annot_level = pd.to_numeric(summary["annot_level"], errors="coerce") if "annot_level" in summary.columns else pd.Series(np.nan, index=summary.index)

        root_unresolved_mask = (
            (annot_level.isna() | annot_level.le(0))
            & (annot_label.eq("unresolved") | annot_status.isin(["unresolved", "stopped_at_parent"]))
        )

        # allow candidate to rescue unknown ones (resovled but weak)
        if rescue_unknown_with_blocked_candidates:
            root_unresolved_mask = root_unresolved_mask | unknown_mask

        # Record candidate information for all root-unresolved rows, but only rewrite
        # export labels when the candidate thresholds are satisfied.
        export_labels["annot_unresolved_candidate_label"] = pd.Series(np.nan, index=export_labels.index, dtype=object)
        export_labels["annot_unresolved_candidate_score"] = np.nan
        export_labels["annot_unresolved_candidate_raw_score"] = np.nan

        if root_unresolved_mask.any():
            export_labels.loc[root_unresolved_mask, "annot_unresolved_candidate_label"] = best_any_label[root_unresolved_mask]
            export_labels.loc[root_unresolved_mask, "annot_unresolved_candidate_score"] = best_any_score[root_unresolved_mask]
            export_labels.loc[root_unresolved_mask, "annot_unresolved_candidate_raw_score"] = best_any_raw[root_unresolved_mask]

        candidate_mask = (
            root_unresolved_mask
            & best_any_label.notna()
            & best_any_label.ne("")
            & best_any_label.str.lower().ne("nan")
            & best_any_score.ge(float(unresolved_candidate_min_score)).fillna(False)
            & best_any_raw.ge(float(unresolved_candidate_min_raw_score)).fillna(False)
        )
        if candidate_mask.any():
            cand_label = best_any_label.map(_sanitize_label_text).str.lower()
            short_label = str(unresolved_candidate_label_prefix) + "_" + cand_label
            export_labels.loc[candidate_mask, "annot_export_label"] = short_label[candidate_mask]
            if "annot_export_label_with_cluster" in export_labels.columns:
                cid = summary["cluster_id"].astype(str) if "cluster_id" in summary.columns else summary.index.astype(str)
                export_labels.loc[candidate_mask, "annot_export_label_with_cluster"] = short_label[candidate_mask] + "_" + cid[candidate_mask]
            export_status.loc[candidate_mask] = "candidate"
            export_reason.loc[candidate_mask] = "blocked_strong_candidate"

    # Build tumor-like export labels after normal export logic so their normal
    # suffix reflects mixed-label and candidate-rescue decisions.
    if malignant_reportable_mask.any():
        suffix = export_labels["annot_export_label"].astype(str).replace({"": unknown_label, "nan": unknown_label, "None": unknown_label, "none": unknown_label})
        tumor_export_labels = _build_tumor_like_export_labels(summary, suffix)
        export_labels.loc[malignant_reportable_mask, "annot_export_tumor_like_label"] = tumor_export_labels[malignant_reportable_mask]
    else:
        tumor_export_labels = pd.Series(np.nan, index=summary.index, dtype=object)

    if use_integrated_source and str(malignant_integration_mode).lower() == "integrate":
        current_label = export_labels["annot_export_label"].astype(str).str.lower()
        current_status = export_status.astype(str).str.lower()
        unknown_value = str(unknown_label).lower()
        unknown_export_mask = current_status.eq("unknown") | current_label.eq(unknown_value)
        resolved_export_mask = current_status.eq("resolved") & ~unknown_export_mask

        if export_label_source == "integrated":
            tumor_export_mask = tumor_like_integrated_mask
        elif export_label_source == "normal_priority":
            tumor_export_mask = tumor_like_integrated_mask & unknown_export_mask
        else:  # auto: integrated for resolved calls, normal safeguards for mixed/candidate, malignant fallback for unknown
            tumor_export_mask = tumor_like_integrated_mask & (resolved_export_mask | unknown_export_mask)

        if tumor_export_mask.any():
            export_labels.loc[tumor_export_mask, "annot_export_label"] = tumor_export_labels[tumor_export_mask]
            if "annot_export_label_with_cluster" in export_labels.columns:
                cid = summary["cluster_id"].astype(str) if "cluster_id" in summary.columns else summary.index.astype(str)
                export_labels.loc[tumor_export_mask, "annot_export_label_with_cluster"] = tumor_export_labels[tumor_export_mask].map(_sanitize_label_text) + "_" + cid[tumor_export_mask]
            export_status.loc[tumor_export_mask] = "tumor_like"
            fallback_mask = tumor_export_mask & unknown_export_mask
            integrated_mask = tumor_export_mask & ~unknown_export_mask
            export_reason.loc[integrated_mask] = "malignant_integrated_label"
            export_reason.loc[fallback_mask] = "malignant_fallback_unknown_normal"

    if not label_with_cluster and "annot_export_label_with_cluster" in export_labels.columns:
        export_labels = export_labels.drop(columns=["annot_export_label_with_cluster"])

    export_labels["annot_export_status"] = export_status
    export_labels["annot_export_reason"] = export_reason.replace("", "resolved")
    overlap = [c for c in export_labels.columns if c in summary.columns]
    if overlap:
        summary = summary.drop(columns=overlap)
    summary = summary.join(export_labels, how="left")
    return _apply_export_view(summary, malignant_integration_mode, export_view)

def expand_cluster_annotation_to_cells(
    cluster_assignments,
    cluster_summary,
    cluster_key: str = "cluster",
    include_columns=None,
    prefix: str | None = None,
    label_with_cluster: bool = True,
    fill_unassigned: bool = True,
    unassigned_label: str = "unassigned",
    unassigned_status: str = "not_analyzed",
):
    """
    Expand a cluster-level annotation summary to a per-cell joinable DataFrame.

    Parameters
    ----------
    cluster_assignments
        Array-like or pandas Series of cluster IDs for each cell/observation.
    cluster_summary
        Complete cluster-level summary table, typically produced by
        `make_cluster_annotation_export_summary()`.
    include_columns
        Optional list of columns from `cluster_summary` to propagate. If omitted,
        uses the common export-label workflow: `annot_export_label` and, when
        `label_with_cluster=True`, `annot_export_label_with_cluster`.
    prefix
        Optional prefix applied to propagated annotation columns in the returned
        table. Leading `annot_` is stripped before adding the custom prefix.
    """

    if not isinstance(cluster_assignments, pd.Series):
        cluster_assignments = pd.Series(cluster_assignments, name=cluster_key)
    else:
        cluster_assignments = cluster_assignments.copy()
        if cluster_assignments.name is None:
            cluster_assignments.name = cluster_key

    summary = cluster_summary.copy()

    if "cluster_id" not in summary.columns:
        summary["cluster_id"] = summary.index

    if include_columns is None:
        include_columns = ["annot_export_label", "annot_export_status"]
    else:
        include_columns = list(include_columns)

    if label_with_cluster and "annot_export_label_with_cluster" in summary.columns:
        if "annot_export_label" in include_columns and "annot_export_label_with_cluster" not in include_columns:
            include_columns.append("annot_export_label_with_cluster")
        elif include_columns == []:
            include_columns = ["annot_export_label_with_cluster", "annot_export_status"]

    include_columns = [c for c in include_columns if c in summary.columns]
    summary_small = summary[["cluster_id"] + include_columns].copy()

    # Merge on a normalized string view of cluster ids to avoid dtype mismatch
    join_df = pd.DataFrame({cluster_key: cluster_assignments.values}, index=cluster_assignments.index)
    join_df["_cluster_merge_key"] = join_df[cluster_key].astype(str)
    summary_small["_cluster_merge_key"] = summary_small["cluster_id"].astype(str)

    expanded = join_df.merge(summary_small, how="left", on="_cluster_merge_key")
    expanded.index = join_df.index

    if fill_unassigned:
        if "annot_export_label" in expanded.columns:
            expanded["annot_export_label"] = expanded["annot_export_label"].fillna(unassigned_label)
        if "annot_export_label_with_cluster" in expanded.columns:
            missing = expanded["annot_export_label_with_cluster"].isna()
            expanded.loc[missing, "annot_export_label_with_cluster"] = (
                unassigned_label + "_" + expanded.loc[missing, cluster_key].astype(str)
            )
        if "annot_status" in expanded.columns:
            expanded["annot_status"] = expanded["annot_status"].fillna(unassigned_status)

    expanded = expanded.drop(columns=["cluster_id", "_cluster_merge_key"], errors="ignore")

    if prefix:
        rename_map = {}
        for c in expanded.columns:
            if c == cluster_key:
                continue
            base = c[len("annot_"):] if c.startswith("annot_") else c
            rename_map[c] = f"{prefix}_{base}"
        expanded = expanded.rename(columns=rename_map)

    return expanded

