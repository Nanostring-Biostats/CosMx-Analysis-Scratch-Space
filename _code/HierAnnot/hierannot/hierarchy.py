
from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from .datamodels import ClusterAnnotation, CompiledMarkerProgram


def pick_best_label(score_df: pd.DataFrame, score_col: str = "decision_score") -> Tuple[Optional[str], float, float, float]:
    valid = score_df.dropna(subset=[score_col]).sort_values(score_col, ascending=False)
    if valid.empty:
        return None, np.nan, np.nan, np.nan
    best = valid.iloc[0]
    second_score = float(valid.iloc[1][score_col]) if len(valid) > 1 else np.nan
    margin = float(best[score_col] - second_score) if not np.isnan(second_score) else np.nan
    raw_score = float(best.get("score", np.nan))
    return str(best["label"]), float(best[score_col]), margin, raw_score


def confidence_from_score(
    score: float,
    margin: float,
    score_threshold: float,
    margin_threshold: float,
    marker_support_fraction: float = 1.0,
) -> str:
    if np.isnan(score):
        return "none"
    if score < score_threshold or marker_support_fraction < 0.2:
        return "low"
    score_high = score >= (score_threshold + max(abs(score_threshold), 0.05))
    margin_ok = np.isnan(margin) or margin >= margin_threshold
    margin_high = np.isnan(margin) or margin >= (margin_threshold + max(abs(margin_threshold), 0.05))
    if marker_support_fraction < 0.5:
        return "medium" if margin_ok else "low"
    if score_high and margin_high:
        return "high"
    if margin_ok:
        return "medium"
    return "low"


def _raw_support_gate(raw: float, weak_high: float = 0.10) -> float:
    if not np.isfinite(raw):
        return 0.0
    if raw < 0:
        return 0.0
    if weak_high <= 0:
        return 1.0
    if raw >= weak_high:
        return 1.0
    return float(raw / weak_high)


def _combine_pair(parent_val: float, child_val: float, parent_weight: float, child_weight: float) -> float:
    if np.isnan(parent_val) and np.isnan(child_val):
        return np.nan
    if np.isnan(child_val):
        return float(parent_val)
    if np.isnan(parent_val):
        return float(child_val)
    return float(parent_weight * parent_val + child_weight * child_val)


def _make_node_lookup(cluster: str, all_scores: pd.DataFrame, score_col: str) -> Dict[str, dict]:
    sub = all_scores[all_scores["cluster"] == cluster].copy()
    out: Dict[str, dict] = {}
    if sub.empty:
        return out
    for _, row in sub.iterrows():
        label = str(row["label"])
        out[label] = {
            "raw": float(row.get("score", np.nan)) if pd.notna(row.get("score", np.nan)) else np.nan,
            "decision": float(row.get(score_col, np.nan)) if pd.notna(row.get(score_col, np.nan)) else np.nan,
            "support": float(row.get("markers_detection_support_fraction", 0.0)) if pd.notna(row.get("markers_detection_support_fraction", 0.0)) else 0.0,
            "level": row.get("level", np.nan),
            "parent_label": row.get("parent_label", np.nan),
        }
    return out


def _compute_recursive_branch_maps(
    current_programs: List[CompiledMarkerProgram],
    node_lookup: Dict[str, dict],
    parent_weight: float = 0.65,
    descendant_weight: float = 0.35,
    raw_weight: float = 0.6,
    weak_raw_score_threshold: float = 0.10,
) -> Dict[str, Dict[str, float]]:
    branch_score_map: Dict[str, float] = {}
    branch_raw_map: Dict[str, float] = {}
    local_score_map: Dict[str, float] = {}
    best_child_branch_map: Dict[str, float] = {}
    best_child_raw_map: Dict[str, float] = {}
    gate_map: Dict[str, float] = {}
    support_scale_map: Dict[str, float] = {}

    def visit(program: CompiledMarkerProgram) -> Tuple[float, float]:
        meta = node_lookup.get(program.name, {})
        parent_raw = float(meta.get("raw", np.nan))
        parent_decision = float(meta.get("decision", np.nan))
        support = float(meta.get("support", 0.0))
        support_scale = 0.5 + 0.5 * max(0.0, min(1.0, support))
        gate = _raw_support_gate(parent_raw, weak_high=weak_raw_score_threshold)

        # Revised local score:
        # keep raw component visible, but gate only the decision/specificity part.
        local_score = np.nan
        raw_part = np.nan
        gated_decision_part = np.nan
        if np.isfinite(parent_raw):
            raw_part = float(raw_weight * parent_raw)
        if np.isfinite(parent_decision):
            gated_decision_part = float((1.0 - raw_weight) * gate * parent_decision)
        local_score = _combine_pair(raw_part, gated_decision_part, 1.0, 1.0)
        if np.isfinite(local_score):
            local_score = float(local_score * support_scale)

        child_branch_scores: List[float] = []
        child_branch_raws: List[float] = []
        for child in program.children:
            child_branch_score, child_branch_raw = visit(child)
            if np.isfinite(child_branch_score):
                child_branch_scores.append(float(max(0.0, child_branch_score)))
            if np.isfinite(child_branch_raw):
                child_branch_raws.append(float(max(0.0, child_branch_raw)))

        best_child_branch = max(child_branch_scores) if child_branch_scores else np.nan
        best_child_raw = max(child_branch_raws) if child_branch_raws else np.nan

        # Recursive propagation uses best child branch-supported score directly
        # (positive-only), without re-gating the propagated term.
        branch_score = _combine_pair(
            parent_val=local_score,
            child_val=max(0.0, best_child_branch) if np.isfinite(best_child_branch) else np.nan,
            parent_weight=parent_weight,
            child_weight=descendant_weight,
        )
        branch_raw = _combine_pair(
            parent_val=parent_raw,
            child_val=max(0.0, best_child_raw) if np.isfinite(best_child_raw) else np.nan,
            parent_weight=parent_weight,
            child_weight=descendant_weight,
        )

        local_score_map[program.name] = local_score
        best_child_branch_map[program.name] = best_child_branch
        best_child_raw_map[program.name] = best_child_raw
        branch_score_map[program.name] = branch_score
        branch_raw_map[program.name] = branch_raw
        gate_map[program.name] = gate
        support_scale_map[program.name] = support_scale
        return branch_score, branch_raw

    for program in current_programs:
        visit(program)

    return {
        "branch_score_map": branch_score_map,
        "branch_raw_map": branch_raw_map,
        "local_score_map": local_score_map,
        "best_child_branch_map": best_child_branch_map,
        "best_child_raw_map": best_child_raw_map,
        "gate_map": gate_map,
        "support_scale_map": support_scale_map,
    }


def compute_cluster_branch_diagnostics(
    cluster: str,
    root_programs: List[CompiledMarkerProgram],
    all_scores: pd.DataFrame,
    score_col: str = "decision_score",
    parent_weight: float = 0.65,
    descendant_weight: float = 0.35,
    raw_weight: float = 0.6,
    weak_raw_score_threshold: float = 0.10,
) -> pd.DataFrame:
    node_lookup = _make_node_lookup(cluster=cluster, all_scores=all_scores, score_col=score_col)
    maps = _compute_recursive_branch_maps(
        current_programs=root_programs,
        node_lookup=node_lookup,
        parent_weight=parent_weight,
        descendant_weight=descendant_weight,
        raw_weight=raw_weight,
        weak_raw_score_threshold=weak_raw_score_threshold,
    )
    rows = []
    for label, meta in node_lookup.items():
        rows.append(
            {
                "cluster": cluster,
                "label": label,
                "local_node_score": maps["local_score_map"].get(label, np.nan),
                "best_child_branch_support_score": maps["best_child_branch_map"].get(label, np.nan),
                "best_child_branch_supported_raw_score": maps["best_child_raw_map"].get(label, np.nan),
                "branch_support_score": maps["branch_score_map"].get(label, np.nan),
                "branch_supported_raw_score": maps["branch_raw_map"].get(label, np.nan),
                "raw_support_gate": maps["gate_map"].get(label, np.nan),
                "support_scale": maps["support_scale_map"].get(label, np.nan),
                "child_rescue_weight_used": float(descendant_weight),
                "parent_weight_used": float(parent_weight),
                "raw_weight_used": float(raw_weight),
            }
        )
    return pd.DataFrame(rows)


def _branch_support_rows(
    cluster: str,
    current_programs: List[CompiledMarkerProgram],
    level_df: pd.DataFrame,
    all_scores: pd.DataFrame,
    score_col: str = "decision_score",
    parent_weight: float = 0.65,
    descendant_weight: float = 0.35,
    raw_weight: float = 0.6,
    weak_raw_score_threshold: float = 0.10,
) -> pd.DataFrame:
    node_lookup = _make_node_lookup(cluster=cluster, all_scores=all_scores, score_col=score_col)
    maps = _compute_recursive_branch_maps(
        current_programs=current_programs,
        node_lookup=node_lookup,
        parent_weight=parent_weight,
        descendant_weight=descendant_weight,
        raw_weight=raw_weight,
        weak_raw_score_threshold=weak_raw_score_threshold,
    )
    rows = []
    for program in current_programs:
        meta = node_lookup.get(program.name, {})
        rows.append(
            {
                "label": program.name,
                "score": float(meta.get("raw", np.nan)),
                "decision_score": float(meta.get("decision", np.nan)),
                "markers_detection_support_fraction": float(meta.get("support", 0.0)),
                "best_child_score": maps["best_child_branch_map"].get(program.name, np.nan),
                "best_child_raw_score": maps["best_child_raw_map"].get(program.name, np.nan),
                "local_node_score": maps["local_score_map"].get(program.name, np.nan),
                "child_rescue_score": maps["best_child_branch_map"].get(program.name, np.nan),
                "branch_support_score": maps["branch_score_map"].get(program.name, np.nan),
                "branch_supported_raw_score": maps["branch_raw_map"].get(program.name, np.nan),
                "raw_support_gate": maps["gate_map"].get(program.name, np.nan),
                "child_rescue_weight_used": float(descendant_weight),
                "parent_weight_used": float(parent_weight),
            }
        )
    return pd.DataFrame(rows)


def _apply_absolute_support_backoff_row(
    row: dict,
    min_branch_supported_raw_score: float = 0.0,
    min_leaf_raw_score: float = 0.0,
) -> dict:
    out = dict(row)
    decision_level = int(out.get("decision_level", 0) or 0)
    out["final_label"] = out.get("decision_label", "Unresolved")
    out["final_path"] = out.get("decision_path", "Unresolved")
    out["final_level"] = decision_level
    out["final_score"] = out.get("decision_score", np.nan)
    out["final_raw_score"] = out.get("decision_raw_score", np.nan)
    out["final_branch_supported_raw_score"] = out.get("decision_branch_supported_raw_score", np.nan)
    out["final_margin"] = out.get(f"level_{decision_level}_margin", np.nan) if decision_level > 0 else np.nan
    out["status"] = out.get("decision_status", out.get("status", "assigned"))
    out["stop_reason"] = out.get("decision_stop_reason", out.get("stop_reason", "leaf_reached"))
    out["confidence"] = out.get(f"level_{decision_level}_confidence", out.get("confidence", "none")) if decision_level > 0 else out.get("confidence", "none")
    out["backoff_applied"] = False
    out["backoff_steps"] = 0
    out["backoff_reason"] = ""

    if decision_level <= 0:
        return out

    branch_raw = float(out.get(f"level_{decision_level}_branch_supported_raw_score", out.get("decision_branch_supported_raw_score", np.nan)))
    leaf_raw = float(out.get(f"level_{decision_level}_raw_score", out.get("decision_raw_score", np.nan)))

    if np.isfinite(branch_raw) and branch_raw >= float(min_branch_supported_raw_score) and np.isfinite(leaf_raw) and leaf_raw >= float(min_leaf_raw_score):
        return out

    fallback_level = None
    for lev in range(decision_level - 1, 0, -1):
        candidate_branch_raw = float(out.get(f"level_{lev}_branch_supported_raw_score", np.nan))
        candidate_raw = float(out.get(f"level_{lev}_raw_score", np.nan))
        if np.isfinite(candidate_branch_raw) and candidate_branch_raw >= float(min_branch_supported_raw_score) and np.isfinite(candidate_raw) and candidate_raw >= float(min_leaf_raw_score):
            fallback_level = lev
            break

    if fallback_level is None:
        return out

    labels = [str(out.get(f"level_{lev}_label")) for lev in range(1, fallback_level + 1) if pd.notna(out.get(f"level_{lev}_label"))]
    out["final_level"] = int(fallback_level)
    out["final_label"] = out.get(f"level_{fallback_level}_label", out["final_label"])
    out["final_path"] = " > ".join(labels) if labels else str(out["final_label"])
    out["final_score"] = float(out.get(f"level_{fallback_level}_score", np.nan))
    out["final_raw_score"] = float(out.get(f"level_{fallback_level}_raw_score", np.nan))
    out["final_branch_supported_raw_score"] = float(out.get(f"level_{fallback_level}_branch_supported_raw_score", np.nan))
    out["final_margin"] = float(out.get(f"level_{fallback_level}_margin", np.nan))
    out["confidence"] = out.get(f"level_{fallback_level}_confidence", out.get("confidence", "none"))
    out["status"] = "assigned"
    out["stop_reason"] = "backed_off_low_absolute_support"
    out["backoff_applied"] = True
    out["backoff_steps"] = int(decision_level - fallback_level)
    out["backoff_reason"] = "low_absolute_support"
    return out


def decide_cluster_annotation_from_scores(
    cluster: str,
    root_programs: List[CompiledMarkerProgram],
    all_scores: pd.DataFrame,
    score_threshold: float = 0.05,
    margin_threshold: float = 0.02,
    branch_support_parent_weight: float = 0.65,
    branch_support_descendant_weight: float = 0.35,
    level1_branch_support_parent_weight: float = 0.60,
    level1_branch_support_descendant_weight: float = 0.40,
    branch_routing_raw_weight: float = 0.60,
    weak_raw_score_threshold: float = 0.10,
    min_branch_supported_raw_score: float = 0.0,
    min_leaf_raw_score: float = 0.0,
) -> ClusterAnnotation:
    current_programs = root_programs
    level = 1
    path: List[str] = []
    level_labels: Dict[str, str] = {}
    level_scores: Dict[str, float] = {}
    level_raw_scores: Dict[str, float] = {}
    level_branch_supported_raw_scores: Dict[str, float] = {}
    level_margins: Dict[str, float] = {}
    level_confidence: Dict[str, str] = {}
    status = "assigned"
    stop_reason = "leaf_reached"
    ambiguity_flag = False

    score_col = "decision_score" if "decision_score" in all_scores.columns else "score"

    while current_programs:
        level_names = {p.name for p in current_programs}
        level_df = all_scores[
            (all_scores["cluster"] == cluster)
            & (all_scores["label"].isin(level_names))
            & (all_scores["level"] == level)
        ].copy()

        current_parent_weight = level1_branch_support_parent_weight if level == 1 else branch_support_parent_weight
        current_descendant_weight = level1_branch_support_descendant_weight if level == 1 else branch_support_descendant_weight

        branch_df = _branch_support_rows(
            cluster=cluster,
            current_programs=current_programs,
            level_df=level_df,
            all_scores=all_scores,
            score_col=score_col,
            parent_weight=current_parent_weight,
            descendant_weight=current_descendant_weight,
            raw_weight=branch_routing_raw_weight,
            weak_raw_score_threshold=weak_raw_score_threshold,
        )
        if branch_df.empty:
            status = "unresolved"
            stop_reason = "no_valid_score"
            ambiguity_flag = True
            break

        best_label, best_score, best_margin, best_raw_score = pick_best_label(branch_df, score_col="branch_support_score")
        level_key = f"level_{level}"
        level_labels[level_key] = best_label if best_label is not None else "Unresolved"
        level_scores[level_key] = best_score
        level_raw_scores[level_key] = best_raw_score
        best_branch_supported_raw = float(branch_df.loc[branch_df["label"] == best_label, "branch_supported_raw_score"].iloc[0]) if best_label is not None and not branch_df.loc[branch_df["label"] == best_label].empty else np.nan
        level_branch_supported_raw_scores[level_key] = best_branch_supported_raw
        level_margins[level_key] = best_margin
        marker_support = float(branch_df.loc[branch_df["label"] == best_label, "markers_detection_support_fraction"].iloc[0]) if best_label is not None and not branch_df.loc[branch_df["label"] == best_label].empty else 0.0
        level_confidence[level_key] = confidence_from_score(best_score, best_margin, score_threshold, margin_threshold, marker_support)

        if best_label is None or np.isnan(best_score):
            status = "unresolved"
            stop_reason = "no_valid_score"
            ambiguity_flag = True
            break
        if best_score < score_threshold:
            status = "stopped_at_parent"
            stop_reason = "below_score_threshold"
            ambiguity_flag = True
            break
        if (np.isfinite(best_branch_supported_raw) and best_branch_supported_raw < min_branch_supported_raw_score) or (np.isfinite(best_raw_score) and best_raw_score < min_leaf_raw_score):
            status = "stopped_at_parent"
            stop_reason = "below_absolute_support"
            ambiguity_flag = True
            break
        if not np.isnan(best_margin) and best_margin < margin_threshold:
            status = "stopped_at_parent"
            stop_reason = "ambiguous_sibling_margin"
            ambiguity_flag = True
            break

        path.append(best_label)
        chosen = next(p for p in current_programs if p.name == best_label)
        if not chosen.children:
            status = "assigned"
            stop_reason = "leaf_reached"
            break
        current_programs = chosen.children
        level += 1

    cluster_df = all_scores[all_scores["cluster"] == cluster].dropna(subset=[score_col])
    if cluster_df.empty:
        best_score_any_level = np.nan
        best_label_any_level = "Unresolved"
    else:
        idx = cluster_df[score_col].idxmax()
        best_score_any_level = float(cluster_df.loc[idx, score_col])
        best_label_any_level = str(cluster_df.loc[idx, "label"])

    decision_label = path[-1] if path else "Unresolved"
    decision_path = " > ".join(path) if path else "Unresolved"
    decision_level = len(path)
    decision_score = level_scores.get(f"level_{decision_level}", np.nan) if decision_level > 0 else np.nan
    decision_raw_score = level_raw_scores.get(f"level_{decision_level}", np.nan) if decision_level > 0 else np.nan
    decision_branch_supported_raw_score = level_branch_supported_raw_scores.get(f"level_{decision_level}", np.nan) if decision_level > 0 else np.nan
    confidence = level_confidence.get(f"level_{decision_level}", "none") if decision_level > 0 else "none"

    base = {
        "cluster": cluster,
        "decision_label": decision_label,
        "decision_path": decision_path,
        "decision_level": decision_level,
        "decision_score": decision_score,
        "decision_raw_score": decision_raw_score,
        "decision_branch_supported_raw_score": decision_branch_supported_raw_score,
        "decision_status": status,
        "decision_stop_reason": stop_reason,
        "confidence": confidence,
        "best_score_any_level": best_score_any_level,
        "best_label_any_level": best_label_any_level,
        "level_labels": level_labels,
        "level_scores": level_scores,
        "level_raw_scores": level_raw_scores,
        "level_branch_supported_raw_scores": level_branch_supported_raw_scores,
        "level_margins": level_margins,
        "level_confidence": level_confidence,
        "status": status,
        "stop_reason": stop_reason,
        "ambiguity_flag": ambiguity_flag,
    }
    final = _apply_absolute_support_backoff_row(
        base,
        min_branch_supported_raw_score=min_branch_supported_raw_score,
        min_leaf_raw_score=min_leaf_raw_score,
    )

    return ClusterAnnotation(
        cluster=cluster,
        final_label=final["final_label"],
        final_path=final["final_path"],
        final_level=final["final_level"],
        status=final["status"],
        stop_reason=final["stop_reason"],
        ambiguity_flag=bool(ambiguity_flag),
        confidence=final["confidence"],
        best_score_any_level=best_score_any_level,
        best_label_any_level=best_label_any_level,
        final_score=final["final_score"],
        final_raw_score=final["final_raw_score"],
        final_branch_supported_raw_score=final["final_branch_supported_raw_score"],
        final_margin=final["final_margin"],
        decision_label=decision_label,
        decision_path=decision_path,
        decision_level=decision_level,
        decision_score=decision_score,
        decision_raw_score=decision_raw_score,
        decision_branch_supported_raw_score=decision_branch_supported_raw_score,
        backoff_applied=bool(final["backoff_applied"]),
        backoff_steps=int(final["backoff_steps"]),
        backoff_reason=str(final["backoff_reason"]),
        level_labels=level_labels,
        level_scores=level_scores,
        level_raw_scores=level_raw_scores,
        level_branch_supported_raw_scores=level_branch_supported_raw_scores,
        level_margins=level_margins,
        level_confidence=level_confidence,
    )
