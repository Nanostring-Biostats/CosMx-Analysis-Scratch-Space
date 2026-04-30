from __future__ import annotations

from typing import Dict, Iterable, List, Optional, Set

import numpy as np
import pandas as pd

from .datamodels import CompiledMarkerProgram, MarkerProgram, ScoreRecord



def collect_all_markers(programs: Iterable[MarkerProgram]) -> Set[str]:
    markers: Set[str] = set()
    for program in programs:
        markers.update(g.upper() for g in program.positive_markers)
        markers.update(g.upper() for g in program.negative_markers)
        if program.children:
            markers.update(collect_all_markers(program.children))
    return markers



def flatten_compiled_programs(programs: Iterable[CompiledMarkerProgram]) -> List[CompiledMarkerProgram]:
    flat: List[CompiledMarkerProgram] = []
    for program in programs:
        flat.extend(list(program.iter_nodes()))
    return flat



def prepare_program_markers(
    program: CompiledMarkerProgram,
    gene_index: pd.Index,
    min_markers_present: int,
    min_marker_fraction: float,
    fallback_to_canonical_if_sparse: bool,
) -> Dict[str, object]:
    gene_set = set(gene_index.astype(str))
    canonical = list(dict.fromkeys(g.upper() for g in program.canonical_markers))
    effective = list(dict.fromkeys(g.upper() for g in program.effective_markers))
    negative = list(dict.fromkeys(g.upper() for g in program.negative_markers))

    canonical_present = [g for g in canonical if g in gene_set]
    canonical_missing = [g for g in canonical if g not in gene_set]
    effective_present = [g for g in effective if g in gene_set]
    effective_missing = [g for g in effective if g not in gene_set]
    negative_present = [g for g in negative if g in gene_set]
    negative_missing = [g for g in negative if g not in gene_set]

    required = program.min_markers_present if program.min_markers_present is not None else min_markers_present
    effective_fraction = len(effective_present) / max(1, len(effective))

    use_fallback = False
    marker_source = "effective"
    scoring_markers = effective_present
    missing_markers = effective_missing
    markers_total = len(effective)

    sparse_effective = len(effective_present) < required or effective_fraction < min_marker_fraction
    if fallback_to_canonical_if_sparse and sparse_effective and len(canonical_present) > 0:
        use_fallback = True
        marker_source = "canonical_fallback"
        scoring_markers = canonical_present
        missing_markers = canonical_missing
        markers_total = len(canonical)

    present_fraction = len(scoring_markers) / max(1, markers_total)
    negative_fraction = len(negative_present) / max(1, len(negative))

    if len(scoring_markers) == 0:
        status = "no_markers_present"
    elif len(scoring_markers) < required:
        status = "insufficient_markers_present"
    elif present_fraction < min_marker_fraction:
        status = "low_marker_fraction"
    else:
        status = "ok"

    return {
        "scoring_markers": scoring_markers,
        "missing_markers": missing_markers,
        "markers_total": markers_total,
        "markers_present": len(scoring_markers),
        "markers_present_fraction": present_fraction,
        "effective_markers_total": len(effective),
        "effective_markers_present": len(effective_present),
        "canonical_markers_total": len(canonical),
        "marker_source": marker_source,
        "used_fallback_marker_set": use_fallback,
        "score_status": status,
        "negative_markers": negative_present,
        "missing_negative_markers": negative_missing,
        "negative_markers_total": len(negative),
        "negative_markers_present": len(negative_present),
        "negative_markers_present_fraction": negative_fraction,
    }



def _clip_deltas(values: List[float], quantile: Optional[float]) -> List[float]:
    if not values or quantile is None:
        return values
    q = float(quantile)
    if not (0.0 < q < 1.0) or len(values) < 3:
        return values
    arr = np.asarray(values, dtype=float)
    limit = float(np.quantile(np.abs(arr), q))
    if not np.isfinite(limit) or limit <= 0:
        return values
    return np.clip(arr, -limit, limit).tolist()



def score_program_across_clusters(
    expr: pd.DataFrame,
    program: CompiledMarkerProgram,
    control_maps: Dict[str, Dict[str, object]],
    min_markers_present: int = 3,
    min_marker_fraction: float = 0.3,
    fallback_to_canonical_if_sparse: bool = True,
    negative_weight: float = 0.5,
    detection_floor: float = 0.0,
    delta_clip_quantile: Optional[float] = None,
) -> List[ScoreRecord]:
    marker_info = prepare_program_markers(
        program=program,
        gene_index=expr.index,
        min_markers_present=min_markers_present,
        min_marker_fraction=min_marker_fraction,
        fallback_to_canonical_if_sparse=fallback_to_canonical_if_sparse,
    )
    scoring_markers = list(marker_info["scoring_markers"])
    negative_markers = list(marker_info["negative_markers"])

    positive_control_map = control_maps.get("positive", {}).get("controls", {})
    negative_control_map = control_maps.get("negative", {}).get("controls", {})
    positive_meta = control_maps.get("positive", {}).get("metadata", {})
    negative_meta = control_maps.get("negative", {}).get("metadata", {})

    records: List[ScoreRecord] = []
    for cluster in expr.columns:
        positive_deltas: List[float] = []
        negative_deltas: List[float] = []
        positive_marker_values: List[float] = []
        positive_control_values: List[float] = []
        negative_marker_values: List[float] = []
        negative_control_values: List[float] = []

        markers_above_detection_floor = 0
        negative_markers_above_detection_floor = 0

        for gene in scoring_markers:
            marker_val = float(expr.at[gene, cluster])
            if marker_val >= detection_floor:
                markers_above_detection_floor += 1
            controls = positive_control_map.get(gene, [])
            if not controls:
                continue
            control_val = float(expr.loc[controls, cluster].mean())
            positive_deltas.append(marker_val - control_val)
            positive_marker_values.append(marker_val)
            positive_control_values.append(control_val)

        for gene in negative_markers:
            marker_val = float(expr.at[gene, cluster])
            if marker_val >= detection_floor:
                negative_markers_above_detection_floor += 1
            controls = negative_control_map.get(gene, [])
            if not controls:
                continue
            control_val = float(expr.loc[controls, cluster].mean())
            negative_deltas.append(marker_val - control_val)
            negative_marker_values.append(marker_val)
            negative_control_values.append(control_val)

        positive_deltas = _clip_deltas(positive_deltas, delta_clip_quantile)
        negative_deltas = _clip_deltas(negative_deltas, delta_clip_quantile)

        if len(scoring_markers) == 0:
            score = np.nan
            positive_score = np.nan
            negative_score = np.nan
            score_status = "no_markers_present"
        elif marker_info["score_status"] != "ok":
            positive_score = float(np.mean(positive_deltas)) if positive_deltas else np.nan
            negative_score = float(np.mean(negative_deltas)) if negative_deltas else 0.0
            score = np.nan
            score_status = str(marker_info["score_status"])
        elif len(positive_deltas) == 0:
            positive_score = np.nan
            negative_score = float(np.mean(negative_deltas)) if negative_deltas else 0.0
            score = np.nan
            score_status = "no_controls_available"
        else:
            positive_score = float(np.mean(positive_deltas))
            negative_score = float(np.mean(negative_deltas)) if negative_deltas else 0.0
            score = positive_score - (negative_weight * negative_score)
            score_status = "ok"
            if positive_meta.get("fallback_used") or negative_meta.get("fallback_used"):
                score_status = "control_count_fallback_used"

        records.append(
            ScoreRecord(
                cluster=str(cluster),
                label=program.name,
                level=program.level,
                parent_label=program.parent_name,
                score=float(score) if not np.isnan(score) else np.nan,
                positive_score=float(positive_score) if not np.isnan(positive_score) else np.nan,
                negative_score=float(negative_score) if not np.isnan(negative_score) else np.nan,
                positive_marker_mean=float(np.mean(positive_marker_values)) if positive_marker_values else np.nan,
                positive_control_mean=float(np.mean(positive_control_values)) if positive_control_values else np.nan,
                negative_marker_mean=float(np.mean(negative_marker_values)) if negative_marker_values else np.nan,
                negative_control_mean=float(np.mean(negative_control_values)) if negative_control_values else np.nan,
                negative_weight=float(negative_weight),
                markers_present=int(marker_info["markers_present"]),
                markers_total=int(marker_info["markers_total"]),
                markers_present_fraction=float(marker_info["markers_present_fraction"]),
                markers_above_detection_floor=int(markers_above_detection_floor),
                markers_detection_support_fraction=float(markers_above_detection_floor / max(1, len(scoring_markers))),
                negative_markers_present=int(marker_info["negative_markers_present"]),
                negative_markers_total=int(marker_info["negative_markers_total"]),
                negative_markers_present_fraction=float(marker_info["negative_markers_present_fraction"]),
                negative_markers_above_detection_floor=int(negative_markers_above_detection_floor),
                negative_markers_detection_support_fraction=float(negative_markers_above_detection_floor / max(1, len(negative_markers))) if negative_markers else 0.0,
                missing_markers=list(marker_info["missing_markers"]),
                missing_negative_markers=list(marker_info["missing_negative_markers"]),
                marker_source=str(marker_info["marker_source"]),
                score_status=score_status,
                effective_markers_total=int(marker_info["effective_markers_total"]),
                effective_markers_present=int(marker_info["effective_markers_present"]),
                canonical_markers_total=int(marker_info["canonical_markers_total"]),
                used_fallback_marker_set=bool(marker_info["used_fallback_marker_set"]),
                positive_controls_total=int(sum(len(positive_control_map.get(g, [])) for g in scoring_markers)),
                negative_controls_total=int(sum(len(negative_control_map.get(g, [])) for g in negative_markers)),
                control_fallback_used=bool(positive_meta.get("fallback_used") or negative_meta.get("fallback_used")),
                median_controls_per_positive_marker=float(positive_meta.get("median_controls_per_marker", 0.0)),
                median_controls_per_negative_marker=float(negative_meta.get("median_controls_per_marker", 0.0)),
            )
        )
    return records



def score_programs_across_clusters(
    expr: pd.DataFrame,
    programs: Iterable[CompiledMarkerProgram],
    control_maps: Dict[str, Dict[str, Dict[str, object]]],
    min_markers_present: int = 3,
    min_marker_fraction: float = 0.3,
    fallback_to_canonical_if_sparse: bool = True,
    negative_weight: float = 0.5,
    detection_floor: float = 0.0,
    delta_clip_quantile: Optional[float] = None,
) -> pd.DataFrame:
    rows = []
    for program in programs:
        records = score_program_across_clusters(
            expr=expr,
            program=program,
            control_maps=control_maps[program.name],
            min_markers_present=min_markers_present,
            min_marker_fraction=min_marker_fraction,
            fallback_to_canonical_if_sparse=fallback_to_canonical_if_sparse,
            negative_weight=negative_weight,
            detection_floor=detection_floor,
            delta_clip_quantile=delta_clip_quantile,
        )
        rows.extend(record.__dict__ for record in records)
    return pd.DataFrame(rows)
