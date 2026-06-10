from __future__ import annotations

from typing import Any, Dict, Iterable, List, Optional, Sequence

import numpy as np
import pandas as pd

from .constants import (
    DEFAULT_STATE_SUPPORT_BONUS_PER_GROUP,
    DEFAULT_STATE_SUPPORT_MAX_BONUS,
    DEFAULT_STATUS_PROGRAM_BONUS_PER_EXTRA,
    DEFAULT_STATUS_PROGRAM_MAX_BONUS,
    DEFAULT_WEAK_STATUS_FRACTION_FOR_STATE_SUPPORT,
    PROGRAM_LABEL_SEPARATOR,
    PROGRAM_LIST_SEPARATOR,
    REPORTING_ROLE_MODIFIER,
    REPORTING_ROLE_STATE,
    REPORTING_ROLE_STATUS,
    REPORTING_ROLES,
)
from .datamodels import CompiledMarkerProgram, MalignantProgram
from .scoring import score_programs_across_clusters


_MALIGNANT_SCORE_COLUMNS = [
    "cluster",
    "label",
    "name",
    "competition_group",
    "reporting_label",
    "reporting_role",
    "report_on_lineages",
    "raw_score",
    "status_score",
    "absolute_support_score",
    "program_robust_zscore",
    "program_specificity_score",
    "specificity_score",
    "decision_score",
    "score",
    "raw_support_gate",
    "support_scale",
    "positive_score",
    "negative_score",
    "positive_marker_mean",
    "positive_control_mean",
    "negative_marker_mean",
    "negative_control_mean",
    "negative_weight",
    "markers_present",
    "markers_total",
    "markers_present_fraction",
    "markers_above_detection_floor",
    "markers_detection_support_fraction",
    "negative_markers_present",
    "negative_markers_total",
    "negative_markers_present_fraction",
    "negative_markers_above_detection_floor",
    "negative_markers_detection_support_fraction",
    "marker_source",
    "score_status",
]

_MALIGNANT_ANNOT_COLUMNS = [
    "cluster_id",
    "annot_malignant_name",
    "annot_malignant_label",
    "annot_malignant_label_concise",
    "annot_malignant_status_program",
    "annot_malignant_status_label",
    "annot_malignant_state_label",
    "annot_malignant_modifier_labels",
    "annot_malignant_status",
    "annot_malignant_tumor_status_pass",
    "annot_malignant_reason",
    "annot_malignant_score",
    "annot_malignant_decision_score",
    "annot_malignant_status_score",
    "annot_malignant_raw_score",
    "annot_malignant_core_status_score",
    "annot_malignant_core_status_raw_score",
    "annot_malignant_combined_status_score",
    "annot_malignant_combined_status_raw_score",
    "annot_malignant_state_support_score",
    "annot_malignant_state_support_programs",
    "annot_malignant_status_program_count",
    "annot_malignant_status_group_count",
    "annot_malignant_status_decision_source",
    "annot_malignant_label_raw_score",
    "annot_malignant_specificity_score",
    "annot_malignant_program_robust_zscore",
    "annot_malignant_margin",
    "annot_malignant_second_name",
    "annot_malignant_second_score",
    "annot_malignant_competition_group",
    "annot_malignant_reporting_role",
    "annot_malignant_program_report_on_lineages",
    "annot_malignant_specific",
    "annot_malignant_mixed",
    "annot_malignant_positive_programs",
    "annot_malignant_strong_programs",
    "annot_malignant_status_programs",
    "annot_malignant_state_programs",
    "annot_malignant_modifier_programs",
    "annot_malignant_positive_groups",
    "annot_malignant_multigroup_positive",
    "annot_malignant_markers_present",
    "annot_malignant_markers_total",
    "annot_malignant_markers_detection_support_fraction",
]


def _normalize_malignant_programs(malignant_programs):
    if malignant_programs is None:
        return []
    out = []
    for prog in malignant_programs:
        if isinstance(prog, MalignantProgram):
            out.append(prog)
        elif isinstance(prog, dict):
            out.append(MalignantProgram.from_dict(prog))
        else:
            raise TypeError("malignant_programs must contain MalignantProgram or dict items")
    return out


def _dedupe_upper(values: Iterable[str]) -> List[str]:
    return list(dict.fromkeys(str(g).strip().upper() for g in values if str(g).strip()))


def _program_metadata(program: MalignantProgram) -> Dict[str, Any]:
    return dict(getattr(program, "metadata", {}) or {})


def _get_competition_group(program: MalignantProgram) -> str:
    meta = _program_metadata(program)
    value = meta.get("competition_group", meta.get("malignant_competition_group", None))
    if value is None:
        return str(program.name)
    text = str(value).strip()
    if text == "" or text.lower() in {"none", "off", "independent", "false"}:
        return str(program.name)
    return text


def _get_reporting_label(program: MalignantProgram) -> str:
    meta = _program_metadata(program)
    value = meta.get("reporting_label", meta.get("label", None))
    if value is None or str(value).strip() == "":
        return str(program.name)
    return str(value).strip()


def _normalize_reporting_role(value) -> str:
    """Normalize program metadata roles used by the tiered reporting layer."""
    text = "" if value is None else str(value).strip().lower()
    text = text.replace("-", "_").replace(" ", "_")
    if text in {"tumor_status", "driver", "identity", "status_program"}:
        return REPORTING_ROLE_STATUS
    if text in {"state", "tumor_state", "substate"}:
        return REPORTING_ROLE_STATE
    if text in {"modifier", "flag", "state_flag", "auxiliary", "aux", "program_flag"}:
        return REPORTING_ROLE_MODIFIER
    return text


def _get_reporting_role(program: MalignantProgram) -> str:
    meta = _program_metadata(program)
    value = meta.get("reporting_role", "")
    return _normalize_reporting_role(value)


def _get_report_on_lineages(program: MalignantProgram) -> str:
    meta = _program_metadata(program)
    value = meta.get("report_on_lineages", meta.get("allowed_lineages", None))
    if value is None:
        return ""
    if isinstance(value, str):
        return value.strip()
    if isinstance(value, (list, tuple, set)):
        return ";".join(str(x).strip() for x in value if str(x).strip())
    return str(value).strip()


def _to_compiled_programs(programs: Sequence[MalignantProgram]) -> List[CompiledMarkerProgram]:
    compiled: List[CompiledMarkerProgram] = []
    for prog in programs:
        pos = _dedupe_upper(prog.positive_markers)
        neg = _dedupe_upper(prog.negative_markers)
        meta = _program_metadata(prog)
        meta.setdefault("competition_group", _get_competition_group(prog))
        meta.setdefault("reporting_label", _get_reporting_label(prog))
        meta.setdefault("reporting_role", _get_reporting_role(prog))
        meta.setdefault("report_on_lineages", _get_report_on_lineages(prog))
        compiled.append(
            CompiledMarkerProgram(
                name=str(prog.name),
                canonical_markers=pos,
                effective_markers=pos,
                negative_markers=neg,
                children=[],
                description=prog.description,
                aliases=[],
                min_markers_present=None,
                level=1,
                parent_name=_get_competition_group(prog),
                cleanup_strategy="none",
                metadata=meta,
            )
        )
    return compiled


def _robust_z(values: pd.Series) -> pd.Series:
    vals = pd.to_numeric(values, errors="coerce").astype(float)
    out = pd.Series(np.nan, index=values.index, dtype=float)
    arr = vals.to_numpy(dtype=float)
    finite = np.isfinite(arr)
    if finite.sum() == 0:
        return out
    if finite.sum() == 1:
        out.iloc[np.where(finite)[0]] = 0.0
        return out
    finite_vals = arr[finite]
    med = float(np.median(finite_vals))
    mad = float(np.median(np.abs(finite_vals - med)))
    scale = 1.4826 * mad if mad > 0 else float(np.std(finite_vals))
    if not np.isfinite(scale) or scale <= 0:
        scale = 1.0
    out.iloc[np.where(finite)[0]] = (finite_vals - med) / scale
    return out


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


def _empty_malignant_scores() -> pd.DataFrame:
    return pd.DataFrame(columns=_MALIGNANT_SCORE_COLUMNS)


def _score_without_controls(
    expr: pd.DataFrame,
    programs: Sequence[MalignantProgram],
    *,
    min_markers_present: int = 1,
    min_marker_fraction: float = 0.0,
    negative_weight: float = 0.5,
    detection_floor: float = 0.0,
) -> pd.DataFrame:
    records: List[Dict[str, Any]] = []
    gene_index = set(expr.index.astype(str))
    for prog in programs:
        name = str(prog.name)
        pos = _dedupe_upper(prog.positive_markers)
        neg = _dedupe_upper(prog.negative_markers)
        pos_present = [g for g in pos if g in gene_index]
        neg_present = [g for g in neg if g in gene_index]
        markers_total = len(pos)
        markers_present = len(pos_present)
        markers_present_fraction = markers_present / max(1, markers_total)
        if len(pos_present) == 0:
            score_status = "no_markers_present"
        elif markers_present < int(min_markers_present):
            score_status = "insufficient_markers_present"
        elif markers_present_fraction < float(min_marker_fraction):
            score_status = "low_marker_fraction"
        else:
            score_status = "ok"

        for cluster in expr.columns:
            if len(pos_present) == 0:
                positive_score = np.nan
                negative_score = np.nan
                raw_score = np.nan
            else:
                positive_score = float(expr.loc[pos_present, cluster].mean())
                negative_score = float(expr.loc[neg_present, cluster].mean()) if neg_present else 0.0
                raw_score = positive_score - float(negative_weight) * negative_score
                if score_status != "ok":
                    raw_score = np.nan

            markers_above = int(sum(float(expr.at[g, cluster]) >= float(detection_floor) for g in pos_present))
            neg_above = int(sum(float(expr.at[g, cluster]) >= float(detection_floor) for g in neg_present))
            records.append(
                {
                    "cluster": str(cluster),
                    "label": name,
                    "name": name,
                    "raw_score": raw_score,
                    "positive_score": positive_score,
                    "negative_score": negative_score,
                    "positive_marker_mean": positive_score,
                    "positive_control_mean": np.nan,
                    "negative_marker_mean": float(expr.loc[neg_present, cluster].mean()) if neg_present else np.nan,
                    "negative_control_mean": np.nan,
                    "negative_weight": float(negative_weight),
                    "markers_present": int(markers_present),
                    "markers_total": int(markers_total),
                    "markers_present_fraction": float(markers_present_fraction),
                    "markers_above_detection_floor": markers_above,
                    "markers_detection_support_fraction": float(markers_above / max(1, len(pos_present))),
                    "negative_markers_present": int(len(neg_present)),
                    "negative_markers_total": int(len(neg)),
                    "negative_markers_present_fraction": float(len(neg_present) / max(1, len(neg))),
                    "negative_markers_above_detection_floor": neg_above,
                    "negative_markers_detection_support_fraction": float(neg_above / max(1, len(neg_present))) if neg_present else 0.0,
                    "marker_source": "malignant_program",
                    "score_status": score_status,
                }
            )
    return pd.DataFrame(records)


def _add_program_metadata(df: pd.DataFrame, programs: Sequence[MalignantProgram]) -> pd.DataFrame:
    out = df.copy()
    program_by_name = {str(p.name): p for p in programs}
    out["name"] = out.get("name", out.get("label", "")).astype(str)
    out["label"] = out.get("label", out["name"]).astype(str)
    out["competition_group"] = out["name"].map(lambda n: _get_competition_group(program_by_name[n]) if n in program_by_name else n)
    out["reporting_label"] = out["name"].map(lambda n: _get_reporting_label(program_by_name[n]) if n in program_by_name else n)
    out["reporting_role"] = out["name"].map(lambda n: _get_reporting_role(program_by_name[n]) if n in program_by_name else "")
    out["report_on_lineages"] = out["name"].map(lambda n: _get_report_on_lineages(program_by_name[n]) if n in program_by_name else "")
    return out


def _add_malignant_decision_metrics(
    scores: pd.DataFrame,
    *,
    raw_weight: float = 0.60,
    weak_raw_score_threshold: float = 0.10,
) -> pd.DataFrame:
    df = scores.copy()
    if df.empty:
        return _empty_malignant_scores()

    df["raw_score"] = pd.to_numeric(df["raw_score"], errors="coerce")
    if "markers_detection_support_fraction" not in df.columns:
        df["markers_detection_support_fraction"] = 1.0
    support = pd.to_numeric(df["markers_detection_support_fraction"], errors="coerce").fillna(0.0).clip(lower=0.0, upper=1.0)
    df["support_scale"] = 0.5 + 0.5 * support

    df["program_robust_zscore"] = df.groupby("name", sort=False)["raw_score"].transform(_robust_z)
    if "competition_group" not in df.columns:
        df["competition_group"] = df["name"].astype(str)
    if "reporting_role" not in df.columns:
        df["reporting_role"] = ""
    if "report_on_lineages" not in df.columns:
        df["report_on_lineages"] = ""
    df["competition_group"] = df["competition_group"].fillna(df["name"]).astype(str)

    df["program_specificity_score"] = np.nan
    for _, idx in df.groupby(["cluster", "competition_group"], dropna=False, sort=False).groups.items():
        sub = pd.to_numeric(df.loc[idx, "program_robust_zscore"], errors="coerce")
        arr = sub.to_numpy(dtype=float)
        finite_mask = np.isfinite(arr)
        if finite_mask.sum() == 0:
            out = np.full(len(arr), np.nan, dtype=float)
        elif len(arr) == 1:
            out = arr.copy()
        else:
            out = np.full(len(arr), np.nan, dtype=float)
            for i in range(len(arr)):
                if not np.isfinite(arr[i]):
                    continue
                others = np.delete(arr, i)
                finite_others = others[np.isfinite(others)]
                out[i] = arr[i] if finite_others.size == 0 else arr[i] - float(np.max(finite_others))
        df.loc[idx, "program_specificity_score"] = out

    df["specificity_score"] = df["program_specificity_score"]
    df["raw_support_gate"] = df["raw_score"].map(lambda x: _raw_support_gate(float(x), weak_high=float(weak_raw_score_threshold)))

    # Status score intentionally remains an absolute, support-weighted tumor-like
    # evidence metric. It should be used for strong/weak malignant status calls.
    df["status_score"] = df["raw_score"] * df["support_scale"]
    df["absolute_support_score"] = df["status_score"]

    rw = float(raw_weight)
    rw = min(1.0, max(0.0, rw))
    raw_part = rw * df["raw_score"]
    specificity_part = (1.0 - rw) * df["raw_support_gate"] * pd.to_numeric(df["program_specificity_score"], errors="coerce")
    df["decision_score"] = (raw_part.fillna(0.0) + specificity_part.fillna(0.0)) * df["support_scale"]
    both_missing = df["raw_score"].isna() & pd.to_numeric(df["program_specificity_score"], errors="coerce").isna()
    df.loc[both_missing, "decision_score"] = np.nan

    # Backward-compatible alias: malignant_scores["score"] is now the label/state
    # decision score, while raw_score/status_score hold absolute evidence.
    df["score"] = df["decision_score"]

    for col in _MALIGNANT_SCORE_COLUMNS:
        if col not in df.columns:
            df[col] = np.nan
    front = [c for c in _MALIGNANT_SCORE_COLUMNS if c in df.columns]
    rest = [c for c in df.columns if c not in front]
    return df[front + rest]


def _score_malignant_programs(
    cluster_means: pd.DataFrame,
    malignant_programs,
    *,
    control_maps: Optional[Dict[str, Dict[str, Dict[str, object]]]] = None,
    min_markers_present: int = 1,
    min_marker_fraction: float = 0.0,
    fallback_to_canonical_if_sparse: bool = True,
    negative_weight: float = 0.5,
    detection_floor: float = 0.0,
    delta_clip_quantile: Optional[float] = None,
    raw_weight: float = 0.60,
    weak_raw_score_threshold: float = 0.10,
) -> pd.DataFrame:
    """Score malignant/transformation programs across clusters.

    ``raw_score`` is the absolute enrichment metric. When ``control_maps`` is
    supplied, it is computed with expression-matched control genes using the same
    scoring helper as the normal hierarchy. ``program_specificity_score`` and
    ``decision_score`` add within-competition-group specificity and should be
    used to describe the malignant state, not to decide whether a cluster is
    tumor-like.
    """
    if malignant_programs is None:
        return _empty_malignant_scores()
    if cluster_means is None or cluster_means.empty:
        return _empty_malignant_scores()

    programs = _normalize_malignant_programs(malignant_programs)
    if not programs:
        return _empty_malignant_scores()

    expr = cluster_means.copy()
    expr.index = expr.index.astype(str).str.upper()

    if control_maps is not None:
        compiled = _to_compiled_programs(programs)
        base = score_programs_across_clusters(
            expr=expr,
            programs=compiled,
            control_maps=control_maps,
            min_markers_present=int(min_markers_present),
            min_marker_fraction=float(min_marker_fraction),
            fallback_to_canonical_if_sparse=bool(fallback_to_canonical_if_sparse),
            negative_weight=float(negative_weight),
            detection_floor=float(detection_floor),
            delta_clip_quantile=delta_clip_quantile,
        )
        if base.empty:
            return _empty_malignant_scores()
        base = base.rename(columns={"score": "raw_score"})
        base["name"] = base["label"].astype(str)
    else:
        base = _score_without_controls(
            expr,
            programs,
            min_markers_present=int(min_markers_present),
            min_marker_fraction=float(min_marker_fraction),
            negative_weight=float(negative_weight),
            detection_floor=float(detection_floor),
        )

    base = _add_program_metadata(base, programs)
    return _add_malignant_decision_metrics(
        base,
        raw_weight=float(raw_weight),
        weak_raw_score_threshold=float(weak_raw_score_threshold),
    )


def _safe_float(value, default=np.nan) -> float:
    try:
        out = float(value)
    except Exception:
        return float(default)
    return out if np.isfinite(out) else float(default)


def _format_program_list(names: Iterable[str]) -> str:
    cleaned = []
    for name in names:
        text = str(name).strip()
        if text and text.lower() != "nan" and text not in cleaned:
            cleaned.append(text)
    return PROGRAM_LIST_SEPARATOR.join(cleaned)


def _format_mixed_label(names: Sequence[str]) -> str:
    cleaned = [str(x).strip() for x in names if str(x).strip() and str(x).strip().lower() != "nan"]
    if not cleaned:
        return "malignant"
    return "mixed_" + PROGRAM_LABEL_SEPARATOR.join(cleaned[:2])


def _empty_annotation_record(cluster) -> Dict[str, Any]:
    return {
        "cluster_id": str(cluster),
        "annot_malignant_name": np.nan,
        "annot_malignant_label": np.nan,
        "annot_malignant_label_concise": np.nan,
        "annot_malignant_status_program": np.nan,
        "annot_malignant_status_label": "",
        "annot_malignant_state_label": "",
        "annot_malignant_modifier_labels": "",
        "annot_malignant_status": "none",
        "annot_malignant_tumor_status_pass": False,
        "annot_malignant_reason": "malignant_support_none",
        "annot_malignant_score": np.nan,
        "annot_malignant_decision_score": np.nan,
        "annot_malignant_status_score": np.nan,
        "annot_malignant_raw_score": np.nan,
        "annot_malignant_core_status_score": np.nan,
        "annot_malignant_core_status_raw_score": np.nan,
        "annot_malignant_combined_status_score": np.nan,
        "annot_malignant_combined_status_raw_score": np.nan,
        "annot_malignant_state_support_score": 0.0,
        "annot_malignant_state_support_programs": "",
        "annot_malignant_status_program_count": 0,
        "annot_malignant_status_group_count": 0,
        "annot_malignant_status_decision_source": "none",
        "annot_malignant_label_raw_score": np.nan,
        "annot_malignant_specificity_score": np.nan,
        "annot_malignant_program_robust_zscore": np.nan,
        "annot_malignant_margin": np.nan,
        "annot_malignant_second_name": np.nan,
        "annot_malignant_second_score": np.nan,
        "annot_malignant_competition_group": np.nan,
        "annot_malignant_reporting_role": np.nan,
        "annot_malignant_program_report_on_lineages": "",
        "annot_malignant_specific": False,
        "annot_malignant_mixed": False,
        "annot_malignant_positive_programs": "",
        "annot_malignant_strong_programs": "",
        "annot_malignant_status_programs": "",
        "annot_malignant_state_programs": "",
        "annot_malignant_modifier_programs": "",
        "annot_malignant_positive_groups": "",
        "annot_malignant_multigroup_positive": False,
        "annot_malignant_markers_present": np.nan,
        "annot_malignant_markers_total": np.nan,
        "annot_malignant_markers_detection_support_fraction": np.nan,
    }


def _row_label(row) -> str:
    name = str(row.get("name", "")).strip()
    label = str(row.get("reporting_label", name)).strip()
    return label if label and label.lower() != "nan" else name


def _sort_programs(df: pd.DataFrame, keys: Sequence[str]) -> pd.DataFrame:
    cols = [c for c in keys if c in df.columns]
    if not cols:
        cols = ["status_score", "raw_score"]
    return df.sort_values(cols, ascending=[False] * len(cols), na_position="last")


def _select_state_label(
    state_candidates: pd.DataFrame,
    *,
    margin_threshold: float,
    specificity_threshold: float,
) -> tuple[str, Optional[pd.Series], object, float, bool, bool]:
    """Return label, row, second_name, margin, specific, mixed for state tier."""
    if state_candidates is None or state_candidates.empty:
        return "", None, np.nan, np.nan, False, False
    ordered = _sort_programs(state_candidates, ["decision_score", "status_score", "raw_score"])
    row = ordered.iloc[0]
    label = _row_label(row)
    best_decision = _safe_float(row.get("decision_score", np.nan))
    best_spec = _safe_float(row.get("program_specificity_score", np.nan))
    best_group = str(row.get("competition_group", ""))
    same_group = ordered[ordered["competition_group"].astype(str) == best_group]
    second_name = np.nan
    second_score = np.nan
    margin = np.nan
    if len(same_group) > 1:
        second = same_group.iloc[1]
        second_name = str(second.get("name", ""))
        second_score = _safe_float(second.get("decision_score", np.nan))
        if np.isfinite(best_decision) and np.isfinite(second_score):
            margin = float(best_decision - second_score)
    margin_ok = (not np.isfinite(margin)) or margin >= float(margin_threshold)
    specificity_ok = (not np.isfinite(best_spec)) or best_spec >= float(specificity_threshold)
    mixed = False
    if len(same_group) > 1 and np.isfinite(margin) and margin < float(margin_threshold):
        mixed = True
        labels = [_row_label(row), _row_label(same_group.iloc[1])]
        label = _format_mixed_label(labels)
    specific = bool(margin_ok and specificity_ok and not mixed)
    return label, row, second_name, margin, specific, mixed


def _format_malignant_annotations(
    *,
    malignant_scores: pd.DataFrame,
    score_threshold: float = 0.35,
    raw_score_threshold: float = 0.15,
    margin_threshold: float = 0.15,
    specificity_threshold: float = 0.0,
) -> pd.DataFrame:
    """Summarize flat tumor/auxiliary program scores into one row per cluster.

    The program track is scored in parallel rather than routed as a hierarchy.
    Reporting is tiered by ``reporting_role`` metadata:

    - ``status`` programs establish tumor-like identity and control whether an
      integrated tumor-like label is allowed. Multiple status programs are
      aggregated with a max-plus-small-support-bonus rule.
    - ``state`` programs decorate the tumor-like label when they have enough
      raw and specificity support. They can add a small gated support bonus only
      when core status evidence is already borderline.
    - ``modifier`` programs are reported as auxiliary flags. They do not
      establish tumor-like identity or contribute to the default tumor-status
      decision.

    Cross-group co-activation is reported through positive-program summary
    columns instead of being treated as a normal-hierarchy-style mixed identity.
    Within-group state ambiguity can still produce a ``mixed_*`` state label.
    If a program set contains no status-role programs, the table behaves as a
    flag-only auxiliary summary: strong support is based on any strong program,
    but ``annot_malignant_tumor_status_pass`` remains ``False``.
    """
    if malignant_scores is None or malignant_scores.empty:
        return pd.DataFrame(columns=_MALIGNANT_ANNOT_COLUMNS)

    df = malignant_scores.copy()
    if "cluster" not in df.columns:
        return pd.DataFrame(columns=_MALIGNANT_ANNOT_COLUMNS)
    df["cluster"] = df["cluster"].astype(str)
    if "name" not in df.columns:
        df["name"] = df.get("label", "").astype(str)
    if "label" not in df.columns:
        df["label"] = df["name"].astype(str)
    if "raw_score" not in df.columns:
        df["raw_score"] = df.get("score", np.nan)
    if "decision_score" not in df.columns:
        df["decision_score"] = df.get("score", np.nan)
    if "status_score" not in df.columns:
        df["status_score"] = df.get("absolute_support_score", df.get("score", df["raw_score"]))
    if "specificity_score" not in df.columns:
        df["specificity_score"] = df.get("program_specificity_score", np.nan)
    if "program_specificity_score" not in df.columns:
        df["program_specificity_score"] = df["specificity_score"]
    if "reporting_label" not in df.columns:
        df["reporting_label"] = df["name"].astype(str)
    if "competition_group" not in df.columns:
        df["competition_group"] = df["name"].astype(str)
    if "reporting_role" not in df.columns:
        df["reporting_role"] = ""
    if "report_on_lineages" not in df.columns:
        df["report_on_lineages"] = ""

    df["reporting_role"] = df["reporting_role"].map(_normalize_reporting_role)
    df.loc[~df["reporting_role"].isin(REPORTING_ROLES), "reporting_role"] = REPORTING_ROLE_STATE

    numeric_cols = [
        "raw_score",
        "score",
        "decision_score",
        "status_score",
        "specificity_score",
        "program_specificity_score",
        "program_robust_zscore",
        "markers_present",
        "markers_total",
        "markers_detection_support_fraction",
    ]
    for col in numeric_cols:
        if col not in df.columns:
            df[col] = np.nan
        df[col] = pd.to_numeric(df[col], errors="coerce")

    records: List[Dict[str, Any]] = []
    for cluster, sub in df.groupby("cluster", sort=False):
        sub = sub.copy()
        finite_any = (
            np.isfinite(sub["raw_score"].to_numpy(dtype=float))
            | np.isfinite(sub["status_score"].to_numpy(dtype=float))
            | np.isfinite(sub["decision_score"].to_numpy(dtype=float))
        )
        if not finite_any.any():
            records.append(_empty_annotation_record(cluster))
            continue

        strong_mask = (
            sub["status_score"].ge(float(score_threshold)).fillna(False)
            & sub["raw_score"].ge(float(raw_score_threshold)).fillna(False)
        )
        weak_status_mask = (
            sub["status_score"].ge(float(score_threshold) * DEFAULT_WEAK_STATUS_FRACTION_FOR_STATE_SUPPORT).fillna(False)
            & sub["raw_score"].ge(float(raw_score_threshold) * DEFAULT_WEAK_STATUS_FRACTION_FOR_STATE_SUPPORT).fillna(False)
        )
        raw_positive_mask = sub["raw_score"].ge(float(raw_score_threshold)).fillna(False)
        has_status_programs = bool(sub["reporting_role"].eq(REPORTING_ROLE_STATUS).any())
        status_mask = sub["reporting_role"].eq(REPORTING_ROLE_STATUS)
        state_mask = sub["reporting_role"].eq(REPORTING_ROLE_STATE)
        modifier_mask = sub["reporting_role"].eq(REPORTING_ROLE_MODIFIER)
        status_strong_mask = strong_mask & status_mask
        state_strong_mask = strong_mask & state_mask
        modifier_strong_mask = strong_mask & modifier_mask
        status_support_mask = weak_status_mask & status_mask

        all_ordered = _sort_programs(sub, ["decision_score", "status_score", "raw_score"])
        aux_or_label_row = all_ordered.iloc[0]

        if has_status_programs:
            status_all_sorted = _sort_programs(sub.loc[status_mask], ["status_score", "raw_score", "decision_score"])
            status_row = status_all_sorted.iloc[0]
        else:
            status_all_sorted = _sort_programs(sub, ["status_score", "raw_score", "decision_score"])
            status_row = status_all_sorted.iloc[0]

        core_status_score = _safe_float(status_row.get("status_score", np.nan)) if has_status_programs else np.nan
        core_status_raw = _safe_float(status_row.get("raw_score", np.nan)) if has_status_programs else np.nan
        n_strong_status = int(status_strong_mask.sum()) if has_status_programs else 0
        if n_strong_status > 0:
            n_strong_status_groups = len(
                [x for x in _format_program_list(sub.loc[status_strong_mask, "competition_group"].astype(str).tolist()).split(";") if x]
            )
        else:
            n_strong_status_groups = 0
        if status_support_mask.any():
            n_support_status_groups = len(
                [x for x in _format_program_list(sub.loc[status_support_mask, "competition_group"].astype(str).tolist()).split(";") if x]
            )
        else:
            n_support_status_groups = 0
        status_bonus = min(
            DEFAULT_STATUS_PROGRAM_MAX_BONUS,
            DEFAULT_STATUS_PROGRAM_BONUS_PER_EXTRA * max(0, n_support_status_groups - 1),
        )
        core_status_score_aggregate = core_status_score + status_bonus if np.isfinite(core_status_score) else np.nan
        core_status_pass = bool(
            has_status_programs
            and np.isfinite(core_status_score_aggregate)
            and core_status_score_aggregate >= float(score_threshold)
            and np.isfinite(core_status_raw)
            and core_status_raw >= float(raw_score_threshold)
        )
        borderline_status = bool(
            has_status_programs
            and np.isfinite(core_status_score)
            and core_status_score >= float(score_threshold) * DEFAULT_WEAK_STATUS_FRACTION_FOR_STATE_SUPPORT
            and np.isfinite(core_status_raw)
            and core_status_raw >= float(raw_score_threshold) * DEFAULT_WEAK_STATUS_FRACTION_FOR_STATE_SUPPORT
        )

        status_candidates = sub.loc[status_strong_mask]
        state_candidates = sub.loc[state_strong_mask]
        modifier_candidates = sub.loc[modifier_strong_mask]

        state_label, state_row, second_name, margin, state_specific, state_mixed = _select_state_label(
            state_candidates,
            margin_threshold=float(margin_threshold),
            specificity_threshold=float(specificity_threshold),
        )

        state_support_programs = ""
        state_support_score = 0.0
        if borderline_status and state_candidates is not None and not state_candidates.empty:
            state_support_programs = _format_program_list(state_candidates["name"].astype(str).tolist())
            n_state_groups = len([x for x in _format_program_list(state_candidates["competition_group"].astype(str).tolist()).split(";") if x])
            state_support_score = min(
                DEFAULT_STATE_SUPPORT_MAX_BONUS,
                DEFAULT_STATE_SUPPORT_BONUS_PER_GROUP * max(1, n_state_groups),
            )
        combined_status_score = (
            core_status_score_aggregate + state_support_score
            if np.isfinite(core_status_score_aggregate)
            else np.nan
        )
        combined_status_raw = core_status_raw
        state_supported_pass = bool(
            has_status_programs
            and not core_status_pass
            and state_support_score > 0
            and np.isfinite(combined_status_score)
            and combined_status_score >= float(score_threshold)
            and np.isfinite(core_status_raw)
            and core_status_raw >= float(raw_score_threshold) * DEFAULT_WEAK_STATUS_FRACTION_FOR_STATE_SUPPORT
        )
        tumor_status_pass = bool(core_status_pass or state_supported_pass)
        any_strong = bool(strong_mask.any())
        track_strong = tumor_status_pass if has_status_programs else any_strong
        status = "strong" if track_strong else ("weak" if finite_any.any() else "none")
        status_decision_source = (
            "status_core"
            if core_status_pass
            else ("status_core_plus_state_support" if state_supported_pass else ("no_status_evidence" if has_status_programs else "no_status_programs"))
        )

        status_label = _row_label(status_row) if tumor_status_pass else ""
        status_program_name = str(status_row.get("name", "")) if tumor_status_pass else ""

        modifier_ordered = _sort_programs(modifier_candidates, ["status_score", "raw_score", "decision_score"])
        modifier_labels = [_row_label(r) for _, r in modifier_ordered.iterrows()]
        modifier_label_text = _format_program_list(modifier_labels)

        if tumor_status_pass:
            if state_label:
                concise = f"{status_label}{PROGRAM_LABEL_SEPARATOR}{state_label}"
                label_row = state_row if state_row is not None else status_row
                specific = bool(state_specific)
                mixed = bool(state_mixed)
                if state_supported_pass and not core_status_pass:
                    reason = "malignant_status_pass_state_supported"
                else:
                    reason = "malignant_status_pass_state_mixed" if mixed else ("malignant_status_pass_state_specific" if specific else "malignant_status_pass_state_unspecified")
            else:
                concise = status_label or "tumor_like"
                label_row = status_row
                specific = True
                mixed = False
                reason = "malignant_status_pass"
        elif status == "strong":
            # Auxiliary flag-only program sets often contain no status-role programs,
            # or tumor program sets can show state/modifier positivity without enough
            # status evidence. Summarize the strongest positive program without
            # allowing tumor integration.
            if state_label:
                label_row = state_row if state_row is not None else _sort_programs(sub.loc[strong_mask], ["decision_score", "status_score", "raw_score"]).iloc[0]
                concise = state_label
                specific = bool(state_specific)
                mixed = bool(state_mixed)
            else:
                label_row = _sort_programs(sub.loc[strong_mask], ["decision_score", "status_score", "raw_score"]).iloc[0]
                concise = _row_label(label_row)
                specific = True
                mixed = False
            reason = "auxiliary_program_support_pass" if not has_status_programs else "malignant_state_positive_without_status"
        elif status == "weak":
            if has_status_programs and bool((state_strong_mask | modifier_strong_mask).any()):
                positive_nonstatus = sub.loc[state_strong_mask | modifier_strong_mask]
                label_row = _sort_programs(positive_nonstatus, ["decision_score", "status_score", "raw_score"]).iloc[0]
                reason = "malignant_state_positive_without_status"
            else:
                label_row = aux_or_label_row
                reason = "malignant_support_weak"
            concise = _row_label(label_row)
            specific = False
            mixed = False
        else:
            label_row = aux_or_label_row
            concise = np.nan
            specific = False
            mixed = False
            reason = "malignant_support_none"

        best_name = str(label_row.get("name", "")) if label_row is not None else ""
        best_decision = _safe_float(label_row.get("decision_score", np.nan)) if label_row is not None else np.nan
        best_specificity = _safe_float(label_row.get("program_specificity_score", np.nan)) if label_row is not None else np.nan
        best_z = _safe_float(label_row.get("program_robust_zscore", np.nan)) if label_row is not None else np.nan
        second_score = np.nan
        if np.isfinite(margin) and state_row is not None and not pd.isna(second_name):
            same_group = _sort_programs(
                state_candidates[state_candidates["competition_group"].astype(str) == str(state_row.get("competition_group", ""))],
                ["decision_score", "status_score", "raw_score"],
            )
            if len(same_group) > 1:
                second_score = _safe_float(same_group.iloc[1].get("decision_score", np.nan))

        positive_programs = _format_program_list(sub.loc[raw_positive_mask, "name"].astype(str).tolist())
        strong_programs = _format_program_list(sub.loc[strong_mask, "name"].astype(str).tolist())
        status_programs = _format_program_list(status_candidates["name"].astype(str).tolist())
        if tumor_status_pass and status_program_name and status_program_name not in [x for x in status_programs.split(";") if x]:
            status_programs = _format_program_list([status_programs, status_program_name])
        state_programs = _format_program_list(state_candidates["name"].astype(str).tolist())
        modifier_programs = _format_program_list(modifier_candidates["name"].astype(str).tolist())
        positive_groups = _format_program_list(sub.loc[strong_mask, "competition_group"].astype(str).tolist())
        multigroup = len([x for x in positive_groups.split(";") if x]) > 1

        # Program-level lineage matching should use the status program when tumor
        # integration is possible; otherwise it follows the reported auxiliary label.
        report_lineage_row = status_row if tumor_status_pass else label_row
        records.append(
            {
                "cluster_id": str(cluster),
                "annot_malignant_name": best_name,
                "annot_malignant_label": concise,
                "annot_malignant_label_concise": concise,
                "annot_malignant_status_program": status_program_name,
                "annot_malignant_status_label": status_label,
                "annot_malignant_state_label": state_label,
                "annot_malignant_modifier_labels": modifier_label_text,
                "annot_malignant_status": status,
                "annot_malignant_tumor_status_pass": bool(tumor_status_pass),
                "annot_malignant_reason": reason,
                "annot_malignant_score": best_decision,
                "annot_malignant_decision_score": best_decision,
                "annot_malignant_status_score": _safe_float(combined_status_score),
                "annot_malignant_raw_score": _safe_float(combined_status_raw),
                "annot_malignant_core_status_score": _safe_float(core_status_score_aggregate),
                "annot_malignant_core_status_raw_score": _safe_float(core_status_raw),
                "annot_malignant_combined_status_score": _safe_float(combined_status_score),
                "annot_malignant_combined_status_raw_score": _safe_float(combined_status_raw),
                "annot_malignant_state_support_score": float(state_support_score),
                "annot_malignant_state_support_programs": state_support_programs,
                "annot_malignant_status_program_count": int(status_support_mask.sum()) if has_status_programs else 0,
                "annot_malignant_status_group_count": int(n_support_status_groups),
                "annot_malignant_status_decision_source": status_decision_source,
                "annot_malignant_label_raw_score": _safe_float(label_row.get("raw_score", np.nan)) if label_row is not None else np.nan,
                "annot_malignant_specificity_score": best_specificity,
                "annot_malignant_program_robust_zscore": best_z,
                "annot_malignant_margin": margin,
                "annot_malignant_second_name": second_name,
                "annot_malignant_second_score": second_score,
                "annot_malignant_competition_group": str(label_row.get("competition_group", "")) if label_row is not None else "",
                "annot_malignant_reporting_role": str(label_row.get("reporting_role", "")) if label_row is not None else "",
                "annot_malignant_program_report_on_lineages": str(report_lineage_row.get("report_on_lineages", "")) if report_lineage_row is not None else "",
                "annot_malignant_specific": bool(specific),
                "annot_malignant_mixed": bool(mixed),
                "annot_malignant_positive_programs": positive_programs,
                "annot_malignant_strong_programs": strong_programs,
                "annot_malignant_status_programs": status_programs,
                "annot_malignant_state_programs": state_programs,
                "annot_malignant_modifier_programs": modifier_programs,
                "annot_malignant_positive_groups": positive_groups,
                "annot_malignant_multigroup_positive": bool(multigroup),
                "annot_malignant_markers_present": _safe_float(label_row.get("markers_present", np.nan)) if label_row is not None else np.nan,
                "annot_malignant_markers_total": _safe_float(label_row.get("markers_total", np.nan)) if label_row is not None else np.nan,
                "annot_malignant_markers_detection_support_fraction": _safe_float(label_row.get("markers_detection_support_fraction", np.nan)) if label_row is not None else np.nan,
            }
        )

    out = pd.DataFrame(records)
    for col in _MALIGNANT_ANNOT_COLUMNS:
        if col not in out.columns:
            out[col] = np.nan
    return out[_MALIGNANT_ANNOT_COLUMNS]
