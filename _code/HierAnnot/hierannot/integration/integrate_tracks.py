from __future__ import annotations

import numpy as np
import pandas as pd

from ..malignant_reporting import resolve_program_report_block_preset


def _as_series(df: pd.DataFrame, col: str, default=np.nan, dtype=None):
    if col in df.columns:
        s = df[col]
    else:
        s = pd.Series(default, index=df.index)
    if isinstance(s, pd.DataFrame):
        s = s.iloc[:, 0]
    if dtype is not None:
        try:
            s = s.astype(dtype)
        except Exception:
            pass
    return s


def _numeric_series_with_fallback(df: pd.DataFrame, columns: list[str], default=np.nan) -> pd.Series:
    out = pd.Series(default, index=df.index, dtype="float64")
    for col in columns:
        if col not in df.columns:
            continue
        vals = pd.to_numeric(_as_series(df, col, np.nan), errors="coerce")
        out = out.where(out.notna(), vals)
    return out


def _malignant_state_from_summary(
    out: pd.DataFrame,
    *,
    malignant_status_score_threshold: float,
    malignant_raw_score_threshold: float,
) -> pd.Series:
    if "annot_malignant_status" in out.columns:
        state = out["annot_malignant_status"].astype(str).str.lower()
        valid = state.isin(["strong", "weak", "none"])
        if valid.any():
            fallback = pd.Series("none", index=out.index, dtype=object)
            fallback.loc[valid] = state.loc[valid]
            missing = ~valid
            if missing.any():
                computed = _malignant_state_from_scores(
                    out.loc[missing],
                    malignant_status_score_threshold=malignant_status_score_threshold,
                    malignant_raw_score_threshold=malignant_raw_score_threshold,
                )
                fallback.loc[missing] = computed
            return fallback
    return _malignant_state_from_scores(
        out,
        malignant_status_score_threshold=malignant_status_score_threshold,
        malignant_raw_score_threshold=malignant_raw_score_threshold,
    )


def _malignant_state_from_scores(
    out: pd.DataFrame,
    *,
    malignant_status_score_threshold: float,
    malignant_raw_score_threshold: float,
) -> pd.Series:
    malignant_status_score = pd.to_numeric(
        _as_series(out, "annot_malignant_status_score", _as_series(out, "annot_malignant_score", np.nan)),
        errors="coerce",
    )
    malignant_raw = pd.to_numeric(_as_series(out, "annot_malignant_raw_score", np.nan), errors="coerce")
    strong = malignant_status_score.ge(float(malignant_status_score_threshold)).fillna(False) & malignant_raw.ge(float(malignant_raw_score_threshold)).fillna(False)
    weak = (malignant_status_score.notna() | malignant_raw.notna()) & ~strong
    return pd.Series(np.where(strong, "strong", np.where(weak, "weak", "none")), index=out.index, dtype=object)


def _normal_state_from_scores(
    out: pd.DataFrame,
    *,
    normal_strong_score_threshold: float,
    normal_strong_raw_threshold: float,
) -> pd.Series:
    normal_score = pd.to_numeric(_as_series(out, "annot_score", np.nan), errors="coerce")
    normal_raw = _numeric_series_with_fallback(out, ["annot_branch_supported_raw_score", "annot_raw_score"])
    normal_conf = _as_series(out, "annot_confidence", "", dtype=str).str.lower()
    strong = (
        normal_score.ge(float(normal_strong_score_threshold)).fillna(False)
        & normal_raw.ge(float(normal_strong_raw_threshold)).fillna(False)
        & ~normal_conf.isin(["low", "none"])
    )
    weak = (normal_score.notna() | normal_raw.notna()) & ~strong
    return pd.Series(np.where(strong, "strong", np.where(weak, "weak", "none")), index=out.index, dtype=object)


def _build_tumor_like_label(prefix: str, malignant_state_label: pd.Series, normal_label: pd.Series) -> pd.Series:
    prefix = str(prefix).strip() or "tumor_like"
    state = malignant_state_label.astype(str).str.strip()
    normal = normal_label.astype(str)
    generic = state.str.lower().isin({"", "nan", "none", "malignant", "unspecified", "tumor_like"})
    out = prefix + "_" + state + "." + normal
    out.loc[generic] = prefix + "." + normal.loc[generic]
    return out


def _raw_delta_pass_mask(
    *,
    malignant_raw: pd.Series,
    normal_raw: pd.Series,
    threshold: float | None,
) -> pd.Series:
    if threshold is None:
        return pd.Series(True, index=malignant_raw.index, dtype=bool)
    # Missing normal raw support is treated as zero support rather than as a hard
    # failure; missing malignant raw support still fails the raw-evidence contrast.
    normal_for_delta = normal_raw.fillna(0.0)
    delta = malignant_raw - normal_for_delta
    return delta.ge(float(threshold)).fillna(False)


def _integrate_normal_and_malignant_annotations(
    normal_summary: pd.DataFrame,
    malignant_summary: pd.DataFrame | None,
    *,
    tumor_like_prefix: str = "tumor_like",
    normal_strong_score_threshold: float = 0.35,
    normal_strong_raw_threshold: float = 0.10,
    malignant_status_score_threshold: float = 0.35,
    malignant_raw_score_threshold: float = 0.15,
    malignant_normal_raw_delta_threshold: float | None = None,
    malignant_integration_mode: str = "flag_only",
    program_report_block_preset=None,
) -> pd.DataFrame:
    """Integrate normal and malignant annotation tracks.

    Malignant status is treated as a tumor-like flagging problem: strong status
    requires absolute malignant evidence. When
    ``malignant_normal_raw_delta_threshold`` is provided, tumor-like integration
    additionally requires malignant raw enrichment to exceed normal-track raw
    evidence by that amount. The default ``None`` relies on malignant absolute
    evidence and the reporting blocklist while still reporting the raw-delta
    diagnostics.
    """
    out = normal_summary.copy()

    if malignant_summary is not None and len(malignant_summary) > 0:
        m = malignant_summary.copy()
        m["cluster_id"] = m["cluster_id"].astype(str)
        if "cluster_id" in out.columns:
            out["cluster_id"] = _as_series(out, "cluster_id", "", dtype=str)
            overlap = [c for c in m.columns if c in out.columns and c != "cluster_id"]
            if overlap:
                out = out.drop(columns=overlap)
            out = out.merge(m, on="cluster_id", how="left")

    normal_score = pd.to_numeric(_as_series(out, "annot_score", np.nan), errors="coerce")
    normal_raw = _numeric_series_with_fallback(out, ["annot_branch_supported_raw_score", "annot_raw_score"])
    malignant_status_score = pd.to_numeric(
        _as_series(out, "annot_malignant_status_score", _as_series(out, "annot_malignant_score", np.nan)),
        errors="coerce",
    )
    malignant_decision_score = pd.to_numeric(
        _as_series(out, "annot_malignant_decision_score", _as_series(out, "annot_malignant_score", np.nan)),
        errors="coerce",
    )
    malignant_raw = pd.to_numeric(_as_series(out, "annot_malignant_raw_score", np.nan), errors="coerce")

    normal_state = _normal_state_from_scores(
        out,
        normal_strong_score_threshold=normal_strong_score_threshold,
        normal_strong_raw_threshold=normal_strong_raw_threshold,
    )
    malignant_state = _malignant_state_from_summary(
        out,
        malignant_status_score_threshold=malignant_status_score_threshold,
        malignant_raw_score_threshold=malignant_raw_score_threshold,
    )
    if "annot_malignant_tumor_status_pass" in out.columns:
        tumor_status_pass = _as_series(out, "annot_malignant_tumor_status_pass", False).map(lambda x: str(x).strip().lower() in {"true", "1", "yes"} if not isinstance(x, (bool, np.bool_)) else bool(x)).astype(bool)
    else:
        # Backward-compatible fallback for cached/pre-tiered summaries.
        tumor_status_pass = malignant_state.eq("strong")

    raw_delta = malignant_raw - normal_raw.fillna(0.0)
    raw_delta_pass = _raw_delta_pass_mask(
        malignant_raw=malignant_raw,
        normal_raw=normal_raw,
        threshold=malignant_normal_raw_delta_threshold,
    )

    out["annot_integrated_normal_state"] = normal_state.to_numpy(dtype=object)
    out["annot_integrated_malignant_state"] = malignant_state.to_numpy(dtype=object)
    out["annot_malignant_normal_raw_delta"] = raw_delta
    out["annot_malignant_normal_raw_delta_pass"] = raw_delta_pass.astype(bool)
    out["annot_malignant_normal_status_score_delta"] = malignant_status_score - normal_score
    out["annot_malignant_normal_decision_score_delta"] = malignant_decision_score - normal_score

    normal_label = _as_series(out, "annot_label", "unknown", dtype=str)
    malignant_name = _as_series(out, "annot_malignant_name", "", dtype=str)
    malignant_label = _as_series(out, "annot_malignant_label_concise", _as_series(out, "annot_malignant_label", "malignant"), dtype=str)
    malignant_label = malignant_label.where(malignant_label.ne(""), _as_series(out, "annot_malignant_label", "malignant", dtype=str))
    malignant_state_label = malignant_label.where(malignant_label.ne(""), malignant_name)

    report_block = resolve_program_report_block_preset(out, program_report_block_preset)
    malignant_block_mask = report_block["blocked"].astype(bool)
    out["annot_malignant_reporting_blocked"] = malignant_block_mask
    out["annot_malignant_reporting_blocked_reason"] = report_block["reason"].to_numpy(dtype=object)
    out["annot_malignant_global_blocked"] = report_block["global_blocked"].astype(bool).to_numpy()
    out["annot_malignant_program_lineage_blocked"] = report_block["program_lineage_blocked"].astype(bool).to_numpy()
    out["annot_malignant_program_lineage_matched"] = report_block["program_lineage_matched"].astype(bool).to_numpy()
    out["annot_malignant_program_report_on_lineages"] = report_block["program_report_on_lineages"].to_numpy(dtype=object)

    out["annot_integrated_label"] = normal_label.copy()
    out["annot_integrated_status"] = pd.Series("resolved", index=out.index, dtype=object)
    out["annot_integrated_reason"] = pd.Series("normal_only", index=out.index, dtype=object)
    out["annot_integrated_source"] = pd.Series("normal", index=out.index, dtype=object)

    weakweak = (malignant_state == "weak") & (normal_state == "weak")
    nonenone = (malignant_state == "none") & (normal_state == "none")
    out.loc[weakweak | nonenone, "annot_integrated_label"] = "unknown"
    out.loc[weakweak, "annot_integrated_status"] = "unknown"
    out.loc[weakweak, "annot_integrated_reason"] = "both_weak"
    out.loc[weakweak, "annot_integrated_source"] = "unknown"
    out.loc[nonenone, "annot_integrated_status"] = "unknown"
    out.loc[nonenone, "annot_integrated_reason"] = "no_support"
    out.loc[nonenone, "annot_integrated_source"] = "unknown"

    mode = str(malignant_integration_mode).lower()
    if mode == "off":
        return out

    normal_strong = normal_state == "strong"
    malignant_strong = malignant_state == "strong"
    malignant_status_ready = malignant_strong & tumor_status_pass
    malignant_reportable = malignant_status_ready & ~malignant_block_mask & raw_delta_pass
    malignant_below_delta = malignant_status_ready & ~malignant_block_mask & ~raw_delta_pass
    malignant_no_status_program = malignant_strong & ~tumor_status_pass

    if mode == "flag_only":
        both = malignant_reportable & normal_strong
        out.loc[both, "annot_integrated_reason"] = "malignant_flag_both_strong"
        mal_only = malignant_reportable & ~normal_strong
        out.loc[mal_only, "annot_integrated_reason"] = "malignant_flag_strong"
        blocked = malignant_status_ready & malignant_block_mask
        out.loc[blocked, "annot_integrated_reason"] = "malignant_blocked_by_blocklist"
        out.loc[malignant_below_delta, "annot_integrated_reason"] = "malignant_below_raw_delta"
        out.loc[malignant_no_status_program, "annot_integrated_reason"] = "auxiliary_program_flag_only"
        return out

    mask = normal_strong & ~malignant_reportable
    out.loc[mask, "annot_integrated_label"] = normal_label[mask]
    out.loc[mask, "annot_integrated_status"] = "resolved"
    out.loc[mask, "annot_integrated_reason"] = "normal_dominant"
    out.loc[mask, "annot_integrated_source"] = "normal"

    blocked_strong = malignant_status_ready & malignant_block_mask
    out.loc[blocked_strong, "annot_integrated_label"] = normal_label[blocked_strong]
    out.loc[blocked_strong, "annot_integrated_status"] = "resolved"
    out.loc[blocked_strong, "annot_integrated_reason"] = "malignant_blocked_by_blocklist"
    out.loc[blocked_strong, "annot_integrated_source"] = "normal"

    out.loc[malignant_below_delta, "annot_integrated_label"] = normal_label[malignant_below_delta]
    out.loc[malignant_below_delta, "annot_integrated_status"] = "resolved"
    out.loc[malignant_below_delta, "annot_integrated_reason"] = "malignant_below_raw_delta"
    out.loc[malignant_below_delta, "annot_integrated_source"] = "normal"

    out.loc[malignant_no_status_program, "annot_integrated_label"] = normal_label[malignant_no_status_program]
    out.loc[malignant_no_status_program, "annot_integrated_status"] = "resolved"
    out.loc[malignant_no_status_program, "annot_integrated_reason"] = "auxiliary_program_not_integrated"
    out.loc[malignant_no_status_program, "annot_integrated_source"] = "normal"

    tumor_mask = malignant_reportable
    out.loc[tumor_mask, "annot_integrated_label"] = _build_tumor_like_label(
        tumor_like_prefix,
        malignant_state_label[tumor_mask],
        normal_label[tumor_mask],
    )
    out.loc[tumor_mask, "annot_integrated_status"] = "tumor_like"
    out.loc[tumor_mask & ~normal_strong, "annot_integrated_reason"] = "malignant_dominant"
    out.loc[tumor_mask & normal_strong, "annot_integrated_reason"] = "malignant_flag_over_normal"
    out.loc[tumor_mask & ~normal_strong, "annot_integrated_source"] = "malignant"
    out.loc[tumor_mask & normal_strong, "annot_integrated_source"] = "combined"

    return out
