from __future__ import annotations

from typing import Iterable, Optional

import pandas as pd
import numpy as np


def _require_matplotlib():
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover - optional dependency
        raise ImportError("Plotting requires matplotlib. Install it separately to use plotting helpers.") from exc
    return plt



def _result_attr(result, name: str):
    if isinstance(result, dict):
        return result.get(name)
    return getattr(result, name)


def _score_heatmap_cmap():
    return "RdBu_r"

def plot_cluster_annotation_heatmap(result, value: str = "score", level: Optional[int] = None, top_n_labels: Optional[int] = None):
    """Plot a cluster-by-label score heatmap from a HierAnnotResult.

    Parameters
    ----------
    result:
        HierAnnotResult from ``HierAnnotPipeline.fit_score``.
    value:
        Which score-like column from ``result.all_scores`` to visualize.
    level:
        Optional hierarchy level filter.
    top_n_labels:
        Optionally keep only the top labels by mean absolute value.
    """
    plt = _require_matplotlib()
    df = _result_attr(result, "all_scores")
    if df is None:
        raise ValueError("result does not contain all_scores")
    df = df.copy()
    if level is not None:
        df = df[df["level"] == level]
    if value not in df.columns:
        raise KeyError(f"'{value}' not found in result.all_scores")
    heatmap = df.pivot_table(index="label", columns="cluster", values=value, aggfunc="mean")
    heatmap = heatmap.fillna(0.0)
    if top_n_labels is not None and top_n_labels < len(heatmap):
        order = heatmap.abs().mean(axis=1).sort_values(ascending=False).index[:top_n_labels]
        heatmap = heatmap.loc[order]
    fig, ax = plt.subplots(figsize=(max(6, 1 + 0.6 * heatmap.shape[1]), max(4, 1 + 0.35 * heatmap.shape[0])))
    vmax = float(max(abs(heatmap.values.min()), abs(heatmap.values.max()))) if heatmap.size else 1.0
    im = ax.imshow(heatmap.values, aspect="auto", cmap=_score_heatmap_cmap(), vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(heatmap.shape[1]))
    ax.set_xticklabels(list(heatmap.columns), rotation=45, ha="right")
    ax.set_yticks(range(heatmap.shape[0]))
    ax.set_yticklabels(list(heatmap.index))
    ax.set_xlabel("Cluster")
    ax.set_ylabel("Label")
    title = f"HierAnnot {value} heatmap"
    if level is not None:
        title += f" (level {level})"
    ax.set_title(title)
    fig.colorbar(im, ax=ax, label=value)
    fig.tight_layout()
    return fig, ax


def plot_cluster_confidence(result):
    """Plot final score versus ambiguity margin, colored by confidence tier."""
    plt = _require_matplotlib()
    df = _result_attr(result, "cluster_annotations")
    if df is None:
        raise ValueError("result does not contain cluster_annotations")
    df = df.copy()

    x_candidates = ["annot_score", "annot_final_score", "annot_decision_score"]
    y_candidates = ["annot_final_call_margin", "annot_top_level_gap", "annot_path_margin", "annot_margin"]
    xcol = next((c for c in x_candidates if c in df.columns), None)
    ycol = next((c for c in y_candidates if c in df.columns), None)
    if xcol is None or ycol is None:
        raise ValueError("cluster_annotations missing confidence plot columns")

    df[xcol] = pd.to_numeric(df[xcol], errors="coerce")
    df[ycol] = pd.to_numeric(df[ycol], errors="coerce")
    sub = df[np.isfinite(df[xcol]) & np.isfinite(df[ycol])].copy()

    fig, ax = plt.subplots(figsize=(7, 5))
    if sub.empty:
        ax.text(0.5, 0.5, "No finite confidence points to plot", ha="center", va="center")
        ax.set_axis_off()
        return fig, ax

    confidence_col = "annot_confidence" if "annot_confidence" in sub.columns else ("confidence" if "confidence" in sub.columns else None)
    if confidence_col is None:
        ax.scatter(sub[xcol], sub[ycol])
    else:
        for conf, grp in sub.groupby(confidence_col, dropna=False):
            ax.scatter(grp[xcol], grp[ycol], label=str(conf))
        ax.legend(title="Confidence")

    if "cluster_id" in sub.columns:
        for _, row in sub.iterrows():
            ax.annotate(str(row["cluster_id"]), (row[xcol], row[ycol]), fontsize=8)

    ax.axhline(0.0, linestyle="--", linewidth=1)
    ax.axvline(0.0, linestyle="--", linewidth=1)
    ax.set_xlabel("Final score")
    ax.set_ylabel(ycol.replace("annot_", "").replace("_", " "))
    ax.set_title("Cluster confidence diagnostics")
    fig.tight_layout()
    return fig, ax


def plot_top_level_compartment_scores(result, clusters: Optional[Iterable[str]] = None):
    """Plot level-1 scores across clusters to inspect mixed-compartment structure."""
    plt = _require_matplotlib()
    df = _result_attr(result, "all_scores")
    if df is None:
        raise ValueError("result does not contain all_scores")
    df = df.copy()
    df = df[df["level"] == 1]
    if clusters is not None:
        cluster_set = {str(c) for c in clusters}
        df = df[df["cluster"].astype(str).isin(cluster_set)]
    table = df.pivot_table(index="cluster", columns="label", values="score", aggfunc="mean").fillna(0.0)
    fig, ax = plt.subplots(figsize=(max(6, 1 + 0.7 * table.shape[1]), max(4, 1 + 0.45 * table.shape[0])))
    table.plot(kind="bar", ax=ax)
    ax.set_ylabel("Score")
    ax.set_title("Top-level compartment scores by cluster")
    ax.legend(title="Label", bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    return fig, ax


def plot_cluster_support_scatter(result, label_points: bool = True):
    """Plot local-vs-propagated cluster support diagnostics."""
    plt = _require_matplotlib()
    df = _result_attr(result, "cluster_annotations")
    if df is None:
        raise ValueError("result does not contain cluster_annotations")
    df = df.copy()
    xcol = "annot_score" if "annot_score" in df.columns else ("annot_decision_score" if "annot_decision_score" in df.columns else None)
    ycol = "annot_branch_supported_raw_score" if "annot_branch_supported_raw_score" in df.columns else None
    if xcol is None or ycol is None:
        raise ValueError("cluster_annotations missing required support columns")
    sub = df[[c for c in ["cluster_id", xcol, ycol, "annot_export_status"] if c in df.columns]].copy()
    sub[xcol] = pd.to_numeric(sub[xcol], errors="coerce")
    sub[ycol] = pd.to_numeric(sub[ycol], errors="coerce")
    sub = sub[np.isfinite(sub[xcol]) | np.isfinite(sub[ycol])].copy()
    fig, ax = plt.subplots(figsize=(7, 5))
    if sub.empty:
        ax.text(0.5, 0.5, "No finite cluster support values to plot", ha="center", va="center")
        ax.set_axis_off()
        return fig, ax
    status_col = "annot_export_status" if "annot_export_status" in sub.columns else None
    if status_col is None:
        ax.scatter(sub[xcol], sub[ycol])
    else:
        for status, grp in sub.groupby(status_col, dropna=False):
            ax.scatter(grp[xcol], grp[ycol], label=str(status))
        ax.legend(title="Export status")
    if label_points and "cluster_id" in sub.columns:
        for _, row in sub.iterrows():
            ax.annotate(str(row["cluster_id"]), (row[xcol], row[ycol]), fontsize=8)
    ax.axhline(0.0, linestyle="--", linewidth=1)
    ax.axvline(0.0, linestyle="--", linewidth=1)
    ax.set_xlabel("Final / export score")
    ax.set_ylabel("Branch-supported raw score")
    ax.set_title("Cluster support diagnostics")
    fig.tight_layout()
    return fig, ax


def plot_backoff_diagnostics(result):
    """Plot clusters where absolute-support backoff was applied."""
    plt = _require_matplotlib()
    df = _result_attr(result, "cluster_annotations")
    if df is None:
        raise ValueError("result does not contain cluster_annotations")
    df = df.copy()
    if "annot_backoff_applied" not in df.columns:
        raise ValueError("result.cluster_annotations does not contain annot_backoff_applied")
    sub = df[df["annot_backoff_applied"].fillna(False).astype(bool)].copy()
    fig, ax = plt.subplots(figsize=(8, max(3, 0.5 * max(len(sub), 1))))
    if sub.empty:
        ax.text(0.5, 0.5, "No backoff-applied clusters", ha="center", va="center")
        ax.set_axis_off()
        return fig, ax
    y = range(len(sub))
    ax.scatter(sub["annot_decision_branch_supported_raw_score"], y, label="Decision node")
    ax.scatter(sub["annot_branch_supported_raw_score"], y, label="Final backed-off node")
    for i, (_, row) in enumerate(sub.iterrows()):
        ax.text(row["annot_branch_supported_raw_score"], i, f'  {row["cluster_id"]}: {row["annot_decision_label"]} -> {row["annot_label"]}', va="center", fontsize=8)
    ax.axvline(0.0, linestyle="--", linewidth=1)
    ax.set_yticks([])
    ax.set_xlabel("Branch-supported raw score")
    ax.set_title("Absolute-support backoff diagnostics")
    ax.legend()
    fig.tight_layout()
    return fig, ax



def plot_hierarchy_fit_summary(result_or_summary):
    """
    Plot a compact two-panel hierarchy-fit diagnostic summary.

    Accepts either:
    - a HierAnnotResult-like object with ``diagnostics_summary``
    - or a precomputed fit-summary DataFrame from ``summarize_hierarchy_fit()``
    """
    plt = _require_matplotlib()

    if isinstance(result_or_summary, pd.DataFrame):
        df = result_or_summary.copy()
    else:
        diag = _result_attr(result_or_summary, "diagnostics_summary")
        if diag is None:
            raise ValueError("result does not contain diagnostics_summary")
        df = diag.copy()

    if "metric" not in df.columns or "value" not in df.columns:
        raise ValueError("Expected a fit-summary style DataFrame with 'metric' and 'value' columns")

    def _metric_value(metric, default=float("nan")):
        sub = df[df["metric"] == metric]
        if sub.empty:
            return default
        return sub["value"].iloc[0]

    status = str(_metric_value("hierarchy_fit_status", "unknown"))
    warning = str(_metric_value("hierarchy_fit_warning", ""))

    frac_metrics = [
        ("fraction_unknown", "Unknown"),
        ("fraction_low_evidence", "Low evidence"),
        ("fraction_low_absolute_support", "Low abs support"),
        ("fraction_backoff_applied", "Backoff"),
        ("fraction_parent_mixing", "Parent mixing"),
        ("fraction_branch_conflict", "Branch conflict"),
    ]
    support_metrics = [
        ("median_final_score", "Median final score"),
        ("median_branch_supported_raw_score", "Median branch raw"),
        ("median_branch_support_score", "Median branch score"),
        ("median_local_node_score", "Median local score"),
        ("median_marker_detection_support_fraction", "Median marker support"),
        ("fraction_clusters_final_score_ge_1", "Frac final score >= 1"),
        ("fraction_clusters_branch_supported_raw_ge_0_25", "Frac branch raw >= 0.25"),
        ("fraction_l1_positive_branch_support", "Frac L1 branch > 0"),
    ]

    frac_vals = [pd.to_numeric(_metric_value(m), errors="coerce") for m, _ in frac_metrics]
    supp_vals = [pd.to_numeric(_metric_value(m), errors="coerce") for m, _ in support_metrics]

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))

    ax = axes[0]
    ax.bar([label for _, label in frac_metrics], frac_vals)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Fraction")
    ax.set_title("Failure-mode fractions")
    ax.tick_params(axis="x", rotation=45)

    ax = axes[1]
    ax.bar([label for _, label in support_metrics], supp_vals)
    ax.axhline(1.0, linestyle="--", linewidth=1)
    ax.axhline(0.25, linestyle="--", linewidth=1)
    ax.axhline(0.60, linestyle="--", linewidth=1)
    ax.set_ylabel("Metric value")
    ax.set_title("Support-strength metrics")
    ax.tick_params(axis="x", rotation=45)

    fig.suptitle(f"Hierarchy fit summary ({status})", y=1.02)
    if warning and warning != "nan":
        fig.text(0.5, 0.01, warning, ha="center", va="bottom", fontsize=9)
    fig.tight_layout()
    return fig, axes

def plot_compare_hierarchy_fit(comparison_df):
    """
    Plot a side-by-side comparison of selected hierarchy-fit metrics.

    Accepts either:
    - the DataFrame returned by `compare_hierarchy_fit()`, or
    - a single `HierAnnotResult`-like object from `pipeline.fit_score()`.
    """
    plt = _require_matplotlib()

    # Accept a direct result object by converting its diagnostics summary
    if isinstance(comparison_df, pd.DataFrame):
        df = comparison_df.copy()
    else:
        diag = _result_attr(comparison_df, "diagnostics_summary")
        if diag is None:
            raise ValueError("Expected either a comparison DataFrame or a result object with diagnostics_summary")
        if not isinstance(diag, pd.DataFrame) or "metric" not in diag.columns or "value" not in diag.columns:
            raise ValueError("diagnostics_summary must be a DataFrame with 'metric' and 'value' columns")
        row = {}
        for metric in [
            "hierarchy_fit_status",
            "fit_rank_score",
            "median_branch_support_score",
            "median_branch_supported_raw_score",
            "median_final_score",
            "fraction_l1_positive_branch_support",
            "median_marker_detection_support_fraction",
            "fraction_unknown",
            "fraction_low_absolute_support",
            "fraction_backoff_applied",
        ]:
            sub = diag[diag["metric"] == metric]
            row[metric] = sub["value"].iloc[0] if not sub.empty else pd.NA
        row["name"] = "current_result"
        df = pd.DataFrame([row])

    if df.empty:
        raise ValueError("comparison_df is empty")

    if "name" not in df.columns:
        df = df.copy()
        df["name"] = [f"candidate_{i+1}" for i in range(len(df))]

    # Backward/forward compatible metric selection
    candidate_metrics = [
        ("median_final_score", "Median final score"),
        ("median_branch_supported_raw_score", "Median branch raw"),
        ("median_marker_detection_support_fraction", "Median marker support"),
        ("fraction_l1_positive_branch_support", "Frac L1 branch > 0"),
        ("fit_rank_score", "Fit rank score"),
        ("median_branch_support_score", "Median branch score"),
    ]
    available = [(col, label) for col, label in candidate_metrics if col in df.columns]
    if not available:
        raise ValueError("comparison_df does not contain any recognized hierarchy-fit metrics to plot")

    # Limit to 3 panels for readability, prioritizing the most informative available metrics
    preferred_order = [
        "median_final_score",
        "median_branch_supported_raw_score",
        "median_marker_detection_support_fraction",
        "fraction_l1_positive_branch_support",
        "fit_rank_score",
        "median_branch_support_score",
    ]
    available_sorted = []
    for col in preferred_order:
        for metric in available:
            if metric[0] == col:
                available_sorted.append(metric)
    metrics = available_sorted[:3]

    n = len(metrics)
    fig, axes = plt.subplots(1, n, figsize=(4 * n, 4))
    if n == 1:
        axes = [axes]

    for ax, (metric, title) in zip(axes, metrics):
        vals = pd.to_numeric(df[metric], errors="coerce")
        sub = pd.DataFrame({"name": df["name"].astype(str), metric: vals}).dropna(subset=[metric]).copy()
        if sub.empty:
            ax.text(0.5, 0.5, f"No finite values for {metric}", ha="center", va="center")
            ax.set_axis_off()
            continue
        ax.bar(sub["name"], sub[metric])
        ax.set_title(title)
        ax.tick_params(axis="x", rotation=45)
        ax.set_ylabel(metric)

    fig.tight_layout()
    return fig, axes





def _annotation_points(ax, df: pd.DataFrame, xcol: str, ycol: str, label_points: bool):
    if label_points and "cluster_id" in df.columns:
        for _, row in df.iterrows():
            if pd.notna(row.get(xcol)) and pd.notna(row.get(ycol)):
                ax.annotate(str(row["cluster_id"]), (row[xcol], row[ycol]), fontsize=8)


def plot_malignant_score_heatmap(result, value: str = "status_score", top_n_programs: Optional[int] = None):
    """Plot a cluster-by-malignant-program heatmap from ``result.malignant_scores``.

    ``status_score`` is the default because it tracks absolute tumor-like
    evidence. Use ``decision_score`` or ``program_specificity_score`` to inspect
    malignant-state specificity instead.
    """
    plt = _require_matplotlib()
    df = _result_attr(result, "malignant_scores")
    if df is None or len(df) == 0:
        raise ValueError("result does not contain malignant_scores")
    df = df.copy()
    requested_value = value
    if value not in df.columns and value == "status_score" and "score" in df.columns:
        value = "score"
    if value not in df.columns:
        raise KeyError(f"{requested_value!r} not found in result.malignant_scores")
    name_col = "name" if "name" in df.columns else "label"
    heatmap = df.pivot_table(index=name_col, columns="cluster", values=value, aggfunc="mean")
    heatmap = heatmap.fillna(0.0)
    if top_n_programs is not None and top_n_programs < len(heatmap):
        order = heatmap.abs().mean(axis=1).sort_values(ascending=False).index[:top_n_programs]
        heatmap = heatmap.loc[order]
    fig, ax = plt.subplots(figsize=(max(6, 1 + 0.6 * heatmap.shape[1]), max(4, 1 + 0.35 * heatmap.shape[0])))
    vmax = float(max(abs(heatmap.values.min()), abs(heatmap.values.max()))) if heatmap.size else 1.0
    im = ax.imshow(heatmap.values, aspect="auto", cmap=_score_heatmap_cmap(), vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(heatmap.shape[1]))
    ax.set_xticklabels(list(heatmap.columns), rotation=45, ha="right")
    ax.set_yticks(range(heatmap.shape[0]))
    ax.set_yticklabels(list(heatmap.index))
    ax.set_xlabel("Cluster")
    ax.set_ylabel("Malignant program")
    title_kind = {
        "raw_score": "Malignant raw enrichment",
        "status_score": "Malignant status evidence",
        "absolute_support_score": "Malignant absolute support",
        "decision_score": "Malignant state decision",
        "score": "Malignant state decision",
        "program_specificity_score": "Malignant state specificity",
        "specificity_score": "Malignant state specificity",
        "program_robust_zscore": "Malignant robust z-score",
    }.get(value, f"Malignant {value}")
    ax.set_title(f"{title_kind} heatmap")
    fig.colorbar(im, ax=ax, label=value)
    fig.tight_layout()
    return fig, ax


def plot_malignant_status_scatter(
    result,
    malignant_status_score_threshold: float = 0.35,
    malignant_raw_score_threshold: float = 0.15,
    label_points: bool = True,
):
    """Plot absolute malignant status evidence for each cluster."""
    plt = _require_matplotlib()
    df = _result_attr(result, "malignant_annotations")
    if df is None or len(df) == 0:
        raise ValueError("result does not contain malignant_annotations")
    df = df.copy()
    xcol = "annot_malignant_raw_score"
    ycol = "annot_malignant_status_score" if "annot_malignant_status_score" in df.columns else "annot_malignant_score"
    if xcol not in df.columns or ycol not in df.columns:
        raise ValueError("malignant_annotations missing malignant raw/status score columns")
    df[xcol] = pd.to_numeric(df[xcol], errors="coerce")
    df[ycol] = pd.to_numeric(df[ycol], errors="coerce")
    sub = df[np.isfinite(df[xcol]) | np.isfinite(df[ycol])].copy()
    fig, ax = plt.subplots(figsize=(7, 5))
    if sub.empty:
        ax.text(0.5, 0.5, "No finite malignant status values to plot", ha="center", va="center")
        ax.set_axis_off()
        return fig, ax
    status_col = "annot_malignant_status" if "annot_malignant_status" in sub.columns else None
    if status_col:
        for status, grp in sub.groupby(status_col, dropna=False):
            ax.scatter(grp[xcol], grp[ycol], label=str(status))
        ax.legend(title="Malignant status")
    else:
        ax.scatter(sub[xcol], sub[ycol])
    _annotation_points(ax, sub, xcol, ycol, label_points)
    ax.axvline(float(malignant_raw_score_threshold), linestyle="--", linewidth=1)
    ax.axhline(float(malignant_status_score_threshold), linestyle="--", linewidth=1)
    ax.set_xlabel("Malignant raw enrichment")
    ax.set_ylabel("Malignant status score")
    ax.set_title("Malignant status evidence")
    fig.tight_layout()
    return fig, ax


def plot_malignant_state_specificity(
    result,
    specificity_threshold: float = 0.0,
    label_points: bool = True,
):
    """Plot tiered malignant-state support and program-level state specificity.

    The tiered malignant/tumor reporting layer no longer relies on a legacy
    per-cluster ``annot_malignant_margin`` column. This diagnostic instead
    shows two complementary views:

    1. Annotation-level state support versus combined tumor-status evidence.
    2. Program-level state specificity from ``result.malignant_scores`` for
       rows with ``reporting_role == 'state'``.

    Returns
    -------
    fig, axes
        Matplotlib figure and the two axes.
    """
    plt = _require_matplotlib()
    ann = _result_attr(result, "malignant_annotations")
    scores = _result_attr(result, "malignant_scores")
    if ann is None or len(ann) == 0:
        raise ValueError("result does not contain malignant_annotations")

    ann = ann.copy()
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    ax0, ax1 = axes

    # Panel 1: per-cluster state support versus combined tumor-status score.
    xcol = "annot_malignant_state_support_score"
    ycol = "annot_malignant_status_score"
    if ycol not in ann.columns:
        ycol = "annot_malignant_combined_status_score" if "annot_malignant_combined_status_score" in ann.columns else None
    if xcol in ann.columns and ycol is not None and ycol in ann.columns:
        ann[xcol] = pd.to_numeric(ann[xcol], errors="coerce")
        ann[ycol] = pd.to_numeric(ann[ycol], errors="coerce")
        sub = ann[np.isfinite(ann[xcol]) | np.isfinite(ann[ycol])].copy()
    else:
        sub = pd.DataFrame()

    if sub.empty:
        ax0.text(0.5, 0.5, "No finite state-support values to plot", ha="center", va="center")
        ax0.set_axis_off()
    else:
        group_col = "annot_malignant_status_decision_source" if "annot_malignant_status_decision_source" in sub.columns else None
        if group_col:
            for group, grp in sub.groupby(group_col, dropna=False):
                ax0.scatter(grp[xcol], grp[ycol], label=str(group))
            ax0.legend(title="Decision source", fontsize=8)
        else:
            ax0.scatter(sub[xcol], sub[ycol])
        _annotation_points(ax0, sub, xcol, ycol, label_points)
        ax0.axvline(0.0, linestyle="--", linewidth=1)
        ax0.set_xlabel("State-support bonus")
        ax0.set_ylabel("Combined tumor-status score")
        ax0.set_title("Cluster-level state support")

    # Panel 2: program-level state specificity from malignant_scores.
    if scores is None or len(scores) == 0:
        state_scores = pd.DataFrame()
    else:
        state_scores = scores.copy()
        if "reporting_role" in state_scores.columns:
            state_scores = state_scores[state_scores["reporting_role"].astype(str).str.lower().eq("state")]
        elif "role" in state_scores.columns:
            state_scores = state_scores[state_scores["role"].astype(str).str.lower().eq("state")]
        # If role annotations are absent, keep rows with a state-like program label/name.
        elif "name" in state_scores.columns:
            state_scores = state_scores[state_scores["name"].astype(str).str.contains("emt|mesenchymal|state", case=False, na=False)]
        else:
            state_scores = pd.DataFrame()

    sx = "status_score" if "status_score" in state_scores.columns else ("raw_score" if "raw_score" in state_scores.columns else None)
    sy = "program_specificity_score" if "program_specificity_score" in state_scores.columns else ("decision_score" if "decision_score" in state_scores.columns else None)
    if sx is not None and sy is not None and len(state_scores) > 0:
        state_scores[sx] = pd.to_numeric(state_scores[sx], errors="coerce")
        state_scores[sy] = pd.to_numeric(state_scores[sy], errors="coerce")
        state_sub = state_scores[np.isfinite(state_scores[sx]) | np.isfinite(state_scores[sy])].copy()
    else:
        state_sub = pd.DataFrame()

    if state_sub.empty:
        ax1.text(0.5, 0.5, "No finite state-program specificity values to plot", ha="center", va="center")
        ax1.set_axis_off()
    else:
        group_col = "competition_group" if "competition_group" in state_sub.columns else None
        if group_col:
            for group, grp in state_sub.groupby(group_col, dropna=False):
                ax1.scatter(grp[sx], grp[sy], label=str(group))
            ax1.legend(title="Competition group", fontsize=8)
        else:
            ax1.scatter(state_sub[sx], state_sub[sy])

        if label_points:
            label_col = None
            for cand in ["cluster_id", "cluster"]:
                if cand in state_sub.columns:
                    label_col = cand
                    break
            name_col = None
            for cand in ["name", "program", "label"]:
                if cand in state_sub.columns:
                    name_col = cand
                    break
            if label_col is not None:
                for _, row in state_sub.iterrows():
                    label = str(row[label_col])
                    if name_col is not None:
                        label = f"{label}:{row[name_col]}"
                    ax1.annotate(label, (row[sx], row[sy]), fontsize=8)

        ax1.axhline(float(specificity_threshold), linestyle="--", linewidth=1)
        ax1.set_xlabel(sx.replace("_", " "))
        ax1.set_ylabel(sy.replace("_", " "))
        ax1.set_title("Program-level state specificity")

    fig.suptitle("Malignant state support and specificity", y=1.02)
    fig.tight_layout()
    return fig, axes


def plot_malignant_annotation_summary(result):
    """Plot malignant annotation status counts and concise state-label counts."""
    plt = _require_matplotlib()
    df = _result_attr(result, "malignant_annotations")
    if df is None or len(df) == 0:
        raise ValueError("result does not contain malignant_annotations")
    df = df.copy()
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    status_counts = df["annot_malignant_status"].astype(str).value_counts()
    axes[0].bar(status_counts.index.tolist(), status_counts.values.tolist())
    axes[0].set_title("Malignant status counts")
    axes[0].set_ylabel("Clusters")
    axes[0].tick_params(axis="x", rotation=30)

    label_col = "annot_malignant_label_concise" if "annot_malignant_label_concise" in df.columns else ("annot_malignant_label" if "annot_malignant_label" in df.columns else "annot_malignant_name")
    prog_counts = df[label_col].astype(str).replace("nan", pd.NA).dropna().value_counts().head(12)
    axes[1].bar(prog_counts.index.tolist(), prog_counts.values.tolist())
    axes[1].set_title("Top malignant state labels")
    axes[1].set_ylabel("Clusters")
    axes[1].tick_params(axis="x", rotation=45)
    fig.tight_layout()
    return fig, axes


def plot_integration_raw_evidence(
    result,
    malignant_normal_raw_delta_threshold: Optional[float] = None,
    label_points: bool = True,
):
    """Plot normal raw evidence against malignant raw evidence used in integration."""
    plt = _require_matplotlib()
    df = _result_attr(result, "integrated_annotations")
    if df is None or len(df) == 0:
        raise ValueError("result does not contain integrated_annotations")
    df = df.copy()
    x_candidates = ["annot_branch_supported_raw_score", "annot_raw_score"]
    xcol = next((c for c in x_candidates if c in df.columns), None)
    ycol = "annot_malignant_raw_score"
    if xcol is None or ycol not in df.columns:
        raise ValueError("integrated_annotations missing normal/malignant raw-score columns")
    df[xcol] = pd.to_numeric(df[xcol], errors="coerce")
    df[ycol] = pd.to_numeric(df[ycol], errors="coerce")
    sub = df[np.isfinite(df[xcol]) | np.isfinite(df[ycol])].copy()
    fig, ax = plt.subplots(figsize=(7, 5))
    if sub.empty:
        ax.text(0.5, 0.5, "No finite integration raw evidence values to plot", ha="center", va="center")
        ax.set_axis_off()
        return fig, ax
    status_col = "annot_integrated_status" if "annot_integrated_status" in sub.columns else None
    if status_col:
        for status, grp in sub.groupby(status_col, dropna=False):
            ax.scatter(grp[xcol], grp[ycol], label=str(status))
        ax.legend(title="Integrated status")
    else:
        ax.scatter(sub[xcol], sub[ycol])
    _annotation_points(ax, sub, xcol, ycol, label_points)
    finite_vals = pd.concat([sub[xcol], sub[ycol]]).dropna()
    if not finite_vals.empty:
        lo = float(finite_vals.min())
        hi = float(finite_vals.max())
        pad = max(0.05, 0.05 * (hi - lo if hi > lo else 1.0))
        xs = np.array([lo - pad, hi + pad])
        ax.plot(xs, xs, linestyle="--", linewidth=1)
        if malignant_normal_raw_delta_threshold is not None:
            ax.plot(xs, xs + float(malignant_normal_raw_delta_threshold), linestyle=":", linewidth=1)
    ax.set_xlabel("Normal branch-supported raw evidence")
    ax.set_ylabel("Malignant raw evidence")
    ax.set_title("Normal vs malignant raw evidence")
    fig.tight_layout()
    return fig, ax


def plot_integration_summary(result):
    """Plot integrated annotation status, reason, and malignant-gate diagnostics."""
    plt = _require_matplotlib()
    df = _result_attr(result, "integrated_annotations")
    if df is None or len(df) == 0:
        raise ValueError("result does not contain integrated_annotations")
    df = df.copy()
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    status_counts = df["annot_integrated_status"].astype(str).value_counts()
    axes[0].bar(status_counts.index.tolist(), status_counts.values.tolist())
    axes[0].set_title("Integrated status counts")
    axes[0].set_ylabel("Clusters")
    axes[0].tick_params(axis="x", rotation=30)

    reason_counts = df["annot_integrated_reason"].astype(str).value_counts().head(12)
    axes[1].bar(reason_counts.index.tolist(), reason_counts.values.tolist())
    axes[1].set_title("Integrated reason counts")
    axes[1].set_ylabel("Clusters")
    axes[1].tick_params(axis="x", rotation=45)

    malignant_state = df["annot_integrated_malignant_state"].astype(str).str.lower() if "annot_integrated_malignant_state" in df.columns else pd.Series("none", index=df.index)
    blocked = df["annot_malignant_reporting_blocked"].fillna(False).astype(bool) if "annot_malignant_reporting_blocked" in df.columns else pd.Series(False, index=df.index)
    raw_pass = df["annot_malignant_normal_raw_delta_pass"].fillna(True).astype(bool) if "annot_malignant_normal_raw_delta_pass" in df.columns else pd.Series(True, index=df.index)
    gate = pd.Series("malignant weak/none", index=df.index, dtype=object)
    gate.loc[malignant_state.eq("strong") & blocked] = "strong + blocked"
    gate.loc[malignant_state.eq("strong") & ~blocked & ~raw_pass] = "strong + raw-delta fail"
    gate.loc[malignant_state.eq("strong") & ~blocked & raw_pass] = "strong + reportable"
    gate_counts = gate.value_counts()
    axes[2].bar(gate_counts.index.tolist(), gate_counts.values.tolist())
    axes[2].set_title("Malignant integration gates")
    axes[2].set_ylabel("Clusters")
    axes[2].tick_params(axis="x", rotation=45)
    fig.tight_layout()
    return fig, axes


def plot_all_available_diagnostics(result):
    """Generate available diagnostic plots in one call.

    Returns
    -------
    dict
        Mapping from plot name to ``(fig, axes_or_ax)``. Plots whose required
        components are missing are skipped safely.
    """
    plots = {}
    candidates = [
        ("cluster_annotation_heatmap", lambda: plot_cluster_annotation_heatmap(result)),
        ("confidence", lambda: plot_cluster_confidence(result)),
        ("support_scatter", lambda: plot_cluster_support_scatter(result)),
        ("backoff", lambda: plot_backoff_diagnostics(result)),
        ("top_level", lambda: plot_top_level_compartment_scores(result)),
        ("hierarchy_fit", lambda: plot_hierarchy_fit_summary(result)),
        ("malignant_score_heatmap", lambda: plot_malignant_score_heatmap(result, value="status_score")),
        ("malignant_status_scatter", lambda: plot_malignant_status_scatter(result)),
        ("malignant_state_specificity", lambda: plot_malignant_state_specificity(result)),
        ("malignant_annotation_summary", lambda: plot_malignant_annotation_summary(result)),
        ("integration_summary", lambda: plot_integration_summary(result)),
        ("integration_raw_evidence", lambda: plot_integration_raw_evidence(result)),
    ]
    for name, fn in candidates:
        try:
            plots[name] = fn()
        except Exception:
            pass
    return plots
