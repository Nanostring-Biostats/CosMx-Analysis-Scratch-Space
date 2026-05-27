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

