from __future__ import annotations

import pandas as pd

def _sanitize_label_text(value: object) -> str:
    text = str(value).strip()
    for old, new in [(" ", "_"), ("/", "_"), ("-", "_")]:
        text = text.replace(old, new)
    while "__" in text:
        text = text.replace("__", "_")
    return text.strip("_")


def summarize_hierarchy_fit(result_or_summary):
    """
    Return a normalized hierarchy-fit summary DataFrame.

    Parameters
    ----------
    result_or_summary
        Either a HierAnnotResult-like object with ``diagnostics_summary`` or a
        precomputed fit-summary DataFrame with ``metric`` and ``value`` columns.
    """

    if isinstance(result_or_summary, pd.DataFrame):
        df = result_or_summary.copy()
    else:
        df = getattr(result_or_summary, "diagnostics_summary", None)
        if df is None:
            raise ValueError("result does not contain diagnostics_summary")
        df = df.copy()

    if "metric" not in df.columns or "value" not in df.columns:
        raise ValueError("Expected a fit-summary DataFrame with 'metric' and 'value' columns")
    return df


def compare_hierarchy_fit(results, names=None):
    """
    Compare hierarchy-fit summaries across multiple results and rank them.

    Parameters
    ----------
    results
        Either:
        - a sequence of HierAnnotResult-like objects or fit-summary DataFrames, or
        - a dictionary mapping candidate names to results / fit summaries.
    names
        Optional names for each candidate when `results` is a sequence.

    Returns
    -------
    pandas.DataFrame
        One row per candidate with raw fit metrics, ``fit_rank_score``, and
        ``fit_rank`` (1 = best).
    """
    metric_order = [
        "hierarchy_fit_status",
        "fit_rank_score",
        "median_branch_supported_raw_score",
        "median_final_score",
        "fraction_l1_positive_branch_support",
        "median_marker_detection_support_fraction",
        "fraction_unknown",
        "fraction_low_absolute_support",
        "fraction_backoff_applied",
        "median_branch_support_score",
    ]

    if isinstance(results, dict):
        items = list(results.items())
    else:
        if names is None:
            names = [f"candidate_{i+1}" for i in range(len(results))]
        items = list(zip(names, results))

    rows = []
    for name, res in items:
        df = summarize_hierarchy_fit(res)
        row = {"name": name}
        for metric in metric_order:
            sub = df[df["metric"] == metric]
            row[metric] = sub["value"].iloc[0] if not sub.empty else pd.NA
        rows.append(row)

    out = pd.DataFrame(rows)
    numeric_cols = [
        "fit_rank_score",
        "median_branch_supported_raw_score",
        "median_final_score",
        "fraction_l1_positive_branch_support",
        "median_marker_detection_support_fraction",
        "fraction_unknown",
        "fraction_low_absolute_support",
        "fraction_backoff_applied",
        "median_branch_support_score",
    ]
    for col in numeric_cols:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")

    if "fit_rank_score" in out.columns:
        out = out.sort_values("fit_rank_score", ascending=False, na_position="last").reset_index(drop=True)
        out["fit_rank"] = range(1, len(out) + 1)
    else:
        out["fit_rank"] = pd.NA
    return out
