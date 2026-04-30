from __future__ import annotations

from typing import Dict, Iterable, List, Optional, Set, Tuple

import numpy as np
import pandas as pd



def compute_bulk_expression(expr: pd.DataFrame) -> pd.Series:
    return expr.mean(axis=1).sort_values()



def assign_expression_bins(
    bulk_expression: pd.Series,
    n_bins: int = 24,
    collapse_low_expr_tail: bool = False,
    low_expr_threshold: Optional[float] = None,
) -> pd.Series:
    n_unique = bulk_expression.nunique()
    bins = min(n_bins, max(1, n_unique))
    if bins == 1:
        return pd.Series(0, index=bulk_expression.index)

    if collapse_low_expr_tail:
        threshold = float(low_expr_threshold) if low_expr_threshold is not None else float(bulk_expression.quantile(0.15))
        low_mask = bulk_expression <= threshold
        out = pd.Series(index=bulk_expression.index, dtype=int)
        out.loc[low_mask] = 0
        remaining = bulk_expression.loc[~low_mask]
        if remaining.empty:
            return out.fillna(0).astype(int)
        rem_bins = min(max(1, bins - 1), max(1, remaining.nunique()))
        if rem_bins == 1:
            out.loc[~low_mask] = 1
        else:
            ranked = remaining.rank(method="average")
            out.loc[~low_mask] = pd.qcut(ranked, q=rem_bins, labels=False, duplicates="drop").astype(int) + 1
        return out.astype(int)

    ranked = bulk_expression.rank(method="average")
    return pd.qcut(ranked, q=bins, labels=False, duplicates="drop").astype(int)



def build_gene_bin_lookup(gene_bins: pd.Series, eligible_genes: Optional[Set[str]] = None) -> Dict[int, List[str]]:
    lookup: Dict[int, List[str]] = {}
    eligible = None if eligible_genes is None else set(eligible_genes)
    for gene, bin_id in gene_bins.items():
        if eligible is not None and gene not in eligible:
            continue
        lookup.setdefault(int(bin_id), []).append(gene)
    return lookup



def sample_control_genes_for_markers(
    marker_genes: Iterable[str],
    gene_bins: pd.Series,
    bin_lookup: Dict[int, List[str]],
    control_size_fallbacks: List[int],
    rng: Optional[np.random.Generator] = None,
    excluded_genes: Optional[Set[str]] = None,
) -> Tuple[Dict[str, List[str]], Dict[str, object]]:
    rng = rng or np.random.default_rng(0)
    excluded_genes = set() if excluded_genes is None else set(excluded_genes)
    marker_genes = [g for g in marker_genes if g in gene_bins.index]
    fallback_sizes = [int(x) for x in control_size_fallbacks if int(x) > 0]
    if not fallback_sizes:
        raise ValueError("control_size_fallbacks must contain at least one positive integer")

    control_map: Dict[str, List[str]] = {}
    controls_per_marker: Dict[str, int] = {}
    fallback_used = False

    for gene in marker_genes:
        bin_id = int(gene_bins[gene])
        candidates = [g for g in bin_lookup.get(bin_id, []) if g != gene and g not in excluded_genes]
        sampled: List[str] = []
        for desired_n in fallback_sizes:
            if len(candidates) >= desired_n:
                sampled = rng.choice(candidates, size=desired_n, replace=False).tolist()
                if desired_n != fallback_sizes[0]:
                    fallback_used = True
                break
        if not sampled and candidates:
            desired_n = min(len(candidates), fallback_sizes[-1])
            if desired_n > 0:
                sampled = rng.choice(candidates, size=desired_n, replace=False).tolist()
                fallback_used = True
        control_map[gene] = sorted(sampled)
        controls_per_marker[gene] = len(sampled)

    counts = list(controls_per_marker.values())
    metadata = {
        "fallback_used": fallback_used,
        "controls_per_marker": controls_per_marker,
        "median_controls_per_marker": float(np.median(counts)) if counts else 0.0,
        "min_controls_per_marker": int(min(counts)) if counts else 0,
        "max_controls_per_marker": int(max(counts)) if counts else 0,
        "markers_with_controls": int(sum(c > 0 for c in counts)),
        "markers_without_controls": int(sum(c == 0 for c in counts)),
        "eligible_candidate_genes": int(sum(len(v) for v in bin_lookup.values())),
    }
    return control_map, metadata
