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


def _as_float_series(df: pd.DataFrame, column: str, default: float = 0.0) -> pd.Series:
    if column not in df.columns:
        return pd.Series([default] * len(df), index=df.index, dtype="float64")
    return pd.to_numeric(df[column], errors="coerce").fillna(default)



def _normalize_confidence_filter(ignore_confidence_levels) -> set[str]:
    if ignore_confidence_levels is None:
        return set()
    if isinstance(ignore_confidence_levels, str):
        ignore_confidence_levels = [ignore_confidence_levels]
    return {str(x).strip().lower() for x in ignore_confidence_levels if str(x).strip()}


def _branch_supported_raw_evidence(ann: pd.DataFrame) -> pd.Series:
    """Return branch-supported raw evidence for the final normal call."""
    if "annot_branch_supported_raw_score" not in ann.columns:
        return pd.Series([pd.NA] * len(ann), index=ann.index, dtype="float64")
    return pd.to_numeric(ann["annot_branch_supported_raw_score"], errors="coerce")


def _valid_normal_label_mask(
    ann: pd.DataFrame,
    *,
    unknown_branch_raw_threshold: float | None = 0.10,
    unknown_min_score: float | None = None,
    ignore_confidence_levels=("low", "none"),
    unknown_on_low_evidence: bool = True,
) -> pd.Series:
    """Export-aligned normal-label validity mask for tissue detection."""
    if len(ann) == 0:
        return pd.Series([], index=ann.index, dtype=bool)

    mask = pd.Series(True, index=ann.index)
    label = ann["annot_label"].astype(str).str.strip().str.lower() if "annot_label" in ann.columns else pd.Series("", index=ann.index)
    status = ann["annot_status"].astype(str).str.strip().str.lower() if "annot_status" in ann.columns else pd.Series("", index=ann.index)
    invalid_label = label.isin(["", "nan", "none", "unknown", "unresolved"])
    invalid_status = status.isin(["unknown", "unresolved"])
    mask = mask & ~invalid_label & ~invalid_status

    ignored = _normalize_confidence_filter(ignore_confidence_levels)
    if ignored and "annot_confidence" in ann.columns:
        conf = ann["annot_confidence"].astype(str).str.strip().str.lower()
        mask = mask & ~conf.isin(ignored)

    if unknown_on_low_evidence and "annot_low_evidence" in ann.columns:
        mask = mask & ~ann["annot_low_evidence"].fillna(False).astype(bool)

    if unknown_branch_raw_threshold is not None:
        branch_raw = _branch_supported_raw_evidence(ann)
        mask = mask & branch_raw.ge(float(unknown_branch_raw_threshold)).fillna(False)

    if unknown_min_score is not None:
        score = _as_float_series(ann, "annot_score")
        mask = mask & score.ge(float(unknown_min_score)).fillna(False)

    return mask.fillna(False).astype(bool)


def _resolve_tissue_origin_hierarchy(hierarchy="tissue_origin_screen"):
    """Resolve a built-in hierarchy name or custom hierarchy for tissue detection."""
    from .builtin import get_builtin_hierarchy
    from .datamodels import MarkerProgram

    if isinstance(hierarchy, str):
        return get_builtin_hierarchy(hierarchy), hierarchy
    if isinstance(hierarchy, MarkerProgram):
        return [hierarchy], None
    if isinstance(hierarchy, (list, tuple)) and all(isinstance(node, MarkerProgram) for node in hierarchy):
        return list(hierarchy), None
    raise TypeError(
        "hierarchy must be a built-in hierarchy name, a MarkerProgram, "
        "or a list/tuple of MarkerProgram roots."
    )




def _compile_tissue_origin_hierarchy_for_panel(root_programs, gene_list, *, hierarchy_name=None):
    """Compile a tissue-origin screen hierarchy against the input panel.

    This mirrors the built-in hierarchy panel compiler but also works for
    custom tissue-origin screens. It additionally prunes sparse tissue-origin
    anchor leaves more strictly than ordinary annotation leaves, because tissue
    detection should be driven by reasonably covered lineage anchors rather than
    partial marker overlap.
    """
    from .builtin.registry import _compile_program_for_panel, _panel_tier_from_gene_count

    gene_set = {str(g).upper() for g in list(gene_list) if pd.notna(g) and str(g).strip() != ""}
    panel_tier = _panel_tier_from_gene_count(len(gene_set))
    report_rows = []
    compiled = []
    for root in root_programs:
        out = _compile_program_for_panel(
            root,
            gene_set=gene_set,
            panel_tier=panel_tier,
            preserve_major_immune_subtypes=True,
            report_rows=report_rows,
            level=1,
        )
        if out is not None:
            compiled.append(out)
    report = pd.DataFrame(report_rows)
    if not report.empty:
        report["builtin_name"] = hierarchy_name if hierarchy_name is not None else "custom_tissue_origin_screen"
        report["panel_name"] = None
        report["panel_gene_count"] = len(gene_set)
        report["tissue_origin_anchor_usable"] = True

        by_name = report.set_index("name", drop=False)
        min_anchor_markers = 3
        min_anchor_fraction = 0.40
        sparse_anchor_names = set()
        for _, row in report.iterrows():
            name = str(row.get("name"))
            present = float(row.get("positive_markers_present", 0) or 0)
            frac = float(row.get("positive_markers_present_fraction", 0) or 0)
            # The metadata role is not part of the panel report in older rows,
            # so the recursive pruning function confirms anchor metadata before
            # removing a node.
            if present < min_anchor_markers or frac < min_anchor_fraction:
                sparse_anchor_names.add(name)
        if sparse_anchor_names:
            report.loc[report["name"].isin(sparse_anchor_names), "tissue_origin_anchor_usable"] = False

        def prune(nodes):
            kept = []
            for node in nodes:
                node.children = prune(list(getattr(node, "children", []) or []))
                meta = getattr(node, "metadata", {}) or {}
                role = str(meta.get("tissue_detection_role", "")).strip()
                if role == "anchor" and node.name in sparse_anchor_names:
                    continue
                kept.append(node)
            return kept
        compiled = prune(compiled)
    return compiled, report

def _tissue_origin_anchor_metadata(hierarchy) -> pd.DataFrame:
    """Return and validate tissue-detection anchor metadata for a hierarchy."""
    from .inspect import describe_hierarchy

    desc = describe_hierarchy(hierarchy)
    rows = []
    for _, row in desc.iterrows():
        meta = row.get("metadata", {})
        if not isinstance(meta, dict):
            meta = {}
        role = str(meta.get("tissue_detection_role", "")).strip()
        if role not in {"anchor", "shared_control"}:
            continue
        tissue_candidate = str(meta.get("tissue_candidate", "")).strip()
        recommended_hierarchy = str(meta.get("recommended_hierarchy", "")).strip()
        if role == "anchor" and not tissue_candidate:
            raise ValueError(
                f"Tissue-detection anchor {row.get('name')!r} is missing "
                "metadata['tissue_candidate']."
            )
        rows.append({
            "label": str(row.get("name")),
            "tissue_detection_role": role,
            "tissue_candidate": tissue_candidate or "generic_tme",
            "recommended_hierarchy": recommended_hierarchy or "tme_core",
            "contributes_to_tissue_detection": bool(meta.get("contributes_to_tissue_detection", role == "anchor")),
            "description": row.get("description"),
        })
    out = pd.DataFrame(rows)
    if out.empty or not (out["tissue_detection_role"] == "anchor").any():
        raise ValueError(
            "detect_tissue_type() requires a tissue-origin screen hierarchy with "
            "at least one node whose metadata includes "
            "tissue_detection_role='anchor', tissue_candidate, and "
            "recommended_hierarchy. Use the built-in 'tissue_origin_screen' "
            "or add this metadata to custom anchor nodes."
        )
    return out


def _cluster_weight_series(ann: pd.DataFrame, cluster_weights=None) -> pd.Series:
    """Align optional cluster weights to result cluster IDs.

    Missing weights are filled with 1.0, so users can pass a Series produced by
    ``adata.obs[cluster_key].value_counts()`` without manually sorting it to the
    cluster-mean matrix.  Non-positive and non-numeric weights are coerced to 0.
    """
    clusters = ann["cluster_id"].astype(str) if "cluster_id" in ann.columns else ann.index.astype(str).to_series(index=ann.index)
    if cluster_weights is None:
        return pd.Series(1.0, index=ann.index, dtype="float64")
    if isinstance(cluster_weights, pd.DataFrame):
        if cluster_weights.shape[1] != 1:
            raise ValueError("cluster_weights DataFrame must have exactly one column")
        weights = cluster_weights.iloc[:, 0]
    elif isinstance(cluster_weights, pd.Series):
        weights = cluster_weights
    else:
        weights = pd.Series(cluster_weights)
    weights.index = weights.index.astype(str)
    aligned = weights.reindex(clusters.values)
    aligned.index = ann.index
    aligned = pd.to_numeric(aligned, errors="coerce").fillna(1.0).astype(float)
    aligned = aligned.clip(lower=0.0)
    return aligned


def _attach_tissue_detection_call(
    out: pd.DataFrame,
    *,
    min_tissue_detection_score: float = 0.05,
    min_tissue_anchor_cluster_count: int = 1,
    min_tissue_anchor_cluster_fraction: float = 0.0,
    ambiguity_score_ratio: float = 0.90,
    generic_tme_min_valid_fraction: float = 0.25,
) -> pd.DataFrame:
    if out.empty:
        return out
    out = out.copy()
    specific = out[~out["is_shared_control"].astype(bool)].copy()
    eligible = (
        pd.to_numeric(specific["tissue_detection_score"], errors="coerce").ge(float(min_tissue_detection_score)).fillna(False)
        & pd.to_numeric(specific["n_valid_anchor_clusters"], errors="coerce").ge(int(min_tissue_anchor_cluster_count)).fillna(False)
        & pd.to_numeric(specific["fraction_valid_anchor_clusters"], errors="coerce").ge(float(min_tissue_anchor_cluster_fraction)).fillna(False)
    )

    call = "insufficient_tissue_specific_evidence"
    detected = "insufficient_tissue_specific_evidence"
    recommended = pd.NA
    reason = "no_tissue_anchor_has_enough_valid_evidence"
    call_conf = "none"

    eligible_rows = specific.loc[eligible].sort_values("tissue_detection_score", ascending=False)
    if not eligible_rows.empty:
        best = eligible_rows.iloc[0]
        call = "specific_tissue"
        detected = best["tissue_candidate"]
        recommended = best["candidate_recommended_hierarchy"]
        call_conf = str(best.get("tissue_detection_confidence", "moderate"))
        reason = "best_tissue_anchor_has_valid_direct_evidence"
        if len(eligible_rows) > 1:
            second = eligible_rows.iloc[1]
            best_score = float(best["tissue_detection_score"])
            second_score = float(second["tissue_detection_score"])
            if best_score > 0 and second_score >= float(ambiguity_score_ratio) * best_score:
                call = "ambiguous"
                detected = "ambiguous"
                call_conf = "ambiguous"
                reason = f"top_two_tissue_anchors_are_close:{best['tissue_candidate']}:{second['tissue_candidate']}"
    else:
        shared = out[out["is_shared_control"].astype(bool)]
        if not shared.empty:
            shared_valid_fraction = float(pd.to_numeric(shared["fraction_valid_anchor_clusters"], errors="coerce").max() or 0.0)
            shared_score = float(pd.to_numeric(shared["tissue_detection_score"], errors="coerce").max() or 0.0)
            if shared_valid_fraction >= float(generic_tme_min_valid_fraction) or shared_score > 0:
                call = "generic_tme"
                detected = "generic_tme"
                recommended = "tme_core"
                call_conf = "generic"
                reason = "shared_tme_controls_have_valid_evidence_but_no_specific_tissue_anchor_passes"

    out["tissue_detection_call"] = call
    out["detected_tissue_type"] = detected
    out["recommended_hierarchy"] = recommended
    out["tissue_detection_call_confidence"] = call_conf
    out["tissue_detection_call_reason"] = reason
    return out


def detect_tissue_type(
    cluster_expr,
    *,
    hierarchy="tissue_origin_screen",
    compile_for_panel: bool = True,
    cluster_weights=None,
    pipeline_kwargs: dict | None = None,
    unknown_branch_raw_threshold: float | None = 0.10,
    unknown_min_score: float | None = None,
    ignore_confidence_levels=("low", "none"),
    unknown_on_low_evidence: bool = True,
    min_tissue_detection_score: float = 0.05,
    min_tissue_anchor_cluster_count: int = 1,
    min_tissue_anchor_cluster_fraction: float = 0.0,
    ambiguity_score_ratio: float = 0.90,
    generic_tme_min_valid_fraction: float = 0.25,
    return_result: bool = False,
):
    """Detect likely tissue type with a tissue-origin screen hierarchy.

    By default this uses the built-in ``"tissue_origin_screen"`` hierarchy, a
    shallow standard :class:`MarkerProgram` hierarchy that contains broad
    tissue-defining anchors in one shared competition space.  Users may pass a
    custom hierarchy with the same tissue-detection metadata on anchor nodes.
    The screen is intended to choose a downstream tissue-specific hierarchy,
    not to replace final cell-type annotation.

    Parameters
    ----------
    cluster_expr
        Cluster-level expression matrix accepted by :class:`HierAnnotPipeline`.
    hierarchy
        Built-in hierarchy name or custom :class:`MarkerProgram` hierarchy used
        as the tissue-origin screen. Anchor nodes must include metadata with
        ``tissue_detection_role="anchor"``, ``tissue_candidate``, and
        ``recommended_hierarchy``. Shared-control nodes may use
        ``tissue_detection_role="shared_control"``. Defaults to the built-in
        ``"tissue_origin_screen"``.
    compile_for_panel
        If True, compile the tissue-origin screen hierarchy against
        ``cluster_expr.index`` before scoring. This is the default because
        tissue-origin detection is sensitive to marker availability, especially
        for targeted RNA panels.
    cluster_weights
        Optional cluster-size weights.  A common choice is
        ``adata.obs[cluster_key].value_counts()``.  Weights are aligned to the
        cluster IDs in the fitted result; missing clusters default to weight 1.
    pipeline_kwargs
        Optional keyword arguments forwarded to :class:`HierAnnotPipeline`.
        Tissue detection uses the normal hierarchy only by default and sets
        ``malignant_integration_mode="off"`` unless the user overrides it.
    unknown_branch_raw_threshold, unknown_min_score, ignore_confidence_levels,
    unknown_on_low_evidence
        Export-aligned controls for deciding whether a final normal label is
        valid enough to contribute tissue-origin evidence.
    min_tissue_detection_score, min_tissue_anchor_cluster_count,
    min_tissue_anchor_cluster_fraction
        Dataset-level evidence requirements before a specific tissue call is
        made.
    ambiguity_score_ratio
        If the second-best tissue score is at least this fraction of the best
        score, the dataset-level call is marked ``"ambiguous"``.
    generic_tme_min_valid_fraction
        Minimum fraction of valid shared-control clusters that can support a
        ``"generic_tme"`` call when no specific tissue anchor passes.
    return_result
        If True, return ``(summary, result)`` where ``result`` is the fitted
        HierAnnot result for ``tissue_origin_screen``.
    """
    from .pipeline import HierAnnotPipeline

    original_root_programs, hierarchy_name = _resolve_tissue_origin_hierarchy(hierarchy)
    root_programs = original_root_programs
    compile_report = None
    if compile_for_panel:
        compiled_roots, compile_report = _compile_tissue_origin_hierarchy_for_panel(
            original_root_programs,
            gene_list=getattr(cluster_expr, "index", []),
            hierarchy_name=hierarchy_name,
        )
        # Use the compiled screen only when it retains at least one valid
        # tissue-detection anchor. Extremely sparse panels may prune all
        # anchors; in that case fall back to the original screen so the result
        # can still report insufficient/generic evidence rather than failing.
        try:
            _tissue_origin_anchor_metadata(compiled_roots)
            root_programs = compiled_roots
        except ValueError:
            root_programs = original_root_programs
    anchors = _tissue_origin_anchor_metadata(root_programs)

    pipeline_kwargs = dict(pipeline_kwargs or {})
    pipeline_kwargs.setdefault("malignant_integration_mode", "off")
    pipeline_kwargs.setdefault("malignant_programs", None)

    pipeline = HierAnnotPipeline(
        root_programs=root_programs,
        **pipeline_kwargs,
    )
    result = pipeline.fit_score(cluster_expr)
    result.metadata = dict(getattr(result, "metadata", None) or {})
    if hierarchy_name is not None:
        result.metadata.setdefault("tissue_origin_screen_hierarchy", hierarchy_name)
    result.metadata.setdefault("tissue_origin_screen_compile_for_panel", bool(compile_for_panel))
    if compile_report is not None:
        result.metadata.setdefault("tissue_origin_screen_panel_gene_count", int(compile_report["panel_gene_count"].iloc[0]) if "panel_gene_count" in compile_report.columns and len(compile_report) else None)
        if "status" in compile_report.columns:
            result.metadata.setdefault("tissue_origin_screen_compile_status_counts", compile_report["status"].value_counts().to_dict())
    ann = getattr(result, "cluster_annotations", None)
    if ann is None or not isinstance(ann, pd.DataFrame):
        raise ValueError("Tissue-origin screen result does not contain cluster_annotations")
    ann = ann.copy()
    n = len(ann)

    branch_raw = _branch_supported_raw_evidence(ann)
    score = _as_float_series(ann, "annot_score")
    margin = _as_float_series(ann, "annot_margin")
    valid_mask = _valid_normal_label_mask(
        ann,
        unknown_branch_raw_threshold=unknown_branch_raw_threshold,
        unknown_min_score=unknown_min_score,
        ignore_confidence_levels=ignore_confidence_levels,
        unknown_on_low_evidence=unknown_on_low_evidence,
    )
    weights = _cluster_weight_series(ann, cluster_weights=cluster_weights)
    total_weight = float(weights.sum()) if float(weights.sum()) > 0 else float(max(n, 1))

    label = ann["annot_label"].astype(str) if "annot_label" in ann.columns else pd.Series([""] * n, index=ann.index)
    rows = []
    for _, anchor in anchors.iterrows():
        anchor_label = str(anchor["label"])
        label_mask = label.eq(anchor_label)
        anchor_mask = label_mask & valid_mask
        weak_mask = label_mask & ~valid_mask & branch_raw.gt(0).fillna(False)
        w = weights[anchor_mask]
        evidence = (branch_raw.clip(lower=0.0) + 0.20 * score.clip(lower=0.0)) * weights
        total_evidence = float(evidence[anchor_mask].sum()) if bool(anchor_mask.any()) else 0.0
        weighted_score = total_evidence / max(total_weight, 1.0)
        n_valid = int(anchor_mask.sum())
        weight_valid = float(w.sum()) if bool(anchor_mask.any()) else 0.0
        fraction_valid = float(n_valid / max(n, 1))
        weight_fraction = float(weight_valid / max(total_weight, 1.0))
        mean_raw = float(branch_raw[anchor_mask].mean()) if bool(anchor_mask.any()) else 0.0
        mean_score = float(score[anchor_mask].mean()) if bool(anchor_mask.any()) else 0.0
        median_margin = float(margin[anchor_mask].median()) if bool(anchor_mask.any()) else 0.0
        detection_score = weighted_score + 0.20 * weight_fraction * max(mean_raw, 0.0) + 0.05 * weight_fraction
        confidence = "no_tissue_specific_evidence"
        if n_valid >= 2 and mean_raw >= max(float(unknown_branch_raw_threshold or 0.0), 0.20):
            confidence = "strong"
        elif n_valid >= 1:
            confidence = "moderate"
        elif int(weak_mask.sum()) > 0:
            confidence = "weak"
        rows.append({
            "tissue_candidate": str(anchor["tissue_candidate"]),
            "anchor_label": anchor_label,
            "is_shared_control": str(anchor["tissue_detection_role"]) == "shared_control",
            "candidate_recommended_hierarchy": str(anchor["recommended_hierarchy"]),
            "n_clusters": int(n),
            "total_cluster_weight": float(total_weight),
            "n_valid_label_clusters": int(valid_mask.sum()),
            "fraction_valid_label_clusters": float(valid_mask.mean()) if n else 0.0,
            "n_valid_anchor_clusters": n_valid,
            "fraction_valid_anchor_clusters": fraction_valid,
            "valid_anchor_cluster_weight": weight_valid,
            "fraction_valid_anchor_cluster_weight": weight_fraction,
            "n_weak_anchor_clusters": int(weak_mask.sum()),
            "mean_anchor_branch_supported_raw_score": mean_raw,
            "mean_anchor_score": mean_score,
            "median_anchor_margin": median_margin,
            "total_anchor_evidence_score": total_evidence,
            "weighted_anchor_evidence_score": weighted_score,
            "tissue_evidence_score": weighted_score,
            "tissue_detection_score": float(detection_score),
            "tissue_detection_confidence": confidence,
            "anchor_description": anchor.get("description"),
        })

    out = pd.DataFrame(rows)
    if not out.empty:
        # Combine multiple anchors for the same tissue, for example skin has both
        # keratinocyte and melanocyte anchors.  Keep anchor labels for audit.
        group_cols = ["tissue_candidate", "is_shared_control", "candidate_recommended_hierarchy"]
        agg = out.groupby(group_cols, dropna=False).agg(
            n_clusters=("n_clusters", "max"),
            total_cluster_weight=("total_cluster_weight", "max"),
            n_valid_label_clusters=("n_valid_label_clusters", "max"),
            fraction_valid_label_clusters=("fraction_valid_label_clusters", "max"),
            n_valid_anchor_clusters=("n_valid_anchor_clusters", "sum"),
            fraction_valid_anchor_clusters=("fraction_valid_anchor_clusters", "sum"),
            valid_anchor_cluster_weight=("valid_anchor_cluster_weight", "sum"),
            fraction_valid_anchor_cluster_weight=("fraction_valid_anchor_cluster_weight", "sum"),
            n_weak_anchor_clusters=("n_weak_anchor_clusters", "sum"),
            mean_anchor_branch_supported_raw_score=("mean_anchor_branch_supported_raw_score", "max"),
            mean_anchor_score=("mean_anchor_score", "max"),
            median_anchor_margin=("median_anchor_margin", "max"),
            total_anchor_evidence_score=("total_anchor_evidence_score", "sum"),
            weighted_anchor_evidence_score=("weighted_anchor_evidence_score", "sum"),
            tissue_evidence_score=("tissue_evidence_score", "sum"),
            tissue_detection_score=("tissue_detection_score", "sum"),
            anchor_label=("anchor_label", lambda x: ";".join([str(v) for v in x if str(v)])),
            tissue_detection_confidence=("tissue_detection_confidence", lambda x: "strong" if "strong" in set(x) else ("moderate" if "moderate" in set(x) else ("weak" if "weak" in set(x) else "no_tissue_specific_evidence"))),
        ).reset_index()
        out = agg.sort_values("tissue_detection_score", ascending=False, na_position="last").reset_index(drop=True)
        out["tissue_detection_rank"] = range(1, len(out) + 1)
        out = _attach_tissue_detection_call(
            out,
            min_tissue_detection_score=min_tissue_detection_score,
            min_tissue_anchor_cluster_count=min_tissue_anchor_cluster_count,
            min_tissue_anchor_cluster_fraction=min_tissue_anchor_cluster_fraction,
            ambiguity_score_ratio=ambiguity_score_ratio,
            generic_tme_min_valid_fraction=generic_tme_min_valid_fraction,
        )
        # Present the detected specific tissue first when a specific/ambiguous
        # tissue call is made. Shared TME controls are useful diagnostics but
        # should not occupy the top row when users consume ``iloc[0]`` to choose
        # a downstream hierarchy. For generic calls, shared controls remain first.
        call = str(out["tissue_detection_call"].iloc[0]) if "tissue_detection_call" in out.columns and len(out) else ""
        detected = str(out["detected_tissue_type"].iloc[0]) if "detected_tissue_type" in out.columns and len(out) else ""
        if call in {"specific_tissue", "ambiguous"}:
            priority = out["tissue_candidate"].astype(str).eq(detected).map({True: 0, False: 1})
            shared_priority = out["is_shared_control"].astype(bool).map({False: 0, True: 1})
            out = (
                out.assign(_display_priority=priority, _shared_display_priority=shared_priority)
                .sort_values(["_display_priority", "_shared_display_priority", "tissue_detection_score"], ascending=[True, True, False])
                .drop(columns=["_display_priority", "_shared_display_priority"])
                .reset_index(drop=True)
            )
            out["tissue_detection_rank"] = range(1, len(out) + 1)
    if return_result:
        return out, result
    return out
