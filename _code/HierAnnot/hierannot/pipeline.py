from __future__ import annotations

from dataclasses import asdict
from typing import Dict, List, Optional, Sequence, Set

import numpy as np
import pandas as pd

from .controls import assign_expression_bins, build_gene_bin_lookup, compute_bulk_expression, sample_control_genes_for_markers
from .datamodels import (
    CompiledMarkerProgram,
    ControlGeneExclusionPolicy,
    HierAnnotResult,
    MarkerCleanupStrategy,
    MarkerProgram,
    PresetName,
)
from .hierarchy import compute_cluster_branch_diagnostics, decide_cluster_annotation_from_scores
from .preprocessing import preprocess_expression
from .scoring import flatten_compiled_programs, score_programs_across_clusters
from .workflows import _attach_final_call_competition


PRESET_DEFAULTS: Dict[str, Dict[str, object]] = {
    "small_panel": {
        "n_bins": 12,
        "control_size_fallbacks": [25, 10, 5],
        "control_gene_exclusion_policy": "current_branch",
        "min_bulk_expr_for_control": None,
        "collapse_low_expr_tail": False,
        "detection_floor_quantile": 0.10,
        "delta_clip_quantile": 0.95,
    },
    "medium_panel": {
        "n_bins": 24,
        "control_size_fallbacks": [50, 25, 10, 5],
        "control_gene_exclusion_policy": "all_pipeline_markers",
        "min_bulk_expr_for_control": None,
        "collapse_low_expr_tail": True,
        "detection_floor_quantile": 0.10,
        "delta_clip_quantile": 0.95,
    },
    "large_panel": {
        "n_bins": 32,
        "control_size_fallbacks": [100, 50, 25, 10],
        "control_gene_exclusion_policy": "all_pipeline_markers",
        "min_bulk_expr_for_control": None,
        "collapse_low_expr_tail": True,
        "detection_floor_quantile": 0.15,
        "delta_clip_quantile": 0.9,
    },
}


class HierAnnotPipeline:
    """Hierarchical marker-scoring pipeline for cluster-level cell type annotation."""

    def __init__(
        self,
        root_programs: List[MarkerProgram],
        input_type: str = "raw_cluster_means",
        scale_factor: float = 1e4,
        preset: PresetName = "auto",
        n_bins: Optional[int] = None,
        ctrl_size: int = 25,
        control_size_fallbacks: Optional[List[int]] = None,
        min_bulk_expr_for_control: Optional[float] = None,
        control_gene_exclusion_policy: Optional[ControlGeneExclusionPolicy] = None,
        min_markers_present: int = 2,
        min_marker_fraction: float = 0.3,
        score_threshold: float = 0.05,
        margin_threshold: float = 0.02,
        branch_child_rescue_weight: float = 0.35,
        level1_branch_child_rescue_weight: float = 0.40,
        negative_weight: float = 0.5,
        uppercase_genes: bool = True,
        exclude_all_marker_genes_from_controls: Optional[bool] = None,
        random_state: Optional[int] = 0,
        child_marker_strategy: MarkerCleanupStrategy = "subtract_direct_parent",
        root_marker_strategy: MarkerCleanupStrategy = "none",
        fallback_to_canonical_if_sparse: bool = True,
        collapse_low_expr_tail: Optional[bool] = None,
        low_expr_threshold: Optional[float] = None,
        detection_floor: Optional[float] = None,
        delta_clip_quantile: Optional[float] = None,
        branch_routing_raw_weight: float = 0.60,
        weak_raw_score_threshold: float = 0.10,
        min_branch_supported_raw_score: float = 0.0,
        min_leaf_raw_score: float = -0.2,
    ) -> None:
        """
        Create a HierAnnot scoring pipeline.

        This is the main public API for fitting a compiled hierarchy to a
        cluster-by-gene expression matrix, scoring all hierarchy nodes, and
        producing final hierarchical annotations plus diagnostics.

        Main public knobs
        -----------------
        score_threshold
            Minimum branch-supported decision score required to continue routing.

        margin_threshold
            Minimum score gap required to treat the best competing label as
            clearly separated from the runner-up.

        branch_child_rescue_weight
            Descendant-rescue weight used below level 1. Higher values allow
            strong child evidence to support weak intermediate nodes more strongly.

        level1_branch_child_rescue_weight
            Descendant-rescue weight used only for top-level routing. This is
            typically set slightly higher because broad level-1 compartments can
            be weaker and more detection-sensitive.

        branch_routing_raw_weight
            Weight of raw marker enrichment relative to sibling-specific decision
            score when computing local node evidence.

        weak_raw_score_threshold
            Raw-score gate threshold. Nodes below this value receive reduced
            contribution from the decision-score term.

        min_branch_supported_raw_score
            Minimum branch-supported raw score required by the decision engine to
            keep a node without backing off. This affects routing/backoff, not
            just export formatting.

        min_leaf_raw_score
            Minimum local raw score required for the selected node itself.
            This is a node-level gate and is distinct from branch-supported rescue.
            A mildly negative value (default ``-0.2``) allows traversal through
            weak intermediate or level-1 nodes when descendant-supported branch
            evidence is strong, while still blocking strongly contradictory nodes.
        """
        self.root_programs = root_programs
        self.input_type = input_type
        self.scale_factor = scale_factor
        self.preset = preset
        self.n_bins = n_bins
        self.ctrl_size = ctrl_size
        self.control_size_fallbacks = control_size_fallbacks
        self.min_bulk_expr_for_control = min_bulk_expr_for_control
        self.control_gene_exclusion_policy = control_gene_exclusion_policy
        self.min_markers_present = min_markers_present
        self.min_marker_fraction = min_marker_fraction
        self.score_threshold = score_threshold
        self.margin_threshold = margin_threshold
        self.branch_child_rescue_weight = float(branch_child_rescue_weight)
        if not (0.0 <= self.branch_child_rescue_weight <= 1.0):
            raise ValueError("branch_child_rescue_weight must be in [0, 1].")
        self.level1_branch_child_rescue_weight = float(level1_branch_child_rescue_weight)
        if not (0.0 <= self.level1_branch_child_rescue_weight <= 1.0):
            raise ValueError("level1_branch_child_rescue_weight must be in [0, 1].")
        self.branch_support_descendant_weight = self.branch_child_rescue_weight
        self.branch_support_parent_weight = 1.0 - self.branch_support_descendant_weight
        self.level1_branch_support_descendant_weight = self.level1_branch_child_rescue_weight
        self.level1_branch_support_parent_weight = 1.0 - self.level1_branch_support_descendant_weight
        self.negative_weight = negative_weight
        self.uppercase_genes = uppercase_genes
        self.exclude_all_marker_genes_from_controls = exclude_all_marker_genes_from_controls
        self.random_state = random_state
        self.child_marker_strategy = child_marker_strategy
        self.root_marker_strategy = root_marker_strategy
        self.fallback_to_canonical_if_sparse = fallback_to_canonical_if_sparse
        self.collapse_low_expr_tail = collapse_low_expr_tail
        self.low_expr_threshold = low_expr_threshold
        self.detection_floor = detection_floor
        self.delta_clip_quantile = delta_clip_quantile
        self.branch_child_rescue_weight = float(branch_child_rescue_weight)
        self.branch_support_descendant_weight = float(self.branch_child_rescue_weight)
        self.branch_support_parent_weight = float(1.0 - self.branch_child_rescue_weight)
        self.branch_routing_raw_weight = float(branch_routing_raw_weight)
        self.weak_raw_score_threshold = float(weak_raw_score_threshold)
        self.min_branch_supported_raw_score = float(min_branch_supported_raw_score)
        self.min_leaf_raw_score = float(min_leaf_raw_score)
        self.compiled_programs = self._compile_programs()
        self.compilation_report = self._build_compilation_report()
        self._resolved_config: Optional[Dict[str, object]] = None

    def _normalize_markers(self, markers: Sequence[str]) -> List[str]:
        cleaned = [str(g).strip() for g in markers if str(g).strip()]
        if self.uppercase_genes:
            cleaned = [g.upper() for g in cleaned]
        return list(dict.fromkeys(cleaned))

    def _apply_cleanup_strategy(self, canonical_markers: List[str], parent_markers: Optional[Set[str]], ancestor_markers: Set[str], sibling_markers: Optional[Set[str]], strategy: MarkerCleanupStrategy) -> List[str]:
        marker_set = set(canonical_markers)
        if strategy == "none":
            effective = canonical_markers
        elif strategy == "subtract_direct_parent":
            effective = [g for g in canonical_markers if g not in (parent_markers or set())]
        elif strategy == "subtract_ancestors":
            effective = [g for g in canonical_markers if g not in ancestor_markers]
        elif strategy == "sibling_unique_only":
            effective = [g for g in canonical_markers if g not in (sibling_markers or set())]
        else:
            raise ValueError(f"Unsupported marker cleanup strategy: {strategy}")
        if len(effective) == 0 and len(marker_set) > 0:
            return canonical_markers
        return list(dict.fromkeys(effective))

    def _compile_children(self, children: List[MarkerProgram], parent_name: str, parent_markers: Set[str], ancestor_markers: Set[str], level: int) -> List[CompiledMarkerProgram]:
        compiled_children: List[CompiledMarkerProgram] = []
        sibling_canonical = {child.name: set(self._normalize_markers(child.positive_markers)) for child in children}
        for child in children:
            canonical = self._normalize_markers(child.positive_markers)
            sibling_other = set().union(*[markers for name, markers in sibling_canonical.items() if name != child.name])
            effective = self._apply_cleanup_strategy(canonical, parent_markers, ancestor_markers, sibling_other, self.child_marker_strategy)
            compiled = CompiledMarkerProgram(
                name=child.name,
                canonical_markers=canonical,
                effective_markers=effective,
                negative_markers=self._normalize_markers(child.negative_markers),
                description=child.description,
                aliases=list(child.aliases),
                min_markers_present=child.min_markers_present,
                metadata=dict(child.metadata),
                level=level,
                parent_name=parent_name,
                cleanup_strategy=self.child_marker_strategy,
            )
            compiled.children = self._compile_children(child.children, child.name, set(canonical), ancestor_markers | set(canonical), level + 1)
            compiled_children.append(compiled)
        return compiled_children

    def _compile_programs(self) -> List[CompiledMarkerProgram]:
        compiled_roots: List[CompiledMarkerProgram] = []
        for root in self.root_programs:
            canonical = self._normalize_markers(root.positive_markers)
            effective = self._apply_cleanup_strategy(canonical, None, set(), None, self.root_marker_strategy)
            compiled = CompiledMarkerProgram(
                name=root.name,
                canonical_markers=canonical,
                effective_markers=effective,
                negative_markers=self._normalize_markers(root.negative_markers),
                description=root.description,
                aliases=list(root.aliases),
                min_markers_present=root.min_markers_present,
                metadata=dict(root.metadata),
                level=1,
                parent_name=None,
                cleanup_strategy=self.root_marker_strategy,
            )
            compiled.children = self._compile_children(root.children, root.name, set(canonical), set(canonical), 2)
            compiled_roots.append(compiled)
        return compiled_roots

    def _build_compilation_report(self) -> pd.DataFrame:
        rows = []
        for program in flatten_compiled_programs(self.compiled_programs):
            rows.append({
                "label": program.name,
                "level": program.level,
                "parent_label": program.parent_name,
                "cleanup_strategy": program.cleanup_strategy,
                "canonical_markers_total": len(program.canonical_markers),
                "effective_markers_total": len(program.effective_markers),
                "negative_markers_total": len(program.negative_markers),
                "effective_marker_loss": len(program.canonical_markers) - len(program.effective_markers),
            })
        return pd.DataFrame(rows)

    def _infer_preset(self, expr: pd.DataFrame) -> Dict[str, object]:
        n_genes = int(expr.shape[0])
        bulk = compute_bulk_expression(expr)
        low_expr_threshold = float(np.quantile(bulk.to_numpy(dtype=float), 0.15)) if n_genes > 0 else 0.0
        low_expr_fraction = float((bulk <= max(low_expr_threshold, 0.05)).mean()) if n_genes > 0 else 0.0
        if n_genes < 2500:
            preset_name = "small_panel"
        elif n_genes > 10000 or (n_genes > 8000 and low_expr_fraction >= 0.4):
            preset_name = "large_panel"
        else:
            preset_name = "medium_panel"
        defaults = dict(PRESET_DEFAULTS[preset_name])
        defaults.update({
            "resolved_preset": preset_name,
            "panel_gene_count": n_genes,
            "low_expr_fraction": low_expr_fraction,
            "low_expr_threshold": low_expr_threshold,
        })
        return defaults

    def _resolve_runtime_config(self, expr: pd.DataFrame) -> Dict[str, object]:
        inferred = self._infer_preset(expr) if self.preset == "auto" else {
            **PRESET_DEFAULTS[self.preset],
            "resolved_preset": self.preset,
            "panel_gene_count": int(expr.shape[0]),
            "low_expr_fraction": float((compute_bulk_expression(expr) <= 0.05).mean()),
            "low_expr_threshold": 0.05,
        }
        resolved = dict(inferred)
        resolved["requested_preset"] = self.preset
        resolved["n_bins"] = self.n_bins if self.n_bins is not None else int(inferred["n_bins"])
        if self.control_size_fallbacks is not None:
            control_size_fallbacks = [int(x) for x in self.control_size_fallbacks if int(x) > 0]
        else:
            default_fs = inferred.get("control_size_fallbacks") or [self.ctrl_size, max(5, self.ctrl_size // 2), 10, 5]
            control_size_fallbacks = [int(x) for x in default_fs if int(x) > 0]
        resolved["control_size_fallbacks"] = list(dict.fromkeys(control_size_fallbacks or [max(1, int(self.ctrl_size))]))
        resolved["control_gene_exclusion_policy"] = self.control_gene_exclusion_policy if self.control_gene_exclusion_policy is not None else inferred["control_gene_exclusion_policy"]
        resolved["min_bulk_expr_for_control"] = self.min_bulk_expr_for_control if self.min_bulk_expr_for_control is not None else inferred.get("min_bulk_expr_for_control")
        if self.exclude_all_marker_genes_from_controls is not None:
            resolved["control_gene_exclusion_policy"] = "all_pipeline_markers" if self.exclude_all_marker_genes_from_controls else "current_program_only"
        resolved["ctrl_size"] = int(self.ctrl_size)
        resolved["negative_weight"] = float(self.negative_weight)
        resolved["child_marker_strategy"] = self.child_marker_strategy
        resolved["root_marker_strategy"] = self.root_marker_strategy
        resolved["fallback_to_canonical_if_sparse"] = bool(self.fallback_to_canonical_if_sparse)
        resolved["collapse_low_expr_tail"] = bool(inferred.get("collapse_low_expr_tail", False) if self.collapse_low_expr_tail is None else self.collapse_low_expr_tail)
        resolved["low_expr_threshold"] = float(inferred.get("low_expr_threshold", 0.0) if self.low_expr_threshold is None else self.low_expr_threshold)
        if self.detection_floor is not None:
            resolved["detection_floor"] = float(self.detection_floor)
        else:
            q = float(inferred.get("detection_floor_quantile", 0.1))
            resolved["detection_floor"] = float(compute_bulk_expression(expr).quantile(q)) if len(expr) else 0.0
        resolved["delta_clip_quantile"] = self.delta_clip_quantile if self.delta_clip_quantile is not None else inferred.get("delta_clip_quantile")
        return resolved

    def _flatten_with_relationships(self):
        flat = flatten_compiled_programs(self.compiled_programs)
        by_name = {p.name: p for p in flat}
        children_lookup = {p.name: [c.name for c in p.children] for p in flat}
        descendants: Dict[str, Set[str]] = {}
        def gather_desc(name: str) -> Set[str]:
            out: Set[str] = set()
            for child_name in children_lookup.get(name, []):
                out.add(child_name)
                out |= gather_desc(child_name)
            return out
        for name in by_name:
            descendants[name] = gather_desc(name)
        ancestors: Dict[str, Set[str]] = {}
        for p in flat:
            lineage: Set[str] = set()
            current = p.parent_name
            while current is not None:
                lineage.add(current)
                current = by_name[current].parent_name if current in by_name else None
            ancestors[p.name] = lineage
        siblings: Dict[str, Set[str]] = {}
        for p in flat:
            if p.parent_name is None:
                siblings[p.name] = {x.name for x in flat if x.parent_name is None and x.name != p.name}
            else:
                siblings[p.name] = {x.name for x in flat if x.parent_name == p.parent_name and x.name != p.name}
        return flat, ancestors, descendants, siblings

    def _program_marker_universe(self, program: CompiledMarkerProgram) -> Set[str]:
        return set(program.canonical_markers) | set(program.effective_markers) | set(program.negative_markers)

    def _build_exclusion_sets(self, policy: ControlGeneExclusionPolicy) -> Dict[str, Set[str]]:
        flat, ancestors, descendants, siblings = self._flatten_with_relationships()
        by_name = {p.name: p for p in flat}
        all_markers = set().union(*(self._program_marker_universe(p) for p in flat)) if flat else set()
        exclusion_sets: Dict[str, Set[str]] = {}
        for program in flat:
            if policy == "all_pipeline_markers":
                exclusion_sets[program.name] = set(all_markers)
            elif policy == "current_program_only":
                exclusion_sets[program.name] = self._program_marker_universe(program)
            elif policy == "current_branch":
                related_names = {program.name} | ancestors[program.name] | descendants[program.name] | siblings[program.name]
                excluded: Set[str] = set()
                for name in related_names:
                    excluded |= self._program_marker_universe(by_name[name])
                exclusion_sets[program.name] = excluded
            else:
                raise ValueError(f"Unsupported control gene exclusion policy: {policy}")
        return exclusion_sets

    def _build_control_maps(self, expr: pd.DataFrame, resolved_config: Dict[str, object]):
        bulk = compute_bulk_expression(expr)
        gene_bins = assign_expression_bins(
            bulk,
            n_bins=int(resolved_config["n_bins"]),
            collapse_low_expr_tail=bool(resolved_config["collapse_low_expr_tail"]),
            low_expr_threshold=float(resolved_config["low_expr_threshold"]),
        )
        min_bulk = resolved_config.get("min_bulk_expr_for_control")
        eligible_control_genes = set(gene_bins.index) if min_bulk is None else set(bulk[bulk >= float(min_bulk)].index)
        bin_lookup = build_gene_bin_lookup(gene_bins, eligible_genes=eligible_control_genes)
        rng = np.random.default_rng(self.random_state)
        flat_programs = flatten_compiled_programs(self.compiled_programs)
        exclusion_sets = self._build_exclusion_sets(resolved_config["control_gene_exclusion_policy"])

        control_maps: Dict[str, Dict[str, Dict[str, object]]] = {}
        for program in flat_programs:
            pos_genes = program.effective_markers if len(program.effective_markers) > 0 else program.canonical_markers
            neg_genes = program.negative_markers
            excluded = exclusion_sets.get(program.name, set())
            pos_controls, pos_meta = sample_control_genes_for_markers(pos_genes, gene_bins, bin_lookup, list(resolved_config["control_size_fallbacks"]), rng=rng, excluded_genes=excluded)
            neg_controls, neg_meta = sample_control_genes_for_markers(neg_genes, gene_bins, bin_lookup, list(resolved_config["control_size_fallbacks"]), rng=rng, excluded_genes=excluded)
            pos_meta["exclusion_policy"] = resolved_config["control_gene_exclusion_policy"]
            neg_meta["exclusion_policy"] = resolved_config["control_gene_exclusion_policy"]
            control_maps[program.name] = {"positive": {"controls": pos_controls, "metadata": pos_meta}, "negative": {"controls": neg_controls, "metadata": neg_meta}}
        return control_maps

    def _add_decision_metrics(self, all_scores: pd.DataFrame) -> pd.DataFrame:
        df = all_scores.copy()
        if df.empty:
            return df

        def robust_z(s: pd.Series) -> pd.Series:
            vals = s.astype(float)
            out = pd.Series(np.nan, index=s.index, dtype=float)
            finite = np.isfinite(vals.to_numpy())
            if finite.sum() == 0:
                return out
            if finite.sum() == 1:
                out.loc[vals.index[finite]] = 0.0
                return out

            finite_vals = vals.to_numpy()[finite]
            med = np.median(finite_vals)
            mad = np.median(np.abs(finite_vals - med))
            scale = 1.4826 * mad if mad > 0 else np.std(finite_vals)
            if not np.isfinite(scale) or scale <= 0:
                scale = 1.0

            out.loc[vals.index[finite]] = (finite_vals - med) / scale
            return out

        df["program_robust_zscore"] = df.groupby("label", sort=False)["score"].transform(robust_z)
        df["sibling_specificity_score"] = np.nan
        group_cols = ["cluster", "level", "parent_label"]
        for _, idx in df.groupby(group_cols, dropna=False, sort=False).groups.items():
            sub = df.loc[idx, "program_robust_zscore"].astype(float)
            arr = sub.to_numpy()
            finite_mask = np.isfinite(arr)

            if finite_mask.sum() == 0:
                df.loc[idx, "sibling_specificity_score"] = np.nan
            elif len(arr) == 1:
                df.loc[idx, "sibling_specificity_score"] = arr
            else:
                out = np.full_like(arr, np.nan, dtype=float)
                for i in range(len(arr)):
                    if not np.isfinite(arr[i]):
                        out[i] = np.nan
                        continue
                    others = np.delete(arr, i)
                    finite_others = others[np.isfinite(others)]
                    if finite_others.size == 0:
                        out[i] = arr[i]
                    else:
                        out[i] = arr[i] - np.max(finite_others)
                df.loc[idx, "sibling_specificity_score"] = out

        support = df["markers_detection_support_fraction"].fillna(0.0).clip(lower=0.0, upper=1.0)
        df["decision_score"] = (
            df["sibling_specificity_score"].fillna(df["program_robust_zscore"]).astype(float)
            * (0.5 + 0.5 * support)
        )

        diag_frames = []
        score_col = "decision_score" if "decision_score" in df.columns else "score"
        for cluster in pd.Index(df["cluster"].astype(str).unique()):
            diag = compute_cluster_branch_diagnostics(
                cluster=str(cluster),
                root_programs=self.compiled_programs,
                all_scores=df,
                score_col=score_col,
                parent_weight=self.branch_support_parent_weight,
                descendant_weight=self.branch_child_rescue_weight,
                raw_weight=self.branch_routing_raw_weight,
                weak_raw_score_threshold=self.weak_raw_score_threshold,
            )
            if diag is not None and not diag.empty:
                diag_frames.append(diag)
        if diag_frames:
            diag_df = pd.concat(diag_frames, ignore_index=True)
            merge_cols = ["cluster", "label"]
            extra_cols = [c for c in diag_df.columns if c not in merge_cols]
            df = df.merge(diag_df[merge_cols + extra_cols], on=merge_cols, how="left")

        return df

    

    def _build_hierarchy_fit_summary(self, cluster_annotations: pd.DataFrame, all_scores: pd.DataFrame) -> pd.DataFrame:
        fit_thresholds = {
            "unknown_frac_severe": 0.50,
            "unknown_frac_moderate": 0.25,
            "low_abs_frac_severe": 0.50,
            "low_abs_frac_moderate": 0.30,
            "backoff_frac_severe": 0.40,
            "backoff_frac_moderate": 0.20,
            "median_branch_raw_severe": 0.10,
            "median_branch_raw_moderate": 0.25,
            "median_final_score_severe": 0.25,
            "median_final_score_moderate": 0.60,
            "median_marker_support_severe": 0.40,
            "median_marker_support_moderate": 0.60,
            "frac_final_score_ge_1_severe": 0.25,
            "frac_final_score_ge_1_moderate": 0.50,
            "frac_branch_raw_ge_025_severe": 0.25,
            "frac_branch_raw_ge_025_moderate": 0.50,
            "frac_l1_positive_branch_support_severe": 0.25,
            "frac_l1_positive_branch_support_moderate": 0.50,
        }
        if cluster_annotations is None or cluster_annotations.empty:
            rows = [
                {"metric": "hierarchy_fit_status", "value": "unknown"},
                {"metric": "hierarchy_fit_warning", "value": "No cluster annotations available."},
            ]
            return pd.DataFrame(rows)

        df = cluster_annotations.copy()
        n_clusters = int(len(df))
        export_status = df["annot_export_status"].astype(str) if "annot_export_status" in df.columns else pd.Series(["resolved"] * n_clusters, index=df.index)
        low_evidence = df["annot_low_evidence"].fillna(False).astype(bool) if "annot_low_evidence" in df.columns else pd.Series(False, index=df.index)
        low_abs = df["annot_low_absolute_support"].fillna(False).astype(bool) if "annot_low_absolute_support" in df.columns else pd.Series(False, index=df.index)
        backoff = df["annot_backoff_applied"].fillna(False).astype(bool) if "annot_backoff_applied" in df.columns else pd.Series(False, index=df.index)
        mixed = df["annot_parent_mixing"].fillna(False).astype(bool) if "annot_parent_mixing" in df.columns else pd.Series(False, index=df.index)
        branch_conflict = df["annot_branch_conflict"].fillna(False).astype(bool) if "annot_branch_conflict" in df.columns else pd.Series(False, index=df.index)

        final_score = pd.to_numeric(df.get("annot_score"), errors="coerce") if "annot_score" in df.columns else pd.Series(dtype=float)
        branch_raw = pd.to_numeric(df.get("annot_branch_supported_raw_score"), errors="coerce") if "annot_branch_supported_raw_score" in df.columns else pd.Series(dtype=float)
        top_gap = pd.to_numeric(df.get("annot_top_level_gap"), errors="coerce") if "annot_top_level_gap" in df.columns else pd.Series(dtype=float)

        unknown_frac = float((export_status == "unknown").mean()) if n_clusters else 0.0
        low_evidence_frac = float(low_evidence.mean()) if n_clusters else 0.0
        low_abs_frac = float(low_abs.mean()) if n_clusters else 0.0
        backoff_frac = float(backoff.mean()) if n_clusters else 0.0
        mixed_frac = float(mixed.mean()) if n_clusters else 0.0
        branch_conflict_frac = float(branch_conflict.mean()) if n_clusters else 0.0
        median_final_score = float(final_score.median()) if len(final_score) else float("nan")
        median_branch_raw = float(branch_raw.median()) if len(branch_raw) else float("nan")
        median_top_gap = float(top_gap.median()) if len(top_gap) else float("nan")
        frac_final_score_ge_1 = float((final_score >= 1.0).mean()) if len(final_score) else 0.0
        frac_branch_raw_ge_025 = float((branch_raw >= 0.25).mean()) if len(branch_raw) else 0.0
        frac_branch_raw_nonneg = float((branch_raw >= 0.0).mean()) if len(branch_raw) else 0.0

        median_marker_support = float("nan")
        top_branch_col = "annot_final_top_branch" if "annot_final_top_branch" in df.columns else None
        n_supported_top_branches = int(df[top_branch_col].dropna().astype(str).nunique()) if top_branch_col is not None else 0
        median_branch_support_score = float("nan")
        median_local_node_score = float("nan")
        frac_l1_positive_branch_support = float("nan")
        positive_raw_programs = 0
        positive_branch_raw_programs = 0

        if all_scores is not None and not all_scores.empty:
            if "markers_detection_support_fraction" in all_scores.columns:
                marker_support = pd.to_numeric(all_scores["markers_detection_support_fraction"], errors="coerce")
                median_marker_support = float(marker_support.median()) if len(marker_support) else float("nan")
            if "score" in all_scores.columns:
                score_max = all_scores.groupby("label")["score"].max()
                positive_raw_programs = int((pd.to_numeric(score_max, errors="coerce") > 0).sum())
            if "branch_supported_raw_score" in all_scores.columns:
                branch_raw_max = all_scores.groupby("label")["branch_supported_raw_score"].max()
                positive_branch_raw_programs = int((pd.to_numeric(branch_raw_max, errors="coerce") > 0).sum())
            if "branch_support_score" in all_scores.columns:
                vals = pd.to_numeric(all_scores["branch_support_score"], errors="coerce")
                median_branch_support_score = float(vals.median()) if len(vals) else float("nan")
                l1 = all_scores[all_scores["level"] == 1].copy() if "level" in all_scores.columns else pd.DataFrame()
                if not l1.empty:
                    top = l1.sort_values(["cluster", "branch_support_score"], ascending=[True, False]).groupby("cluster", sort=False).head(1)
                    top_vals = pd.to_numeric(top["branch_support_score"], errors="coerce")
                    frac_l1_positive_branch_support = float((top_vals > 0).mean()) if len(top_vals) else float("nan")
            if "local_node_score" in all_scores.columns:
                vals = pd.to_numeric(all_scores["local_node_score"], errors="coerce")
                median_local_node_score = float(vals.median()) if len(vals) else float("nan")

        severe_flags = 0
        moderate_flags = 0
        if unknown_frac >= fit_thresholds["unknown_frac_severe"]:
            severe_flags += 1
        elif unknown_frac >= fit_thresholds["unknown_frac_moderate"]:
            moderate_flags += 1
        if low_abs_frac >= fit_thresholds["low_abs_frac_severe"]:
            severe_flags += 1
        elif low_abs_frac >= fit_thresholds["low_abs_frac_moderate"]:
            moderate_flags += 1
        if backoff_frac >= fit_thresholds["backoff_frac_severe"]:
            severe_flags += 1
        elif backoff_frac >= fit_thresholds["backoff_frac_moderate"]:
            moderate_flags += 1
        if pd.notna(median_branch_raw) and median_branch_raw < fit_thresholds["median_branch_raw_severe"]:
            severe_flags += 1
        elif pd.notna(median_branch_raw) and median_branch_raw < fit_thresholds["median_branch_raw_moderate"]:
            moderate_flags += 1
        if pd.notna(median_final_score) and median_final_score < fit_thresholds["median_final_score_severe"]:
            severe_flags += 1
        elif pd.notna(median_final_score) and median_final_score < fit_thresholds["median_final_score_moderate"]:
            moderate_flags += 1
        if pd.notna(median_marker_support) and median_marker_support < fit_thresholds["median_marker_support_severe"]:
            severe_flags += 1
        elif pd.notna(median_marker_support) and median_marker_support < fit_thresholds["median_marker_support_moderate"]:
            moderate_flags += 1
        if frac_final_score_ge_1 < fit_thresholds["frac_final_score_ge_1_severe"]:
            severe_flags += 1
        elif frac_final_score_ge_1 < fit_thresholds["frac_final_score_ge_1_moderate"]:
            moderate_flags += 1
        if frac_branch_raw_ge_025 < fit_thresholds["frac_branch_raw_ge_025_severe"]:
            severe_flags += 1
        elif frac_branch_raw_ge_025 < fit_thresholds["frac_branch_raw_ge_025_moderate"]:
            moderate_flags += 1
        if pd.notna(frac_l1_positive_branch_support) and frac_l1_positive_branch_support < fit_thresholds["frac_l1_positive_branch_support_severe"]:
            severe_flags += 1
        elif pd.notna(frac_l1_positive_branch_support) and frac_l1_positive_branch_support < fit_thresholds["frac_l1_positive_branch_support_moderate"]:
            moderate_flags += 1

        if severe_flags >= 2:
            fit_status = "poor"
            warning = "Hierarchy fit appears poor; absolute branch support is frequently weak or unresolved."
        elif severe_flags >= 1 or moderate_flags >= 3:
            fit_status = "questionable"
            warning = "Hierarchy fit is questionable; review branch-supported raw support and unresolved fractions."
        else:
            fit_status = "good"
            warning = "Hierarchy fit looks acceptable."

        fit_rank_score = (
            1.4 * (median_branch_raw if pd.notna(median_branch_raw) else 0.0)
            + 1.0 * (median_final_score if pd.notna(median_final_score) else 0.0)
            + 0.8 * (frac_l1_positive_branch_support if pd.notna(frac_l1_positive_branch_support) else 0.0)
            + 0.5 * (median_marker_support if pd.notna(median_marker_support) else 0.0)
            - 1.0 * unknown_frac
            - 0.8 * low_abs_frac
            - 0.6 * backoff_frac
        )

        rows = [
            {"metric": "hierarchy_fit_status", "value": fit_status},
            {"metric": "hierarchy_fit_warning", "value": warning},
            {"metric": "fit_rank_score", "value": fit_rank_score},
            {"metric": "n_clusters", "value": n_clusters},
            {"metric": "fraction_unknown", "value": unknown_frac},
            {"metric": "fraction_low_evidence", "value": low_evidence_frac},
            {"metric": "fraction_low_absolute_support", "value": low_abs_frac},
            {"metric": "fraction_backoff_applied", "value": backoff_frac},
            {"metric": "fraction_parent_mixing", "value": mixed_frac},
            {"metric": "fraction_branch_conflict", "value": branch_conflict_frac},
            {"metric": "median_final_score", "value": median_final_score},
            {"metric": "median_branch_supported_raw_score", "value": median_branch_raw},
            {"metric": "median_marker_detection_support_fraction", "value": median_marker_support},
            {"metric": "fraction_l1_positive_branch_support", "value": frac_l1_positive_branch_support},
            {"metric": "median_branch_support_score", "value": median_branch_support_score},
            {"metric": "median_local_node_score", "value": median_local_node_score},
            {"metric": "fraction_clusters_final_score_ge_1", "value": frac_final_score_ge_1},
            {"metric": "fraction_clusters_branch_supported_raw_ge_0_25", "value": frac_branch_raw_ge_025},
            {"metric": "fraction_clusters_branch_supported_raw_nonnegative", "value": frac_branch_raw_nonneg},
            {"metric": "median_top_level_gap", "value": median_top_gap},
            {"metric": "n_supported_top_level_branches", "value": n_supported_top_branches},
            {"metric": "n_programs_with_positive_raw_support", "value": positive_raw_programs},
            {"metric": "n_programs_with_positive_branch_supported_raw_support", "value": positive_branch_raw_programs},
        ]
        return pd.DataFrame(rows)

    def _build_diagnostics_summary(self, all_scores: pd.DataFrame, resolved_config: Dict[str, object]) -> pd.DataFrame:
        rows = [
            {"metric": "resolved_preset", "value": str(resolved_config["resolved_preset"])},
            {"metric": "panel_gene_count", "value": int(resolved_config["panel_gene_count"])},
            {"metric": "low_expr_fraction", "value": float(resolved_config["low_expr_fraction"])},
            {"metric": "low_expr_threshold", "value": float(resolved_config["low_expr_threshold"])},
            {"metric": "detection_floor", "value": float(resolved_config["detection_floor"])},
            {"metric": "collapse_low_expr_tail", "value": bool(resolved_config["collapse_low_expr_tail"])},
            {"metric": "n_bins", "value": int(resolved_config["n_bins"])},
            {"metric": "control_gene_exclusion_policy", "value": str(resolved_config["control_gene_exclusion_policy"])},
            {"metric": "control_size_fallbacks", "value": list(resolved_config["control_size_fallbacks"])},
            {"metric": "branch_support_parent_weight", "value": float(self.branch_support_parent_weight)},
            {"metric": "branch_support_descendant_weight", "value": float(self.branch_child_rescue_weight)},
            {"metric": "branch_support_weights_sum_to_one", "value": bool(abs((self.branch_support_parent_weight + self.branch_child_rescue_weight) - 1.0) < 1e-12)},
        ]
        if all_scores.empty:
            return pd.DataFrame(rows)

        valid_mask = all_scores["score"].notna()
        score_status = all_scores["score_status"].astype(str)
        summary_rows = [
            {"metric": "n_scored_rows", "value": int(len(all_scores))},
            {"metric": "n_valid_scores", "value": int(valid_mask.sum())},
            {"metric": "n_unscorable_rows", "value": int((~valid_mask).sum())},
            {"metric": "n_scored_clusters", "value": int(all_scores.loc[valid_mask, "cluster"].nunique())},
            {"metric": "n_clusters_with_any_unscorable_program", "value": int(all_scores.loc[~valid_mask, "cluster"].nunique())},
            {"metric": "median_valid_programs_per_cluster", "value": float(all_scores.loc[valid_mask].groupby("cluster").size().median()) if valid_mask.any() else 0.0},
            {"metric": "median_unscorable_programs_per_cluster", "value": float(all_scores.loc[~valid_mask].groupby("cluster").size().median()) if (~valid_mask).any() else 0.0},
            {"metric": "n_control_fallback_rows", "value": int(all_scores["control_fallback_used"].sum())},
            {"metric": "n_no_controls_rows", "value": int((score_status == "no_controls_available").sum())},
            {"metric": "n_insufficient_markers_rows", "value": int(score_status.isin(["insufficient_markers_present", "no_markers_present"]).sum())},
            {"metric": "median_positive_controls_per_marker", "value": float(all_scores["median_controls_per_positive_marker"].median())},
            {"metric": "median_marker_detection_support_fraction", "value": float(all_scores["markers_detection_support_fraction"].median())},
        ]
        return pd.DataFrame(rows + summary_rows)
    def _expand_annotation_dict_columns(self, annotations_df: pd.DataFrame) -> pd.DataFrame:
        df = annotations_df.copy()
        dict_cols = [
            "level_labels",
            "level_scores",
            "level_raw_scores",
            "level_branch_supported_raw_scores",
            "level_margins",
            "level_confidence",
        ]
        for col in dict_cols:
            if col not in df.columns:
                continue
            series = df[col]
            if len(series) == 0:
                continue
            keys = set()
            for val in series:
                if isinstance(val, dict):
                    keys.update(val.keys())
            for key in sorted(keys):
                if col == "level_labels":
                    out_col = f"{key}_label"
                elif col == "level_scores":
                    out_col = f"{key}_score"
                elif col == "level_raw_scores":
                    out_col = f"{key}_raw_score"
                elif col == "level_branch_supported_raw_scores":
                    out_col = f"{key}_branch_supported_raw_score"
                elif col == "level_margins":
                    out_col = f"{key}_margin"
                elif col == "level_confidence":
                    out_col = f"{key}_confidence"
                else:
                    continue
                df[out_col] = series.map(lambda d: d.get(key) if isinstance(d, dict) else None)
        return df


    def _apply_absolute_support_backoff(self, annotations_df: pd.DataFrame) -> pd.DataFrame:
        df = annotations_df.copy()
        if df.empty:
            return df

        min_raw = 0.0
        for idx, row in df.iterrows():
            decision_level = int(row.get("final_level", 0) or 0)
            decision_label = row.get("final_label", "Unresolved")
            decision_path = row.get("final_path", "Unresolved")
            decision_score = float(row.get("final_score", np.nan))
            decision_raw_score = float(row.get("final_raw_score", np.nan))
            decision_branch_raw = float(row.get("final_branch_supported_raw_score", np.nan))

            df.at[idx, "decision_label"] = decision_label
            df.at[idx, "decision_path"] = decision_path
            df.at[idx, "decision_level"] = decision_level
            df.at[idx, "decision_score"] = decision_score
            df.at[idx, "decision_raw_score"] = decision_raw_score
            df.at[idx, "decision_branch_supported_raw_score"] = decision_branch_raw
            df.at[idx, "backoff_applied"] = False
            df.at[idx, "backoff_steps"] = 0
            df.at[idx, "backoff_reason"] = ""

            if decision_level <= 0:
                continue

            branch_raw = float(row.get(f"level_{decision_level}_branch_supported_raw_score", decision_branch_raw))
            if np.isfinite(branch_raw) and branch_raw >= min_raw:
                continue

            fallback_level = None
            for lev in range(decision_level - 1, 0, -1):
                candidate_raw = float(row.get(f"level_{lev}_branch_supported_raw_score", np.nan))
                if np.isfinite(candidate_raw) and candidate_raw >= min_raw:
                    fallback_level = lev
                    break

            if fallback_level is None:
                continue

            df.at[idx, "final_level"] = int(fallback_level)
            df.at[idx, "final_label"] = row.get(f"level_{fallback_level}_label", decision_label)
            labels = [str(row.get(f"level_{lev}_label")) for lev in range(1, fallback_level + 1) if pd.notna(row.get(f"level_{lev}_label"))]
            df.at[idx, "final_path"] = " > ".join(labels) if labels else str(df.at[idx, "final_label"])
            df.at[idx, "final_score"] = float(row.get(f"level_{fallback_level}_score", np.nan))
            df.at[idx, "final_raw_score"] = float(row.get(f"level_{fallback_level}_raw_score", np.nan))
            df.at[idx, "final_branch_supported_raw_score"] = float(row.get(f"level_{fallback_level}_branch_supported_raw_score", np.nan))
            df.at[idx, "final_margin"] = float(row.get(f"level_{fallback_level}_margin", np.nan))
            df.at[idx, "confidence"] = row.get(f"level_{fallback_level}_confidence", row.get("confidence", "none"))
            df.at[idx, "status"] = "assigned"
            df.at[idx, "stop_reason"] = "backed_off_low_absolute_support"
            df.at[idx, "backoff_applied"] = True
            df.at[idx, "backoff_steps"] = int(decision_level - fallback_level)
            df.at[idx, "backoff_reason"] = "low_absolute_support"
        return df


    def _format_cluster_annotations(self, cluster_annotations: pd.DataFrame, all_scores: pd.DataFrame) -> pd.DataFrame:
        if cluster_annotations.empty:
            return cluster_annotations.copy()

        deeper_unresolved = cluster_annotations["stop_reason"].astype(str).isin(
            ["below_score_threshold", "ambiguous_sibling_margin"]
        ) & cluster_annotations["final_level"].fillna(0).astype(int).gt(0)

        low_evidence = (
            cluster_annotations["confidence"].astype(str).isin(["low", "none"])
            | cluster_annotations["final_score"].isna()
            | cluster_annotations["final_level"].fillna(0).astype(int).eq(0)
        )

        summary = pd.DataFrame({
            "cluster_id": cluster_annotations["cluster"].astype(str),
            "annot_label": cluster_annotations["final_label"],
            "annot_label_with_cluster": cluster_annotations["final_label"].astype(str).map(lambda x: str(x).strip().replace(" ", "_").replace("/", "_").replace("-", "_")) + "_" + cluster_annotations["cluster"].astype(str),
            "annot_path": cluster_annotations["final_path"],
            "annot_level": cluster_annotations["final_level"],
            "annot_status": cluster_annotations["status"],
            "annot_stop_reason": cluster_annotations["stop_reason"],
            "annot_confidence": cluster_annotations["confidence"],
            "annot_score": cluster_annotations["final_score"],
            "annot_raw_score": cluster_annotations["final_raw_score"],
            "annot_branch_supported_raw_score": cluster_annotations["final_branch_supported_raw_score"],
            "annot_margin": cluster_annotations["final_margin"],
            "annot_best_any_level_label": cluster_annotations["best_label_any_level"],
            "annot_best_any_level_score": cluster_annotations["best_score_any_level"],
            "annot_low_evidence": low_evidence,
            "annot_low_absolute_support": pd.to_numeric(cluster_annotations["final_branch_supported_raw_score"], errors="coerce").lt(self.min_branch_supported_raw_score).fillna(True),
            "annot_deeper_level_unresolved": deeper_unresolved,
            "annot_branch_conflict": False,
            "annot_decision_label": cluster_annotations["decision_label"],
            "annot_decision_path": cluster_annotations["decision_path"],
            "annot_decision_level": cluster_annotations["decision_level"],
            "annot_decision_score": cluster_annotations["decision_score"],
            "annot_decision_raw_score": cluster_annotations["decision_raw_score"],
            "annot_decision_branch_supported_raw_score": cluster_annotations["decision_branch_supported_raw_score"],
            "annot_backoff_applied": cluster_annotations["backoff_applied"],
            "annot_backoff_steps": cluster_annotations["backoff_steps"],
            "annot_backoff_reason": cluster_annotations["backoff_reason"],
        })

        # Attach raw score for the best-any-level label so downstream export
        # logic can evaluate blocked strong candidates consistently.
        if "annot_best_any_level_raw_score" not in summary.columns and not all_scores.empty:
            lookup = all_scores[["cluster", "label", "score"]].copy()
            lookup["cluster"] = lookup["cluster"].astype(str)
            lookup["label"] = lookup["label"].astype(str)
            lookup = lookup.rename(columns={
                "cluster": "cluster_id",
                "label": "annot_best_any_level_label",
                "score": "annot_best_any_level_raw_score",
            })
            lookup = lookup.drop_duplicates(subset=["cluster_id", "annot_best_any_level_label"], keep="first")
            summary = summary.merge(
                lookup,
                on=["cluster_id", "annot_best_any_level_label"],
                how="left",
            )
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
                    "annot_mixed_compartment": bool(np.isfinite(gap) and gap < self.margin_threshold),
                    "annot_parent_mixing": False,
                    "annot_mixed_primary": None,
                    "annot_mixed_secondary": None,
                })
            purity_df = pd.DataFrame(purity_rows)
            summary = summary.merge(purity_df, on="cluster_id", how="left")

        primary_score = pd.to_numeric(summary.get("annot_primary_compartment_score"), errors="coerce")
        secondary_score = pd.to_numeric(summary.get("annot_secondary_compartment_score"), errors="coerce")
        top_gap = pd.to_numeric(summary.get("annot_top_level_gap"), errors="coerce")
        mix_min_score = max(float(self.score_threshold), 0.2)
        mix_max_gap = max(float(self.margin_threshold), 0.15)
        parent_mixing = (primary_score.ge(mix_min_score) & secondary_score.ge(mix_min_score) & top_gap.le(mix_max_gap)).fillna(False)
        summary["annot_parent_mixing"] = parent_mixing
        summary["annot_mixed_primary"] = summary["annot_mixed_primary"].astype(object)
        summary["annot_mixed_secondary"] = summary["annot_mixed_secondary"].astype(object)
        summary.loc[parent_mixing, "annot_mixed_primary"] = summary.loc[parent_mixing, "annot_primary_compartment"].astype(object)
        summary.loc[parent_mixing, "annot_mixed_secondary"] = summary.loc[parent_mixing, "annot_secondary_compartment"].astype(object)

        if "annot_best_any_level_label" in summary.columns and not all_scores.empty:
            rel = all_scores[["label", "parent_label"]].drop_duplicates()
            parent_map = {str(r["label"]): (None if pd.isna(r["parent_label"]) else str(r["parent_label"])) for _, r in rel.iterrows()}
            def top_ancestor(label):
                cur = str(label)
                seen = set()
                while cur in parent_map and parent_map[cur] is not None and cur not in seen:
                    seen.add(cur)
                    cur = parent_map[cur]
                return cur
            summary["annot_best_any_level_top_branch"] = summary["annot_best_any_level_label"].astype(str).map(top_ancestor)
            summary["annot_final_top_branch"] = summary["annot_path"].astype(str).str.split(" > ").str[0]
            summary["annot_branch_conflict"] = (
                summary["annot_best_any_level_top_branch"].notna()
                & summary["annot_final_top_branch"].notna()
                & (summary["annot_best_any_level_top_branch"] != summary["annot_final_top_branch"])
                & (pd.to_numeric(summary["annot_best_any_level_score"], errors="coerce") > pd.to_numeric(summary["annot_score"], errors="coerce"))
            )

            summary = _attach_final_call_competition(
                summary=summary,
                all_scores=all_scores,
                compiled_programs=self.compiled_programs,
                cluster_id_column="cluster_id",
                final_label_column="annot_label",
                final_score_column="annot_score",
                path_margin_column="annot_margin",
                overwrite=True,
            )
        return summary
    def fit_score(self, expr: pd.DataFrame) -> HierAnnotResult:
        norm_expr = preprocess_expression(expr, input_type=self.input_type, scale_factor=self.scale_factor, uppercase_genes=self.uppercase_genes)
        resolved_config = self._resolve_runtime_config(norm_expr)
        resolved_config.update({
            "score_threshold": float(self.score_threshold),
            "margin_threshold": float(self.margin_threshold),
            "branch_child_rescue_weight": float(self.branch_child_rescue_weight),
            "level1_branch_child_rescue_weight": float(self.level1_branch_child_rescue_weight),
            "branch_routing_raw_weight": float(self.branch_routing_raw_weight),
            "weak_raw_score_threshold": float(self.weak_raw_score_threshold),
            "min_branch_supported_raw_score": float(self.min_branch_supported_raw_score),
            "min_leaf_raw_score": float(self.min_leaf_raw_score),
        })
        self._resolved_config = resolved_config
        control_maps = self._build_control_maps(norm_expr, resolved_config)
        flat_programs = flatten_compiled_programs(self.compiled_programs)
        all_scores = score_programs_across_clusters(
            expr=norm_expr,
            programs=flat_programs,
            control_maps=control_maps,
            min_markers_present=self.min_markers_present,
            min_marker_fraction=self.min_marker_fraction,
            fallback_to_canonical_if_sparse=self.fallback_to_canonical_if_sparse,
            negative_weight=self.negative_weight,
            detection_floor=float(resolved_config["detection_floor"]),
            delta_clip_quantile=resolved_config.get("delta_clip_quantile"),
        )
        all_scores = self._add_decision_metrics(all_scores)
        annotations = [asdict(decide_cluster_annotation_from_scores(
            cluster=str(cluster),
            root_programs=self.compiled_programs,
            all_scores=all_scores,
            score_threshold=self.score_threshold,
            margin_threshold=self.margin_threshold,
            branch_support_parent_weight=self.branch_support_parent_weight,
            branch_support_descendant_weight=self.branch_child_rescue_weight,
            level1_branch_support_parent_weight=self.level1_branch_support_parent_weight,
            level1_branch_support_descendant_weight=self.level1_branch_child_rescue_weight,
            branch_routing_raw_weight=self.branch_routing_raw_weight,
            weak_raw_score_threshold=self.weak_raw_score_threshold,
            min_branch_supported_raw_score=self.min_branch_supported_raw_score,
            min_leaf_raw_score=self.min_leaf_raw_score,
        )) for cluster in norm_expr.columns]
        cluster_annotations_raw = pd.DataFrame(annotations)
        cluster_annotations_raw = self._expand_annotation_dict_columns(cluster_annotations_raw)
        cluster_annotations = self._format_cluster_annotations(cluster_annotations_raw, all_scores)
        level_scores = all_scores[[c for c in [
            "cluster", "label", "level", "parent_label", "score", "decision_score", "program_robust_zscore", "sibling_specificity_score",
            "positive_score", "negative_score", "positive_marker_mean", "positive_control_mean", "negative_marker_mean", "negative_control_mean",
            "markers_present", "markers_total", "markers_present_fraction", "markers_above_detection_floor", "markers_detection_support_fraction",
            "negative_markers_present", "negative_markers_total", "negative_markers_present_fraction", "negative_markers_above_detection_floor",
            "negative_markers_detection_support_fraction", "marker_source", "score_status", "used_fallback_marker_set", "control_fallback_used",
            "median_controls_per_positive_marker", "median_controls_per_negative_marker",
            "local_node_score", "best_child_branch_support_score", "best_child_branch_supported_raw_score",
            "branch_support_score", "branch_supported_raw_score", "raw_support_gate",
            "child_rescue_weight_used", "parent_weight_used", "raw_weight_used"
        ] if c in all_scores.columns]].copy()
        diagnostics_summary = self._build_diagnostics_summary(all_scores, resolved_config)
        hierarchy_fit_summary = self._build_hierarchy_fit_summary(cluster_annotations, all_scores)
        diagnostics_summary = pd.concat([diagnostics_summary, hierarchy_fit_summary], ignore_index=True)
        return HierAnnotResult(
            cluster_annotations=cluster_annotations,
            level_scores=level_scores,
            all_scores=all_scores,
            normalized_matrix=norm_expr,
            control_gene_map=control_maps,
            compiled_programs=self.compiled_programs,
            compilation_report=self.compilation_report,
            diagnostics_summary=diagnostics_summary,
            resolved_config=resolved_config,
        )
