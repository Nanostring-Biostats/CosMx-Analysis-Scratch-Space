from ._version import __version__
from .builtin import (
    compose_hierarchy,
    compile_builtin_hierarchy_for_panel,
    describe_builtin_hierarchy,
    get_builtin_hierarchy,
    get_builtin_programs,
    match_builtin_hierarchy,
    suggest_builtin_hierarchies,
    list_builtin_hierarchies,
    list_builtin_programs,
    list_builtin_marker_program_sets,
    get_builtin_marker_program_set,
    describe_builtin_marker_program_set,
    is_immune_like_lineage,
)
from .datamodels import (
    ClusterAnnotation,
    CompiledMarkerProgram,
    ControlGeneExclusionPolicy,
    HierAnnotResult,
    MarkerProgram,
    MalignantProgram,
    PresetName,
    ScoreRecord,
)
from .inspect import describe_hierarchy, format_hierarchy_tree, print_hierarchy_tree, summarize_marker_program_set, format_marker_program_set, print_marker_program_set
from .io import hierarchy_from_dict, hierarchy_to_dict, load_hierarchy, save_hierarchy
from .pipeline import HierAnnotPipeline
from .plotting import (
    plot_cluster_annotation_heatmap,
    plot_cluster_confidence,
    plot_cluster_support_scatter,
    plot_backoff_diagnostics,
    plot_top_level_compartment_scores,
    plot_hierarchy_fit_summary,
    plot_compare_hierarchy_fit,
    plot_malignant_score_heatmap,
    plot_malignant_status_scatter,
    plot_malignant_state_specificity,
    plot_malignant_annotation_summary,
    plot_integration_raw_evidence,
    plot_integration_summary,
    plot_all_available_diagnostics,
)
from .result_io import load_result_bundle, save_result_bundle
from .utils import summarize_hierarchy_fit, compare_hierarchy_fit
from .markers import (
    collect_positive_marker_genes,
    validate_marker_programs,
    compile_marker_programs_for_panel,
    compile_builtin_marker_program_set_for_panel,
)
from .aggregate import (
    aggregate_anndata_to_cluster_means,
    aggregate_expression_to_cluster_means,
    cluster_anndata_on_representation,
    cluster_and_aggregate_anndata,
)
from .workflows import (
    make_cluster_annotation_export_summary,
    rerun_cluster_annotation_result_from_scores,
    expand_cluster_annotation_to_cells,
)
from .malignant_reporting import resolve_program_report_block_preset

__all__ = [
    "MarkerProgram",
    "MalignantProgram",
    "CompiledMarkerProgram",
    "ScoreRecord",
    "ClusterAnnotation",
    "HierAnnotResult",
    "HierAnnotPipeline",

    "describe_hierarchy",
    "format_hierarchy_tree",
    "print_hierarchy_tree",
    "summarize_marker_program_set",
    "format_marker_program_set",
    "print_marker_program_set",
    "hierarchy_to_dict",
    "hierarchy_from_dict",
    "save_hierarchy",
    "load_hierarchy",
    "save_result_bundle",
    "load_result_bundle",

    "ControlGeneExclusionPolicy",
    "PresetName",

    "aggregate_expression_to_cluster_means",
    "aggregate_anndata_to_cluster_means",
    "cluster_anndata_on_representation",
    "cluster_and_aggregate_anndata",
    "make_cluster_annotation_export_summary",
    "rerun_cluster_annotation_result_from_scores",
    "expand_cluster_annotation_to_cells",

    "plot_all_available_diagnostics",
    "plot_cluster_annotation_heatmap",
    "plot_cluster_confidence",
    "plot_cluster_support_scatter",
    "plot_backoff_diagnostics",
    "plot_top_level_compartment_scores",
    "plot_hierarchy_fit_summary",
    "plot_compare_hierarchy_fit",
    "summarize_hierarchy_fit",
    "compare_hierarchy_fit",
    "collect_positive_marker_genes",
    "validate_marker_programs",
    "compile_marker_programs_for_panel",
    "compile_builtin_marker_program_set_for_panel",

    "list_builtin_hierarchies",
    "get_builtin_hierarchy",
    "describe_builtin_hierarchy",
    "list_builtin_programs",
    "get_builtin_programs",
    "match_builtin_hierarchy",
    "suggest_builtin_hierarchies",
    "compose_hierarchy",
    "compile_builtin_hierarchy_for_panel",

    "list_builtin_marker_program_sets",
    "get_builtin_marker_program_set",
    "describe_builtin_marker_program_set",
    "is_immune_like_lineage",
    "resolve_program_report_block_preset",
    "plot_malignant_score_heatmap",
    "plot_malignant_status_scatter",
    "plot_malignant_state_specificity",
    "plot_malignant_annotation_summary",
    "plot_integration_raw_evidence",
    "plot_integration_summary",
]
