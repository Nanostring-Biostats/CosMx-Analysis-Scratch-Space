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
)
from .datamodels import (
    ClusterAnnotation,
    CompiledMarkerProgram,
    ControlGeneExclusionPolicy,
    HierAnnotResult,
    MarkerProgram,
    PresetName,
    ScoreRecord,
)
from .inspect import describe_hierarchy, format_hierarchy_tree, print_hierarchy_tree
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
)
from .result_io import load_result_bundle, save_result_bundle
from .utils import summarize_hierarchy_fit, compare_hierarchy_fit
from .workflows import (
    aggregate_anndata_to_cluster_means,
    aggregate_expression_to_cluster_means,
    make_cluster_annotation_export_summary,
    expand_cluster_annotation_to_cells,
    cluster_anndata_on_representation,
    cluster_and_aggregate_anndata,
)

__all__ = [
    "MarkerProgram",
    "CompiledMarkerProgram",
    "ScoreRecord",
    "ClusterAnnotation",
    "HierAnnotResult",
    "HierAnnotPipeline",
    "describe_hierarchy",
    "format_hierarchy_tree",
    "print_hierarchy_tree",
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
    "expand_cluster_annotation_to_cells",
    "plot_cluster_annotation_heatmap",
    "plot_cluster_confidence",
    "plot_cluster_support_scatter",
    "plot_backoff_diagnostics",
    "plot_top_level_compartment_scores",
    "plot_hierarchy_fit_summary",
    "plot_compare_hierarchy_fit",
    "summarize_hierarchy_fit",
    "compare_hierarchy_fit",
    "list_builtin_hierarchies",
    "get_builtin_hierarchy",
    "describe_builtin_hierarchy",
    "list_builtin_programs",
    "get_builtin_programs",
    "match_builtin_hierarchy",
    "suggest_builtin_hierarchies",
    "compose_hierarchy",
    "compile_builtin_hierarchy_for_panel",
]

