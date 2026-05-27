"""Example: cluster on a custom embedding, score on raw expression, and build a join-ready annotation table."""

from hierannot import HierAnnotPipeline, get_builtin_hierarchy
from hierannot.plotting import plot_cluster_annotation_heatmap, plot_cluster_confidence
from hierannot.workflows import (
    cluster_and_aggregate_anndata,
    make_cluster_annotation_export_summary,
    expand_cluster_annotation_to_cells,
)


# Assume `adata` already exists and includes:
# - raw counts in `adata.X` / `adata.layers[...]` or `adata.raw`
# - a custom cell embedding in `adata.obsm["X_custom"]`

adata, cluster_means = cluster_and_aggregate_anndata(
    adata,
    clustering_source="obsm",
    clustering_source_key="X_custom",
    cluster_key="hierannot_leiden",
    aggregation_source="layer",
    expression_layer="counts",
    n_pcs=30,
    n_neighbors=15,
    leiden_resolution=0.8,
    compute_umap=True,
)

pipe = HierAnnotPipeline(
    root_programs=get_builtin_hierarchy("tme_core"),
    preset="auto",
    input_type="raw_cluster_means",
)
result = pipe.fit_score(cluster_means)

# Cluster-level summary is already available from fit_score.
cluster_summary = result.cluster_annotations.copy()
print(cluster_summary)

# Retrieve more per-level information and decide rules for "unknown" labels.
cluster_export = make_cluster_annotation_export_summary(
    result,
    unknown_min_score=0.15,
    unknown_max_margin=0.03,
    label_with_cluster=True,
)

# Build a per-cell table that can be joined into adata.obs when desired.
cell_join = expand_cluster_annotation_to_cells(
    adata.obs["hierannot_leiden"],
    cluster_export,
    cluster_key="hierannot_leiden",
    include_columns=[
        "annot_export_label",
        "annot_export_label_with_cluster",
        "annot_export_status",
        "annot_confidence",
        "annot_stop_reason",
    ],
    prefix="tme",
)

print(cluster_export.head())
print(cell_join.head())

adata.obs = adata.obs.join(cell_join.drop(columns=["hierannot_leiden"], errors="ignore"))
fig1, ax1 = plot_cluster_annotation_heatmap(result, level=1)
fig2, ax2 = plot_cluster_confidence(result)


