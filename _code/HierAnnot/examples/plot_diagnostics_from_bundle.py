from hierannot import (
    load_result_bundle,
    plot_cluster_annotation_heatmap,
    plot_cluster_confidence,
    plot_cluster_support_scatter,
    plot_backoff_diagnostics,
    plot_top_level_compartment_scores,
    plot_hierarchy_fit_summary,
)

# Load a previously saved result bundle. By default this returns a lightweight
# result-like object with attributes compatible with plotting helpers.
result = load_result_bundle("hierannot_result_bundle")

# 1) Broad score heatmap
fig, ax = plot_cluster_annotation_heatmap(result, value="decision_score", level=1, top_n_labels=10)
fig.savefig("diagnostic_heatmap.png", dpi=150, bbox_inches="tight")

# 2) Cluster confidence scatter
fig, ax = plot_cluster_confidence(result)
fig.savefig("diagnostic_confidence.png", dpi=150, bbox_inches="tight")

# 3) Decision score vs branch-supported raw score
fig, ax = plot_cluster_support_scatter(result)
fig.savefig("diagnostic_support_scatter.png", dpi=150, bbox_inches="tight")

# 4) Backoff diagnostics
fig, ax = plot_backoff_diagnostics(result)
fig.savefig("diagnostic_backoff.png", dpi=150, bbox_inches="tight")

# 5) Top-level branch competition
fig, ax = plot_top_level_compartment_scores(result)
fig.savefig("diagnostic_top_level.png", dpi=150, bbox_inches="tight")

# 6) Ontology/hierarchy fit diagnostics
fig, ax = plot_hierarchy_fit_summary(result)
fig.savefig("diagnostic_hierarchy_fit.png", dpi=150, bbox_inches="tight")

