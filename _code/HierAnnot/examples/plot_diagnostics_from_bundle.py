from pathlib import Path
from hierannot import (
    load_result_bundle,
    plot_all_available_diagnostics,
)

# Load a previously saved result bundle. By default this returns a lightweight
# result-like object with attributes compatible with plotting helpers.
result = load_result_bundle("hierannot_result_bundle")

outdir = Path("diagnostics")
outdir.mkdir(exist_ok=True)

# Use wrapper to generate available diagnostic plots in one call
all_plots = plot_all_available_diagnostics(result)
for name, (fig, _axes) in all_plots.items():
    fig.savefig(outdir / f"diagnostic_{name}.png", dpi=150, bbox_inches="tight")

## Plot specific diagnostics separately if you want more control over figure saving or further customization.

# Part I: Standard diagnostics plots for normal hierarhcy annotation
from hierannot import (
    plot_cluster_annotation_heatmap,
    plot_cluster_confidence,
    plot_cluster_support_scatter,
    plot_backoff_diagnostics,
    plot_top_level_compartment_scores,
    plot_hierarchy_fit_summary,
)

# 1) Broad score heatmap
fig, ax = plot_cluster_annotation_heatmap(result, value="decision_score", level=1, top_n_labels=10)
fig.savefig(outdir / "diagnostic_heatmap.png", dpi=150, bbox_inches="tight")

# 2) Cluster confidence scatter
fig, ax = plot_cluster_confidence(result)
fig.savefig(outdir / "diagnostic_confidence.png", dpi=150, bbox_inches="tight")

# 3) Decision score vs branch-supported raw score
fig, ax = plot_cluster_support_scatter(result)
fig.savefig(outdir / "diagnostic_support_scatter.png", dpi=150, bbox_inches="tight")

# 4) Backoff diagnostics
fig, ax = plot_backoff_diagnostics(result)
fig.savefig(outdir / "diagnostic_backoff.png", dpi=150, bbox_inches="tight")

# 5) Top-level branch competition
fig, ax = plot_top_level_compartment_scores(result)
fig.savefig(outdir / "diagnostic_top_level.png", dpi=150, bbox_inches="tight")

# 6) Ontology/hierarchy fit diagnostics
fig, ax = plot_hierarchy_fit_summary(result)
fig.savefig(outdir / "diagnostic_hierarchy_fit.png", dpi=150, bbox_inches="tight")


# Part II: Malignant scoring and integration diagnostics

# If the results bundle were generated with `malignant_integration_mode="integrate" when running HierAnnotPipeline and export summary, 
# malignant scoring diagnostics plots are also available to evaluate the malignant scoring and integration performance.

from hierannot import (
    plot_malignant_score_heatmap,
    plot_malignant_status_scatter, 
    plot_malignant_state_specificity, 
    plot_malignant_annotation_summary,
    plot_integration_summary,
    plot_integration_raw_evidence, 
    plot_all_available_diagnostics,
)

# 7) cluster-by-malignant-program heatmap
fig, ax = plot_malignant_score_heatmap(result, value="status_score")
fig.savefig(outdir / "diagnostic_malignant_score_heatmap.png", dpi=150, bbox_inches="tight")

# 8) absolute malignant status evidence for each cluster
fig, ax = plot_malignant_status_scatter(result)
fig.savefig(outdir / "diagnostic_malignant_status_scatter.png", dpi=150, bbox_inches="tight")

# 9) per-cluster malignant-state support and program-level malignant-state specificity 
fig, ax = plot_malignant_state_specificity(result)
fig.savefig(outdir / "diagnostic_malignant_state_specificity.png", dpi=150, bbox_inches="tight")

# 10) malignant annotation status counts and concise state-label counts
fig, ax = plot_malignant_annotation_summary(result)
fig.savefig(outdir / "diagnostic_malignant_annotation_summary.png", dpi=150, bbox_inches="tight")

# 11) normal raw evidence against malignant raw evidence used in integration
fig, ax = plot_integration_raw_evidence(result)
fig.savefig(outdir / "diagnostic_integration_raw_evidence.png", dpi=150, bbox_inches="tight")

# 12) integrated annotation status, reason, and malignant-gate diagnostics
fig, ax = plot_integration_summary(result)
fig.savefig(outdir / "diagnostic_integration_summary.png", dpi=150, bbox_inches="tight")
