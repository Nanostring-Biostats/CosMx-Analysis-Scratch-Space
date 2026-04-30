# HierAnnot

HierAnnot annotates cluster-level expression profiles with hierarchical marker programs. It scores every program with a bulk-expression-matched-control enrichment after LogNormalize-style preprocessing, routes labels with branch-support logic, and backs off weak subtypes to the deepest ancestor on the chosen path that still has sufficient absolute support.

## Core ideas

HierAnnot implements a hierarchical competitive marker-program workflow adapted for **cluster-level expression profiles**:

- Match control gene sampling by expression bin.
- Score all marker programs first, so competitive sibling-aware decision logic is separated from raw score calculation.
- Branch-aware routing combines current-node evidence with the strongest direct-child evidence in the same branch.
- Absolute-support backoff prevents a weak subtype from forcing an `unknown` label when its parent branch is still supported.
- Export logic separates:
  - resolved labels
  - low-evidence `unknown`
  - strong dual-parent `mixed_<primary>.<secondary>` patterns


## Installation

```bash
pip install -e .
```

## Quick start

```python
import pandas as pd
from hierannot import (
    get_builtin_hierarchy,
    HierAnnotPipeline,
    make_cluster_annotation_export_summary,
)

expr = pd.DataFrame(
    {
        "cluster_1": [8, 7, 0, 0, 1, 0, 0],
        "cluster_2": [0, 0, 5, 6, 0, 8, 7],
        "cluster_3": [1, 0, 0, 1, 8, 0, 0],
    },
    index=["CD3D", "TRBC1", "LYZ", "S100A8", "EPCAM", "KRT18", "KRT19"],
)

roots = get_builtin_hierarchy("tme_core")
pipeline = HierAnnotPipeline(
    root_programs=roots,
    input_type="raw_cluster_means",
    score_threshold=0.01,
    margin_threshold=0.0,
    random_state=0,
)

results = pipeline.fit_score(expr)
print(results.cluster_annotations)
print(results.level_scores.head())

cluster_report = make_cluster_annotation_export_summary(result)
print(cluster_report.head())
```

Diagnostic plots are intended to answer different questions:

- **cluster support scatter**: final decision score vs branch-supported raw score
- **top-level branch scores**: broad branch competition only
- **backoff diagnostics**: clusters that fell back from a weak subtype to a supported parent
- **hierarchy fit summary**: whether the selected ontology appears appropriate for the dataset

A basic plotting example on results is included in:
- `examples/plot_diagnostics_from_bundle.py`

### Input expectations

The main API expects a pandas `DataFrame` with:
- rows = genes
- columns = clusters

Supported input modes:
- `raw_cluster_means`: pseudo-count-like cluster means; LogNormalize is applied
- `normalized`: already normalized/log-scaled values; normalization is skipped

For starting with single-cell expression matrix in `AnnData`, see [End-to-end workflow](#end-to-end-workflow ).

### Branch-support routing

All marker programs are scored first. Hierarchical routing then uses a
branch-support score for each parent node:

- parent node raw contribution
- parent node decision contribution
- plus the best direct-child rescue contribution within that branch

This helps top-level routing when a broad parent program is weaker than a
highly specific child program, such as epithelial subtypes in breast tissue.

You can tune the weights in `HierAnnotPipeline`:

```python
pipe = HierAnnotPipeline(
    root_programs=roots,
    branch_child_rescue_weight=0.35, 
    level1_branch_child_rescue_weight=0.40, # slightly stronger rescue for L1 given their broad program
    branch_routing_raw_weight = 0.60,

    # raw-score gating on decision contribution
    weak_raw_score_threshold = 0.10,
    min_branch_supported_raw_score = 0.0,
    min_leaf_raw_score = 0.0,
)
```


## Hierarchies and marker programs 

### Defining custom marker programs

```python
from hierannot import MarkerProgram

cd4 = MarkerProgram(name="CD4 T cell", positive_markers=["IL7R", "LTB", "MAL"])
cd8 = MarkerProgram(name="CD8 T cell", positive_markers=["NKG7", "CCL5", "PRF1"])

root = MarkerProgram(
    name="T cell",
    positive_markers=["CD3D", "CD3E", "TRBC1"],
    children=[cd4, cd8],
)
```

### Built-in hierarchies

HierAnnot ships with modular built-in packs that follow explicit canonical design:

- **shared backbones**
  - `immune_core`
  - `stromal_core`
  - `vascular_core`
  - `mural_core`
  - `solid_tissue_core`

- **tissue-aware parenchymal modules**
  - `Pancreatic parenchymal`
  - `Mammary epithelial`
  - `Renal parenchymal`
  - `Hepatic parenchymal`
  - `Pulmonary parenchymal`
  - `Cutaneous epithelial`
  - `Intestinal epithelial`

- **pre-assembled tissue-aware microenvironment hierarchies**
  - `brain_core`: core neuro hierarchy for brain-like samples.
  - `tme_core`: combined solid tissue plus immune hierarchy for mixed tissue microenvironments.
  - `breast_tme`
  - `colon_tme`
  - `kidney_tme`
  - `pancreas_tme`
  - `tonsil_tme`
  - `liver_tme`
  - `lung_tme`
  - `skin_tme`

```python
from hierannot import list_builtin_hierarchies, describe_builtin_hierarchy, get_builtin_hierarchy

print(list_builtin_hierarchies())
print(describe_builtin_hierarchy("tme_core"))

# Modify safely without mutating package defaults
roots = get_builtin_hierarchy("tme_core")
roots[0].positive_markers.append("KRT17")
```

Helper functions are provided to recommend the most relevant built-in hierarchies based on free-text descriptions of tissue type. 

```python
from hierannot import match_builtin_hierarchy, suggest_builtin_hierarchies

queries = [
    "breast_tumor",
    "pancreatic cancer",
    "hepatic tumor",
    "pulmonary lesion",
    "epidermal tumor",
    "solid tumor",
]

for query in queries:
    best = match_builtin_hierarchy(query)
    print(query, "->", best)
    ranked = suggest_builtin_hierarchies(query, top_k=3)
    print(ranked)
    print()
```

Built-in hierarchies can be compiled automatically against an input RNA panel gene list to prune weakly covered leaves if desired.

```python
from hierannot import compile_builtin_hierarchy_for_panel

compiled_roots, report = compile_builtin_hierarchy_for_panel(
    gene_list=panel_genes,
    builtin_name="pancreas_tme",
    panel_name="medium_6K_panel",
    preserve_major_immune_subtypes=True,
)
```

### Extending built-ins for a new tissue hierarchy

A practical pattern is to start from an existing built-in tissue microenvironment hierarchy and then replace or extend only the branch you need. This is often cleaner than starting from scratch because you can keep the shared immune, stromal, endothelial, and mural programs unchanged.

```python
from hierannot import MarkerProgram, get_builtin_hierarchy, describe_hierarchy, format_hierarchy_tree

# Start with an existing tissue hierarchy and replace only the epithelial branch.
roots = get_builtin_hierarchy("breast_tme")

breast_epithelial = MarkerProgram(
    name="breast epithelial",
    positive_markers=["EPCAM", "KRT8", "KRT18", "KRT19"],
    children=[
        MarkerProgram(name="luminal epithelial", positive_markers=["EPCAM", "KRT8", "KRT18", "ESR1", "PGR"]),
        MarkerProgram(name="basal epithelial", positive_markers=["KRT5", "KRT14", "KRT17", "TP63"]),
    ],
)

breast_tme = []
for node in roots:
    if node.name == "Epithelial":
        breast_tme.append(breast_epithelial)
    else:
        breast_tme.append(node)

print(format_hierarchy_tree(breast_tme))
print(describe_hierarchy(breast_tme)[["name", "level", "parent", "n_positive_markers", "n_children"]])

pipeline = HierAnnotPipeline(root_programs=breast_tme, preset="auto", input_type="raw_cluster_means")
```

A good rule of thumb is to keep broad identity programs at the parent level and put subtype or state-like refinements as children.

### Saving and reloading custom hierarchies

Custom hierarchies can be saved as JSON and reloaded later for reuse.

```python
from hierannot import save_hierarchy, load_hierarchy

save_hierarchy(breast_tme, "breast_tme.json")
reloaded = load_hierarchy("breast_tme.json")
```

This makes it easy to version a custom hierarchy alongside your analysis code or share it with collaborators.

### Hierarchy-fit guardrails

HierAnnot adds dataset-level hierarchy-fit diagnostics to `result.diagnostics_summary`, which allows evaluate whether the selected ontology appears appropriate for the dataset.

```python
from hierannot import summarize_hierarchy_fit, plot_hierarchy_fit_summary

fit_summary = summarize_hierarchy_fit(result)
print(fit_summary)

plot_hierarchy_fit_summary(result)
```

One could also compare multiple candidate hierarchies side by side:

```python
from hierannot import compare_hierarchy_fit, plot_compare_hierarchy_fit

comparison = compare_hierarchy_fit({
    "breast_tme": breast_result,
    "brain_core": brain_result,
    "tme_core": tme_result,
})

print(comparison)
plot_compare_hierarchy_fit(comparison)
```

## End-to-end workflow 

### Practical workflow with AnnData and custom embeddings

A common workflow is:

1. cluster cells using a custom embedding in ``adata.obsm``
2. aggregate raw expression to gene x cluster means
3. annotate the cluster means with ``HierAnnotPipeline``

```python
from hierannot import HierAnnotPipeline, get_builtin_hierarchy
from hierannot.workflows import cluster_and_aggregate_anndata

adata, cluster_means = cluster_and_aggregate_anndata(
    adata,
    clustering_source="obsm",
    clustering_source_key="X_custom",
    cluster_key="hierannot_leiden",
    aggregation_source="layer",
    aggregation_source_key="counts",
    n_pcs=30,
    n_neighbors=15,
    leiden_resolution=0.8,
)

pipeline = HierAnnotPipeline(
    root_programs=get_builtin_hierarchy("tme_core"),
    preset="auto",
    input_type="raw_cluster_means",
)
result = pipeline.fit_score(cluster_means)
```

For smaller utility steps, the package also exposes:
- ``aggregate_expression_to_cluster_means``: Aggregate a feature x cell matrix into feature x cluster profiles.
- ``aggregate_anndata_to_cluster_means``: Wrapper for aggregate an AnnData object to feature x cluster profiles.
- ``cluster_anndata_on_representation``: Cluster cells from a chosen AnnData representation.
    - Source-aware PCA preprocessing defaults: `obsm` avoids extra centering/scaling; `X` and `layer` zero-center by default.
    - To clustering on representation directly without zero-center or PCA preprocessing, set `do_pca = False`.
- ``cluster_and_aggregate_anndata``: Cluster cells from one AnnData source, then aggregate another source by cluster.


### Efficient cluster-mean aggregation

The ``aggregate_*_to_cluster_means`` functions support multiple aggregation methods and are equipped with automatic backend handling for large-scale sparse datasets and graph-based central-cell selection. 

Supported aggregation methods include:
- ``mean`` / ``median``: standard aggregation over all included cells.
- ``trimmed_mean``: feature-wise trimmed mean over all included cells.
- ``central_cells_mean``: keep only the cells closest to each cluster centroid in ``centrality_matrix``, then average expression.
- ``central_cells_trimmed_mean``: select central cells first, then apply the feature-wise trimmed mean.

Automatic backend handling:
- `mean` and `central_cells_mean` use sparse-aware or chunked backends without densifying the full matrix.
- `trimmed_mean` and `median` fall back to dense/order-statistic reducers and warn on very large inputs.

Graph-based central-cell selection is for `central_cells_*` methods:
- `centrality_mode="centroid"` uses distance to cluster centroid in a chosen representation.
- `centrality_mode="graph_degree"` reuses a precomputed neighbor graph and keeps the cells with highest within-cluster weighted degree (i.e. the sum of weights to other cells in the same cluster).

```python
from hierannot.workflows import cluster_anndata_on_representation, aggregate_anndata_to_cluster_means

adata = cluster_anndata_on_representation(
    adata,
    cluster_key="hierannot_leiden",
    pca_key="X_pca_custom",
    neighbors_key="X_neighbors_custom",
    n_pcs=50,
    n_neighbors=30,
    leiden_resolution=0.8,
    compute_umap=False,

    source="obsm",
    source_key="X_custom",
    pca_zero_center=False,
    pca_scale=False,
)

cluster_means = aggregate_anndata_to_cluster_means(
    adata,
    cluster_key="hierannot_leiden",
    source="layer",
    source_key="counts",
    method="central_cells_mean",
    centrality_mode="graph_degree",
    neighbors_key="X_neighbors_custom",
)

```

### AnnData round-trip and plotting

```python
from hierannot import HierAnnotPipeline, get_builtin_hierarchy
from hierannot.workflows import (
    cluster_and_aggregate_anndata,
    make_cluster_annotation_export_summary,
    expand_cluster_annotation_to_cells,
)
from hierannot.plotting import plot_cluster_annotation_heatmap, plot_cluster_confidence

adata, cluster_means = cluster_and_aggregate_anndata(adata, obsm_key="X_custom")
pipe = HierAnnotPipeline(root_programs=get_builtin_hierarchy("tme_core"), preset="auto")
result = pipe.fit_score(cluster_means)
cluster_export = make_cluster_annotation_export_summary(result)
cell_join = expand_cluster_annotation_to_cells(adata.obs["hierannot_leiden"], cluster_export, cluster_key="hierannot_leiden")
adata.obs = adata.obs.join(cell_join.drop(columns=["hierannot_leiden"]))
fig, ax = plot_cluster_annotation_heatmap(result, level=1)
```

You can restrict cluster-profile aggregation to a subset of cells using `include_obs_mask`, for example `include_obs_mask="include"` to use `adata.obs["include"]`. Clusters excluded from analysis can later be filled as `unassigned` in the join table, while analyzed clusters with weak or ambiguous support can be exported as `unknown` via `make_cell_annotation_join_table(...)`.

Use `make_cluster_annotation_export_summary()` to tune unknown/mixed/export labeling at the cluster level before expanding back to cells. One could also run it with `rerun_decision=True` and new weights/thresholds to tune the decision routing in addition to the export logic if requested. Then use `expand_cluster_annotation_to_cells()` to build a join-ready per-cell table.

### Saving and reloading a result bundle

You can persist the main score tables, hierarchy, and resolved config without
pickling the entire result object:

```python
from hierannot import save_result_bundle, load_result_bundle

save_result_bundle(
    result,
    "hierannot_result_bundle",
    hierarchy=roots,
    resolved_config=getattr(result, "resolved_config", None),
    metadata={"dataset": "sample_01"},
    table_format="csv",  # default; use "parquet" optionally
)

bundle = load_result_bundle("hierannot_result_bundle")
print(bundle["cluster_annotations"].head())
print(bundle.get("resolved_config"))
```

## Design guidance

- Hierarchical assignment is the primary decision path.
- The best score across all levels is also reported as a diagnostic.
- Scores should be interpreted together with marker coverage and sibling margin.
- Use parent nodes for broad identity programs and child nodes for more specific or state-like refinements.
- ``annot_score`` is the decision score used for label assignment; ``annot_raw_score`` remains available for diagnostics.
- Top-level purity diagnostics such as ``annot_primary_compartment`` and ``annot_mixed_compartment`` help flag clusters that may be mixed or boundary-adjacent.

### Cluster-specific expression profiles for spatial data 

When minor neighborhood contamination is expected, you can aggregate cluster profiles from the cells closest to the cluster centroid in the clustering representation:

```python
from hierannot.workflows import cluster_and_aggregate_anndata

adata, cluster_means = cluster_and_aggregate_anndata(
    adata,
    clustering_source="obsm",
    clustering_source_key="X_custom",
    aggregation_source="layer",
    aggregation_source_key="counts",
    aggregation_method="central_cells_mean",
    central_fraction=0.8,
)
```

This computes centroid distances in the clustering representation, keeps the most central cells in each cluster, and then averages raw expression over those cells.

One could also use `aggregation_method="trimmed_mean"` for feature-wise trimmed averaging, which excludes oultiers on per-feature basis within each cluster separately. To combine both central-cell selection with feature-wise trimmed averaging, use `aggregation_method="central_cells_trimmed_mean"`.

