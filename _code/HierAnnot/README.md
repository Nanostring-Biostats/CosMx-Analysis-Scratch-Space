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
  - Optional flagging of clusters based on auxiliary cell-state tracks, e.g. the built-in malignant/transformation programs for tumor-like status.

## Table of Contents

- [Installation](#installation)
- [Quick start](#quick-start)
    - [Input expectations](#nput-expectations)
    - [Branch-support routing](#branch-support-routing)
    - [Examples and diagnostic plots](#examples-and-diagnostic-plots)
- [Hierarchies and marker programs](#hierarchies-and-marker-programs)
    - [Hierarchy design guidelines](#hierarchy-design-guidelines)
    - [Defining custom marker programs](#defining-custom-marker-programs)
    - [Built-in hierarchies](#built-in-hierarchies)
    - [Extending built-ins for a new tissue hierarchy](#extending-built-ins-for-a-new-tissue-hierarchy)
    - [Saving and reloading custom hierarchies](#saving-and-reloading-custom-hierarchies)
    - [Hierarchy-fit guardrails](#hierarchy-fit-guardrails)
- [Auxillary track annotations](#auxillary-track-annotations)
    - [Malignant scoring and tier-based reporting](#malignant-scoring-and-tier-based-reporting)
    - [Built-in and custom malignant programs](#built-in--and-custom-malignant-programs)
    - [Malignant reporting controls](#malignant-reporting-controls)
    - [Export labels with tumor-like reporting](#export-labels-with-tumor-like-reporting)
- [End-to-end workflow](#end-to-end-workflow)
    - [Practical workflow with AnnData and custom embeddings](#practical-workflow-with-anndata-and-custom-embeddings)
    - [Efficient cluster-mean aggregation](#efficient-cluster-mean-aggregation)
    - [AnnData round-trip](#anndata-round-trip)
    - [Saving and reloading a result bundle](#saving-and-reloading-a-result-bundle)
- [Notes and Tips](#notes-and-tips)
    - [Cluster-specific expression profiles for spatial data](#cluster-specific-expression-profiles-for-spatial-data)
- [License](LICENSE)

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
    malignant_integration_mode="flag_only",
    program_report_block_preset="tumor_reportable",
)

result = pipeline.fit_score(expr)
print(result.cluster_annotations)
print(result.level_scores.head())

# The export helper uses the fitted pipeline decisions by default. Set
# rerun_decision=True only when you intentionally want to recompute normal
# routing or malignant integration with new thresholds.
cluster_report = make_cluster_annotation_export_summary(result)
print(cluster_report.head())
```

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

### Examples and diagnostic plots

You can find expanded example scripts under the `examples` subfolder.
- `examples/basic_usage.py`: normal hierarchy annotation along with optional maligant integration.
- `examples/end_to_end_anndata.py`: end-to-end workflow with `AnnData` and custom embeddings.
- `examples/plot_diagnostics_from_bundle.py`: all individual diagnosit plots

Various diagnostic plots could be generated from the results and are intended to answer different questions:

- **cluster support scatter**: final decision score vs branch-supported raw score
- **top-level branch scores**: broad branch competition only
- **backoff diagnostics**: clusters that fell back from a weak subtype to a supported parent
- **hierarchy fit summary**: whether the selected ontology appears appropriate for the dataset
- **malignant status scatter / heatmap / summary**: absolute tumor-like evidence and malignant-state labels 
- **integration summary / raw-evidence scatter**: normal-vs-malignant integration status, reasons, blocklist and raw-delta diagnostics

The wrapper `plot_all_available_diagnostics(result)` generates all available diagnostic plots if one would like a complete view of the annotation quality on your clusters. 

## Hierarchies and marker programs 

### Hierarchy design guidelines

Design principles of current built-in hierarchies:

1. **Stable broad top-level branches**
   - Keep top branches broad and robust across panel sizes (for example immune, endothelial, mural, stromal/fibroblast, epithelial/parenchymal).
2. **Concentrate tissue specificity in epithelial/parenchymal branches**
   - Most tissue-specific variation should enter through epithelial or other parenchymal branches.
   - Stromal refinements can be tissue-specific when marker support is strong and clearly distinct.
   - Immune, endothelial, and mural branches should remain conservative unless a subtype is reliably separable on realistic panels.
3. **Add new nodes only when sibling separation is clear**
   - A new node should have markers that are distinct from nearby siblings. Otherwise it will compete poorly and destabilize routing.
4. **Prefer robust parents over fragile subtypes**
   - If a subtype cannot be supported across realistic panel sizes, keep the stronger parent instead of adding a brittle leaf.
5. **Keep lineage and state concepts separate**
   - State programs such as cycling, stress/IFN, EMT, or transformation are often better handled as auxiliary state tracks rather than as lineage nodes.
   - `hierannot>=0.8.2` introduces built-in `MalignantProgram` sets that could be used to flag cluster for tumor status in additional to hierarhical lineage-based cell types and allow incorporate the tumor-focused labels into the final export label. 
   - When use with `malignant_integration_mode="flag_only"`, the `malignant_programs` of `HierAnnotPipeline` can also be used as an auxiliary program-status layer for proliferation, stress, fibrosis, inflammation, or other disease-state programs without changing normal hierarchy export labels.

These same principles are useful for user-defined hierarchies so custom branches remain compatible with HierAnnot routing, panel compilation, and export logic.

### Defining custom marker programs

```python
from hierannot import MarkerProgram

cd4 = MarkerProgram(name="CD4 T cell", positive_markers=["IL7R", "LTB", "MAL"])
cd8 = MarkerProgram(name="CD8 T cell", positive_markers=["NKG7", "CCL5", "PRF1"])

roots = MarkerProgram(
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

Built-in hierarchies can be compiled against an input RNA panel gene list before scoring if desired. The hierarchy compiler returns a report so you can see which nodes remain fully supported, which fall back to a parent-level interpretation, and which are disabled.

```python
from hierannot import compile_builtin_hierarchy_for_panel

compiled_roots, hierarchy_report = compile_builtin_hierarchy_for_panel(
    gene_list=panel_genes,
    builtin_name="pancreas_tme",
    panel_name="medium_6K_panel",
    preserve_major_immune_subtypes=True,
)

hierarchy_report[[
    "name",
    "level",
    "positive_markers_present",
    "positive_markers_present_fraction",
    "status",
]].sort_values("positive_markers_present_fraction").head()
```

The `status` column is the main field to inspect:

- `full`: enough marker coverage for the node to remain in the compiled hierarchy.
- `fallback_to_parent`: the node has limited marker coverage and should be interpreted conservatively through its parent context.
- `disabled`: the node was not retained in the compiled hierarchy because marker support was too sparse.

Low marker coverage does not always mean the hierarchy is unusable. Broad parent labels may still be supported even when some fine-grained leaves are under-covered. If important L1 or L2 branches are disabled or fall back, use a broader hierarchy, expand the panel, or interpret fine-grained labels more conservatively.

Flat tumor/auxiliary marker-program sets can also be compiled against a panel. Unlike the normal hierarchy, flat programs can optionally be dropped or merged when coverage is sparse or within-group overlap is high.

```python
from hierannot import compile_builtin_marker_program_set_for_panel

compiled_programs, program_report = compile_builtin_marker_program_set_for_panel(
    program_set="tumor_colon",
    genes=panel_genes,
    min_markers=3,
    merge_high_overlap=True,
)

program_report[[
    "name",
    "reporting_role",
    "competition_group",
    "positive_markers_available",
    "marker_coverage_fraction",
    "high_overlap_within_group",
]]
```

Use `collect_positive_marker_genes()` to collect canonical positive markers from a hierarchy and/or built-in flat program sets for quick marker heatmaps.

```python
from hierannot import collect_positive_marker_genes

heatmap_genes = collect_positive_marker_genes(
    hierarchy="breast_tme",
    program_sets=["tumor_breast", "cell_cycle"],
    return_format="list", # either "list" or "dataframe"
)

marker_table = collect_positive_marker_genes(
    hierarchy=roots, # custom hierarchy
    programs=[program_1, program_2], # custom list of MalignantProgram
    return_format="dataframe",
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

## Auxillary track annotations

HierAnnot uses the normal hierarchy to capture **lineage / cell-type identity** and a separate flat program track to capture **tumor/transformation evidence** or auxiliary disease-state programs. These tracks solve different problems. The normal hierarchy is routed through a tree, while the flat program track is scored in parallel and then summarized by a tiered reporting layer.

This keeps four questions separate:

- normal hierarchy result: what lineage or cell type does the cluster resemble?
- status-tier program result: is aggregated tumor-like evidence strong enough to allow integration?
- state/modifier result: which tumor-associated or auxiliary programs are active and specific enough to report?
- integrated/export report: concise labels for downstream use, such as `tumor_like_transformed_emt.Epithelial`, while flag-only auxiliary programs leave `annot_export_label` unchanged.

`HierAnnotResult.malignant_annotations` stores the flat program-track details, while `HierAnnotResult.integrated_annotations` stores the combined normal/tumor interpretation for each cluster. The export summary is intentionally concise so `annot_export_label` can be joined back to single-cell metadata.

Use tumor-oriented program sets with integration when you want tumor-like labels in the final export:

```python
pipeline = HierAnnotPipeline(
    root_programs=roots,
    malignant_programs=get_builtin_marker_program_set("tumor_breast"),
    malignant_integration_mode="integrate",
    malignant_control_gene_exclusion_policy="current_program_only",  # default
)
```

Use non-tumor auxiliary program sets in flag-only mode when normal hierarchy labels should remain the final export labels:

```python
pipeline = HierAnnotPipeline(
    root_programs=roots,
    malignant_programs=get_builtin_marker_program_set("cell_cycle"),
    malignant_integration_mode="flag_only",
    program_report_block_preset="off", # disable blocking to see auxiliary states of all clusters
)

result = pipeline.fit_score(cluster_means)

# Normal hierarchy labels for joining back to cells.
# export["annot_export_label"] remains normal-hierarchy based.
export = make_cluster_annotation_export_summary(result)

# result.malignant_annotations stores the auxiliary program status.
aux_status = result.malignant_annotations[[
    "cluster_id",
    "annot_malignant_status",
    "annot_malignant_label_concise",
    "annot_malignant_status_label",
    "annot_malignant_state_label",
    "annot_malignant_modifier_labels",
    "annot_malignant_tumor_status_pass",
    "annot_malignant_positive_programs",
    "annot_malignant_raw_score",
    "annot_malignant_status_score",
]]
```

### Malignant scoring and tier-based reporting

Malignant scoring separates absolute evidence from state specificity:

- `raw_score`: control-matched malignant-program enrichment, analogous to normal hierarchy raw enrichment.
- `status_score`: support-weighted absolute evidence for each program. In `result.malignant_annotations`, `annot_malignant_status_score` is the combined tumor-status score aggregated from status programs, with a small gated state-support bonus when core status evidence is borderline.
- `program_specificity_score`: within-competition-group separation, useful for exact state labeling.
- `decision_score` / `score`: state-label decision score, not the primary tumor-like status metric.

The malignant reporting layer is tier-based and uses `metadata["reporting_role"]` as active tier metadata:

- `"status"`: establishes tumor-like identity and controls whether a tumor-like label can be integrated with the normal label. Multiple status programs are aggregated by taking the strongest status evidence and adding a small bonus for additional independent strong status groups.
- `"state"`: decorates a tumor-like label when there is sufficient raw evidence and within-group specificity, for example EMT or mesenchymal shift. State programs can provide a small support bonus only when core status evidence is already borderline.
- `"modifier"`: reports auxiliary context such as cycling, IFN-high, hypoxia, or fibrosis; modifiers do not establish tumor-like identity or contribute to the default tumor-status score.

Cross-group co-activation is expected and is reported through fields such as `annot_malignant_positive_programs`, `annot_malignant_positive_groups`, and `annot_malignant_multigroup_positive`. Note that the state and modifier programs can be strong and visible in `result.malignant_annotations`, but they do not establish tumor-like identity for integration unless a status-tier program is also strong. This prevents labels such as `tumor_like_cycling.T cell` or `tumor_like_emt.Fibroblast` when no tumor-status program supports the call.

### Built-in and custom malignant programs

The package provides a neutral built-in registry of flat marker-program sets for direct usage. Tumor sets are intended for tumor-like reporting/integration. Non-tumor sets such as cell-cycle, stress, inflammation, and fibrosis are intended for `malignant_integration_mode="flag_only"`, where they behave as auxiliary status programs and do not alter `annot_export_label`.

```python
from hierannot.builtin import (
    list_builtin_marker_program_sets,
    get_builtin_marker_program_set,
    describe_builtin_marker_program_set,
)

# Summary metadata for available built-in marker program sets 
program_sets = list_builtin_marker_program_sets()
print(program_sets[["name", "program_family", "tissue_scope", "recommended_integration_mode", "recommended_hierarchy"]])

# Details of one particular built-in marker program set
print(describe_builtin_marker_program_set("cell_cycle"))

# Default tumor-focused shared set.
tumor_general = get_builtin_marker_program_set("tumor_general")

# Tissue-type specifc tumor-focused add-ons with or without the shared set.
full_breast_tumor_set = get_builtin_marker_program_set("tumor_breast")
breast_only = get_builtin_marker_program_set("tumor_breast", include_general=False)

# Non-tumor auxiliary sets do not include tumor_general by default.
cell_cycle_programs = get_builtin_marker_program_set("cell_cycle")
stress_programs = get_builtin_marker_program_set("stress")
```

Note that tumor-specific program sets are designed to avoid duplicating normal lineage calls where possible. For example, the `skin_tme` hierarchy should identify whether a cluster is squamous epithelial, keratinocyte-like, or melanocyte. The `tumor_skin` maligant program set then adds tumor-state context, such as squamous activation or melanocytic dedifferentiation, rather than using pure lineage markers as standalone malignant-status evidence. Similarly, lymphoma-context marker sets are provided as `immune_lymphoid_states` for subtype/state annotation in `flag_only` mode; marker enrichment alone is usually not enough to prove lymphoma malignancy without clonality, genotype, copy-number, pathology, or sample-design context.

You can inspect a flat marker-program set in a tree-like text view grouped by reporting role and competition group:

```python
from hierannot import format_marker_program_set, summarize_marker_program_set

programs = get_builtin_marker_program_set("tumor_skin", include_general=False)
print(format_marker_program_set(programs, title="tumor_skin"))

program_summary = summarize_marker_program_set(programs)
program_summary[[
    "name",
    "reporting_role",
    "competition_group",
    "n_positive_markers",
    "report_on_lineages",
]]
```

Custom malignant programs can be authored with the `MalignantProgram` data class and supplied to `HierAnnotPipeline(malignant_programs=...)`. The class name is tumor-focused because tumor reporting is the default supported integration workflow; in `flag_only` mode the same scorer can be used for auxiliary disease/state programs.

```python
from hierannot import MalignantProgram

auxiliary_programs = [
    # Minimal valid custom program: only name + positive_markers are required.
    MalignantProgram(
        name="ifn_response",
        positive_markers=["STAT1", "ISG15", "IFIT1", "IFIT3"],
        description="Interferon-response auxiliary program",
        metadata={"reporting_label": "ifn_high", "reporting_role": "modifier"},
    ),
    MalignantProgram(
        name="proliferation",
        positive_markers=["MKI67", "TOP2A", "UBE2C", "BIRC5"],
        metadata={"reporting_label": "proliferating", "reporting_role": "modifier"},
    ),
]
```
Only `name` and `positive_markers` are required. `negative_markers`, `description`, `tags`, and `metadata` are optional. When `metadata` is omitted, HierAnnot treats the program as an independent competition group and uses `name` as the concise reporting label.

The metadata fields that matter most for custom programs are:

- `competition_group`: affects program specificity scoring. Programs in the same group compete against one another when computing `program_specificity_score`; unrelated programs should usually be left in separate groups. Missing, empty, `"none"`, `"off"`, or `"independent"` means the program forms its own group.
- `reporting_label`: affects concise program labels such as `annot_malignant_label_concise`, mixed-state labels, and tumor-like export labels. It does not change raw/status scoring. Missing values default to the program `name`.
- `reporting_role`: active reporting-tier metadata. Use `"status"` for programs that can establish tumor-like identity, `"state"` for specific tumor-state descriptors, and `"modifier"` for auxiliary flags. 

In the example below, two related EMT/mesenchymal tumor-state programs can share a competition group for state specificity, while an IFN/stress program can remain independent:

```python
from hierannot import MalignantProgram

custom_malignant_programs = [
    # two related EMT/mesenchymal programs share a competition group 
    MalignantProgram(
        name="emt_program_v1",
        positive_markers=["VIM", "FN1", "ZEB1", "SNAI2"],
        negative_markers=["EPCAM", "KRT8"],
        metadata={
            "competition_group": "emt_mesenchymal",
            "reporting_label": "emt",
            "reporting_role": "state",

            # "report_on_lineages" field is used only when program_report_block_preset="lineage_aware".
            # Score is still computed everywhere, but tumor-like reporting is allowed only when the 
            # final normal hierarchy path matches one of these lineages.
            "report_on_lineages": ["Epithelial", "Mammary epithelial", "Basal epithelial"],
        },
    ),
    MalignantProgram(
        name="mesenchymal_shift_v1",
        positive_markers=["COL1A1", "COL1A2", "SPARC", "TAGLN"],
        metadata={
            "competition_group": "emt_mesenchymal",
            "reporting_label": "mesenchymal_shift",
            "reporting_role": "state",
            "report_on_lineages": ["Epithelial", "Mammary epithelial", "Basal epithelial"],
        },
    ),

    # one independent IFN/stress program
    MalignantProgram(
        name="stress_ifn",
        positive_markers=["STAT1", "ISG15", "IFIT1"],
        metadata={"reporting_label": "stress_ifn", "reporting_role": "modifier"},
    ),
]
```

### Malignant reporting controls

Tumor-like reporting is conservative but normal-label preserving. A cluster receives a tumor-like integrated flag only when all of the following are true:

1. aggregated tumor-status evidence passes the configured `malignant_status_score_threshold` and `malignant_raw_score_threshold`. The aggregate is driven by status-role programs; state-role programs can only add a small gated support bonus when core status evidence is already borderline;
2. the normal-track label is not blocked by `program_report_block_preset`;
3. when `malignant_normal_raw_delta_threshold` is not `None`, malignant status-program raw enrichment exceeds normal-track raw evidence by at least that amount.
    - The default `malignant_normal_raw_delta_threshold=None` keeps the raw-delta value as a diagnostic instead of a hard gate. Use `malignant_normal_raw_delta_threshold=0.05` or `0.10` for high-specificity export behavior.

The most relevant audit columns are `annot_malignant_core_status_score`, `annot_malignant_state_support_score`, `annot_malignant_combined_status_score`, `annot_malignant_status_programs`, `annot_malignant_state_support_programs`, and `annot_malignant_status_decision_source`.

The tumor-like reporting control is intentionally blocklist-oriented:

- `program_report_block_preset=None` or `"off"`: no report blocking.
- `program_report_block_preset="tumor_reportable"` (default): report tumor-like labels only on curated tumor-reportable built-in brachnes of normal hierarchies. 
    - The reportable branches in this preset include epithelial/parenchymal branches, unresolved normal calls, and selected non-epithelial tumor lineages such as melanocyte for skin. It also applies immune, stromal, endothelial, and mural guards to avoid common false-positive tumor reporting.
    - This is useful for epithelial tumor workflows where fibroblast, mural, endothelial, and immune labels should not become `tumor_like_emt.*`. 
    - If you use a custom hierarchy with different branch names, pass an explicit list of blocked final normal nodes, such as `["Immune", "Fibroblast", "Mural", "Endothelial"]`, or use `"lineage_aware"` with program-level `metadata["report_on_lineages"]`.
- `program_report_block_preset="immune_like"`: block tumor-like relabeling on immune-like normal calls. 
    - This preset includes immune, stromal, endothelial, and mural guards and is intended to work out-of-the-box for built-in epithelial/parenchymal tumor workflows. 
- `program_report_block_preset="lineage_aware"`: use the selected program's `metadata["report_on_lineages"]` when present. Programs without this metadata remain broadly reportable, which is useful for flag-only auxiliary sets such as cell cycle.
- `program_report_block_preset=[...]`: block exact normal node names or combine preset names, for example `["immune_like", "Fibroblast", "Mural"]`. A cluster is blocked if the final normal label or any exact node in the final normal annotation path matches one of the supplied names (i.e. itself and all its descendants).

Use `malignant_integration_mode="flag_only"` to keep normal labels primary while carrying malignant diagnostics, or `malignant_integration_mode="integrate"` to compute tumor-like integrated labels such as `tumor_like_emt.Epithelial`. 

### Export labels with tumor-like reporting

`make_cluster_annotation_export_summary()` builds an export-ready cluster annotation summary from a HierAnnot result. It inherits and respects the fitted pipeline decisions by default: it uses `result.malignant_annotations` and `result.integrated_annotations` instead of recomputing malignant status or normal/malignant integration with export-time thresholds. Pass `rerun_decision=True` only when you intentionally want to recompute normal routing, malignant status formatting, and integrated labels from the stored score tables under new controls.

The export helper also exposes `export_label_source`:

- `export_label_source="auto"` (default): hybrid behavior. In integrate mode, resolved normal labels become tumor-like integrated labels; mixed and rescued-candidate normal labels stay primary; unknown normal labels can fall back to `tumor_like_<malignant_label>.unknown`.
- `export_label_source="integrated"`: force reportable tumor-like calls to become the primary `annot_export_label`, including mixed/candidate rows.
- `export_label_source="normal_priority"`: keep normal, mixed, or rescued candidate labels as the primary `annot_export_label` whenever available and carry the tumor-like label in `annot_export_tumor_like_label` as separate column.

Normal-side export logic is applied before tumor-like relabeling. For sparse spatial data, the default therefore preserves useful safeguards such as `mixed_*` and `candidate_*` labels while still reporting tumor-like interpretation in `annot_export_tumor_like_label`. Inspect `HierAnnotResult.integrated_annotations` and the export columns `annot_export_malignant_flag` / `annot_export_tumor_like_label` for the full malignant-track interpretation.

## End-to-end workflow 

### Practical workflow with AnnData and custom embeddings

A common workflow is:

1. cluster cells using a custom embedding in ``adata.obsm``
2. aggregate raw expression to gene x cluster means
3. annotate the cluster means with ``HierAnnotPipeline``

```python
from hierannot import HierAnnotPipeline, get_builtin_hierarchy
from hierannot.aggregate import cluster_and_aggregate_anndata

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
from hierannot.aggregate import cluster_anndata_on_representation, aggregate_anndata_to_cluster_means

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

### AnnData round-trip

```python
from hierannot import HierAnnotPipeline, get_builtin_hierarchy
from hierannot.aggregate import cluster_and_aggregate_anndata
from hierannot.workflows import (
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

Use `make_cluster_annotation_export_summary()` to tune unknown/mixed/export labeling at the cluster level before expanding back to cells. The default `export_view="compact"` returns the join-ready label columns plus the most useful normal, malignant, and integrated diagnostics; use `export_view="diagnostic"` to return the full derived table. One could also run it with `rerun_decision=True` and new weights/thresholds to tune the decision routing in addition to the export logic if requested. Then use `expand_cluster_annotation_to_cells()` to build a join-ready per-cell table.

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

## Notes and Tips

- Hierarchical assignment is the primary decision path.
- The best score across all levels is also reported as a diagnostic.
- Scores should be interpreted together with marker coverage and sibling margin.
- Use parent nodes for broad identity programs and child nodes for more specific or state-like refinements.
- ``annot_score`` is the decision score used for label assignment; ``annot_raw_score`` remains available for diagnostics.
- Top-level purity diagnostics such as ``annot_primary_compartment`` and ``annot_mixed_compartment`` help flag clusters that may be mixed or boundary-adjacent.
- When designing custom hierarchies, recommend to keep lineage identity in the hierarchy, while keep transformation-state / stress programs outside the lineage tree when possible. You can put cell state relevant programs as `MalignantProgram` sets and integrate it to the export-ready report using relevant control knobs (e.g. `malignant_integration_mode`).

### Cluster-specific expression profiles for spatial data 

When minor neighborhood contamination is expected, you can aggregate cluster profiles from the cells closest to the cluster centroid in the clustering representation:

```python
from hierannot.aggregate import cluster_and_aggregate_anndata

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

