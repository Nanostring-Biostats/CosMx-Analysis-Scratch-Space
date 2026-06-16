"""
For tissue-type detection, use `detect_tissue_type()`. This helper runs the built-in `tissue_origin_screen` hierarchy, 
which is a shallow standard `MarkerProgram` hierarchy designed specifically for tissue-origin screening. It contains 
broad tissue-defining anchors such as mammary epithelial, intestinal epithelial, renal parenchymal, cutaneous 
squamous/keratinocyte, melanocyte, neural/glial, and other major tissue-specific cell types in one shared competition 
space. Shared immune, stromal, endothelial, and mural controls absorb generic TME evidence but do not count as a 
specific tissue call.
"""

from hierannot import (
    list_builtin_hierarchies, 
    get_builtin_hierarchy, 
    format_hierarchy_tree,
)

# description of `tissue_origin_screen` hierarchy,
df = list_builtin_hierarchies()
df.query("name == 'tissue_origin_screen'")

root_programs = get_builtin_hierarchy('tissue_origin_screen')
print(format_hierarchy_tree(root_programs))

# get cluster-mean profiles and cell number 
# cluster_means is a gene x cluster expression matrix.

from hierannot import (
    aggregate_anndata_to_cluster_means, 
    detect_tissue_type,
)

cluster_key = "leiden"
cluster_means = aggregate_anndata_to_cluster_means(
    adata,
    cluster_key=cluster_key,
    source="layer",
    source_key="counts",
    method="mean",
    uppercase_genes=True,
)

cluster_sizes = adata.obs[cluster_key].value_counts()

# Detect tissue type with built-in 'tissue_origin_screen' hierarchy 
tissue_summary = detect_tissue_type(
    cluster_means,
    
    hierarchy="tissue_origin_screen", 
    # default to use built-in hierarchy "tissue_origin_screen", 
    # optional to pass in a list of custom `MarkerProgram` that has tissue-detection metadata 

    compile_for_panel=True, 
    # whether to compile the tissue-origin screen hierarchy against input gene panels, 
    # default to True and will drop tissue-specific anchor cell types with too fewer requested
    # marker genes (< 3 markers or < 0.4x of all ) from tissue-type candidates. 
    
    cluster_weights=cluster_sizes, 
    # default None, not to weight by cluster size; 
    # or to weight tissue detection score by the cluster size, pass in pd.Series or 1-column 
    # pd.DataFrame for cell number in each cluster with cluster_id as index. 
)

"""
# example custom tissue-origin screen hierarchy 
from hierannot import MarkerProgram, detect_tissue_type

# recommend not to include children to anchor nodes
custom_origin_screen = [
    MarkerProgram(
        name="Custom epithelial anchor",
        positive_markers=["EPCAM", "KRT8", "KRT18"],
        negative_markers=["PTPRC", "COL1A1"],
        metadata={
            "tissue_detection_role": "anchor",
            "tissue_candidate": "custom_epithelial_tissue",
            "recommended_hierarchy": "tme_core",  # or your custom downstream hierarchy name
        },
    ),
    MarkerProgram(
        name="Shared immune control",
        positive_markers=["PTPRC", "CD3D", "CD3E"],
        metadata={
            "tissue_detection_role": "shared_control",
            "tissue_candidate": "generic_tme",
            "recommended_hierarchy": "tme_core",
        },
    ),
]

tissue_summary = detect_tissue_type(
    cluster_means,
    hierarchy=custom_origin_screen,
)

"""

# dataset-level call 
tissue_summary[[
    "detected_tissue_type",
    "tissue_detection_call",
    "recommended_hierarchy",
]].head(1)

"""
detected_tissue_type
    The inferred tissue type, or a non-specific call such as
    "generic_tme", "ambiguous", or "insufficient_tissue_specific_evidence".

tissue_detection_call
    "specific_tissue", "generic_tme", "ambiguous", or
    "insufficient_tissue_specific_evidence".

recommended_hierarchy
    The hierarchy to use downstream when the dataset-level call is specific or generic.


Meaning of `tissue_detection_call` values:
- "specific_tissue": sucessful identified a specific tissue type. 
- "ambiguous": the second-best tissue score is at least `ambiguity_score_ratio` (default = 0.90) of the best call. 
- "generic_tme": no specific tissue anchor has enough direct evidence but the shared TME controls are supported.
- "insufficient_tissue_specific_evidence": neither specific anchors nor generic TME controls have enough evidence

"""

# anchor-level results for each major cell types 
tissue_summary[[
    "tissue_candidate",
    "anchor_label",
    "candidate_recommended_hierarchy",
    "tissue_detection_rank",
    "tissue_detection_score",
    "n_valid_anchor_clusters",
    "mean_anchor_branch_supported_raw_score",
    
    # extra details 
    "fraction_valid_anchor_clusters",
    "valid_anchor_cluster_weight",
    "fraction_valid_anchor_cluster_weight",
    "mean_anchor_score", 
    "median_anchor_margin"
]]


"""
By default, a cluster contributes tissue-origin evidence only when its final screen label passes export-like validity filtering. 
The evidence gate uses the same core columns as the export unknown mask: `annot_branch_supported_raw_score` is compared with 
`unknown_branch_raw_threshold`, and `annot_score` is compared with `unknown_min_score` when that optional threshold is provided. 
Medium- and high-confidence labels can contribute; low and none confidence labels are ignored. Margin is reported as a diagnostic 
but is not used as a hard default gate, because closely related tissue anchors can have low margins while still supporting the 
same tissue source.
"""

# Choosing a downstream hierarchy after tissue detection
import pandas as pd
from hierannot import (
    HierAnnotPipeline,
    list_builtin_marker_program_sets, 
    get_builtin_marker_program_set,
    save_result_bundle, 
)

# convinent wrapper to get best matched hierarhcy and tumor names from tissue_summary outcome
def get_best_setup(tissue_summary):
    call = tissue_summary.iloc[0]["tissue_detection_call"]
    recommended_hierarchy = tissue_summary.iloc[0]["recommended_hierarchy"]
    detected_tissue = tissue_summary.iloc[0]["detected_tissue_type"]

    # get normal hierarhcy 
    if call == "generic_tme":
        recommended_hierarchy="tme_core"
    elif call != "specific_tissue":
        print(
            "No confident built-in tissue hierarchy was detected. "
            "Fall back to tme_core. "
        )
        recommended_hierarchy="tme_core"
    
    # Note that squamous_epithelial would support both tonsil_tme and skin_tme the same, 
    # internally prefer skin_tme due to higher coverage in diverse cell types

    # auto-pick tissue-specific tumor program set 
    program_sets = list_builtin_marker_program_sets()
    all_aux_program_names = program_sets['name'].values
    tumor_programs_name = "tumor_" + recommended_hierarchy.removesuffix("_tme")
    if not tumor_programs_name in all_aux_program_names:
        tumor_programs_name = "tumor_general" 
    
    print(f"{call}, recommend use {recommended_hierarchy} and {tumor_programs_name}")
    df = pd.DataFrame([{
        "tissue_detection_call": call, 
        "detected_tissue_type": detected_tissue,
        "recommended_hierarchy": recommended_hierarchy,
        "tumor_programs_name": tumor_programs_name
    }])
    return df 

df = get_best_setup(tissue_summary)
detected_tissue= df.iloc[0]["detected_tissue_type"]
tme_name = df.iloc[0]["recommended_hierarchy"]
tumor_programs_name = df.iloc[0]["tumor_programs_name"]
print(f"Tissue-type detected = {detected_tissue}, use hierarhcy = {tme_name}, tumor = {tumor_programs_name}")

# get normal hierarhcy 
root_programs = get_builtin_hierarchy(tme_name)

# also pick a matching tumor program set
tumor_programs= get_builtin_marker_program_set(tumor_programs_name)

# setup tissue-specific pipeline with integrated tumor annotation 
pipeline = HierAnnotPipeline(
    root_programs=root_programs,
    malignant_integration_mode="integrate",
    malignant_programs=tumor_programs,
    program_report_block_preset="tumor_reportable",
)

result = pipeline.fit_score(cluster_means)

save_result_bundle(
    result, 
    path="hierannot_result_bundle", 
    metadata={
        "tissue_type": detected_tissue, 
        "hierarchy_name": tme_name, 
        "tumor_programs_name": tumor_programs_name
    },
)
