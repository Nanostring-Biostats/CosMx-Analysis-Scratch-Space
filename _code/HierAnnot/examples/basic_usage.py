import pandas as pd

from hierannot import (
    HierAnnotPipeline,
    get_builtin_hierarchy,
    get_builtin_marker_program_set,
    make_cluster_annotation_export_summary,
)

# cluster-mean expression: genes x clusters pandas.DataFrame
expr = pd.DataFrame(
    {
        "cluster_T": [9, 8, 0, 0, 0, 0, 0, 1, 1],
        "cluster_Myeloid": [0, 0, 8, 7, 0, 0, 0, 6, 5],
        "cluster_Epi": [0, 0, 0, 0, 8, 7, 6, 0, 0],
    },
    index=["CD3D", "TRBC1", "LYZ", "S100A8", "EPCAM", "KRT18", "KRT19", "FCER1G", "CTSS"],
)

pipeline = HierAnnotPipeline(
    root_programs=get_builtin_hierarchy("tme_core"),
    input_type="raw_cluster_means",
    score_threshold=0.01, # default 0.05
    margin_threshold=0.0, # default 0.02
    random_state=0,

    # Tumor-oriented program sets can be integrated into export labels.
    malignant_integration_mode="integrate", # default "flag_only"

    malignant_programs=get_builtin_marker_program_set("tumor_general"),
    # default None to automatically use "general" built-in malignant programs if malignant_integration_mode is not "off"
    
    program_report_block_preset="tumor_reportable",
    # This controls how malignant annotations are reported in the final integrated annotations. 
    # Default = "tumor_reportable" to report tumor-like labels only on epithelial-like ans skin melanocyte final normal calls;
    # Set to "off" to disable any gating and report malignant labels for all clusters regardless of the top-level lineage/compartment;
    # or, set to "immune_like" to only report malignant labels for non-immune clusters in the integrated labels;
    # or, set to "lineage_aware" to apply preset blocklists based on the "report_on_lineages" metadata of each malignant program;
    # or, pass in a list of exact hierarchy node names that you want to block malignant annotation for their descendants (e.g. ["Immune", "Fibroblast", "Mural"]) 

    malignant_control_gene_exclusion_policy="current_program_only",
    # Control-gene exclusion policy used only for malignant programs.
    # Default `"current_program_only"` treats flat malignant programs as independent evidence tracks when sampling matched controls. 
    # Use `"current_branch"` to exclude all markers in the same malignant competition group, or `"all_pipeline_markers"` for the most conservative exclusion on large panels.

)

result = pipeline.fit_score(expr)

print("\n=== Normal hierarchy annotations ===")
print(result.cluster_annotations)

print("\n=== Normal hierarchy scores at all levels ===")
print(result.level_scores)

print("\n=== Malignant scoring summary ===")
print(result.malignant_annotations[[
    # overall malignant detection
    "cluster_id",
    "annot_malignant_status",
    "annot_malignant_label_concise",

    # strong programs under different roles 
    "annot_malignant_status_label",
    "annot_malignant_state_label",
    "annot_malignant_modifier_labels",

    # tumor-like evidence metrics 
    "annot_malignant_tumor_status_pass",
    "annot_malignant_positive_programs",
    "annot_malignant_raw_score",
    "annot_malignant_status_score",
]])

print("\n=== Integrated annotations before applying other export logics ===")
print(result.integrated_annotations[["cluster_id", "annot_integrated_label", "annot_integrated_reason"]])

print("\n=== Export summary ===")
cluster_export = make_cluster_annotation_export_summary(
    result,
    export_view="compact", 
    # default "compact", use "diagnostic" to include all columns of the intermediate results, or "minimal" for the core export columns.

    unknown_on_low_confidence=True, 
    # default True to override low-confidence annotations to "unknown"

    mixed_on_parent_mixing=True, 
    # default True to call "mixed" if top candidates belong to different lineages/compartments and are both strongly supported

    rescue_unknown_with_blocked_candidates = True, 
    # default True to attempt rescuing "unknown" annotations with blocked candidates that have high raw scores

    rerun_decision = False, # default False, if True would rerun the annotation decision and integration logic with new threshold provided to this function.

    # Malignant mode, thresholds, raw-delta gate, and blocklist decisions are
    # inherited from the fitted result unless rerun_decision=True is requested.
    malignant_integration_mode="integrate", # default "flag_only"

    # Other thresholds and controls relevant to malignant integration would take effect when `malignant_integration_mode!="off"` and `rerun_decision=True`.
    normal_strong_score_threshold = 0.35,
    normal_strong_raw_threshold= 0.10,
    malignant_status_score_threshold= 0.35,
    malignant_raw_score_threshold= 0.15,
    malignant_normal_raw_delta_threshold = None,
    program_report_block_preset="tumor_reportable",
)

print(cluster_export[["cluster_id", "annot_export_label", "annot_export_status", "annot_export_malignant_flag"]])

print("\n=== Create a new result object under alternate routing and malignant integration cutoffs ===")
from hierannot import rerun_cluster_annotation_result_from_scores
rerun_result =  rerun_cluster_annotation_result_from_scores(
    result,
    weak_raw_score_threshold= 0.10,
    min_branch_supported_raw_score= 0.0,
    min_leaf_raw_score = -0.2,
    malignant_integration_mode= "flag_only",
    malignant_status_score_threshold= 0.35,
    malignant_raw_score_threshold= 0.15,
    malignant_normal_raw_delta_threshold= 0.05, 
    # default to None to disable the reporting criteria on raw score delta between malignant and normal tracks.

    program_report_block_preset=["Immune", "Fibroblast", "Mural", "Endothelial"], 
    # default "tumor_reportable" to report tumor-like labels only on epithelial-like and selected skin melanocyte final normal calls;
    # here, pass in a list of exact hierarchy node names that you want to block malignant annotation for their descendants.
)


print("\n=== Flag clusters with non-tumor programs ===")
# One could also pass in non-tumor program sets and run in "flag_only" mode 
# to get the program status reported in the malignant_annotations without 
# integrating them into the export labels.
cell_cycle_pipeline = HierAnnotPipeline(
    root_programs=get_builtin_hierarchy("tme_core"),
    input_type="raw_cluster_means",
    score_threshold=0.01, # default 0.05
    margin_threshold=0.0, # default 0.02
    random_state=0,
    malignant_integration_mode="flag_only",
    malignant_programs=get_builtin_marker_program_set("cell_cycle"),
    program_report_block_preset="off", # keep reporting cell-cycle flags for all clusters regardless of their lineage
)

cell_cycle_result = cell_cycle_pipeline.fit_score(expr)
cell_cycle_export = make_cluster_annotation_export_summary(
    cell_cycle_result,
    unknown_on_low_confidence=True,
    mixed_on_parent_mixing=True,
    rescue_unknown_with_blocked_candidates = True, 
)

print("\nCell-cycle auxiliary flag-only scoring")
print(cell_cycle_result.malignant_annotations[[
    "cluster_id",
    "annot_malignant_status",
    "annot_malignant_tumor_status_pass",
    "annot_malignant_modifier_labels",
    "annot_malignant_label_concise",
]])

print(cell_cycle_export[["cluster_id", "annot_export_label", "annot_export_status", "annot_export_malignant_flag"]])

