from __future__ import annotations

import pandas as pd

try:
    from ..datamodels import MalignantProgram
except Exception:  # pragma: no cover - import-time safety for docs/builds
    MalignantProgram = None

_TUMOR_REPORT_ON_LINEAGES = [
    "Epithelial",
    "Mammary epithelial",
    "Basal epithelial",
    "Luminal epithelial",
    "Secretory/alveolar epithelial",
    "Intestinal epithelial",
    "Colon epithelial",
    "Absorptive-like",
    "Enterocyte",
    "Goblet-like",
    "Goblet cell",
    "Stem/TA-like",
    "Renal parenchymal",
    "Proximal tubule",
    "Distal/TAL-like",
    "Collecting duct",
    "Podocyte",
    "Hepatic parenchymal",
    "Hepatocyte-like",
    "Cholangiocyte-like",
    "Pancreatic parenchymal",
    "Ductal",
    "Acinar",
    "Endocrine",
    "Pulmonary parenchymal",
    "Alveolar epithelial",
    "Airway epithelial",
    "Cutaneous epithelial",
    "Basal keratinocyte",
    "Suprabasal keratinocyte",
    "Granular/cornified keratinocyte",
    "Differentiated keratinocyte",
    "Squamous epithelial",
    "Basal squamous epithelial",
    "Suprabasal squamous epithelial",
    "Keratinizing squamous epithelial",
    "Activated/stress squamous epithelial",
    "Hair follicle epithelial",
    "Sebaceous/ductal epithelial",
    "Melanocyte",
    "Tumor epithelial",
]


_TUMOR_GENERAL = [
    {
        "name": "core_transformation",
        "positive_markers": ["SOX9", "KRT17", "TACSTD2", "CLDN4", "MSLN", "S100A14", "LCN2", "CEACAM5", "PROM1"],
        "negative_markers": ["PTPRC"],
        "description": "Core epithelial transformation status module. This avoids using broad epithelial lineage markers as the main tumor-status evidence.",
        "metadata": {"competition_group": "core_transformation", "reporting_label": "transformed", "reporting_role": "status", "program_family": "tumor", "merged_from": ["transformed", "dedifferentiated"], "report_on_lineages": _TUMOR_REPORT_ON_LINEAGES},
    },
    {
        "name": "epithelial_tumor_stress",
        "positive_markers": ["KRT17", "KRT6A", "KRT6B", "LCN2", "S100A8", "S100A9", "SERPINB3", "SERPINB4", "MMP7", "CLDN4"],
        "negative_markers": ["PTPRC"],
        "description": "Epithelial tumor/stress status module that supports tumor-like detection when core transformation markers are borderline.",
        "metadata": {"competition_group": "tumor_epithelial_stress", "reporting_label": "tumor_stress", "reporting_role": "status", "program_family": "tumor", "report_on_lineages": _TUMOR_REPORT_ON_LINEAGES},
    },
    {
        "name": "stemlike_transformation",
        "positive_markers": ["PROM1", "LGR5", "ASCL2", "OLFM4", "SOX9", "CD44", "ALCAM", "MSI1", "SMOC2", "MYC"],
        "negative_markers": ["PTPRC"],
        "description": "Stem-like transformation status module used as additional tumor-status evidence in epithelial/parenchymal contexts.",
        "metadata": {"competition_group": "stemlike_transformation", "reporting_label": "stemlike", "reporting_role": "status", "program_family": "tumor", "report_on_lineages": _TUMOR_REPORT_ON_LINEAGES},
    },
    {
        "name": "emt_mesenchymal_transition",
        "positive_markers": ["VIM", "FN1", "ITGA5", "ITGB1", "ITGA6", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "LAMB3", "MMP7", "MMP14", "KRT17"],
        "negative_markers": ["PTPRC"],
        "description": "EMT/mesenchymal-transition tumor state module. ",
        "metadata": {"competition_group": "emt_mesenchymal", "reporting_label": "emt_mesenchymal", "reporting_role": "state", "program_family": "tumor", "merged_from": ["emt", "mesenchymal_shift"], "report_on_lineages": _TUMOR_REPORT_ON_LINEAGES},
    },
    {
        "name": "cycling",
        "positive_markers": ["MKI67", "TOP2A", "UBE2C", "PCNA", "TYMS", "BIRC5"],
        "negative_markers": [],
        "description": "Cycling/proliferative tumor-associated modifier. It is reported as context and does not establish tumor-like status by itself.",
        "metadata": {"competition_group": "proliferation", "reporting_label": "cycling", "reporting_role": "modifier", "program_family": "tumor", "report_on_lineages": _TUMOR_REPORT_ON_LINEAGES},
    },
    {
        "name": "stress_ifn",
        "positive_markers": ["STAT1", "ISG15", "IFIT1", "IFIT3", "MX1", "OAS1"],
        "negative_markers": [],
        "description": "Interferon/stress-associated tumor modifier. It is reported as context and does not establish tumor-like status by itself.",
        "metadata": {"competition_group": "stress_ifn", "reporting_label": "stress_ifn", "reporting_role": "modifier", "program_family": "tumor", "report_on_lineages": _TUMOR_REPORT_ON_LINEAGES},
    },
]

_TUMOR_BREAST = [
    {
        "name": "basal_program",
        "positive_markers": ["KRT5", "KRT14", "KRT17", "TP63", "EGFR"],
        "negative_markers": ["PTPRC"],
        "description": "Breast tumor basal-like state descriptor.",
        "metadata": {"competition_group": "breast_subtype", "reporting_label": "basal", "reporting_role": "state", "program_family": "tumor", "report_on_lineages": _TUMOR_REPORT_ON_LINEAGES},
    },
    {
        "name": "luminal_program",
        "positive_markers": ["EPCAM", "KRT8", "KRT18", "KRT19", "MUC1"],
        "negative_markers": ["PTPRC"],
        "description": "Breast tumor luminal-like state descriptor.",
        "metadata": {"competition_group": "breast_subtype", "reporting_label": "luminal", "reporting_role": "state", "program_family": "tumor", "report_on_lineages": _TUMOR_REPORT_ON_LINEAGES},
    },
]

_TUMOR_COLON = [
    {
        "name": "wnt_stem_like",
        "positive_markers": ["LGR5", "ASCL2", "SOX9", "AXIN2", "MYC"],
        "negative_markers": ["PTPRC"],
        "description": "Colon tumor WNT/stem-like state descriptor.",
        "metadata": {"competition_group": "colon_subtype", "reporting_label": "wnt_stem_like", "reporting_role": "state", "program_family": "tumor", "report_on_lineages": _TUMOR_REPORT_ON_LINEAGES},
    },
    {
        "name": "goblet_mucinous_like",
        "positive_markers": ["MUC2", "SPINK4", "TFF3", "AGR2", "CLCA1"],
        "negative_markers": ["PTPRC"],
        "description": "Colon tumor goblet/mucinous-like state descriptor.",
        "metadata": {"competition_group": "colon_subtype", "reporting_label": "goblet_mucinous_like", "reporting_role": "state", "program_family": "tumor", "report_on_lineages": _TUMOR_REPORT_ON_LINEAGES},
    },
]

_SQUAMOUS_REPORT_ON_LINEAGES = [
    "Cutaneous epithelial",
    "Squamous epithelial",
    "Basal keratinocyte",
    "Suprabasal keratinocyte",
    "Granular/cornified keratinocyte",
    "Basal squamous epithelial",
    "Suprabasal squamous epithelial",
    "Keratinizing squamous epithelial",
    "Activated/stress squamous epithelial",
    "Hair follicle epithelial",
    "Sebaceous/ductal epithelial",
]

_MELANOCYTIC_REPORT_ON_LINEAGES = ["Melanocyte"]

_LYMPHOID_REPORT_ON_LINEAGES = [
    "Immune",
    "B cell",
    "Plasma cell",
    "T cell",
    "CD4 T cell",
    "CD8 T cell",
    "Treg",
    "NK cell",
]

_TUMOR_SQUAMOUS = [
    {
        "name": "squamous_activation_state",
        "positive_markers": ["KRT6A", "KRT16", "KRT17", "S100A7", "SERPINB3", "SERPINB4", "LAMC2", "MMP10"],
        "negative_markers": ["PTPRC", "MLANA"],
        "description": "Squamous epithelial activation/invasive-state descriptor. It is a state program used after squamous lineage is established by the normal hierarchy, not a standalone malignant-status call.",
        "metadata": {"competition_group": "squamous_state", "reporting_label": "squamous_activation", "reporting_role": "state", "program_family": "tumor", "tissue_scope": "squamous", "report_on_lineages": _SQUAMOUS_REPORT_ON_LINEAGES},
    },
    {
        "name": "basal_squamous_state",
        "positive_markers": ["TP63", "KRT5", "KRT14", "ITGA6", "COL17A1", "DST"],
        "negative_markers": ["MLANA", "PTPRC"],
        "description": "Basal squamous/keratinocyte state descriptor. These are lineage-rich markers and should decorate a tumor-like call only when normal hierarchy and status evidence support the context.",
        "metadata": {"competition_group": "squamous_state", "reporting_label": "basal_squamous", "reporting_role": "state", "program_family": "tumor", "tissue_scope": "squamous", "lineage_rich": True, "report_on_lineages": _SQUAMOUS_REPORT_ON_LINEAGES},
    },
    {
        "name": "keratinizing_squamous_state",
        "positive_markers": ["KRT1", "KRT10", "IVL", "DSG1", "DSC1", "FLG", "LOR", "TGM1"],
        "negative_markers": ["KRT14", "MLANA", "PTPRC"],
        "description": "Keratinizing/differentiated squamous state descriptor for skin and tonsil/oropharyngeal squamous contexts.",
        "metadata": {"competition_group": "squamous_state", "reporting_label": "keratinizing_squamous", "reporting_role": "state", "program_family": "tumor", "tissue_scope": "squamous", "lineage_rich": True, "report_on_lineages": _SQUAMOUS_REPORT_ON_LINEAGES},
    },
]

_TUMOR_SKIN = _TUMOR_SQUAMOUS + [
    {
        "name": "melanoma_dedifferentiation_state",
        "positive_markers": ["AXL", "WNT5A", "NGFR", "VIM", "FN1", "ZEB1", "JUN", "FOSL1"],
        "negative_markers": ["PTPRC", "KRT14", "COL1A1"],
        "description": "Melanoma/neural-crest-like dedifferentiation or invasive-state descriptor. Melanocyte identity is handled by the normal skin hierarchy; this program does not establish tumor-like status by itself.",
        "metadata": {"competition_group": "melanocytic_state", "reporting_label": "melanoma_dedifferentiation", "reporting_role": "state", "program_family": "tumor", "tissue_scope": "skin", "report_on_lineages": _MELANOCYTIC_REPORT_ON_LINEAGES},
    },

]

_TUMOR_TONSIL = list(_TUMOR_SQUAMOUS)

_IMMUNE_LYMPHOID_STATES = [
    {
        "name": "germinal_center_bcell_state",
        "positive_markers": ["BCL6", "AICDA", "MEF2B", "LMO2", "CD83", "CD79A", "CD79B"],
        "negative_markers": ["CD3D", "LYZ"],
        "description": "Germinal-center B-cell state program for lymphoid/lymphoma-context flag-only annotation; it is not a marker-only malignant-status detector.",
        "metadata": {"competition_group": "bcell_state", "reporting_label": "germinal_center_bcell", "reporting_role": "state", "program_family": "immune_state", "tissue_scope": "lymphoid", "report_on_lineages": ["Immune", "B cell"]},
    },
    {
        "name": "activated_bcell_state",
        "positive_markers": ["IRF4", "PRDM1", "XBP1", "CD44", "NFKBIA", "TNFAIP3"],
        "negative_markers": ["CD3D", "LYZ"],
        "description": "Activated B-cell/plasmacytic state program for lymphoid-state annotation.",
        "metadata": {"competition_group": "bcell_state", "reporting_label": "activated_bcell", "reporting_role": "state", "program_family": "immune_state", "tissue_scope": "lymphoid", "report_on_lineages": ["Immune", "B cell", "Plasma cell"]},
    },
    {
        "name": "plasmablast_plasma_cell_state",
        "positive_markers": ["XBP1", "PRDM1", "SDC1", "MZB1", "JCHAIN", "IGHG1", "IGKC"],
        "negative_markers": ["CD3D", "LYZ"],
        "description": "Plasmablast/plasma-cell state program for lymphoid-state annotation.",
        "metadata": {"competition_group": "bcell_state", "reporting_label": "plasmablast_plasma", "reporting_role": "state", "program_family": "immune_state", "tissue_scope": "lymphoid", "report_on_lineages": ["Immune", "B cell", "Plasma cell"]},
    },
    {
        "name": "cytotoxic_t_nk_state",
        "positive_markers": ["NKG7", "GNLY", "GZMB", "PRF1", "CTSW", "KLRD1"],
        "negative_markers": ["MS4A1", "CD79A"],
        "description": "Cytotoxic T/NK lymphoid state program.",
        "metadata": {"competition_group": "t_nk_state", "reporting_label": "cytotoxic_t_nk", "reporting_role": "state", "program_family": "immune_state", "tissue_scope": "lymphoid", "report_on_lineages": ["Immune", "T cell", "CD8 T cell", "NK cell"]},
    },
    {
        "name": "exhausted_tcell_state",
        "positive_markers": ["PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX", "CTLA4"],
        "negative_markers": ["MS4A1", "CD79A"],
        "description": "Exhausted T-cell-like state program.",
        "metadata": {"competition_group": "t_nk_state", "reporting_label": "exhausted_tcell", "reporting_role": "state", "program_family": "immune_state", "tissue_scope": "lymphoid", "report_on_lineages": ["Immune", "T cell", "CD4 T cell", "CD8 T cell"]},
    },
    {
        "name": "treg_state",
        "positive_markers": ["FOXP3", "IL2RA", "CTLA4", "IKZF2", "TIGIT"],
        "negative_markers": ["MS4A1", "NKG7"],
        "description": "Regulatory T-cell-like state program.",
        "metadata": {"competition_group": "t_nk_state", "reporting_label": "treg_like", "reporting_role": "state", "program_family": "immune_state", "tissue_scope": "lymphoid", "report_on_lineages": ["Immune", "T cell", "CD4 T cell", "Treg"]},
    },
    {
        "name": "cycling_lymphoid_modifier",
        "positive_markers": ["MKI67", "TOP2A", "UBE2C", "BIRC5", "CENPF", "CCNB1"],
        "negative_markers": [],
        "description": "Cycling/proliferative lymphoid modifier. It is useful in lymphoma-context studies but does not establish malignancy by itself.",
        "metadata": {"competition_group": "lymphoid_proliferation", "reporting_label": "cycling_lymphoid", "reporting_role": "modifier", "program_family": "immune_state", "tissue_scope": "lymphoid", "report_on_lineages": _LYMPHOID_REPORT_ON_LINEAGES},
    },
    {
        "name": "ifn_response_lymphoid_modifier",
        "positive_markers": ["STAT1", "ISG15", "IFIT1", "IFIT3", "MX1", "OAS1", "IRF7"],
        "negative_markers": [],
        "description": "Interferon-response lymphoid modifier for flag-only immune-state annotation.",
        "metadata": {"competition_group": "lymphoid_ifn_response", "reporting_label": "ifn_high", "reporting_role": "modifier", "program_family": "immune_state", "tissue_scope": "lymphoid", "report_on_lineages": _LYMPHOID_REPORT_ON_LINEAGES},
    },
]


_CELL_CYCLE = [
    {
        "name": "proliferation",
        "positive_markers": ["MKI67", "TOP2A", "PCNA", "TYMS", "MCM2", "MCM5", "STMN1"],
        "negative_markers": [],
        "description": "General proliferation/cell-cycle activity program.",
        "metadata": {"competition_group": "cell_cycle_activity", "reporting_label": "proliferating", "reporting_role": "modifier", "program_family": "cell_cycle"},
    },
    {
        "name": "g1_s",
        "positive_markers": ["MCM2", "MCM5", "PCNA", "TYMS", "RRM2"],
        "negative_markers": [],
        "description": "G1/S-phase cell-cycle program.",
        "metadata": {"competition_group": "cell_cycle_phase", "reporting_label": "g1_s", "reporting_role": "modifier", "program_family": "cell_cycle"},
    },
    {
        "name": "g2_m",
        "positive_markers": ["TOP2A", "UBE2C", "BIRC5", "CDC20", "CCNB1", "CENPF"],
        "negative_markers": [],
        "description": "G2/M and mitotic cell-cycle program.",
        "metadata": {"competition_group": "cell_cycle_phase", "reporting_label": "g2_m", "reporting_role": "modifier", "program_family": "cell_cycle"},
    },
]

_STRESS = [
    {
        "name": "hypoxia",
        "positive_markers": ["CA9", "VEGFA", "SLC2A1", "LDHA", "PGK1", "ENO1"],
        "negative_markers": [],
        "description": "Hypoxia/glycolytic stress program.",
        "metadata": {"competition_group": "hypoxia", "reporting_label": "hypoxic", "reporting_role": "modifier", "program_family": "stress"},
    },
    {
        "name": "ifn_response",
        "positive_markers": ["STAT1", "ISG15", "IFIT1", "IFIT3", "MX1", "OAS1", "IRF7"],
        "negative_markers": [],
        "description": "Interferon-response program.",
        "metadata": {"competition_group": "ifn_response", "reporting_label": "ifn_high", "reporting_role": "modifier", "program_family": "stress"},
    },
    {
        "name": "unfolded_protein_response",
        "positive_markers": ["XBP1", "HSPA5", "DDIT3", "ATF3", "ATF4"],
        "negative_markers": [],
        "description": "Unfolded-protein/stress-response program.",
        "metadata": {"competition_group": "upr", "reporting_label": "upr_high", "reporting_role": "modifier", "program_family": "stress"},
    },
]

_INFLAMMATION = [
    {
        "name": "antigen_presentation",
        "positive_markers": ["HLA-DRA", "HLA-DRB1", "CD74", "B2M", "TAP1"],
        "negative_markers": [],
        "description": "Antigen-presentation activity program.",
        "metadata": {"competition_group": "antigen_presentation", "reporting_label": "antigen_presenting", "reporting_role": "modifier", "program_family": "inflammation"},
    },
    {
        "name": "nfkb_inflammation",
        "positive_markers": ["NFKBIA", "TNFAIP3", "IL1B", "CXCL8", "CCL2"],
        "negative_markers": [],
        "description": "NF-kB/inflammatory activation program.",
        "metadata": {"competition_group": "inflammatory_activation", "reporting_label": "inflammatory", "reporting_role": "modifier", "program_family": "inflammation"},
    },
]

_FIBROSIS = [
    {
        "name": "fibroblast_activation",
        "positive_markers": ["COL1A1", "COL1A2", "COL3A1", "FN1", "POSTN", "THY1"],
        "negative_markers": [],
        "description": "Fibroblast/stromal activation program.",
        "metadata": {"competition_group": "fibrosis", "reporting_label": "fibroblast_activation", "reporting_role": "modifier", "program_family": "fibrosis", "report_on_lineages": ["Fibroblast", "Mural", "Stromal", "Pericyte", "Smooth muscle"]},
    },
    {
        "name": "myofibroblast",
        "positive_markers": ["ACTA2", "TAGLN", "MYL9", "COL1A1", "FN1"],
        "negative_markers": [],
        "description": "Myofibroblast-like contractile stromal program.",
        "metadata": {"competition_group": "fibrosis", "reporting_label": "myofibroblast", "reporting_role": "modifier", "program_family": "fibrosis", "report_on_lineages": ["Fibroblast", "Mural", "Stromal", "Pericyte", "Smooth muscle"]},
    },
    {
        "name": "ecm_remodeling",
        "positive_markers": ["MMP2", "MMP11", "LOX", "LUM", "DCN", "POSTN"],
        "negative_markers": [],
        "description": "Extracellular-matrix remodeling program.",
        "metadata": {"competition_group": "fibrosis", "reporting_label": "ecm_remodeling", "reporting_role": "modifier", "program_family": "fibrosis", "report_on_lineages": ["Fibroblast", "Mural", "Stromal", "Pericyte", "Smooth muscle"]},
    },
]

_SETS = {
    "tumor_general": _TUMOR_GENERAL,
    "tumor_breast": _TUMOR_BREAST,
    "tumor_colon": _TUMOR_COLON,
    "tumor_squamous": _TUMOR_SQUAMOUS,
    "tumor_skin": _TUMOR_SKIN,
    "tumor_tonsil": _TUMOR_TONSIL,
    "immune_lymphoid_states": _IMMUNE_LYMPHOID_STATES,
    "cell_cycle": _CELL_CYCLE,
    "stress": _STRESS,
    "inflammation": _INFLAMMATION,
    "fibrosis": _FIBROSIS,
}


_SET_METADATA = {
    "tumor_general": {"program_family": "tumor", "recommended_integration_mode": "integrate", "tissue_scope": "pan-solid-tumor", "recommended_hierarchy": None, "recommended_report_block_preset": "tumor_reportable", "description": "Shared solid-tumor status, state, and modifier programs. Status programs establish tumor-like support; state/modifier programs add context."},
    "tumor_breast": {"program_family": "tumor", "recommended_integration_mode": "integrate", "tissue_scope": "breast", "recommended_hierarchy": "breast_tme", "recommended_report_block_preset": "tumor_reportable", "description": "Breast-focused tumor-state marker programs, plus tumor_general by default."},
    "tumor_colon": {"program_family": "tumor", "recommended_integration_mode": "integrate", "tissue_scope": "colon", "recommended_hierarchy": "colon_tme", "recommended_report_block_preset": "tumor_reportable", "description": "Colon-focused tumor-state marker programs, plus tumor_general by default."},
    "tumor_squamous": {"program_family": "tumor", "recommended_integration_mode": "integrate", "tissue_scope": "squamous", "recommended_hierarchy": "skin_tme or tonsil_tme", "recommended_report_block_preset": "tumor_reportable", "description": "Reusable squamous tumor-state programs. These decorate solid-tumor calls after squamous lineage is established by the normal hierarchy."},
    "tumor_skin": {"program_family": "tumor", "recommended_integration_mode": "integrate", "tissue_scope": "skin", "recommended_hierarchy": "skin_tme", "recommended_report_block_preset": "tumor_reportable", "description": "Skin-focused tumor-state programs: reusable squamous states plus melanocytic state descriptors, with tumor_general added by default."},
    "tumor_tonsil": {"program_family": "tumor", "recommended_integration_mode": "integrate", "tissue_scope": "tonsil/oropharyngeal squamous", "recommended_hierarchy": "tonsil_tme", "recommended_report_block_preset": "tumor_reportable", "description": "Tonsil/oropharynx squamous tumor-state programs, plus tumor_general by default."},
    "immune_lymphoid_states": {"program_family": "immune_state", "recommended_integration_mode": "flag_only", "tissue_scope": "lymphoid/lymphoma-context", "recommended_hierarchy": "immune_core or tme_core", "recommended_report_block_preset": "lineage_aware", "description": "Lymphoid immune subtype/state programs for lymphoma-context flag-only annotation. These do not establish malignant status by themselves."},
    "cell_cycle": {"program_family": "cell_cycle", "recommended_integration_mode": "flag_only", "tissue_scope": "pan-tissue", "recommended_hierarchy": None, "recommended_report_block_preset": "off", "description": "Cell-cycle/proliferation auxiliary programs for flag-only scoring."},
    "stress": {"program_family": "stress", "recommended_integration_mode": "flag_only", "tissue_scope": "pan-tissue", "recommended_hierarchy": None, "recommended_report_block_preset": "off", "description": "Stress, hypoxia, interferon, and UPR auxiliary programs for flag-only scoring."},
    "inflammation": {"program_family": "inflammation", "recommended_integration_mode": "flag_only", "tissue_scope": "pan-tissue", "recommended_hierarchy": None, "recommended_report_block_preset": "off", "description": "Inflammatory and antigen-presentation auxiliary programs for flag-only scoring."},
    "fibrosis": {"program_family": "fibrosis", "recommended_integration_mode": "flag_only", "tissue_scope": "stromal", "recommended_hierarchy": None, "recommended_report_block_preset": "lineage_aware", "description": "Fibrosis and stromal remodeling auxiliary programs."},
}



def _canonical_set_name(name: str) -> str:
    key = str(name).strip()
    if key not in _SETS:
        available = ", ".join(sorted(_SETS))
        raise KeyError(f"Unknown marker program set: {name!r}. Available sets: {available}")
    return key


def _normalize_program_entry(item):
    if MalignantProgram is not None and isinstance(item, MalignantProgram):
        return {
            "name": item.name,
            "positive_markers": list(item.positive_markers),
            "negative_markers": list(item.negative_markers),
            "description": getattr(item, "description", None),
            "tags": list(getattr(item, "tags", [])),
            "metadata": dict(getattr(item, "metadata", {}) or {}),
        }
    if isinstance(item, dict):
        return {
            "name": item.get("name", item.get("label")),
            "positive_markers": list(item.get("positive_markers", [])),
            "negative_markers": list(item.get("negative_markers", [])),
            "description": item.get("description"),
            "tags": list(item.get("tags", [])),
            "metadata": dict(item.get("metadata", {}) or {}),
        }
    raise TypeError(f"Unsupported marker program entry: {type(item)!r}")


def _program_roles_for_set(name: str) -> dict[str, int]:
    roles = {"status": 0, "state": 0, "modifier": 0}
    for item in _SETS[name]:
        role = str((item.get("metadata") or {}).get("reporting_role", "")).lower()
        if role in roles:
            roles[role] += 1
    return roles


def list_builtin_marker_program_sets() -> pd.DataFrame:
    """Return available flat marker-program sets as a metadata table.

    Tumor sets use explicit ``tumor_*`` names such as ``tumor_general``,
    ``tumor_breast``, ``tumor_colon`` and ``tumor_skin``. Non-tumor sets are
    intended for flag-only auxiliary scoring so normal hierarchy export labels
    remain unchanged. The returned DataFrame is designed to mirror
    ``list_builtin_hierarchies()`` and help users choose an appropriate program
    set.
    """
    rows = []
    for name in sorted(_SETS):
        meta = dict(_SET_METADATA.get(name, {}))
        roles = _program_roles_for_set(name)
        rows.append({
            "name": name,
            "program_family": meta.get("program_family"),
            "tissue_scope": meta.get("tissue_scope"),
            "recommended_hierarchy": meta.get("recommended_hierarchy"),
            "recommended_integration_mode": meta.get("recommended_integration_mode"),
            "recommended_report_block_preset": meta.get("recommended_report_block_preset"),
            "description": meta.get("description"),
            "n_programs": len(_SETS[name]),
            "n_status_programs": roles["status"],
            "n_state_programs": roles["state"],
            "n_modifier_programs": roles["modifier"],
        })
    return pd.DataFrame(rows).sort_values(["program_family", "name"]).reset_index(drop=True)


def describe_builtin_marker_program_set(name: str):
    """Return metadata for a built-in flat marker-program set."""
    canonical = _canonical_set_name(name)
    out = dict(_SET_METADATA.get(canonical, {}))
    out["name"] = canonical
    out["n_programs"] = len(_SETS[canonical])
    return out


def get_builtin_marker_program_set(name: str, include_general: bool | None = None):
    """Return a built-in flat marker-program set as dictionaries.

    For tumor-specific sets such as ``tumor_breast`` and ``tumor_colon``,
    ``include_general=None`` prepends the shared ``tumor_general`` programs.
    For non-tumor auxiliary sets such as
    ``cell_cycle`` or ``stress``, ``include_general=None`` does not add tumor
    programs. Explicitly pass ``include_general=True`` to prepend
    ``tumor_general`` to any set.
    """
    canonical = _canonical_set_name(name)
    if canonical == "tumor_general":
        return [_normalize_program_entry(x) for x in _TUMOR_GENERAL]

    if include_general is None:
        include_general = canonical.startswith("tumor_")

    specific = [_normalize_program_entry(x) for x in _SETS[canonical]]
    if not include_general:
        return specific

    out = [_normalize_program_entry(x) for x in _TUMOR_GENERAL]
    seen = {x["name"] for x in out}
    for item in specific:
        if item["name"] not in seen:
            out.append(item)
            seen.add(item["name"])
    return out

