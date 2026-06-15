ANNOT_PATH_SEPARATOR = " > "

# Program-report blocking presets for the flat tumor/auxiliary program track.
PROGRAM_REPORT_BLOCK_OFF_TOKENS = {"none", "off", "false", "no_block", "no_blocklist"}
PROGRAM_REPORT_BLOCK_PRESETS = {"immune_like", "tumor_reportable", "lineage_aware"}

# Curated built-in hierarchy lineages where tumor-like reporting is allowed by
# program_report_block_preset="tumor_reportable". This is intentionally a
# package-level reporting allowlist for built-in hierarchies, not a general
# ontology resolver. Users with custom hierarchies should pass an explicit node
# block list or use lineage_aware with program-level report_on_lineages metadata.
TUMOR_REPORTABLE_BUILTIN_LINEAGES = {
    "epithelial",
    "epithelium",
    "epithelial cell",
    "epithelial cells",
    "unresolved",
    "parenchymal",
    "parenchyma",
    "parenchymal cell",
    "parenchymal cells",
    "mammary epithelial",
    "colon epithelial",
    "intestinal epithelial",
    "renal parenchymal",
    "kidney parenchymal",
    "hepatic parenchymal",
    "liver parenchymal",
    "pancreatic parenchymal",
    "pulmonary parenchymal",
    "cutaneous epithelial",
    "basal keratinocyte",
    "differentiated keratinocyte/squamous epithelial",
    "adnexal epithelial",
    "squamous epithelial",
    "melanocyte",
    "melanocytic",
}

# Short non-immune guard list used by the tumor_reportable preset. Immune labels
# are handled by is_immune_like_lineage().
TUMOR_REPORTABLE_GUARD_LINEAGES = {
    "fibroblast",
    "myofibroblast",
    "stromal",
    "stroma",
    "endothelial",
    "endothelium",
    "vascular",
    "mural",
    "pericyte",
    "smooth muscle",
}

# Reporting roles used by the tiered tumor/auxiliary program reporting layer.
REPORTING_ROLE_STATUS = "status"
REPORTING_ROLE_STATE = "state"
REPORTING_ROLE_MODIFIER = "modifier"
REPORTING_ROLE_UNKNOWN = "unknown"
REPORTING_ROLES = {REPORTING_ROLE_STATUS, REPORTING_ROLE_STATE, REPORTING_ROLE_MODIFIER}

# Aggregate tumor-status reporting defaults. Status-role programs establish
# tumor-like identity. State-role programs can provide a small, gated support
# bonus only when core status evidence is already borderline. Modifiers remain
# diagnostic and do not contribute to tumor-status decisions by default.
DEFAULT_STATUS_PROGRAM_BONUS_PER_EXTRA = 0.05
DEFAULT_STATUS_PROGRAM_MAX_BONUS = 0.15
DEFAULT_STATE_SUPPORT_BONUS_PER_GROUP = 0.05
DEFAULT_STATE_SUPPORT_MAX_BONUS = 0.10
DEFAULT_WEAK_STATUS_FRACTION_FOR_STATE_SUPPORT = 0.50

PROGRAM_LABEL_SEPARATOR = "_"
PROGRAM_LIST_SEPARATOR = ";"

