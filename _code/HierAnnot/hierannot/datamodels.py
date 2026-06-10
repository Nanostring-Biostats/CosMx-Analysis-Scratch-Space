from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
from typing import Any, Dict, List, Literal, Optional


MarkerCleanupStrategy = Literal[
    "none",
    "subtract_direct_parent",
    "subtract_ancestors",
    "sibling_unique_only",
]

ControlGeneExclusionPolicy = Literal[
    "all_pipeline_markers",
    "current_program_only",
    "current_branch",
]

PresetName = Literal["auto", "small_panel", "medium_panel", "large_panel"]


@dataclass
class MarkerProgram:
    """Hierarchical marker program for cell type annotation.

    A practical pattern is to place broad identity markers at a parent node and
    optional state-like refinements (for example activated, cycling, interferon-high)
    as children beneath the identity branch. This keeps the main annotation path
    biologically interpretable while still allowing progressively finer labels.
    """

    name: str
    positive_markers: List[str]
    negative_markers: List[str] = field(default_factory=list)
    children: List["MarkerProgram"] = field(default_factory=list)
    description: Optional[str] = None
    aliases: List[str] = field(default_factory=list)
    min_markers_present: Optional[int] = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    def add_child(self, child: "MarkerProgram") -> None:
        self.children.append(child)

    def clone(self) -> "MarkerProgram":
        return deepcopy(self)


    def to_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "positive_markers": list(self.positive_markers),
            "negative_markers": list(self.negative_markers),
            "children": [child.to_dict() for child in self.children],
            "description": self.description,
            "aliases": list(self.aliases),
            "min_markers_present": self.min_markers_present,
            "metadata": dict(self.metadata),
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "MarkerProgram":
        return cls(
            name=data["name"],
            positive_markers=list(data.get("positive_markers", [])),
            negative_markers=list(data.get("negative_markers", [])),
            children=[cls.from_dict(child) for child in data.get("children", [])],
            description=data.get("description"),
            aliases=list(data.get("aliases", [])),
            min_markers_present=data.get("min_markers_present"),
            metadata=dict(data.get("metadata", {})),
        )

    def iter_nodes(self, level: int = 1, parent: Optional[str] = None):
        yield {
            "name": self.name,
            "level": level,
            "parent": parent,
            "n_positive_markers": len(self.positive_markers),
            "n_negative_markers": len(self.negative_markers),
            "n_children": len(self.children),
            "metadata": dict(self.metadata),
        }
        for child in self.children:
            yield from child.iter_nodes(level=level + 1, parent=self.name)


@dataclass
class MalignantProgram:
    """Flat marker program scored by the tumor/auxiliary program track.

    Unlike :class:`MarkerProgram`, these programs are not routed as a hierarchy.
    They are scored in parallel to capture tumor/transformation signals by
    default, and can also be used in ``malignant_integration_mode="flag_only"``
    to score auxiliary disease or state programs such as proliferation, hypoxia,
    fibrosis, or inflammatory activation without changing normal hierarchy
    export labels.

    Required fields are ``name`` and ``positive_markers``. ``negative_markers``
    are optional. ``metadata`` may include ``competition_group`` to make related
    programs compete for state specificity, ``reporting_label`` to control
    concise program labels, and ``report_on_lineages`` to restrict reportability
    under ``program_report_block_preset="lineage_aware"``.
    ``reporting_role`` is active reporting metadata with three recommended
    values: ``"status"`` programs establish tumor-like identity for integration,
    ``"state"`` programs decorate a tumor-like label when specific, and
    ``"modifier"`` programs are reported as auxiliary flags.
    """

    name: str
    positive_markers: List[str]
    negative_markers: List[str] = field(default_factory=list)
    description: Optional[str] = None
    tags: List[str] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "positive_markers": list(self.positive_markers),
            "negative_markers": list(self.negative_markers),
            "description": self.description,
            "tags": list(self.tags),
            "metadata": dict(self.metadata),
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "MalignantProgram":
        return cls(
            name=data.get("name", data.get("label", "transformed")),
            positive_markers=list(data.get("positive_markers", [])),
            negative_markers=list(data.get("negative_markers", [])),
            description=data.get("description"),
            tags=list(data.get("tags", [])),
            metadata=dict(data.get("metadata", {})),
        )


@dataclass
class CompiledMarkerProgram:
    """Internal representation with canonical and effective scoring markers."""

    name: str
    canonical_markers: List[str]
    effective_markers: List[str]
    negative_markers: List[str] = field(default_factory=list)
    children: List["CompiledMarkerProgram"] = field(default_factory=list)
    description: Optional[str] = None
    aliases: List[str] = field(default_factory=list)
    min_markers_present: Optional[int] = None
    level: int = 1
    parent_name: Optional[str] = None
    cleanup_strategy: MarkerCleanupStrategy = "none"
    metadata: Dict[str, Any] = field(default_factory=dict)

    def iter_nodes(self):
        yield self
        for child in self.children:
            yield from child.iter_nodes()


@dataclass
class ScoreRecord:
    cluster: str
    label: str
    level: int
    parent_label: Optional[str]
    score: float
    positive_score: float
    negative_score: float
    positive_marker_mean: float
    positive_control_mean: float
    negative_marker_mean: float
    negative_control_mean: float
    negative_weight: float
    markers_present: int
    markers_total: int
    markers_present_fraction: float
    markers_above_detection_floor: int
    markers_detection_support_fraction: float
    negative_markers_present: int
    negative_markers_total: int
    negative_markers_present_fraction: float
    negative_markers_above_detection_floor: int
    negative_markers_detection_support_fraction: float
    missing_markers: List[str]
    missing_negative_markers: List[str]
    marker_source: str
    score_status: str
    effective_markers_total: int
    effective_markers_present: int
    canonical_markers_total: int
    used_fallback_marker_set: bool
    positive_controls_total: int
    negative_controls_total: int
    control_fallback_used: bool
    median_controls_per_positive_marker: float
    median_controls_per_negative_marker: float


@dataclass
class ClusterAnnotation:
    cluster: str
    final_label: str
    final_path: str
    final_level: int
    status: str
    stop_reason: str
    ambiguity_flag: bool
    confidence: str
    best_score_any_level: float
    best_label_any_level: str
    final_score: float
    final_raw_score: float
    final_branch_supported_raw_score: float
    final_margin: float
    decision_label: str
    decision_path: str
    decision_level: int
    decision_score: float
    decision_raw_score: float
    decision_branch_supported_raw_score: float
    backoff_applied: bool
    backoff_steps: int
    backoff_reason: str
    level_labels: Dict[str, str]
    level_scores: Dict[str, float]
    level_raw_scores: Dict[str, float]
    level_branch_supported_raw_scores: Dict[str, float]
    level_margins: Dict[str, float]
    level_confidence: Dict[str, str]


@dataclass
class HierAnnotResult:
    """Container returned by :meth:`HierAnnotPipeline.fit_score`.

    Key tables
    ----------
    cluster_annotations
        Normal hierarchy annotation table. This is the primary lineage/cell-type
        call per cluster.
    malignant_scores
        Per-cluster, per-program flat marker-program score table. Although the
        column names are tumor-focused, the same table is also used for
        auxiliary flag-only program sets such as cell cycle, stress,
        inflammation, or fibrosis.
    malignant_annotations
        One-row-per-cluster summary of the flat tumor/auxiliary program track.
        In ``malignant_integration_mode="flag_only"``, normal export labels
        are unchanged and this is the main table to inspect auxiliary status via
        columns such as ``annot_malignant_status``,
        ``annot_malignant_label_concise``, ``annot_malignant_status_label``,
        ``annot_malignant_state_label``, ``annot_malignant_modifier_labels``,
        ``annot_malignant_positive_programs``, ``annot_malignant_raw_score`` and
        ``annot_malignant_status_score``. For tumor integration, aggregated
        status-role evidence sets ``annot_malignant_tumor_status_pass=True``;
        state-role programs can provide a small gated support bonus only when
        core status evidence is already borderline, while modifier programs
        remain diagnostic by default.
    integrated_annotations
        Combined normal + tumor-track interpretation. In integrate mode this
        table records tumor-like labels and report-blocking diagnostics such as
        ``annot_malignant_reporting_blocked`` and
        ``annot_malignant_reporting_blocked_reason``.
    malignant_programs
        Flat tumor/auxiliary marker programs used to compute malignant_scores
        and malignant_annotations. These definitions are persisted in result
        bundles as ``malignant_programs.json`` when available so rerun/export
        workflows remain auditable and reproducible across package versions.
    """

    cluster_annotations: "object"
    level_scores: "object"
    all_scores: "object"
    malignant_scores: "object" = None
    malignant_annotations: "object" = None
    integrated_annotations: "object" = None
    normalized_matrix: "object" = None
    control_gene_map: Dict[str, Dict[str, Dict[str, List[str]]]] = field(default_factory=dict)
    compiled_programs: List[CompiledMarkerProgram] = field(default_factory=list)
    compilation_report: "object" = None
    diagnostics_summary: "object" = None
    resolved_config: Dict[str, object] = field(default_factory=dict)
    metadata: Optional[Dict[str, Any]] = None
    hierarchy: Optional[List[MarkerProgram]] = None
    malignant_programs: Optional[List[MalignantProgram]] = None
