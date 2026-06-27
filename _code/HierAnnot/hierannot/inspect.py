from __future__ import annotations

from typing import Iterable, List

import pandas as pd

from .datamodels import MarkerProgram


def _iter_all_nodes(programs: Iterable[MarkerProgram]):
    for root in programs:
        yield from root.iter_nodes()


def describe_hierarchy(programs: List[MarkerProgram]) -> pd.DataFrame:
    """Return a flat summary table for any hierarchy.

    Parameters
    ----------
    programs:
        Root marker programs defining a hierarchy.

    Returns
    -------
    pandas.DataFrame
        One row per node with level, parent, marker counts, and metadata.
    """
    rows = list(_iter_all_nodes(programs))
    if not rows:
        return pd.DataFrame(
            columns=[
                "name",
                "level",
                "parent",
                "n_positive_markers",
                "n_negative_markers",
                "n_children",
                "metadata",
            ]
        )
    return pd.DataFrame(rows)


def format_hierarchy_tree(programs: List[MarkerProgram], indent: str = "  ") -> str:
    """Format any hierarchy as a readable text tree.

    Each node is rendered with compact marker-count information using the form
    ``(+<n_positive>, -<n_negative>, children=<n_children>)``.
    """
    lines: List[str] = []

    def visit(node: MarkerProgram, depth: int = 0) -> None:
        lines.append(
            f"{indent * depth}- {node.name} (+{len(node.positive_markers)}, -{len(node.negative_markers)}, children={len(node.children)})"
        )
        for child in node.children:
            visit(child, depth + 1)

    for root in programs:
        visit(root)
    return "\n".join(lines)


def print_hierarchy_tree(programs: List[MarkerProgram], indent: str = "  ") -> None:
    """Print a readable text tree for any hierarchy."""
    print(format_hierarchy_tree(programs, indent=indent))

def _marker_program_to_dict(program):
    if isinstance(program, dict):
        return {
            "name": program.get("name", program.get("label", "<unnamed>")),
            "positive_markers": list(program.get("positive_markers", [])),
            "negative_markers": list(program.get("negative_markers", [])),
            "description": program.get("description"),
            "metadata": dict(program.get("metadata", {}) or {}),
        }
    return {
        "name": getattr(program, "name", "<unnamed>"),
        "positive_markers": list(getattr(program, "positive_markers", []) or []),
        "negative_markers": list(getattr(program, "negative_markers", []) or []),
        "description": getattr(program, "description", None),
        "metadata": dict(getattr(program, "metadata", {}) or {}),
    }


def summarize_marker_program_set(programs) -> pd.DataFrame:
    """Return a flat summary table for a marker-program set.

    Marker-program sets are flat rather than hierarchical. This summary groups
    each program by reporting role and competition group so users can quickly
    inspect how a tumor/auxiliary program set will be scored and reported.
    """
    rows = []
    for item in programs:
        rec = _marker_program_to_dict(item)
        meta = rec["metadata"]
        rows.append({
            "name": rec["name"],
            "reporting_label": meta.get("reporting_label", rec["name"]),
            "reporting_role": meta.get("reporting_role", "unknown"),
            "competition_group": meta.get("competition_group", rec["name"]),
            "program_family": meta.get("program_family"),
            "tissue_scope": meta.get("tissue_scope"),
            "n_positive_markers": len(rec["positive_markers"]),
            "n_negative_markers": len(rec["negative_markers"]),
            "positive_markers": list(rec["positive_markers"]),
            "negative_markers": list(rec["negative_markers"]),
            "report_on_lineages": meta.get("report_on_lineages"),
            "description": rec.get("description"),
            "metadata": meta,
        })
    if not rows:
        return pd.DataFrame(columns=[
            "name", "reporting_label", "reporting_role", "competition_group",
            "program_family", "tissue_scope", "n_positive_markers",
            "n_negative_markers", "positive_markers", "negative_markers",
            "report_on_lineages", "description", "metadata",
        ])
    role_order = {"status": 0, "state": 1, "modifier": 2, "unknown": 3}
    df = pd.DataFrame(rows)
    df["_role_order"] = df["reporting_role"].map(lambda x: role_order.get(str(x).lower(), 3))
    df = df.sort_values(["_role_order", "competition_group", "name"]).drop(columns=["_role_order"]).reset_index(drop=True)
    return df


def _format_marker_list(markers: List[str], max_markers: int) -> str:
    markers = list(markers)
    if max_markers is not None and len(markers) > max_markers:
        shown = markers[:max_markers]
        return ", ".join(shown) + f", ... (+{len(markers) - max_markers} more)"
    return ", ".join(markers)


def format_marker_program_set(
    programs,
    title: str | None = None,
    indent: str = "  ",
    show_markers: bool = True,
    max_markers: int = 8,
    show_negative_markers: bool = False,
    show_lineages: bool = True,
    show_descriptions: bool = False,
) -> str:
    """Format a flat marker-program set as grouped readable text.

    Programs are grouped by ``reporting_role`` and then ``competition_group``.
    This mirrors :func:`format_hierarchy_tree` for normal hierarchies, but uses
    role/group blocks instead of tree indentation because marker programs are
    scored in parallel rather than routed as parent/child nodes.
    """
    df = summarize_marker_program_set(programs)
    if df.empty:
        return "Marker program set: <empty>"

    lines: List[str] = []
    header = title or "Marker program set"
    lines.append(f"{header}: {len(df)} programs")

    for role, role_df in df.groupby("reporting_role", sort=False, dropna=False):
        role_label = str(role) if str(role) and str(role) != "nan" else "unknown"
        lines.append(f"\n[{role_label}]")
        for group, group_df in role_df.groupby("competition_group", sort=False, dropna=False):
            group_label = str(group) if str(group) and str(group) != "nan" else "independent"
            lines.append(f"{indent}competition_group: {group_label}")
            for _, row in group_df.iterrows():
                label = row.get("reporting_label") or row["name"]
                lines.append(
                    f"{indent * 2}- {row['name']} (label={label}; +{row['n_positive_markers']}, -{row['n_negative_markers']})"
                )
                if show_markers:
                    pos = _format_marker_list(row["positive_markers"], max_markers=max_markers)
                    if pos:
                        lines.append(f"{indent * 3}+ markers: {pos}")
                    if show_negative_markers:
                        neg = _format_marker_list(row["negative_markers"], max_markers=max_markers)
                        if neg:
                            lines.append(f"{indent * 3}- markers: {neg}")
                if show_lineages and row.get("report_on_lineages"):
                    lineages = row["report_on_lineages"]
                    if isinstance(lineages, (list, tuple, set)):
                        lineages = ", ".join(map(str, lineages))
                    lines.append(f"{indent * 3}report_on_lineages: {lineages}")
                if show_descriptions and row.get("description"):
                    lines.append(f"{indent * 3}description: {row['description']}")
    return "\n".join(lines)


def print_marker_program_set(programs, **kwargs) -> None:
    """Print :func:`format_marker_program_set` output."""
    print(format_marker_program_set(programs, **kwargs))

