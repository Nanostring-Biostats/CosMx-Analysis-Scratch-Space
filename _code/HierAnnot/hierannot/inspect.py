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
