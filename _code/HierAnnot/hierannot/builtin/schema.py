from __future__ import annotations

from typing import List, Optional

from ..datamodels import MarkerProgram
from .._version import __version__ as PACKAGE_VERSION


def _meta(
    category: str,
    tissue_scope: str,
    pack: str,
    version: str = PACKAGE_VERSION,
    role: str = "major",
    lineage_module: Optional[str] = None,
    requires_ancestor_support: bool = True,
    fallback_to_parent: bool = True,
) -> dict:
    return {
        "category": category,
        "tissue_scope": tissue_scope,
        "pack": pack,
        "version": version,
        "role": role,
        "lineage_module": lineage_module or category,
        "requires_ancestor_support": requires_ancestor_support,
        "fallback_to_parent": fallback_to_parent,
    }


def _mp(
    name: str,
    pos: List[str],
    neg: Optional[List[str]] = None,
    children: Optional[List[MarkerProgram]] = None,
    description: Optional[str] = None,
    aliases: Optional[List[str]] = None,
    metadata: Optional[dict] = None,
) -> MarkerProgram:
    return MarkerProgram(
        name=name,
        positive_markers=pos,
        negative_markers=neg or [],
        children=children or [],
        description=description,
        aliases=aliases or [],
        metadata=metadata or {},
    )
