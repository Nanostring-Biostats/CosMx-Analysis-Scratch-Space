from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Union

from .datamodels import MarkerProgram

PathLike = Union[str, Path]


def _program_to_dict(program: MarkerProgram) -> Dict[str, Any]:
    return {
        "name": program.name,
        "positive_markers": list(program.positive_markers),
        "negative_markers": list(program.negative_markers),
        "children": [_program_to_dict(child) for child in program.children],
        "description": program.description,
        "aliases": list(program.aliases),
        "min_markers_present": program.min_markers_present,
        "metadata": dict(program.metadata),
    }



def _program_from_dict(data: Dict[str, Any]) -> MarkerProgram:
    return MarkerProgram(
        name=data["name"],
        positive_markers=list(data.get("positive_markers", [])),
        negative_markers=list(data.get("negative_markers", [])),
        children=[_program_from_dict(child) for child in data.get("children", [])],
        description=data.get("description"),
        aliases=list(data.get("aliases", [])),
        min_markers_present=data.get("min_markers_present"),
        metadata=dict(data.get("metadata", {})),
    )



def hierarchy_to_dict(programs: List[MarkerProgram]) -> List[Dict[str, Any]]:
    """Serialize a hierarchy into plain Python dictionaries."""
    return [_program_to_dict(program) for program in programs]



def hierarchy_from_dict(data: List[Dict[str, Any]]) -> List[MarkerProgram]:
    """Reconstruct a hierarchy from dictionary records."""
    return [_program_from_dict(item) for item in data]



def save_hierarchy(programs: List[MarkerProgram], path: PathLike, indent: int = 2) -> Path:
    """Save a hierarchy to a JSON file for reuse.

    Parameters
    ----------
    programs:
        Root marker programs to serialize.
    path:
        Output JSON path.
    indent:
        JSON indentation level.
    """
    out_path = Path(path)
    out_path.write_text(json.dumps(hierarchy_to_dict(programs), indent=indent), encoding="utf-8")
    return out_path



def load_hierarchy(path: PathLike) -> List[MarkerProgram]:
    """Load a hierarchy previously saved with :func:`save_hierarchy`."""
    in_path = Path(path)
    data = json.loads(in_path.read_text(encoding="utf-8"))
    if not isinstance(data, list):
        raise ValueError("Serialized hierarchy JSON must contain a top-level list of root programs.")
    return hierarchy_from_dict(data)
