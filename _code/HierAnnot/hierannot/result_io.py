from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from .datamodels import HierAnnotResult
from .io import hierarchy_from_dict, hierarchy_to_dict


def _ensure_path(path: str | Path) -> Path:
    return path if isinstance(path, Path) else Path(path)


def _table_filename(name: str, table_format: str) -> str:
    if table_format == "csv":
        return f"{name}.csv"
    if table_format == "parquet":
        return f"{name}.parquet"
    raise ValueError("table_format must be 'csv' or 'parquet'")


def _write_table(df: pd.DataFrame | None, path: Path, table_format: str) -> None:
    if df is None:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    if table_format == "csv":
        df.to_csv(path, index=True)
    elif table_format == "parquet":
        df.to_parquet(path, index=True)
    else:
        raise ValueError("table_format must be 'csv' or 'parquet'")


def _read_table(path: Path) -> pd.DataFrame | None:
    if not path.exists():
        return None
    if path.suffix == ".csv":
        return pd.read_csv(path, index_col=0)
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    raise ValueError(f"Unsupported table file format for {path}")


def save_result_bundle(
    result: HierAnnotResult,
    path: str | Path,
    *,
    hierarchy=None,
    resolved_config: dict[str, Any] | None = None,
    metadata: dict[str, Any] | None = None,
    table_format: str = "csv",
) -> Path:
    outdir = _ensure_path(path)
    outdir.mkdir(parents=True, exist_ok=True)

    _write_table(getattr(result, "cluster_annotations", None), outdir / _table_filename("cluster_annotations", table_format), table_format)
    _write_table(getattr(result, "level_scores", None), outdir / _table_filename("level_scores", table_format), table_format)
    _write_table(getattr(result, "all_scores", None), outdir / _table_filename("all_scores", table_format), table_format)
    _write_table(getattr(result, "diagnostics_summary", None), outdir / _table_filename("diagnostics_summary", table_format), table_format)

    if hierarchy is None:
        hierarchy = getattr(result, "hierarchy", None)
    if resolved_config is None:
        resolved_config = getattr(result, "resolved_config", None)
    if metadata is None:
        metadata = getattr(result, "metadata", None)

    manifest = {
        "bundle_version": 1,
        "table_format": table_format,
        "has_hierarchy": hierarchy is not None,
        "has_resolved_config": resolved_config is not None,
        "metadata": metadata or {},
    }
    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2))

    if hierarchy is not None:
        (outdir / "hierarchy.json").write_text(json.dumps(hierarchy_to_dict(hierarchy), indent=2))

    if resolved_config is not None:
        (outdir / "resolved_config.json").write_text(json.dumps(resolved_config, indent=2, default=str))

    if metadata is not None:
        (outdir / "metadata.json").write_text(json.dumps(metadata, indent=2, default=str))

    return outdir



def load_result_bundle(path: str | Path, as_dict: bool = False):
    """
    Load a structured HierAnnot result bundle from a directory.

    By default this returns a `HierAnnotResult`-compatible object with attributes such as
    `.cluster_annotations`, `.all_scores`, `.diagnostics_summary`, optional `.metadata`,
    and optional `.hierarchy`, so it works naturally with plotting helpers.

    Set `as_dict=True` to get the raw dictionary representation instead.
    """
    indir = _ensure_path(path)
    if not indir.exists():
        raise FileNotFoundError(f"Result bundle path does not exist: {indir}")

    manifest = {}
    manifest_path = indir / "manifest.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())

    bundle: dict[str, Any] = {"manifest": manifest} if manifest else {}

    table_names = [
        "cluster_annotations",
        "level_scores",
        "all_scores",
        "diagnostics_summary",
    ]
    for name in table_names:
        csv_path = indir / f"{name}.csv"
        parquet_path = indir / f"{name}.parquet"
        if csv_path.exists():
            bundle[name] = _read_table(csv_path)
        elif parquet_path.exists():
            bundle[name] = _read_table(parquet_path)
        else:
            bundle[name] = None

    # Normalize cluster identifier columns to string so loaded bundles behave the
    # same way as in-memory results during downstream merges.
    if bundle.get("cluster_annotations") is not None and "cluster_id" in bundle["cluster_annotations"].columns:
        bundle["cluster_annotations"]["cluster_id"] = bundle["cluster_annotations"]["cluster_id"].astype(str)
    if bundle.get("all_scores") is not None and "cluster" in bundle["all_scores"].columns:
        bundle["all_scores"]["cluster"] = bundle["all_scores"]["cluster"].astype(str)
    if bundle.get("level_scores") is not None and "cluster" in bundle["level_scores"].columns:
        bundle["level_scores"]["cluster"] = bundle["level_scores"]["cluster"].astype(str)

    hierarchy_path = indir / "hierarchy.json"
    if hierarchy_path.exists():
        bundle["hierarchy"] = hierarchy_from_dict(json.loads(hierarchy_path.read_text()))

    resolved_config_path = indir / "resolved_config.json"
    if resolved_config_path.exists():
        bundle["resolved_config"] = json.loads(resolved_config_path.read_text())

    metadata_path = indir / "metadata.json"
    if metadata_path.exists():
        bundle["metadata"] = json.loads(metadata_path.read_text())

    if as_dict:
        return bundle

    return HierAnnotResult(
        cluster_annotations=bundle.get("cluster_annotations"),
        level_scores=bundle.get("level_scores"),
        all_scores=bundle.get("all_scores"),
        normalized_matrix=None,
        control_gene_map={},
        compiled_programs=[],
        compilation_report=None,
        diagnostics_summary=bundle.get("diagnostics_summary"),
        resolved_config=bundle.get("resolved_config", {}) or {},
        metadata=bundle.get("metadata"),
        hierarchy=bundle.get("hierarchy"),
    )
