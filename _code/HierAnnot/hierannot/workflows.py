from __future__ import annotations

from .hierarchy import decide_cluster_annotation_from_scores

from typing import Iterable, Mapping, Optional, Sequence, List, Dict, Set

import numpy as np
import re
import pandas as pd
import warnings

try:
    from scipy import sparse as sp
except Exception:  # pragma: no cover
    sp = None



def aggregate_expression_to_cluster_means(
    expression,
    cluster_labels: Sequence,
    feature_names: Optional[Sequence[str]] = None,
    cell_names: Optional[Sequence[str]] = None,
    method: str = "mean",
    trim_fraction: float = 0.1,
    include_mask: Optional[Sequence[bool] | pd.Series | np.ndarray] = None,
    min_cells_per_cluster: int = 1,
    centrality_matrix=None,
    central_fraction: float = 0.8,
    centrality_mode: str = "centroid",
    neighbor_connectivities=None,
    backend: str = "auto",
    chunk_size: int = 50000,
    warn_on_expensive_dense: bool = True,
    uppercase_features: bool = False,
) -> pd.DataFrame:
    """Aggregate an expression matrix into feature x cluster profiles.

    Orientation contract
    --------------------
    `expression` is interpreted as **cells x features** for all supported input
    types, including pandas DataFrame, NumPy array, and SciPy sparse matrix.

    Parameters
    ----------
    expression
        Cells x features expression matrix.
    cluster_labels
        Cluster label per cell (one label per row of `expression`).
    feature_names, cell_names
        Optional names used when `expression` is not a DataFrame.
    method
        One of: mean, median, trimmed_mean, central_cells_mean, central_cells_trimmed_mean
    centrality_mode
        For central-cells methods:
        - ``centroid``: distance to cluster centroid in `centrality_matrix`
        - ``graph_degree``: within-cluster weighted degree in `neighbor_connectivities`
    """
    method = str(method)
    centrality_mode = str(centrality_mode).lower()
    if backend not in {"auto", "sparse", "chunked", "dense"}:
        raise ValueError("backend must be one of: 'auto', 'sparse', 'chunked', 'dense'")
    if centrality_mode not in {"centroid", "graph_degree"}:
        raise ValueError("centrality_mode must be one of: 'centroid', 'graph_degree'")

    if isinstance(expression, pd.DataFrame):
        matrix = expression.to_numpy(dtype=float)
        if cell_names is None:
            cell_names = expression.index.astype(str)
        else:
            cell_names = pd.Index(pd.Index(cell_names).astype(str))
        if feature_names is None:
            feature_names = expression.columns.astype(str)
        else:
            feature_names = pd.Index(pd.Index(feature_names).astype(str))
    else:
        matrix = expression
        if getattr(matrix, "ndim", 2) != 2:
            raise ValueError("expression must be 2D")
        n_cells, n_features = matrix.shape
        if feature_names is None:
            feature_names = pd.Index([f"feature_{i}" for i in range(n_features)])
        else:
            feature_names = pd.Index(pd.Index(feature_names).astype(str))
        if cell_names is None:
            cell_names = pd.Index([str(i) for i in range(n_cells)])
        else:
            cell_names = pd.Index(pd.Index(cell_names).astype(str))

    feature_names = pd.Index(pd.Index(feature_names).astype(str))
    if uppercase_features:
        feature_names = feature_names.str.upper()
        if feature_names.has_duplicates:
            raise ValueError("uppercase_features=True produced duplicate feature names; merge duplicates before aggregation.")

    n_cells = matrix.shape[0]
    if len(cluster_labels) != n_cells:
        raise ValueError("Length of cluster_labels must match the number of cells/rows in expression.")

    cluster_series = pd.Series(list(cluster_labels), index=np.arange(n_cells), name="cluster").astype(str)

    if include_mask is None:
        include_mask_arr = np.ones(n_cells, dtype=bool)
    elif isinstance(include_mask, pd.Series):
        if len(include_mask) == n_cells:
            include_mask_arr = include_mask.astype(bool).to_numpy()
        else:
            include_mask_arr = include_mask.reindex(cell_names).astype(bool).to_numpy()
    else:
        include_mask_arr = np.asarray(include_mask, dtype=bool)
        if include_mask_arr.shape[0] != n_cells:
            raise ValueError("Length of include_mask must match the number of cells in expression.")

    positions = np.flatnonzero(include_mask_arr)
    if len(positions) == 0:
        return pd.DataFrame(index=feature_names)

    matrix = _subset_matrix_rows(matrix, positions)
    cluster_series = cluster_series.loc[positions]

    cluster_series, categories, _, keep_mask = _factorize_valid_clusters(cluster_series, min_cells_per_cluster)
    if len(categories) == 0:
        return pd.DataFrame(index=feature_names)

    keep_positions = np.flatnonzero(keep_mask.to_numpy())
    matrix = _subset_matrix_rows(matrix, keep_positions)
    cluster_series = cluster_series.reset_index(drop=True)
    cluster_codes = pd.Categorical(cluster_series, categories=categories, ordered=True).codes

    if method in {"central_cells_mean", "central_cells_trimmed_mean"}:
        if centrality_mode == "centroid":
            if centrality_matrix is None:
                raise ValueError("centrality_matrix is required for central_cells_* when centrality_mode='centroid'.")
            cent = _subset_matrix_rows(centrality_matrix, positions)
            cent = _subset_matrix_rows(cent, keep_positions)
            selected_pos = _select_central_cell_positions(cent, cluster_series.reset_index(drop=True), central_fraction)
        else:
            if neighbor_connectivities is None:
                raise ValueError("neighbor_connectivities is required for central_cells_* when centrality_mode='graph_degree'.")
            conn = neighbor_connectivities
            if sp is not None and sp.issparse(conn):
                conn = conn.tocsr()[positions][:, positions]
                conn = conn[keep_positions][:, keep_positions]
            else:
                conn = np.asarray(conn)[np.ix_(positions, positions)]
                conn = conn[np.ix_(keep_positions, keep_positions)]
            selected_pos = _select_central_cell_positions_graph(conn, cluster_series.reset_index(drop=True), central_fraction)

        matrix = _subset_matrix_rows(matrix, selected_pos)
        cluster_series = cluster_series.reset_index(drop=True).iloc[selected_pos]
        cluster_codes = pd.Categorical(cluster_series, categories=categories, ordered=True).codes

    n_cells2, n_features2 = matrix.shape
    if warn_on_expensive_dense:
        _warn_if_expensive_dense_method(method, n_cells=n_cells2, n_features=n_features2, source="expression")

    if method in {"mean", "central_cells_mean"}:
        out, _ = _aggregate_mean_matrix(matrix, cluster_codes=cluster_codes, n_clusters=len(categories), chunk_size=chunk_size)
        return pd.DataFrame(out.T.astype(float), index=feature_names, columns=categories.astype(str))

    values = _to_dense_rows(matrix)
    cluster_series = cluster_series.reset_index(drop=True)

    def _trimmed_mean_nd(arr: np.ndarray, frac: float) -> np.ndarray:
        if arr.shape[0] <= 2 or frac <= 0:
            return arr.mean(axis=0)
        k = int(np.floor(arr.shape[0] * frac))
        if k <= 0 or (2 * k) >= arr.shape[0]:
            return arr.mean(axis=0)
        arr = np.sort(arr, axis=0)[k:arr.shape[0] - k, :]
        return arr.mean(axis=0)

    out_cols = []
    out_data = []
    for cluster_id, idx in cluster_series.groupby(cluster_series, sort=False).groups.items():
        pos = np.asarray(idx, dtype=int)
        arr = np.asarray(values[pos], dtype=float)
        if method == "median":
            agg = np.median(arr, axis=0)
        elif method in {"trimmed_mean", "central_cells_trimmed_mean"}:
            agg = _trimmed_mean_nd(arr, trim_fraction)
        else:
            raise ValueError("method must be one of: mean, median, trimmed_mean, central_cells_mean, central_cells_trimmed_mean")
        out_cols.append(str(cluster_id))
        out_data.append(agg)
    out = np.vstack(out_data).T if out_data else np.zeros((len(feature_names), 0))
    return pd.DataFrame(out, index=feature_names, columns=out_cols)
def _to_dense_2d(x):
    if hasattr(x, "toarray"):
        return x.toarray()
    return np.asarray(x)





def _resolve_adata_matrix_source(
    adata,
    source: str,
    key: Optional[str] = None,
    uppercase_features: bool = False,
):
    """Return (cells x features matrix, feature names, cell names) without densifying by default."""
    source_norm = str(source).lower()
    cell_names = pd.Index(adata.obs_names.astype(str))

    if source_norm == "raw":
        if getattr(adata, "raw", None) is None:
            raise ValueError("source='raw' but adata.raw is None.")
        matrix = adata.raw.X
        features = pd.Index(adata.raw.var_names.astype(str))
    elif source_norm == "x":
        matrix = adata.X
        features = pd.Index(adata.var_names.astype(str))
    elif source_norm == "layer":
        if key is None:
            raise ValueError("source='layer' requires key to be set to a layer name.")
        if key not in adata.layers:
            raise KeyError(f"layer '{key}' not found in adata.layers")
        matrix = adata.layers[key]
        features = pd.Index(adata.var_names.astype(str))
    elif source_norm == "obsm":
        if key is None:
            raise ValueError("source='obsm' requires key to be set to an obsm key.")
        if key not in adata.obsm:
            raise KeyError(f"obsm key '{key}' not found in adata.obsm")
        matrix = adata.obsm[key]
        if hasattr(matrix, "shape") and len(matrix.shape) != 2:
            raise ValueError("adata.obsm[key] must be a 2D matrix.")
        if isinstance(matrix, pd.DataFrame):
            features = pd.Index(matrix.columns.astype(str))
            matrix = matrix.to_numpy()
        elif f"{key}_feature_names" in adata.uns:
            features = pd.Index(pd.Index(adata.uns[f"{key}_feature_names"]).astype(str))
        else:
            n_features = int(matrix.shape[1])
            features = pd.Index([f"{_sanitize_label_text(key) or 'feature'}_{i}" for i in range(n_features)])
    else:
        raise ValueError("source must be one of: 'X', 'raw', 'layer', 'obsm'")

    features = features.astype(str)
    if uppercase_features:
        features = features.str.upper()
        if features.has_duplicates:
            raise ValueError("uppercase_features=True produced duplicate feature names for AnnData source; use a preprocessing step to merge duplicates first.")
    return matrix, features, cell_names


def _subset_matrix_rows(matrix, row_idx):
    if sp is not None and sp.issparse(matrix):
        return matrix[row_idx]
    return np.asarray(matrix)[row_idx]


def _to_dense_rows(matrix):
    if sp is not None and sp.issparse(matrix):
        return matrix.toarray()
    return np.asarray(matrix)


def _factorize_valid_clusters(cluster_series: pd.Series, min_cells_per_cluster: int):
    counts = cluster_series.astype(str).value_counts(sort=False)
    valid = counts[counts >= int(min_cells_per_cluster)].index.astype(str)
    keep = cluster_series.astype(str).isin(valid)
    cluster_series = cluster_series.astype(str).loc[keep]
    categories = pd.Index(cluster_series.unique().astype(str))
    codes = pd.Categorical(cluster_series, categories=categories, ordered=True).codes
    return cluster_series, categories, codes, keep


def _select_central_cell_positions(centrality_matrix, cluster_series: pd.Series, central_fraction: float):
    if not (0 < float(central_fraction) <= 1.0):
        raise ValueError("central_fraction must be in the interval (0, 1].")
    selected_positions = []
    for cluster_id, idx in cluster_series.groupby(cluster_series, sort=False).groups.items():
        pos = np.asarray(idx, dtype=int)
        rep = _subset_matrix_rows(centrality_matrix, pos)
        rep = _to_dense_rows(rep)
        if rep.shape[0] == 0:
            continue
        if rep.shape[0] == 1:
            keep_pos = pos
        else:
            centroid = rep.mean(axis=0)
            distances = np.linalg.norm(rep - centroid[None, :], axis=1)
            n_keep = max(1, int(np.ceil(len(pos) * float(central_fraction))))
            order = np.argsort(distances, kind="mergesort")[:n_keep]
            keep_pos = pos[order]
        selected_positions.extend(keep_pos.tolist())
    return np.asarray(selected_positions, dtype=int)




def _resolve_neighbors_connectivities(adata, neighbors_key: Optional[str] = None):
    """Return a sparse connectivity matrix aligned to adata.obs_names."""
    if sp is None:
        raise ImportError("Graph-based central cell selection requires scipy.sparse.")
    candidates = []
    if neighbors_key:
        candidates.extend([
            f"{neighbors_key}_connectivities",
            neighbors_key,
        ])
    candidates.append("connectivities")
    for key in candidates:
        if key in getattr(adata, "obsp", {}):
            conn = adata.obsp[key]
            if conn.shape[0] != adata.n_obs or conn.shape[1] != adata.n_obs:
                raise ValueError(f"Neighbor graph '{key}' shape does not match adata.n_obs.")
            return conn
    raise KeyError(
        "Could not find a neighbor connectivity matrix in adata.obsp. "
        "Provide neighbors_key matching a stored connectivities key."
    )


def _select_central_cell_positions_graph(connectivities, cluster_series: pd.Series, central_fraction: float):
    """Select central cells using within-cluster weighted degree on a precomputed neighbor graph."""
    if sp is None or not sp.issparse(connectivities):
        if sp is None:
            raise ImportError("Graph-based central cell selection requires scipy.sparse.")
        connectivities = sp.csr_matrix(connectivities)
    else:
        connectivities = connectivities.tocsr()

    if not (0 < float(central_fraction) <= 1.0):
        raise ValueError("central_fraction must be in the interval (0, 1].")

    selected_positions = []
    for _, idx in cluster_series.groupby(cluster_series, sort=False).groups.items():
        pos = np.asarray(idx, dtype=int)
        if pos.size == 0:
            continue
        if pos.size == 1:
            selected_positions.extend(pos.tolist())
            continue
        sub = connectivities[pos][:, pos]
        degree = np.asarray(sub.sum(axis=1)).ravel()
        n_keep = max(1, int(np.ceil(len(pos) * float(central_fraction))))
        order = np.argsort(-degree, kind="mergesort")[:n_keep]
        keep_pos = pos[order]
        selected_positions.extend(keep_pos.tolist())
    return np.asarray(selected_positions, dtype=int)

def _aggregate_mean_matrix(matrix, cluster_codes: np.ndarray, n_clusters: int, chunk_size: int = 50000):
    n_cells, n_features = matrix.shape
    counts = np.bincount(cluster_codes, minlength=n_clusters).astype(float)

    if sp is not None and sp.issparse(matrix):
        matrix = matrix.tocsr()
        membership = sp.csr_matrix(
            (np.ones(n_cells, dtype=np.float32), (cluster_codes, np.arange(n_cells))),
            shape=(n_clusters, n_cells),
        )
        sums = membership @ matrix
        out = sums.toarray().astype(float, copy=False)
    else:
        out = np.zeros((n_clusters, n_features), dtype=float)
        for start in range(0, n_cells, int(chunk_size)):
            stop = min(n_cells, start + int(chunk_size))
            block = np.asarray(matrix[start:stop], dtype=float)
            block_codes = cluster_codes[start:stop]
            np.add.at(out, block_codes, block)
    nonzero = counts > 0
    out[nonzero] = out[nonzero] / counts[nonzero, None]
    return out, counts


def _warn_if_expensive_dense_method(method: str, n_cells: int, n_features: int, source: str):
    method = str(method)
    if method not in {"median", "trimmed_mean", "central_cells_trimmed_mean"}:
        return
    size = int(n_cells) * int(n_features)
    if size >= 100_000_000:
        warnings.warn(
            f"{method} aggregation on a large matrix ({n_cells} cells x {n_features} features from source={source!r}) "
            "may be slow and memory-intensive because it cannot use the sparse mean fast path. "
            "Consider method='mean' or 'central_cells_mean' for very large sparse datasets.",
            RuntimeWarning,
            stacklevel=2,
        )

def _sanitize_label_text(value: object) -> str:
    text = str(value).strip()
    for old, new in [(" ", "_"), ("/", "_"), ("-", "_")]:
        text = text.replace(old, new)
    while "__" in text:
        text = text.replace("__", "_")
    return text.strip("_")




def aggregate_anndata_to_cluster_means(
    adata,
    cluster_key: str,
    source: str = "X",
    source_key: Optional[str] = None,
    uppercase_genes: bool = True,
    method: str = "mean",
    trim_fraction: float = 0.1,
    include_obs_mask: Optional[str | pd.Series | Sequence[bool] | np.ndarray] = None,
    min_cells_per_cluster: int = 1,
    centrality_source: Optional[str] = None,
    centrality_source_key: Optional[str] = None,
    central_fraction: float = 0.8,
    centrality_mode: str = "centroid",
    neighbors_key: Optional[str] = None,
    backend: str = "auto",
    chunk_size: int = 50000,
    warn_on_expensive_dense: bool = True,
) -> pd.DataFrame:
    """Aggregate an AnnData object to feature x cluster profiles.

    This is a convenience wrapper around `aggregate_expression_to_cluster_means()`
    that resolves matrices, masks, and optional neighbor graphs from AnnData.
    """
    if cluster_key not in adata.obs:
        raise KeyError(f"cluster_key '{cluster_key}' not found in adata.obs")

    source_norm = str(source).lower()
    matrix, features, cell_names = _resolve_adata_matrix_source(
        adata,
        source=source,
        key=source_key,
        uppercase_features=uppercase_genes if source_norm != "obsm" else False,
    )

    if include_obs_mask is None:
        resolved_include_mask = None
    elif isinstance(include_obs_mask, str):
        if include_obs_mask not in adata.obs:
            raise KeyError(f"include_obs_mask column '{include_obs_mask}' not found in adata.obs")
        resolved_include_mask = adata.obs[include_obs_mask].astype(bool).to_numpy()
    else:
        resolved_include_mask = include_obs_mask

    resolved_centrality = None
    resolved_graph = None
    if method in {"central_cells_mean", "central_cells_trimmed_mean"}:
        if centrality_mode == "centroid":
            if centrality_source is None:
                raise ValueError("centrality_source is required for central_cells_* when centrality_mode='centroid'.")
            resolved_centrality, _, _ = _resolve_adata_matrix_source(
                adata,
                source=centrality_source,
                key=centrality_source_key,
                uppercase_features=False,
            )
        else:
            resolved_graph = _resolve_neighbors_connectivities(adata, neighbors_key=neighbors_key)

    return aggregate_expression_to_cluster_means(
        matrix,
        adata.obs[cluster_key].astype(str).tolist(),
        feature_names=features,
        cell_names=cell_names,
        method=method,
        trim_fraction=trim_fraction,
        include_mask=resolved_include_mask,
        min_cells_per_cluster=min_cells_per_cluster,
        centrality_matrix=resolved_centrality,
        central_fraction=central_fraction,
        centrality_mode=centrality_mode,
        neighbor_connectivities=resolved_graph,
        backend=backend,
        chunk_size=chunk_size,
        warn_on_expensive_dense=warn_on_expensive_dense,
        uppercase_features=False,
    )

def _prepare_clustering_representation(
    adata,
    source: str = "obsm",
    source_key: Optional[str] = None,
    n_components: int = 30,
    random_state: int = 0,
    pca_zero_center: Optional[bool] = None,
    pca_scale: bool = False,
    do_pca: bool = True,
):
    """Prepare a clustering representation from AnnData.

    Returns a tuple of ``(representation, source_label, resolved_zero_center)`` where
    ``representation`` is a cell x feature numpy array suitable for ``scanpy.pp.neighbors``
    via ``use_rep``.

    PCA/SVD preprocessing is controlled by ``do_pca``:
    - ``do_pca=True``: attempt PCA/SVD, but skip it when ``n_components`` is not meaningful
    - ``do_pca=False``: use the selected representation directly

    Default centering behavior remains source-aware when PCA is used:
    - ``obsm`` inputs: no extra centering before dimensionality reduction
    - ``X``/``layer`` inputs: zero-center before PCA by default
    """
    source_norm = str(source).lower()
    rep_matrix, _, _ = _resolve_adata_matrix_source(adata, source=source_norm, key=source_key, uppercase_features=False)
    embedding = _to_dense_rows(rep_matrix).astype(float)
    if embedding.ndim != 2:
        raise ValueError("Selected clustering source must resolve to a 2D cell x feature matrix.")

    resolved_zero_center = pca_zero_center
    if resolved_zero_center is None:
        resolved_zero_center = source_norm != "obsm"

    source_label = source_norm if source_norm == "x" else f"{source_norm}_{source_key}"

    if embedding.shape[1] == 0:
        raise ValueError("Selected clustering source has zero features.")

    requested_components = int(n_components)
    n_components = max(1, min(requested_components, embedding.shape[1], embedding.shape[0]))

    if pca_scale:
        try:
            from sklearn.preprocessing import StandardScaler
        except Exception as exc:  # pragma: no cover - optional dependency
            raise ImportError("PCA scaling requires scikit-learn. Install HierAnnot with preprocessing extras.") from exc
        scaler = StandardScaler(with_mean=bool(resolved_zero_center), with_std=True)
        embedding = scaler.fit_transform(embedding)

    if not bool(do_pca):
        return embedding, source_label, bool(resolved_zero_center)

    if requested_components >= embedding.shape[1]:
        # No real dimensionality reduction to perform; keep the representation as-is.
        return embedding, source_label, bool(resolved_zero_center)

    if resolved_zero_center:
        try:
            from sklearn.decomposition import PCA
        except Exception as exc:  # pragma: no cover - optional dependency
            raise ImportError("PCA-based clustering requires scikit-learn. Install HierAnnot with preprocessing extras.") from exc
        reducer = PCA(n_components=n_components, random_state=random_state)
        representation = reducer.fit_transform(embedding)
    else:
        try:
            from sklearn.decomposition import TruncatedSVD
        except Exception as exc:  # pragma: no cover - optional dependency
            raise ImportError("SVD-based clustering requires scikit-learn. Install HierAnnot with preprocessing extras.") from exc
        reducer = TruncatedSVD(n_components=n_components, random_state=random_state)
        representation = reducer.fit_transform(embedding)

    return representation, source_label, bool(resolved_zero_center)


def cluster_anndata_on_representation(
    adata,
    source: str = "obsm",
    source_key: Optional[str] = None,
    cluster_key: str = "hierannot_leiden",
    pca_key: Optional[str] = None,
    neighbors_key: Optional[str] = None,
    umap_key: Optional[str] = None,
    n_pcs: int = 30,
    n_neighbors: int = 15,
    leiden_resolution: float = 1.0,
    random_state: int = 0,
    compute_umap: bool = True,
    copy: bool = False,
    pca_zero_center: Optional[bool] = None,
    pca_scale: bool = False,
    do_pca: bool = True,
):
    """Cluster cells from a chosen AnnData representation.

    Supports ``source='obsm'``, ``'layer'``, and ``'X'``. The helper is written
    to work with Scanpy 1.10.3 and 1.11.x by relying on stable ``neighbors_key``
    / ``key_added`` behavior and by storing custom UMAP output under a separate
    ``obsm`` key without overwriting existing embeddings by default.

    PCA preprocessing defaults are source-aware:
    - ``obsm``: no extra zero-centering or scaling before SVD/PCA
    - ``X``/``layer``: zero-center before PCA, no scaling by default
    """
    try:
        import scanpy as sc
    except Exception as exc:  # pragma: no cover - optional dependency
        raise ImportError("cluster_anndata_on_representation requires scanpy. Install HierAnnot with preprocessing extras.") from exc

    target = adata.copy() if copy else adata
    representation, source_label, resolved_zero_center = _prepare_clustering_representation(
        target,
        source=source,
        source_key=source_key,
        n_components=n_pcs,
        random_state=random_state,
        pca_zero_center=pca_zero_center,
        pca_scale=pca_scale,
        do_pca=do_pca,
    )

    safe_stem = _sanitize_label_text(source_label) or "representation"
    rep_key = pca_key or f"X_rep_{safe_stem}"
    neighbors_key = neighbors_key or f"neighbors_{safe_stem}"
    umap_key = umap_key or f"X_umap_{safe_stem}"

    target.obsm[rep_key] = representation
    if not hasattr(target, "uns") or target.uns is None:
        target.uns = {}
    target.uns[f"{rep_key}_params"] = {
        "source": str(source).lower(),
        "source_key": source_key,
        "n_components": int(max(1, min(int(n_pcs), representation.shape[1], representation.shape[0]))),
        "pca_zero_center": bool(resolved_zero_center),
        "pca_scale": bool(pca_scale),
    }

    sc.pp.neighbors(
        target,
        use_rep=rep_key,
        n_neighbors=n_neighbors,
        random_state=random_state,
        key_added=neighbors_key,
    )
    sc.tl.leiden(
        target,
        key_added=cluster_key,
        resolution=leiden_resolution,
        random_state=random_state,
        neighbors_key=neighbors_key,
    )
    if compute_umap:
        old_x_umap = target.obsm.pop("X_umap", None)
        sc.tl.umap(target, neighbors_key=neighbors_key)
        if umap_key != "X_umap":
            target.obsm[umap_key] = target.obsm.pop("X_umap")
            if old_x_umap is not None:
                target.obsm["X_umap"] = old_x_umap
    return target


def cluster_and_aggregate_anndata(
    adata,
    cluster_key: str = "hierannot_leiden",
    clustering_source: str = "obsm",
    clustering_source_key: Optional[str] = None,
    aggregation_source: str = "X",
    aggregation_source_key: Optional[str] = None,
    pca_key: Optional[str] = None,
    neighbors_key: Optional[str] = None,
    umap_key: Optional[str] = None,
    n_pcs: int = 30,
    n_neighbors: int = 15,
    leiden_resolution: float = 1.0,
    random_state: int = 0,
    compute_umap: bool = True,
    copy: bool = False,
    aggregation_method: str = "mean",
    trim_fraction: float = 0.1,
    pca_zero_center: Optional[bool] = None,
    pca_scale: bool = False,
    do_pca: bool = True,
    include_obs_mask: Optional[str | pd.Series | Sequence[bool] | np.ndarray] = None,
    min_cells_per_cluster: int = 1,
    centrality_source: Optional[str] = None,
    centrality_source_key: Optional[str] = None,
    central_fraction: float = 0.8,
):
    """Cluster cells from one AnnData source, then aggregate another source by cluster.

    Parameters
    ----------
    clustering_source, clustering_source_key
        Define the matrix used for clustering. Supported values are ``"obsm"``,
        ``"layer"``, and ``"X"``. ``clustering_source_key`` is required for
        ``"obsm"`` and ``"layer"``.
    aggregation_source, aggregation_source_key
        Define the matrix used for cluster-level aggregation. Supported values are
        ``"X"``, ``"raw"``, ``"layer"``, and ``"obsm"``. ``aggregation_source_key``
        is required for ``"layer"`` and ``"obsm"``.
    centrality_source, centrality_source_key
        Optional representation used only for ``central_cells_*`` aggregation.
        If omitted for those methods, the clustering representation is used.
    """
    resolved_rep_key = pca_key
    if resolved_rep_key is None:
        source_label = str(clustering_source).lower() if str(clustering_source).lower() == "x" else f"{str(clustering_source).lower()}_{clustering_source_key}"
        safe_stem = _sanitize_label_text(source_label) or "representation"
        resolved_rep_key = f"X_rep_{safe_stem}"

    clustered = cluster_anndata_on_representation(
        adata=adata,
        cluster_key=cluster_key,
        pca_key=pca_key,
        neighbors_key=neighbors_key,
        umap_key=umap_key,
        n_pcs=n_pcs,
        n_neighbors=n_neighbors,
        leiden_resolution=leiden_resolution,
        random_state=random_state,
        compute_umap=compute_umap,
        copy=copy,
        source=clustering_source,
        source_key=clustering_source_key,
        pca_zero_center=pca_zero_center,
        pca_scale=pca_scale,
        do_pca=do_pca,
    )

    if centrality_source is None and aggregation_method in {"central_cells_mean", "central_cells_trimmed_mean"}:
        centrality_source = "obsm"
        centrality_source_key = resolved_rep_key

    cluster_means = aggregate_anndata_to_cluster_means(
        clustered,
        cluster_key=cluster_key,
        source=aggregation_source,
        source_key=aggregation_source_key,
        method=aggregation_method,
        trim_fraction=trim_fraction,
        include_obs_mask=include_obs_mask,
        min_cells_per_cluster=min_cells_per_cluster,
        centrality_source=centrality_source,
        centrality_source_key=centrality_source_key,
        central_fraction=central_fraction,
    )
    return clustered, cluster_means




def _build_per_level_summary(result) -> pd.DataFrame:
    """Build a wide per-cluster summary from result.level_scores."""
    level_scores = getattr(result, "level_scores", None)
    if level_scores is None or len(level_scores) == 0:
        return pd.DataFrame(columns=["cluster_id"])

    df = level_scores.copy()
    df["cluster"] = df["cluster"].astype(str)
    score_col = "decision_score" if "decision_score" in df.columns else "score"
    best_idx = df.groupby(["cluster", "level"], sort=False)[score_col].idxmax()
    best = df.loc[best_idx].copy()

    wide = pd.DataFrame({"cluster_id": sorted(df["cluster"].unique(), key=lambda x: x)})
    level_values = sorted(best["level"].dropna().unique())
    for level in level_values:
        sub = best[best["level"] == level].copy()
        rename_map = {
            "label": f"annot_l{int(level)}_label",
            score_col: f"annot_l{int(level)}_score",
            "parent_label": f"annot_l{int(level)}_parent_label",
            "positive_score": f"annot_l{int(level)}_positive_score",
            "negative_score": f"annot_l{int(level)}_negative_score",
            "marker_source": f"annot_l{int(level)}_marker_source",
            "score_status": f"annot_l{int(level)}_score_status",
            "branch_supported_raw_score": f"annot_l{int(level)}_branch_supported_raw_score",
            "markers_present": f"annot_l{int(level)}_markers_present",
            "markers_total": f"annot_l{int(level)}_markers_total",
            "markers_present_fraction": f"annot_l{int(level)}_markers_present_fraction",
            "markers_detection_support_fraction": f"annot_l{int(level)}_markers_detection_support_fraction",
            "negative_markers_present": f"annot_l{int(level)}_negative_markers_present",
            "negative_markers_total": f"annot_l{int(level)}_negative_markers_total",
            "negative_markers_present_fraction": f"annot_l{int(level)}_negative_markers_present_fraction",
            "used_fallback_marker_set": f"annot_l{int(level)}_used_fallback_marker_set",
            "control_fallback_used": f"annot_l{int(level)}_control_fallback_used",
            "median_controls_per_positive_marker": f"annot_l{int(level)}_median_controls_per_positive_marker",
            "median_controls_per_negative_marker": f"annot_l{int(level)}_median_controls_per_negative_marker",
            "sibling_specificity_score": f"annot_l{int(level)}_sibling_specificity_score",
            "program_robust_zscore": f"annot_l{int(level)}_program_robust_zscore",
        }
        if score_col != "score" and "score" in sub.columns:
            rename_map["score"] = f"annot_l{int(level)}_raw_score"
        keep_cols = ["cluster", *[c for c in rename_map if c in sub.columns]]
        sub = sub[keep_cols].rename(columns={k: v for k, v in rename_map.items() if k in keep_cols})
        sub = sub.rename(columns={"cluster": "cluster_id"})
        wide = wide.merge(sub, on="cluster_id", how="left")
    return wide





def _normalize_confidence_levels(levels: Optional[Iterable[str]]) -> set[str]:
    if levels is None:
        return {"low", "none"}
    return {str(x).strip().lower() for x in levels}




def _rerun_decision_summary_from_scores(
    result,
    score_threshold: float = 0.05,
    margin_threshold: float = 0.02,
    branch_support_parent_weight: float = 0.65,
    branch_support_descendant_weight: float = 0.35,
    level1_branch_support_parent_weight: float = 0.60,
    level1_branch_support_descendant_weight: float = 0.40,
    branch_routing_raw_weight: float = 0.60,
    weak_raw_score_threshold: float = 0.10,
    min_branch_supported_raw_score: float = 0.0,
    min_leaf_raw_score: float = -0.2,
):
    compiled = getattr(result, "compiled_programs", None)
    all_scores = getattr(result, "all_scores", None)
    if compiled is None or all_scores is None or len(compiled) == 0 or len(all_scores) == 0:
        return None
    rows = []
    clusters = pd.Index(all_scores["cluster"].astype(str).unique())
    for cluster in clusters:
        ann = decide_cluster_annotation_from_scores(
            cluster=str(cluster),
            root_programs=compiled,
            all_scores=all_scores,
            score_threshold=score_threshold,
            margin_threshold=margin_threshold,
            branch_support_parent_weight=branch_support_parent_weight,
            branch_support_descendant_weight=branch_support_descendant_weight,
            level1_branch_support_parent_weight=level1_branch_support_parent_weight,
            level1_branch_support_descendant_weight=level1_branch_support_descendant_weight,
            branch_routing_raw_weight=branch_routing_raw_weight,
            weak_raw_score_threshold=weak_raw_score_threshold,
            min_branch_supported_raw_score=float(min_branch_supported_raw_score),
            min_leaf_raw_score=min_leaf_raw_score,
        )
        rows.append(ann.__dict__)
    if not rows:
        return None
    raw = pd.DataFrame(rows)
    expanded = pd.DataFrame(raw)
    # expand dict columns same style as pipeline
    for col in ["level_labels", "level_scores", "level_raw_scores", "level_branch_supported_raw_scores", "level_margins", "level_confidence"]:
        if col in expanded.columns:
            mapped = expanded[col].apply(lambda x: x if isinstance(x, dict) else {})
            wide = pd.DataFrame(list(mapped))
            rename_prefix = {
                "level_labels": "_label",
                "level_scores": "_score",
                "level_raw_scores": "_raw_score",
                "level_branch_supported_raw_scores": "_branch_supported_raw_score",
                "level_margins": "_margin",
                "level_confidence": "_confidence",
            }[col]
            wide.columns = [str(c) + rename_prefix for c in wide.columns]
            expanded = pd.concat([expanded.drop(columns=[col]), wide], axis=1)
    out = pd.DataFrame({
        "cluster_id": expanded["cluster"].astype(str),
        "annot_label": expanded["final_label"],
        "annot_label_with_cluster": expanded["final_label"].astype(str).map(_sanitize_label_text) + "_" + expanded["cluster"].astype(str),
        "annot_path": expanded["final_path"],
        "annot_level": expanded["final_level"],
        "annot_status": expanded["status"],
        "annot_stop_reason": expanded["stop_reason"],
        "annot_confidence": expanded["confidence"],
        "annot_score": expanded["final_score"],
        "annot_raw_score": expanded["final_raw_score"],
        "annot_branch_supported_raw_score": expanded["final_branch_supported_raw_score"],
        "annot_margin": expanded["final_margin"],
        "annot_best_any_level_label": expanded["best_label_any_level"],
        "annot_best_any_level_score": expanded["best_score_any_level"],
        "annot_low_evidence": expanded["confidence"].astype(str).isin(["low", "none"]) | expanded["final_score"].isna() | expanded["final_level"].fillna(0).astype(int).eq(0),
        "annot_low_absolute_support": pd.to_numeric(expanded["final_branch_supported_raw_score"], errors="coerce").lt(float(min_branch_supported_raw_score)).fillna(True),
        "annot_deeper_level_unresolved": expanded["stop_reason"].astype(str).isin(["below_score_threshold", "ambiguous_sibling_margin", "below_absolute_support"]) & expanded["final_level"].fillna(0).astype(int).gt(0),
        "annot_branch_conflict": False,
        "annot_decision_label": expanded["decision_label"],
        "annot_decision_path": expanded["decision_path"],
        "annot_decision_level": expanded["decision_level"],
        "annot_decision_score": expanded["decision_score"],
        "annot_decision_raw_score": expanded["decision_raw_score"],
        "annot_decision_branch_supported_raw_score": expanded["decision_branch_supported_raw_score"],
        "annot_backoff_applied": expanded["backoff_applied"],
        "annot_backoff_steps": expanded["backoff_steps"],
        "annot_backoff_reason": expanded["backoff_reason"],
    })
    keep_cols = [c for c in expanded.columns if c.startswith("level_")]
    if keep_cols:
        rename = {}
        for c in keep_cols:
            m = re.match(r"level_(\d+)_(.+)", c)
            if m:
                lvl, suffix = m.groups()
                rename[c] = f"annot_l{int(lvl)}_{suffix}"
        out = out.join(expanded[keep_cols].rename(columns=rename))
    # derive top-level competition from all_scores
    l1 = all_scores[all_scores["level"] == 1].copy()
    if not l1.empty:
        score_col = "decision_score" if "decision_score" in l1.columns else "score"
        purity_rows = []
        for cluster, sub in l1.groupby("cluster", sort=False):
            sub = sub.dropna(subset=[score_col]).sort_values(score_col, ascending=False)
            top1 = sub.iloc[0] if len(sub) >= 1 else None
            top2 = sub.iloc[1] if len(sub) >= 2 else None
            gap = float(top1[score_col] - top2[score_col]) if top1 is not None and top2 is not None else np.nan
            purity_rows.append({
                "cluster_id": str(cluster),
                "annot_primary_compartment": str(top1["label"]) if top1 is not None else np.nan,
                "annot_secondary_compartment": str(top2["label"]) if top2 is not None else np.nan,
                "annot_primary_compartment_score": float(top1[score_col]) if top1 is not None else np.nan,
                "annot_secondary_compartment_score": float(top2[score_col]) if top2 is not None else np.nan,
                "annot_top_level_gap": gap,
            })
        purity_df = pd.DataFrame(purity_rows)
        out = out.merge(purity_df, on="cluster_id", how="left")
    return out


def _flatten_compiled_programs_for_export(programs) -> List:
    out = []
    def visit(node):
        out.append(node)
        for child in getattr(node, "children", []) or []:
            visit(child)
    for root in programs or []:
        visit(root)
    return out


def _relationship_maps_for_export(programs) -> tuple[Dict[str, Set[str]], Dict[str, Set[str]]]:
    flat = _flatten_compiled_programs_for_export(programs)
    by_name = {p.name: p for p in flat}
    children_lookup = {p.name: [c.name for c in getattr(p, "children", []) or []] for p in flat}

    descendants: Dict[str, Set[str]] = {}
    def gather_desc(name: str) -> Set[str]:
        out: Set[str] = set()
        for child_name in children_lookup.get(name, []):
            out.add(child_name)
            out |= gather_desc(child_name)
        return out

    for name in by_name:
        descendants[name] = gather_desc(name)

    ancestors: Dict[str, Set[str]] = {name: set() for name in by_name}
    for parent_name, child_names in children_lookup.items():
        for child_name in child_names:
            ancestors.setdefault(child_name, set()).add(parent_name)
    changed = True
    while changed:
        changed = False
        for name in list(ancestors):
            expanded = set(ancestors[name])
            for anc in list(ancestors[name]):
                expanded |= ancestors.get(anc, set())
            if expanded != ancestors[name]:
                ancestors[name] = expanded
                changed = True
    return ancestors, descendants



def _attach_final_call_competition(
    summary: pd.DataFrame,
    all_scores: pd.DataFrame,
    compiled_programs,
    *,
    cluster_id_column: str = "cluster_id",
    final_label_column: str = "annot_label",
    final_score_column: str = "annot_score",
    path_margin_column: str = "annot_margin",
    overwrite: bool = True,
) -> pd.DataFrame:
    if summary is None or summary.empty:
        return summary
    out = summary.copy()
    if all_scores is None or len(all_scores) == 0:
        if "annot_path_margin" not in out.columns:
            out["annot_path_margin"] = out.get(path_margin_column, np.nan)
        if "annot_second_best_call" not in out.columns:
            out["annot_second_best_call"] = np.nan
        if "annot_second_best_call_score" not in out.columns:
            out["annot_second_best_call_score"] = np.nan
        if "annot_final_call_margin" not in out.columns:
            out["annot_final_call_margin"] = np.nan
        if "annot_final_call_margin" in out.columns:
            out["annot_margin"] = out["annot_final_call_margin"]
        return out

    ancestors, descendants = _relationship_maps_for_export(compiled_programs)
    score_col = "branch_support_score" if "branch_support_score" in all_scores.columns else ("decision_score" if "decision_score" in all_scores.columns else "score")

    records = []
    for _, row in out.iterrows():
        cluster = str(row[cluster_id_column]) if cluster_id_column in out.columns else str(row.get("cluster"))
        final_label = str(row.get(final_label_column, row.get("final_label", "Unresolved")))
        final_score = pd.to_numeric(pd.Series([row.get(final_score_column, np.nan)]), errors="coerce").iloc[0]

        sub = all_scores[all_scores["cluster"].astype(str) == cluster].copy()
        if sub.empty:
            records.append({
                cluster_id_column: cluster,
                "annot_second_best_call": np.nan,
                "annot_second_best_call_score": np.nan,
                "annot_final_call_margin": np.nan,
            })
            continue

        sub[score_col] = pd.to_numeric(sub[score_col], errors="coerce")
        sub = sub[np.isfinite(sub[score_col])].copy()
        if sub.empty:
            records.append({
                cluster_id_column: cluster,
                "annot_second_best_call": np.nan,
                "annot_second_best_call_score": np.nan,
                "annot_final_call_margin": np.nan,
            })
            continue

        excluded = {final_label} | ancestors.get(final_label, set()) | descendants.get(final_label, set())
        alt = sub[~sub["label"].astype(str).isin(excluded)].sort_values(score_col, ascending=False)
        if alt.empty:
            second_label = np.nan
            second_score = np.nan
            final_margin = np.nan
        else:
            second = alt.iloc[0]
            second_label = str(second["label"])
            second_score = float(second[score_col])
            final_margin = float(final_score - second_score) if np.isfinite(final_score) and np.isfinite(second_score) else np.nan

        records.append({
            cluster_id_column: cluster,
            "annot_second_best_call": second_label,
            "annot_second_best_call_score": second_score,
            "annot_final_call_margin": final_margin,
        })

    comp = pd.DataFrame(records)
    refresh_cols = ["annot_second_best_call", "annot_second_best_call_score", "annot_final_call_margin"]
    out = out.copy()
    if cluster_id_column in out.columns:
        out[cluster_id_column] = out[cluster_id_column].astype(str)
    if cluster_id_column in comp.columns:
        comp[cluster_id_column] = comp[cluster_id_column].astype(str)
    if overwrite:
        existing_refresh = [c for c in refresh_cols if c in out.columns]
        if existing_refresh:
            out = out.drop(columns=existing_refresh)
    out = out.merge(comp, on=cluster_id_column, how="left")
    if "annot_path_margin" not in out.columns:
        out["annot_path_margin"] = out.get(path_margin_column, np.nan)
    if "annot_final_call_margin" in out.columns:
        out["annot_margin"] = out["annot_final_call_margin"]
    else:
        out["annot_final_call_margin"] = np.nan
        out["annot_margin"] = np.nan
    return out


def _compute_final_call_competition(summary: pd.DataFrame, result) -> pd.DataFrame:
    all_scores = getattr(result, "all_scores", None)
    compiled = getattr(result, "compiled_programs", None)
    return _attach_final_call_competition(
        summary=summary,
        all_scores=all_scores,
        compiled_programs=compiled,
        cluster_id_column="cluster_id",
        final_label_column="annot_label",
        final_score_column="annot_score",
        path_margin_column="annot_margin",
        overwrite=True,
    )


def _attach_best_any_level_raw_score(summary: pd.DataFrame, result) -> pd.DataFrame:
    """Attach annot_best_any_level_raw_score from result.all_scores when possible."""
    if summary is None or summary.empty:
        return summary
    if "annot_best_any_level_raw_score" in summary.columns:
        return summary
    required = {"cluster_id", "annot_best_any_level_label"}
    if not required.issubset(set(summary.columns)):
        return summary

    all_scores = getattr(result, "all_scores", None)
    if all_scores is None or len(all_scores) == 0:
        out = summary.copy()
        out["annot_best_any_level_raw_score"] = np.nan
        return out
    if not {"cluster", "label", "score"}.issubset(set(all_scores.columns)):
        out = summary.copy()
        out["annot_best_any_level_raw_score"] = np.nan
        return out

    out = summary.copy()
    lookup = all_scores[["cluster", "label", "score"]].copy()
    lookup["cluster"] = lookup["cluster"].astype(str)
    lookup["label"] = lookup["label"].astype(str)
    lookup = lookup.rename(columns={
        "cluster": "cluster_id",
        "label": "annot_best_any_level_label",
        "score": "annot_best_any_level_raw_score",
    })
    lookup = lookup.drop_duplicates(subset=["cluster_id", "annot_best_any_level_label"], keep="first")

    out["cluster_id"] = out["cluster_id"].astype(str)
    out["annot_best_any_level_label"] = out["annot_best_any_level_label"].astype(str)

    overlap = [c for c in ["annot_best_any_level_raw_score"] if c in out.columns]
    if overlap:
        out = out.drop(columns=overlap)
    out = out.merge(lookup, on=["cluster_id", "annot_best_any_level_label"], how="left")
    return out



def _build_export_labels(
    summary: pd.DataFrame,
    label_source_column: str,
    label_with_cluster_source_column: str,
    label_sep: str = "_",
    sanitize_label: bool = True,
    unknown_on_low_confidence: bool = True,
    unknown_confidence_levels: Optional[Iterable[str]] = None,
    unknown_label: str = "unknown",
    unknown_min_score: Optional[float] = None,
    unknown_max_margin: Optional[float] = None,
    unknown_score_column: str = "annot_score",
    unknown_margin_column: str = "annot_final_call_margin",
) -> pd.DataFrame:
    out = summary.copy()
    if label_source_column not in out.columns:
        raise KeyError(f"label_source_column '{label_source_column}' not found in cluster summary")
    levels = _normalize_confidence_levels(unknown_confidence_levels)
    export_label = out[label_source_column].astype(str).copy()
    if unknown_on_low_confidence and "annot_confidence" in out.columns:
        mask = out["annot_confidence"].astype(str).str.lower().isin(levels)
        export_label.loc[mask] = unknown_label
    if unknown_min_score is not None and unknown_score_column in out.columns:
        score_mask = pd.to_numeric(out[unknown_score_column], errors="coerce") < float(unknown_min_score)
        export_label.loc[score_mask.fillna(False)] = unknown_label
    if unknown_max_margin is not None and unknown_margin_column in out.columns:
        margin_mask = pd.to_numeric(out[unknown_margin_column], errors="coerce") < float(unknown_max_margin)
        export_label.loc[margin_mask.fillna(False)] = unknown_label
    out["annot_export_label"] = export_label
    if label_with_cluster_source_column in out.columns:
        cluster_ids = out["cluster_id"].astype(str)
        base = export_label.copy()
        if sanitize_label:
            base = base.map(_sanitize_label_text)
        out["annot_export_label_with_cluster"] = base + label_sep + cluster_ids
    return out


def make_cluster_annotation_export_summary(
    result,
    label_with_cluster: bool = True,
    unknown_on_low_confidence: bool = True,
    unknown_confidence_levels=("low", "none"),
    unknown_on_branch_conflict: bool = False,
    unknown_on_low_absolute_support: bool = False,
    unknown_branch_raw_threshold: Optional[float] = None,
    mixed_on_parent_mixing: bool = True,
    mixed_label_prefix: str = "mixed",
    mixed_branch_separator: str = ".",
    mixed_min_score: float = 0.35,
    unresolved_candidate_min_score: float = 1.0,
    unresolved_candidate_min_raw_score: float = 0.2,
    unresolved_candidate_label_prefix: str = "candidate",
    unknown_label: str = "unknown",
    unknown_min_score=None,
    unknown_max_margin=None,
    unknown_score_column: str = "annot_score",
    unknown_margin_column: str = "annot_final_call_margin",
    rescue_unknown_with_blocked_candidates: bool = True,
    rerun_decision: bool = False,
    score_threshold: float = 0.05,
    margin_threshold: float = 0.02,
    branch_child_rescue_weight: float = 0.35,
    level1_branch_child_rescue_weight: float = 0.40,
    branch_routing_raw_weight: float = 0.60,
    weak_raw_score_threshold: float = 0.10,
    min_branch_supported_raw_score: float = 0.0,
    min_leaf_raw_score: float = -0.2,
):
    """
    Build an export-ready cluster annotation summary from a HierAnnot result.

    This is the main public workflow helper for turning `result.cluster_annotations`
    into a flat cluster-level table with export labels, export status, optional
    rerun-based decision refresh, and extra diagnostics.

    Parameters controlling export masking
    ------------------------------------
    label_with_cluster
        If True, also create `annot_export_label_with_cluster` by appending the
        cluster id to the export label.

    unknown_on_low_confidence
        If True, mark clusters as unknown when `annot_low_evidence` is True.

    unknown_confidence_levels
        Reserved compatibility argument. The current export helper relies on the
        precomputed low-evidence flag in `cluster_annotations`.

    unknown_on_branch_conflict
        If True, mark clusters as unknown when `annot_branch_conflict` is True.

    unknown_on_low_absolute_support
        If True, mark clusters as unknown when the precomputed
        `annot_low_absolute_support` flag is True.

    unknown_branch_raw_threshold
        Optional export-time threshold on `annot_branch_supported_raw_score`.
        When provided, clusters below this threshold are exported as unknown
        regardless of `unknown_on_low_absolute_support`.

    unknown_label
        Export label used for unknown clusters.

    unknown_min_score
        Optional export-time lower bound on `unknown_score_column`. Clusters
        below this threshold are exported as unknown.

    unknown_max_margin
        Optional export-time upper bound on `unknown_margin_column`. Clusters
        below this ambiguity margin are exported as unknown.

    unknown_score_column
        Column used with `unknown_min_score`. Default is `annot_score`.

    unknown_margin_column
        Column used with `unknown_max_margin`. Default is
        `annot_final_call_margin`.

    rescue_unknown_with_blocked_candidates
        If True, allow blocked strong candidates to rescue unknown clusters.

    Parameters controlling mixed export labels
    -----------------------------------------
    mixed_on_parent_mixing
        If True, enable mixed-label export based on final-call competition.

    mixed_label_prefix, mixed_branch_separator
        Formatting controls for mixed export labels.

    mixed_min_score
        Minimum score required for both the final call and second-best call
        before a mixed label is emitted.

    Parameters controlling blocked-candidate hints
    ---------------------------------------------
    unresolved_candidate_min_score
        Minimum best-any-level score required to surface a blocked strong
        candidate label. This hint is only applied to unresolved-at-root cases
        and only when `annot_best_any_level_label`,
        `annot_best_any_level_score`, and `annot_best_any_level_raw_score`
        are all present.

    unresolved_candidate_min_raw_score
        Minimum raw score required for the blocked strong candidate.

    unresolved_candidate_label_prefix
        Prefix used for short blocked-candidate labels such as
        `candidate_plasma_cell`.

    Parameters controlling rerun
    ---------------------------
    rerun_decision
        If False, use the existing decision stored in `result.cluster_annotations`.
        If True, rerun hierarchical decision logic using the thresholds and
        routing weights passed to this function.

    score_threshold, margin_threshold
        Decision thresholds used only when `rerun_decision=True`.

    branch_child_rescue_weight, level1_branch_child_rescue_weight
        Descendant-rescue weights used only when `rerun_decision=True`.

    branch_routing_raw_weight, weak_raw_score_threshold
        Local evidence and gating controls used only when `rerun_decision=True`.

    min_branch_supported_raw_score
        Minimum branch-supported raw score required by the rerun decision logic.

    min_leaf_raw_score
        Minimum local raw score required for the selected node when
        `rerun_decision=True`. This is a node-level gate, distinct from
        branch-supported rescue and `min_branch_supported_raw_score`. The
        default `-0.2` allows routing to continue through mildly weak parent
        nodes when child-supported branch evidence is strong.
    """
    
    summary = None
    if rerun_decision:
        summary = _rerun_decision_summary_from_scores(
            result,
            score_threshold=score_threshold,
            margin_threshold=margin_threshold,
            branch_support_parent_weight=(1.0 - float(branch_child_rescue_weight)),
            branch_support_descendant_weight=float(branch_child_rescue_weight),
            level1_branch_support_parent_weight=(1.0 - float(level1_branch_child_rescue_weight)),
            level1_branch_support_descendant_weight=float(level1_branch_child_rescue_weight),
            branch_routing_raw_weight=branch_routing_raw_weight,
            weak_raw_score_threshold=weak_raw_score_threshold,
            min_branch_supported_raw_score=float(min_branch_supported_raw_score),
            min_leaf_raw_score=float(min_leaf_raw_score),
        )
    if summary is None:
        cluster_summary = result.cluster_annotations.copy()
        per_level_summary = _build_per_level_summary(result)
        if per_level_summary is not None and len(per_level_summary) > 0:
            if "cluster_id" in cluster_summary.columns and "cluster_id" in per_level_summary.columns:
                per_level_summary = per_level_summary.drop(columns=["cluster_id"])
            summary = cluster_summary.join(per_level_summary, how="left")
        else:
            summary = cluster_summary

    summary = _compute_final_call_competition(summary, result)
    summary = _attach_best_any_level_raw_score(summary, result)

    export_labels = _build_export_labels(
        summary,
        label_source_column="annot_label",
        label_with_cluster_source_column="annot_label_with_cluster",
        unknown_on_low_confidence=unknown_on_low_confidence,
        unknown_confidence_levels=unknown_confidence_levels,
        unknown_label=unknown_label,
        unknown_min_score=unknown_min_score,
        unknown_max_margin=unknown_max_margin,
        unknown_score_column=unknown_score_column,
        unknown_margin_column=unknown_margin_column,
    )

    export_status = pd.Series("resolved", index=summary.index, dtype=object)
    export_reason = pd.Series("", index=summary.index, dtype=object)

    unknown_mask = pd.Series(False, index=summary.index)
    if unknown_on_low_confidence and "annot_low_evidence" in summary.columns:
        unknown_mask = unknown_mask | summary["annot_low_evidence"].fillna(False).astype(bool)
    if unknown_on_branch_conflict and "annot_branch_conflict" in summary.columns:
        unknown_mask = unknown_mask | summary["annot_branch_conflict"].fillna(False).astype(bool)
    if unknown_on_low_absolute_support and "annot_low_absolute_support" in summary.columns:
        unknown_mask = unknown_mask | summary["annot_low_absolute_support"].fillna(False).astype(bool)
    if unknown_branch_raw_threshold is not None and "annot_branch_supported_raw_score" in summary.columns:
        unknown_mask = unknown_mask | (
            pd.to_numeric(summary["annot_branch_supported_raw_score"], errors="coerce") < float(unknown_branch_raw_threshold)
        ).fillna(False)
    if unknown_min_score is not None and unknown_score_column in summary.columns:
        unknown_mask = unknown_mask | (
            pd.to_numeric(summary[unknown_score_column], errors="coerce") < float(unknown_min_score)
        ).fillna(False)
    if unknown_max_margin is not None and unknown_margin_column in summary.columns:
        unknown_mask = unknown_mask | (
            pd.to_numeric(summary[unknown_margin_column], errors="coerce") < float(unknown_max_margin)
        ).fillna(False)

    mixed_mask = pd.Series(False, index=summary.index)
    required_mixed_cols = {
        "annot_final_call_margin",
        "annot_score",
        "annot_second_best_call_score",
        "annot_second_best_call",
    }
    if mixed_on_parent_mixing and required_mixed_cols.issubset(set(summary.columns)):
        final_margin = pd.to_numeric(summary["annot_final_call_margin"], errors="coerce")
        final_score = pd.to_numeric(summary["annot_score"], errors="coerce")
        second_score = pd.to_numeric(summary["annot_second_best_call_score"], errors="coerce")
        final_top = summary["annot_final_top_branch"].astype(str) if "annot_final_top_branch" in summary.columns else pd.Series("", index=summary.index)
        second_top = summary["annot_second_best_call"].astype(str).str.split(" > ").str[0]
        effective_margin_threshold = float(margin_threshold)
        mixed_mask = (
            final_margin.notna()
            & (final_margin <= effective_margin_threshold)
            & final_score.notna()
            & second_score.notna()
            & (final_score >= float(mixed_min_score))
            & (second_score >= float(mixed_min_score))
            & (final_top != second_top)
        )
        if "annot_parent_mixing" in summary.columns:
            mixed_mask = mixed_mask | summary["annot_parent_mixing"].fillna(False).astype(bool)

    # resolved by default; mixed can apply, but unknown explicitly overrides mixed at the end
    if mixed_mask.any():
        primary = summary.get("annot_label").astype(str).map(_sanitize_label_text).str.lower()
        secondary = summary.get("annot_second_best_call").astype(str).map(_sanitize_label_text).str.lower()
        mixed_label = mixed_label_prefix + "_" + primary + mixed_branch_separator + secondary
        export_status.loc[mixed_mask] = "mixed"
        export_reason.loc[mixed_mask] = "final_call_competition"
        export_labels.loc[mixed_mask, "annot_export_label"] = mixed_label[mixed_mask]
        if "annot_export_label_with_cluster" in export_labels.columns:
            cid = summary["cluster_id"].astype(str) if "cluster_id" in summary.columns else summary.index.astype(str)
            export_labels.loc[mixed_mask, "annot_export_label_with_cluster"] = mixed_label[mixed_mask] + "_" + cid[mixed_mask]

    export_status.loc[unknown_mask] = "unknown"
    export_reason.loc[unknown_mask] = "unknown_rule"
    export_labels.loc[unknown_mask, "annot_export_label"] = unknown_label

    if "annot_export_label_with_cluster" in export_labels.columns:
        cid = summary["cluster_id"].astype(str) if "cluster_id" in summary.columns else summary.index.astype(str)
        export_labels.loc[unknown_mask, "annot_export_label_with_cluster"] = unknown_label + "_" + cid[unknown_mask]

    # Surface blocked strong candidates only for true unresolved-at-root cases.
    required_candidate_cols = {
        "annot_best_any_level_label",
        "annot_best_any_level_score",
        "annot_best_any_level_raw_score",
    }
    if required_candidate_cols.issubset(set(summary.columns)):
        best_any_label = summary["annot_best_any_level_label"].astype(str)
        best_any_score = pd.to_numeric(summary["annot_best_any_level_score"], errors="coerce")
        best_any_raw = pd.to_numeric(summary["annot_best_any_level_raw_score"], errors="coerce")
        annot_status = summary["annot_status"].astype(str).str.lower() if "annot_status" in summary.columns else pd.Series("", index=summary.index)
        annot_label = summary["annot_label"].astype(str).str.lower() if "annot_label" in summary.columns else pd.Series("", index=summary.index)
        annot_level = pd.to_numeric(summary["annot_level"], errors="coerce") if "annot_level" in summary.columns else pd.Series(np.nan, index=summary.index)

        root_unresolved_mask = (
            (annot_level.isna() | annot_level.le(0))
            & (
                annot_label.eq("unresolved")
                | annot_status.isin(["unresolved", "stopped_at_parent"])
            )
        ) 

        # v0.7.55 allow candidate to rescue unknown ones (resovled but weak)
        if rescue_unknown_with_blocked_candidates:
            root_unresolved_mask = root_unresolved_mask | unknown_mask

        # Record candidate information for all root-unresolved rows, but only rewrite
        # export labels when the candidate thresholds are satisfied.
        export_labels["annot_unresolved_candidate_label"] = np.nan
        export_labels["annot_unresolved_candidate_score"] = np.nan
        export_labels["annot_unresolved_candidate_raw_score"] = np.nan

        if root_unresolved_mask.any():
            export_labels.loc[root_unresolved_mask, "annot_unresolved_candidate_label"] = best_any_label[root_unresolved_mask]
            export_labels.loc[root_unresolved_mask, "annot_unresolved_candidate_score"] = best_any_score[root_unresolved_mask]
            export_labels.loc[root_unresolved_mask, "annot_unresolved_candidate_raw_score"] = best_any_raw[root_unresolved_mask]

        candidate_mask = (
            root_unresolved_mask
            & best_any_label.notna()
            & best_any_label.ne("")
            & best_any_label.str.lower().ne("nan")
            & best_any_score.ge(float(unresolved_candidate_min_score)).fillna(False)
            & best_any_raw.ge(float(unresolved_candidate_min_raw_score)).fillna(False)
        )

        if candidate_mask.any():
            cand_label = best_any_label.map(_sanitize_label_text).str.lower()
            short_label = str(unresolved_candidate_label_prefix) + "_" + cand_label
            export_labels.loc[candidate_mask, "annot_export_label"] = short_label[candidate_mask]
            if "annot_export_label_with_cluster" in export_labels.columns:
                cid = summary["cluster_id"].astype(str) if "cluster_id" in summary.columns else summary.index.astype(str)
                export_labels.loc[candidate_mask, "annot_export_label_with_cluster"] = short_label[candidate_mask] + "_" + cid[candidate_mask]
            export_reason.loc[candidate_mask] = "blocked_strong_candidate"
           
    if not label_with_cluster and "annot_export_label_with_cluster" in export_labels.columns:
        export_labels = export_labels.drop(columns=["annot_export_label_with_cluster"])

    export_labels["annot_export_status"] = export_status
    export_labels["annot_export_reason"] = export_reason.replace("", "resolved")
    overlap = [c for c in export_labels.columns if c in summary.columns]
    if overlap:
        summary = summary.drop(columns=overlap)
    summary = summary.join(export_labels, how="left")
    return summary


def expand_cluster_annotation_to_cells(
    cluster_assignments,
    cluster_summary,
    cluster_key: str = "cluster",
    include_columns=None,
    prefix: str | None = None,
    label_with_cluster: bool = True,
    fill_unassigned: bool = True,
    unassigned_label: str = "unassigned",
    unassigned_status: str = "not_analyzed",
):
    """
    Expand a cluster-level annotation summary to a per-cell joinable DataFrame.

    Parameters
    ----------
    cluster_assignments
        Array-like or pandas Series of cluster IDs for each cell/observation.
    cluster_summary
        Complete cluster-level summary table, typically produced by
        `make_cluster_annotation_export_summary()`.
    include_columns
        Optional list of columns from `cluster_summary` to propagate. If omitted,
        uses the common export-label workflow: `annot_export_label` and, when
        `label_with_cluster=True`, `annot_export_label_with_cluster`.
    prefix
        Optional prefix applied to propagated annotation columns in the returned
        table. Leading `annot_` is stripped before adding the custom prefix.
    """

    if not isinstance(cluster_assignments, pd.Series):
        cluster_assignments = pd.Series(cluster_assignments, name=cluster_key)
    else:
        cluster_assignments = cluster_assignments.copy()
        if cluster_assignments.name is None:
            cluster_assignments.name = cluster_key

    summary = cluster_summary.copy()

    if "cluster_id" not in summary.columns:
        summary["cluster_id"] = summary.index

    if include_columns is None:
        include_columns = ["annot_export_label", "annot_export_status"]
    else:
        include_columns = list(include_columns)

    if label_with_cluster and "annot_export_label_with_cluster" in summary.columns:
        if "annot_export_label" in include_columns and "annot_export_label_with_cluster" not in include_columns:
            include_columns.append("annot_export_label_with_cluster")
        elif include_columns == []:
            include_columns = ["annot_export_label_with_cluster", "annot_export_status"]

    include_columns = [c for c in include_columns if c in summary.columns]
    summary_small = summary[["cluster_id"] + include_columns].copy()

    # Merge on a normalized string view of cluster ids to avoid dtype mismatch
    join_df = pd.DataFrame({cluster_key: cluster_assignments.values}, index=cluster_assignments.index)
    join_df["_cluster_merge_key"] = join_df[cluster_key].astype(str)
    summary_small["_cluster_merge_key"] = summary_small["cluster_id"].astype(str)

    expanded = join_df.merge(summary_small, how="left", on="_cluster_merge_key")
    expanded.index = join_df.index

    if fill_unassigned:
        if "annot_export_label" in expanded.columns:
            expanded["annot_export_label"] = expanded["annot_export_label"].fillna(unassigned_label)
        if "annot_export_label_with_cluster" in expanded.columns:
            missing = expanded["annot_export_label_with_cluster"].isna()
            expanded.loc[missing, "annot_export_label_with_cluster"] = (
                unassigned_label + "_" + expanded.loc[missing, cluster_key].astype(str)
            )
        if "annot_status" in expanded.columns:
            expanded["annot_status"] = expanded["annot_status"].fillna(unassigned_status)

    expanded = expanded.drop(columns=["cluster_id", "_cluster_merge_key"], errors="ignore")

    if prefix:
        rename_map = {}
        for c in expanded.columns:
            if c == cluster_key:
                continue
            base = c[len("annot_"):] if c.startswith("annot_") else c
            rename_map[c] = f"{prefix}_{base}"
        expanded = expanded.rename(columns=rename_map)

    return expanded

