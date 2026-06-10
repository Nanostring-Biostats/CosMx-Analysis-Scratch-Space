from __future__ import annotations

from typing import Optional, Sequence
from .utils import _sanitize_label_text

import numpy as np
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
