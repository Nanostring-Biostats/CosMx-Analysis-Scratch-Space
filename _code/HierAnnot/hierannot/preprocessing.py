from __future__ import annotations

import numpy as np
import pandas as pd



def standardize_gene_index(expr: pd.DataFrame, uppercase: bool = True) -> pd.DataFrame:
    df = expr.copy()
    genes = df.index.astype(str).str.strip()
    if uppercase:
        genes = genes.str.upper()
    df.index = genes
    if df.index.has_duplicates:
        df = df.groupby(df.index, sort=False).mean()
    return df



def validate_expression_matrix(expr: pd.DataFrame) -> None:
    if not isinstance(expr, pd.DataFrame):
        raise TypeError("Expression matrix must be a pandas DataFrame.")
    if expr.empty:
        raise ValueError("Expression matrix is empty.")
    if expr.index.empty or expr.columns.empty:
        raise ValueError("Expression matrix must have non-empty gene index and cluster columns.")
    if not np.issubdtype(expr.to_numpy().dtype, np.number):
        raise TypeError("Expression matrix must be numeric.")



def log_normalize(expr: pd.DataFrame, scale_factor: float = 1e4) -> pd.DataFrame:
    arr = expr.to_numpy(dtype=float)
    libsize = arr.sum(axis=0)
    libsize[libsize == 0] = 1.0
    norm = np.log1p((arr / libsize) * scale_factor)
    return pd.DataFrame(norm, index=expr.index, columns=expr.columns)



def preprocess_expression(
    expr: pd.DataFrame,
    input_type: str = "raw_cluster_means",
    scale_factor: float = 1e4,
    uppercase_genes: bool = True,
) -> pd.DataFrame:
    validate_expression_matrix(expr)
    df = standardize_gene_index(expr, uppercase=uppercase_genes)
    df = df.loc[df.sum(axis=1) > 0]

    if input_type == "raw_cluster_means":
        return log_normalize(df, scale_factor=scale_factor)
    if input_type == "normalized":
        return df.astype(float)
    raise ValueError("input_type must be either 'raw_cluster_means' or 'normalized'.")
