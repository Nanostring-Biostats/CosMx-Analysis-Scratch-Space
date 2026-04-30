import pandas as pd

from hierannot import HierAnnotPipeline, get_builtin_hierarchy

# cluster-mean expression: genes x clusters pandas.DataFrame
expr = pd.DataFrame(
    {
        "cluster_T": [9, 8, 0, 0, 0, 0, 0, 1, 1],
        "cluster_Myeloid": [0, 0, 8, 7, 0, 0, 0, 6, 5],
        "cluster_Epi": [0, 0, 0, 0, 8, 7, 6, 0, 0],
    },
    index=["CD3D", "TRBC1", "LYZ", "S100A8", "EPCAM", "KRT18", "KRT19", "FCER1G", "CTSS"],
)

pipeline = HierAnnotPipeline(
    root_programs=get_builtin_hierarchy("tme_core"),
    input_type="raw_cluster_means",
    score_threshold=0.01,
    margin_threshold=0.0,
    random_state=0,
)

result = pipeline.fit_score(expr)
print(result.cluster_annotations)
print(result.level_scores)
