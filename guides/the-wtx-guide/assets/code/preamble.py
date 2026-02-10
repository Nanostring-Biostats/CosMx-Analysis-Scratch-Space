
import os
import pandas as pd
import scipy.sparse as sparse
import anndata as ad
from tqdm import tqdm
import hdf5plugin
import scipy as sp
import numpy as np
from scipy.sparse import issparse
from sklearn.neighbors import NearestNeighbors
import warnings
import scanpy as sc
import json


def load_jsons(rlf):
  if os.path.exists(rlf):
    with open(rlf, 'r') as f:
      rl = json.load(f)
  else:
      print("rlf does not exist. Creating one now.")
      rl = {}
  return rl

def save_jsons(rl, rlf):
  print("Saving results")
  with open(rlf, "w") as outfile:
    json.dump(rl, outfile) 

# formats ranked genes from marker analysis
def format_ranked_genes(Ann: ad.AnnData, n_genes=50, 
  key_use='rank_genes_groups'):
    # Ensure the key exists before trying to access it
    if 'rank_genes_groups' not in Ann.uns:
        print("Error: 'rank_genes_groups' not found in Ann.uns. Please run sc.tl.rank_genes_groups first.")
        return None, None
    ranked_genes_info = Ann.uns[key_use]
    group_names = ranked_genes_info['names'].dtype.names
    top_markers_with_scores_by_group = {}
    print(f"Extracting top {n_genes} marker genes for each group...")
    # This loop gathers the top markers for each cluster
    for group in group_names:
        # Extract all ranked gene names and scores for the current group
        all_genes_for_group = ranked_genes_info['names'][group]
        all_scores_for_group = ranked_genes_info['scores'][group]
        
        # Determine how many genes to take (in case a cluster has fewer than n_genes)
        num_to_take = min(n_genes, len(all_genes_for_group))
        
        # Slice to get top N genes and scores
        top_n_genes = all_genes_for_group[:num_to_take].tolist()
        top_n_scores = all_scores_for_group[:num_to_take].tolist()
        
        # Pair them up into (gene_name, score) tuples
        paired_markers = list(zip(top_n_genes, top_n_scores))
        top_markers_with_scores_by_group[group] = paired_markers
        
        # Print a preview to the console
        print_markers = [f"('{gene}', {score:.2f})" for gene, score in paired_markers[:min(5, len(paired_markers))]]
        print(f"  Group '{group}': [{', '.join(print_markers)}{'...' if len(paired_markers) > 5 else ''}]")
    
    print('Making DataFrame...')
    data_for_df = []
    
    # This loop unnests the dictionary into a long-format list for the DataFrame
    for group_id, markers in top_markers_with_scores_by_group.items():
        for marker_name, score_value in markers:
            data_for_df.append({
                'cluster': group_id,
                'marker': marker_name,
                'score': score_value
            })
            
    df_marker_genes = pd.DataFrame(data_for_df)
    
    return top_markers_with_scores_by_group, df_marker_genes
    

