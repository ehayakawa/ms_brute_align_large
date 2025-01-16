import pandas as pd
import networkx as nx

def find_cliques(G):
    return list(nx.find_cliques(G))

def generate_clique_tables(G, cliques, all_list_features):
    aligned_features = {}
    aligned_mz = {}
    aligned_intensity = {}

    for clique_index, clique in enumerate(cliques):
        aligned_features[clique_index] = {}
        aligned_mz[clique_index] = []
        aligned_intensity[clique_index] = []

        for node in clique:
            # Parse node ID string (format: "file_idx_feat_idx")
            file_index, feature_index = map(int, node.split('_'))
            feature = all_list_features[file_index][1][feature_index]
            
            aligned_features[clique_index][file_index] = feature_index
            aligned_mz[clique_index].append(feature['precursor_mz'])
            
            # Use feature.get() method with a default value of float('nan')
            intensity = feature.get('signal_intensity', float('nan'))
            aligned_intensity[clique_index].append(intensity)

    return aligned_features, aligned_mz, aligned_intensity
