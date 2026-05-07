"""
Module for detecting cliques in feature similarity graphs.

This module implements clique detection algorithms to find groups of mass features
where every feature is similar to every other feature in the group. This provides
a 'hard' grouping approach with strict all-to-all connectivity requirements.

Main functions/classes:
    - find_cliques: Finds all maximal cliques in the graph
    - generate_clique_tables: Creates structured output tables from clique results
    - filter_cliques_by_dataset: Ensures cliques span multiple datasets

Inputs:
    - NetworkX graph with features as nodes and similarities as edges
    - Minimum clique size and dataset requirements

Outputs:
    - List of cliques (sets of fully connected nodes)
    - Structured tables with clique membership and feature information

Important arguments:
    - G: NetworkX graph from graph_construction module
    - min_size: Minimum number of nodes required for valid clique
    - min_datasets: Minimum number of datasets required in each clique
"""
import pandas as pd
import networkx as nx
import numpy as np
from collections import defaultdict

def find_cliques(G):
    """
    Find maximal cliques in the graph.
    
    Parameters:
    -----------
    G : networkx.Graph
        Graph with nodes representing features and edges representing similarity
        
    Returns:
    --------
    filtered_cliques : list
        List of lists of node IDs for cliques with at least 3 nodes
    G : networkx.Graph
        The original graph (passed through for use in other functions)
    """
    print("Finding maximal cliques...")
    
    # Find all maximal cliques
    all_cliques = list(nx.find_cliques(G))
    
    # Filter cliques to only include those with at least 3 nodes
    filtered_cliques = [clique for clique in all_cliques if len(clique) >= 3]
    
    print(f"Found {len(filtered_cliques)} cliques with at least 3 nodes")
    
    # Print top 5 cliques by size
    top_cliques = sorted(filtered_cliques, key=len, reverse=True)[:5]
    for i, clique in enumerate(top_cliques):
        print(f"Clique {i}: {len(clique)} nodes")
    
    return filtered_cliques, G

def group_features_by_clique(cliques, G):
    """
    Group features by clique.
    
    Parameters:
    -----------
    cliques : list
        List of lists of node IDs for cliques
    G : networkx.Graph
        Graph with nodes representing features
        
    Returns:
    --------
    aligned_features : dict
        Dictionary mapping clique IDs to aligned feature groups
    """
    print("Grouping features by clique...")
    
    # Create aligned features dictionary
    aligned_features = {}
    
    for clique_id, nodes in enumerate(cliques):
        # Group nodes by dataset
        datasets = defaultdict(list)
        for node in nodes:
            # Get node attributes from the graph
            node_attrs = G.nodes[node]
            dataset_id = node_attrs.get('dataset_id')
            feature_id = node_attrs.get('feature_id')
            intensity = node_attrs.get('intensity', 0)
            
            datasets[dataset_id].append((feature_id, intensity))
        
        # Only keep one feature per dataset (the one with highest intensity)
        aligned_group = {}
        for dataset_id, feature_info in datasets.items():
            if len(feature_info) == 1:
                # If only one feature, use it
                aligned_group[dataset_id] = feature_info[0][0]
            else:
                # If multiple features from the same dataset, choose the one with highest intensity
                # Sort by intensity (descending) and take the highest
                sorted_features = sorted(feature_info, key=lambda x: x[1], reverse=True)
                aligned_group[dataset_id] = sorted_features[0][0]
        
        # Only add groups with at least 2 features
        if len(aligned_group) >= 2:
            aligned_features[clique_id] = aligned_group
    
    print(f"Created {len(aligned_features)} aligned feature groups from cliques")
    return aligned_features

def generate_clique_tables(G, cliques, all_list_features):
    """
    Generate tables of aligned features by clique.
    
    Parameters:
    -----------
    G : networkx.Graph
        Graph with nodes representing features and edges representing similarity
    cliques : list
        List of lists of node IDs for cliques
    all_list_features : list of tuples
        List of (filename, features) tuples
        
    Returns:
    --------
    aligned_features : dict
        Dictionary mapping clique IDs to aligned feature groups
    feature_mzs : dict
        Dictionary mapping (clique_id, dataset_id, feature_id) to m/z values
    """
    # Group features by clique
    aligned_features = group_features_by_clique(cliques, G)
    
    # Get m/z values for each feature
    feature_mzs = {}
    for clique_id, aligned_group in aligned_features.items():
        for dataset_id, feature_id in aligned_group.items():
            # Get feature data
            _, features = all_list_features[dataset_id]
            feature = features[feature_id]
            
            # Store m/z value
            feature_mzs[(clique_id, dataset_id, feature_id)] = feature.get('mz', 0)
    
    return aligned_features, feature_mzs

def get_top_cliques(cliques, top_n=10):
    """
    Get the top N cliques by size.
    
    Parameters:
    -----------
    cliques : list
        List of lists of node IDs for cliques
    top_n : int
        Number of top cliques to return
        
    Returns:
    --------
    top_cliques : list
        List of lists of node IDs for the top N cliques
    """
    # Sort cliques by size and return the top N
    top_cliques = sorted(cliques, key=len, reverse=True)[:top_n]
    return top_cliques
