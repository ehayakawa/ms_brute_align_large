"""
Module for detecting communities in feature similarity graphs.

This module implements community detection algorithms to group mass features
that likely represent the same chemical compound across different datasets.
It uses the Louvain algorithm for community detection, providing a 'soft'
grouping approach that allows flexible feature associations.

Main functions/classes:
    - detect_communities: Applies Louvain algorithm to find feature communities
    - group_features_by_community: Groups features based on detected communities
    - detect_cliques: Finds maximal cliques for stricter grouping (deprecated, use clique_detection.py)
    - group_features_by_clique: Groups features based on clique membership

Inputs:
    - NetworkX graph with features as nodes and similarities as edges
    - Optional parameters for community detection algorithm

Outputs:
    - Dictionary mapping community IDs to lists of features
    - Community membership information for each node

Important arguments:
    - G: NetworkX graph from graph_construction module
    - resolution: Resolution parameter for Louvain algorithm (higher = smaller communities)
"""
import networkx as nx
import community as community_louvain
import pandas as pd
import numpy as np
from collections import defaultdict
import random
import logging

# Configure logger for this module
logger = logging.getLogger(__name__)

def detect_communities(G, resolution=1.0, hard_separation=False):
    """
    Detect communities in the graph using the Louvain algorithm.
    
    Parameters:
    -----------
    G : networkx.Graph
        Graph to detect communities in
    resolution : float
        Resolution parameter for the Louvain algorithm. Higher values lead to smaller communities.
    hard_separation : bool
        If True, use a higher resolution and post-process communities to ensure hard separation
        
    Returns:
    --------
    partition : dict
        Dictionary mapping node IDs to community IDs
    """
    if G.number_of_nodes() == 0:
        logger.warning("Empty graph, no communities to detect")
        return {}
    
    # Use higher resolution if hard separation is requested
    if hard_separation:
        resolution = 1.5  # Higher resolution for smaller, more distinct communities
    
    logger.info(f"Detecting communities with resolution={resolution}...")
    
    # Apply Louvain algorithm
    partition = community_louvain.best_partition(G, resolution=resolution)
    
    # Post-process communities if hard separation is requested
    if hard_separation:
        logger.info("Applying hard separation to communities...")
        partition = refine_communities_by_mz_rt(G, partition)
    
    # Count communities
    communities = set(partition.values())
    logger.info(f"Detected {len(communities)} communities")
    
    # Print sizes of largest communities
    community_sizes = {}
    for comm_id in partition.values():
        community_sizes[comm_id] = community_sizes.get(comm_id, 0) + 1
    
    sorted_sizes = sorted(community_sizes.items(), key=lambda x: x[1], reverse=True)
    logger.info("Largest communities:")
    for i, (comm_id, size) in enumerate(sorted_sizes[:5]):
        logger.info(f"  Community {comm_id}: {size} nodes")
    
    return partition

def refine_communities_by_mz_rt(G, partition, mz_tolerance=0.01, rt_tolerance=0.5):
    """
    Refine communities based on m/z and RT values to ensure they represent distinct chemical entities.
    
    Parameters:
    -----------
    G : networkx.Graph
        Graph with nodes representing features
    partition : dict
        Dictionary mapping node IDs to community IDs
    mz_tolerance : float
        m/z tolerance for considering features as the same entity
    rt_tolerance : float
        RT tolerance for considering features as the same entity
        
    Returns:
    --------
    refined_partition : dict
        Refined partition with communities split if necessary
    """
    logger.info("Refining communities based on m/z and RT values...")
    
    # Group nodes by community
    communities = {}
    for node, comm_id in partition.items():
        if comm_id not in communities:
            communities[comm_id] = []
        communities[comm_id].append(node)
    
    # Check each community for m/z and RT consistency
    refined_partition = partition.copy()
    next_comm_id = max(partition.values()) + 1 if partition else 0
    
    for comm_id, nodes in communities.items():
        if len(nodes) <= 3:  # Skip very small communities
            continue
            
        # Extract m/z and RT values
        node_data = []
        for node in nodes:
            if 'precursor_mz' in G.nodes[node] and 'retention_time' in G.nodes[node]:
                node_data.append({
                    'node': node,
                    'mz': G.nodes[node]['precursor_mz'],
                    'rt': G.nodes[node]['retention_time']
                })
        
        if len(node_data) <= 3:  # Skip if not enough data
            continue
            
        # Calculate m/z and RT ranges
        mz_values = [data['mz'] for data in node_data]
        rt_values = [data['rt'] for data in node_data]
        
        mz_range = max(mz_values) - min(mz_values)
        rt_range = max(rt_values) - min(rt_values)
        
        # If ranges exceed tolerances, split the community
        if mz_range > 3 * mz_tolerance or rt_range > 3 * rt_tolerance:
            logger.debug(f"Splitting community {comm_id} (mz_range={mz_range:.4f}, rt_range={rt_range:.2f})")
            
            # Use a simple clustering approach based on m/z values
            # Sort nodes by m/z
            node_data.sort(key=lambda x: x['mz'])
            
            # Find gaps in m/z values
            current_cluster = [node_data[0]]
            clusters = [current_cluster]
            
            for i in range(1, len(node_data)):
                if node_data[i]['mz'] - node_data[i-1]['mz'] > mz_tolerance:
                    # Start a new cluster
                    current_cluster = [node_data[i]]
                    clusters.append(current_cluster)
                else:
                    current_cluster.append(node_data[i])
            
            # Assign new community IDs to clusters (except the first one)
            for i, cluster in enumerate(clusters):
                if i == 0:
                    continue  # Keep the first cluster in the original community
                
                # Assign a new community ID to this cluster
                for data in cluster:
                    refined_partition[data['node']] = next_comm_id
                
                logger.debug(f"  Created new community {next_comm_id} with {len(cluster)} nodes")
                next_comm_id += 1
    
    # Count communities after refinement
    communities = set(refined_partition.values())
    logger.info(f"After refinement: {len(communities)} communities")
    
    return refined_partition

def group_features_by_community(G, partition):
    """
    Group features by community.
    
    Parameters:
    -----------
    G : networkx.Graph
        Graph with features as nodes
    partition : dict
        Dictionary mapping node IDs to community IDs
        
    Returns:
    --------
    grouped_features : dict
        Dictionary mapping community IDs to lists of features
    """
    logger.info("Grouping features by community...")
    
    # Group features by community
    grouped_features = defaultdict(list)
    for node, comm_id in partition.items():
        # Get feature data from node attributes
        feature_data = G.nodes[node]
        
        # Create a feature dictionary with all relevant information
        feature = {
            'node_id': node,
            'dataset_id': feature_data.get('dataset_id'),
            'feature_id': feature_data.get('feature_id'),
            'mz': feature_data.get('mz'),
            'rt': feature_data.get('rt'),
            'intensity': feature_data.get('intensity'),
            'filename': feature_data.get('filename')
        }
        
        grouped_features[comm_id].append(feature)
    
    # Filter communities to ensure at most one feature per dataset per community
    filtered_grouped_features = {}
    for comm_id, features in grouped_features.items():
        # Group features by dataset_id
        features_by_dataset = defaultdict(list)
        for feature in features:
            dataset_id = feature.get('dataset_id')
            features_by_dataset[dataset_id].append(feature)
        
        # Select the best feature from each dataset (highest intensity)
        filtered_features = []
        for dataset_id, dataset_features in features_by_dataset.items():
            if dataset_features:
                # Sort by intensity (descending) and take the highest
                best_feature = sorted(dataset_features, key=lambda x: x.get('intensity', 0), reverse=True)[0]
                filtered_features.append(best_feature)
        
        # Only keep the community if it has at least 2 features
        if len(filtered_features) >= 2:
            filtered_grouped_features[comm_id] = filtered_features
    
    # Sort communities by size
    sorted_communities = sorted(filtered_grouped_features.items(), key=lambda x: len(x[1]), reverse=True)
    
    # Create a new dictionary with sorted communities
    result = {}
    for i, (comm_id, features) in enumerate(sorted_communities):
        result[i] = features
    
    logger.info(f"Grouped features into {len(result)} communities after filtering")
    logger.info(f"Removed {len(grouped_features) - len(filtered_grouped_features)} communities that didn't meet criteria")
    
    return result

def generate_community_tables(G, partition, all_list_features):
    """
    Generate tables of aligned features by community.
    
    Parameters:
    -----------
    G : networkx.Graph
        Graph with nodes representing features and edges representing similarity
    partition : dict
        Dictionary mapping node IDs to community IDs
    all_list_features : list of tuples
        List of (filename, features) tuples
        
    Returns:
    --------
    aligned_features : dict
        Dictionary mapping community IDs to aligned feature groups
    feature_mzs : dict
        Dictionary mapping (community_id, dataset_id, feature_id) to m/z values
    """
    # Group features by community
    aligned_features = group_features_by_community(G, partition)
    
    # Get m/z values for each feature
    feature_mzs = {}
    for comm_id, features in aligned_features.items():
        for feature in features:
            # Get feature data
            _, features = all_list_features[feature['dataset_id']]
            feature_data = features[feature['feature_id']]
            
            # Store m/z value
            feature_mzs[(comm_id, feature['dataset_id'], feature['feature_id'])] = feature_data.get('mz', 0)
    
    return aligned_features, feature_mzs

def get_top_communities(partition, top_n=10):
    """
    Get the top N communities by size.
    
    Parameters:
    -----------
    partition : dict
        Dictionary mapping node IDs to community IDs
    top_n : int
        Number of top communities to return
        
    Returns:
    --------
    top_communities : list
        List of lists of node IDs for the top N communities
    """
    # Group nodes by community
    communities = {}
    for node, comm_id in partition.items():
        if comm_id not in communities:
            communities[comm_id] = []
        communities[comm_id].append(node)
    
    # Sort communities by size and return the top N
    top_communities = sorted(communities.items(), key=lambda x: len(x[1]), reverse=True)[:top_n]
    return [nodes for _, nodes in top_communities]

def detect_cliques(G, min_size=3):
    """
    Detect cliques in the graph.
    
    Parameters:
    -----------
    G : networkx.Graph
        Graph to detect cliques in
    min_size : int
        Minimum size of cliques to detect
        
    Returns:
    --------
    cliques : list
        List of lists of node IDs for cliques
    """
    logger.info("Detecting cliques...")
    
    # Check if the graph is empty
    if G.number_of_nodes() == 0:
        logger.warning("Empty graph, no cliques to detect")
        return []
    
    # Find all cliques of size at least min_size
    cliques = list(nx.find_cliques(G))
    
    # Filter cliques by size
    cliques = [clique for clique in cliques if len(clique) >= min_size]
    
    # Sort cliques by size (largest first)
    cliques.sort(key=len, reverse=True)
    
    logger.info(f"Detected {len(cliques)} cliques of size at least {min_size}")
    
    # Print the sizes of the largest cliques
    if cliques:
        logger.info("Largest cliques:")
        for i, clique in enumerate(cliques[:10]):
            logger.info(f"  Clique {i}: {len(clique)} nodes")
    
    return cliques

def group_features_by_clique(G, cliques):
    """
    Group features by clique.
    
    Parameters:
    -----------
    G : networkx.Graph
        Graph with features as nodes
    cliques : list
        List of lists of node IDs for cliques
        
    Returns:
    --------
    grouped_features : dict
        Dictionary mapping clique IDs to lists of features
    """
    logger.info("Grouping features by clique...")
    
    # Group features by clique
    grouped_features = {}
    
    # Create a mapping from node to cliques it belongs to
    node_to_cliques = defaultdict(list)
    for i, clique in enumerate(cliques):
        for node in clique:
            node_to_cliques[node].append(i)
    
    # For each node, assign it to the largest clique it belongs to
    node_to_largest_clique = {}
    for node, clique_ids in node_to_cliques.items():
        # Sort cliques by size (using the original cliques list)
        sorted_cliques = sorted(clique_ids, key=lambda x: len(cliques[x]), reverse=True)
        node_to_largest_clique[node] = sorted_cliques[0]
    
    # Group features by their assigned clique
    for node, clique_id in node_to_largest_clique.items():
        if clique_id not in grouped_features:
            grouped_features[clique_id] = []
        
        # Get feature data from node attributes
        feature_data = G.nodes[node]
        
        # Create a feature dictionary with all relevant information
        feature = {
            'node_id': node,
            'dataset_id': feature_data.get('dataset_id'),
            'feature_id': feature_data.get('feature_id'),
            'mz': feature_data.get('mz'),
            'rt': feature_data.get('rt'),
            'intensity': feature_data.get('intensity'),
            'filename': feature_data.get('filename')
        }
        
        grouped_features[clique_id].append(feature)
    
    # Sort cliques by size
    sorted_cliques = sorted(grouped_features.items(), key=lambda x: len(x[1]), reverse=True)
    
    # Create a new dictionary with sorted cliques
    result = {}
    for i, (_, features) in enumerate(sorted_cliques):
        result[i] = features
    
    logger.info(f"Grouped features into {len(result)} cliques")
    
    return result
