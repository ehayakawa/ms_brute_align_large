"""
Module for constructing similarity graphs from mass spectrometry features.

This module builds a graph where nodes represent mass features and edges
represent similarity relationships between features from different datasets.
It uses KD-trees for efficient similarity searching based on m/z and retention
time tolerances.

Main functions/classes:
    - GraphBuilder: Main class for constructing feature similarity graphs
    - build_graph: Creates graph with nodes and edges based on feature similarity
    - clean_multiple_connections: Resolves ambiguous connections between datasets

Inputs:
    - List of features from multiple datasets (with mz, rt, intensity values)
    - mz_tolerance: Tolerance for m/z matching in Daltons
    - rt_tolerance: Tolerance for retention time matching in minutes

Outputs:
    - NetworkX graph with features as nodes and similarities as weighted edges
    - Edge weights based on combined m/z and RT similarity scores

Important arguments:
    - features: List of feature dictionaries from read_files module
    - mz_tolerance: Maximum allowed m/z difference (default: 0.01 Da)
    - rt_tolerance: Maximum allowed RT difference (default: 0.5 min)
"""
from typing import List, Dict
import networkx as nx
import numpy as np
from scipy.spatial import KDTree
import pandas as pd
from tqdm import tqdm
from collections import defaultdict
import random
import os
import logging
from spectral_similarity import calculate_spectral_similarity, has_msms_data

# Configure logger for this module
logger = logging.getLogger(__name__)

class GraphBuilder:
    """
    Class for building a graph from mass spectrometry features.
    
    Nodes represent features, and edges represent similarity between features.
    """
    
    def __init__(self, mz_tolerance=0.01, rt_tolerance=0.5, cosine_threshold=0.5, min_shared_peaks=3):
        """
        Initialize the GraphBuilder with tolerance parameters and MS/MS similarity settings.
        
        Parameters:
        -----------
        mz_tolerance : float
            Tolerance for m/z values (in Da)
        rt_tolerance : float
            Tolerance for retention time values (in minutes)
        cosine_threshold : float
            Minimum cosine similarity for MS/MS-based edges (default: 0.5)
        min_shared_peaks : int
            Minimum number of shared peaks required for MS/MS similarity (default: 3)
        """
        self.mz_tolerance = mz_tolerance
        self.rt_tolerance = rt_tolerance
        self.cosine_threshold = cosine_threshold
        self.min_shared_peaks = min_shared_peaks
        self.G = nx.Graph()
        
        logger.info(f"GraphBuilder initialized - mz_tol: {mz_tolerance}, rt_tol: {rt_tolerance}, "
                   f"cosine_threshold: {cosine_threshold}, min_shared_peaks: {min_shared_peaks}")
    
    def build_graph(self, all_list_features):
        """
        Build a graph from a list of features using two-case matching logic.
        
        Case 1: No MS/MS data - use m/z and RT similarity for weight
        Case 2: MS/MS data available for both - use cosine similarity for weight
        
        Parameters:
        -----------
        all_list_features : list
            List of tuples (filename, features) where features is a list of dictionaries
            
        Returns:
        --------
        G : networkx.Graph
            Graph with nodes representing features and edges representing similarity
        """
        logger.info("Building graph with two-case matching logic...")
        self.G = nx.Graph()
        
        # Statistics tracking
        mz_rt_edges = 0
        msms_edges = 0
        msms_rejected = 0
        total_features = 0
        features_with_msms = 0
        
        # Add nodes
        node_id = 0
        for dataset_id, (filename, features) in enumerate(all_list_features):
            logger.info(f"Adding nodes for dataset {dataset_id}: {os.path.basename(filename)}")
            
            for feature_id, feature in enumerate(features):
                # Create a unique ID for the node
                node_id_str = f"{dataset_id}_{feature_id}"
                
                # Add node with attributes
                self.G.add_node(node_id_str, 
                               dataset_id=dataset_id,
                               feature_id=feature_id,
                               mz=feature.get('mz', 0),
                               rt=feature.get('rt', 0),
                               intensity=feature.get('intensity', 0),
                               has_msms=feature.get('has_msms', False),
                               filename=os.path.basename(filename))
                
                total_features += 1
                if feature.get('has_msms', False):
                    features_with_msms += 1
                node_id += 1
        
        logger.info(f"Added {self.G.number_of_nodes()} nodes to the graph")
        logger.info(f"Features with MS/MS data: {features_with_msms}/{total_features} ({100*features_with_msms/total_features:.1f}%)")
        
        # Add edges using two-case logic
        logger.info("Adding edges with two-case matching logic...")
        edge_count = 0
        
        # Compare features across different datasets
        for i, (filename_i, features_i) in enumerate(all_list_features):
            for j, (filename_j, features_j) in enumerate(all_list_features):
                if i >= j:  # Skip same dataset and avoid duplicate comparisons
                    continue
                
                logger.info(f"Comparing dataset {i} and {j}...")
                
                # Compare features between datasets
                for feature_id_i, feature_i in enumerate(features_i):
                    for feature_id_j, feature_j in enumerate(features_j):
                        # Step 1: Apply m/z and RT gate (same for both cases)
                        mz_diff = abs(feature_i.get('mz', 0) - feature_j.get('mz', 0))
                        rt_diff = abs(feature_i.get('rt', 0) - feature_j.get('rt', 0))
                        
                        # Only proceed if within tolerance
                        if mz_diff <= self.mz_tolerance and rt_diff <= self.rt_tolerance:
                            
                            # Step 2: Determine which case applies
                            both_have_msms = (has_msms_data(feature_i) and has_msms_data(feature_j))
                            
                            if both_have_msms:
                                # Case 2: MS/MS available - use cosine similarity as weight
                                weight = calculate_spectral_similarity(
                                    feature_i, feature_j, 
                                    min_shared_peaks=self.min_shared_peaks,
                                    cosine_threshold=self.cosine_threshold
                                )
                                
                                if weight > 0:  # Only add edge if above threshold
                                    # Calculate detailed spectral information for storage
                                    from spectral_similarity import fast_cosine_similarity
                                    peaks1 = feature_i['msms_peaks']
                                    intensities1 = feature_i['msms_intensities']
                                    peaks2 = feature_j['msms_peaks']
                                    intensities2 = feature_j['msms_intensities']
                                    cosine_score, shared_peaks_count = fast_cosine_similarity(
                                        peaks1, intensities1, peaks2, intensities2, self.min_shared_peaks
                                    )
                                    
                                    node_i = f"{i}_{feature_id_i}"
                                    node_j = f"{j}_{feature_id_j}"
                                    self.G.add_edge(node_i, node_j, 
                                                  weight=weight, 
                                                  edge_type='msms',
                                                  cosine_similarity=cosine_score,
                                                  shared_peaks=shared_peaks_count)
                                    edge_count += 1
                                    msms_edges += 1
                                else:
                                    msms_rejected += 1
                                    
                            else:
                                # Case 1: No MS/MS - use m/z/RT weight (current behavior)
                                weight = 1.0 - (mz_diff / self.mz_tolerance + rt_diff / self.rt_tolerance) / 2.0
                                node_i = f"{i}_{feature_id_i}"
                                node_j = f"{j}_{feature_id_j}"
                                self.G.add_edge(node_i, node_j, weight=weight, edge_type='mz_rt')
                                edge_count += 1
                                mz_rt_edges += 1
        
        # Log comprehensive statistics
        logger.info(f"Edge creation completed:")
        logger.info(f"  Total edges added: {edge_count}")
        logger.info(f"  Case 1 (m/z/RT) edges: {mz_rt_edges} ({100*mz_rt_edges/edge_count:.1f}%)")
        logger.info(f"  Case 2 (MS/MS) edges: {msms_edges} ({100*msms_edges/edge_count:.1f}%)")
        logger.info(f"  MS/MS edges rejected (cosine < {self.cosine_threshold}): {msms_rejected}")
        
        # Remove isolated nodes
        isolated_nodes = list(nx.isolates(self.G))
        self.G.remove_nodes_from(isolated_nodes)
        logger.info(f"Removed {len(isolated_nodes)} isolated nodes")
        logger.info(f"Final graph has {self.G.number_of_nodes()} nodes and {self.G.number_of_edges()} edges")
        
        return self.G

    def get_feature_data(self, node_id):
        """
        Get feature data for a node.
        
        Parameters:
        -----------
        node_id : str
            Node ID in the graph
            
        Returns:
        --------
        feature_data : dict
            Dictionary with feature data
        """
        if node_id not in self.G:
            return None
        
        return self.G.nodes[node_id]

    def calculate_similarity_score(self, feature1: Dict, feature2: Dict) -> float:
        """
        Calculate similarity score between two features
        Returns score between 0 and 1
        """
        mz_diff = abs(feature1.get('mz', 0) - feature2.get('mz', 0))
        rt_diff = abs(feature1.get('rt', 0) - feature2.get('rt', 0))
        
        # Normalize differences by tolerances
        mz_score = 1 - (mz_diff / self.mz_tolerance)
        rt_score = 1 - (rt_diff / self.rt_tolerance)
        
        # Weight m/z more heavily than RT (0.7 vs 0.3)
        return 0.7 * mz_score + 0.3 * rt_score if mz_score > 0 and rt_score > 0 else 0

    def clean_multiple_connections(self) -> nx.Graph:
        """
        Pre-process graph to resolve multiple connections with edge type priority.
        
        MS/MS edges are prioritized over m/z/RT edges, then by weight within each type.
        This ensures structurally similar features (MS/MS) are preferred over 
        proximity-based matches (m/z/RT).
        
        Returns:
        --------
        nx.Graph: Cleaned graph with resolved multiple connections
        """
        logger.info("Cleaning multiple connections with edge type prioritization...")
        
        # Statistics tracking
        total_multiple_connections = 0
        msms_preferred = 0
        mz_rt_kept = 0
        edges_removed = 0
        
        for node in list(self.G.nodes()):
            neighbors = list(self.G.neighbors(node))
            by_dataset = {}
            
            # Group neighbors by dataset
            for neighbor in neighbors:
                dataset_id = neighbor.split('_')[0]  # Assuming ID format: "dataset_featurenum"
                by_dataset.setdefault(dataset_id, []).append(neighbor)
            
            # Resolve multiple connections within each dataset
            for dataset, features in by_dataset.items():
                if len(features) > 1:
                    total_multiple_connections += 1
                    
                    # Get edge information: (feature, weight, edge_type)
                    edges_info = []
                    for feature in features:
                        edge_data = self.G[node][feature]
                        weight = edge_data.get('weight', 0.0)
                        edge_type = edge_data.get('edge_type', 'mz_rt')
                        edges_info.append((feature, weight, edge_type))
                    
                    # Prioritize by edge type, then by weight
                    msms_edges = [(f, w) for f, w, t in edges_info if t == 'msms']
                    mz_rt_edges = [(f, w) for f, w, t in edges_info if t == 'mz_rt']
                    
                    # Select best connection
                    if msms_edges:
                        # Prefer MS/MS edges - select highest weight MS/MS edge
                        best_feature = max(msms_edges, key=lambda x: x[1])[0]
                        msms_preferred += 1
                        logger.debug(f"Node {node}: MS/MS edge preferred to {best_feature} "
                                   f"over {len(features)-1} alternatives")
                    else:
                        # No MS/MS edges available - select highest weight m/z/RT edge
                        best_feature = max(mz_rt_edges, key=lambda x: x[1])[0]
                        mz_rt_kept += 1
                        logger.debug(f"Node {node}: Best m/z/RT edge kept to {best_feature}")
                    
                    # Remove other edges
                    for feature in features:
                        if feature != best_feature:
                            self.G.remove_edge(node, feature)
                            edges_removed += 1
        
        # Log statistics
        logger.info(f"Multiple connection resolution completed:")
        logger.info(f"  Cases with multiple connections: {total_multiple_connections}")
        logger.info(f"  MS/MS edges preferred: {msms_preferred}")
        logger.info(f"  m/z/RT edges kept: {mz_rt_kept}")
        logger.info(f"  Total edges removed: {edges_removed}")
        
        return self.G

    def get_graph_stats(self) -> Dict:
        """
        Return basic statistics about the graph
        """
        return {
            'num_nodes': self.G.number_of_nodes(),
            'num_edges': self.G.number_of_edges(),
            'avg_degree': sum(dict(self.G.degree()).values()) / self.G.number_of_nodes(),
            'connected_components': nx.number_connected_components(self.G)
        }
