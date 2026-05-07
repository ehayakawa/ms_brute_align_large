"""
Module for generating aligned feature output files.

This module takes the grouped features from community or clique detection
and creates structured TSV output files. It handles feature merging,
average m/z calculation, and formatting of the final alignment results.

Main functions/classes:
    - write_aligned_features_tsv: Writes aligned features to TSV format
    - filter_aligned_features: Filters features based on dataset requirements
    - calculate_average_mz: Computes representative m/z for feature groups
    - merge_similar_groups: Combines groups with overlapping features

Inputs:
    - Aligned feature dictionaries from community/clique detection
    - Original feature data with intensities and metadata
    - Output file paths and filtering parameters

Outputs:
    - TSV files with aligned features across datasets
    - Average m/z values, intensities per dataset, and feature IDs
    - Summary statistics for alignment quality

Important arguments:
    - aligned_features: Dictionary of grouped features
    - min_datasets: Minimum datasets required for valid alignment
    - output_file: Path to output TSV file
"""
import os
import pandas as pd
import numpy as np
from collections import defaultdict
from typing import Dict, List, Tuple, Any

def get_msms_matching_info(features_in_group, graph):
    """
    Extract MS/MS matching information for features in a group.
    
    Parameters:
    -----------
    features_in_group : list or dict
        Features in the current group
    graph : networkx.Graph
        Graph containing edge information
        
    Returns:
    --------
    msms_matches : list
        List of tuples (node1, node2, cosine_similarity, shared_peaks)
    """
    msms_matches = []
    
    if graph is None:
        return msms_matches
    
    # Get node IDs for features in this group
    if isinstance(features_in_group, list):
        # Community detection format
        node_ids = [f"{feature['dataset_id']}_{feature['feature_id']}" for feature in features_in_group]
    else:
        # Clique detection format
        node_ids = [f"{dataset_id}_{feature_id}" for dataset_id, feature_id in features_in_group.items()]
    
    # Check all pairs of features in the group for MS/MS edges
    for i, node1 in enumerate(node_ids):
        for j, node2 in enumerate(node_ids):
            if i < j and graph.has_edge(node1, node2):
                edge_data = graph[node1][node2]
                if edge_data.get('edge_type') == 'msms':
                    cosine_sim = edge_data.get('cosine_similarity', edge_data.get('weight', 0))
                    shared_peaks = edge_data.get('shared_peaks', 0)
                    msms_matches.append((node1, node2, cosine_sim, shared_peaks))
    
    return msms_matches

def write_aligned_features_tsv(aligned_features, feature_mzs, all_list_features, output_file, graph=None):
    """
    Write aligned features to a TSV file with MS/MS matching information.
    
    Parameters:
    -----------
    aligned_features : dict
        Dictionary mapping group IDs to either:
        - lists of feature dictionaries (from community detection)
        - dictionaries mapping dataset_id to feature_id (from clique detection)
    feature_mzs : dict
        Dictionary mapping (group_id, dataset_id, feature_id) to m/z values
    all_list_features : list
        List of tuples (filename, features) where features is a list of dictionaries
    output_file : str
        Path to the output TSV file
    graph : networkx.Graph, optional
        Graph containing edge information for MS/MS similarity data
    """
    print(f"Writing aligned features to {output_file}...")
    
    # Create a list to store rows for the TSV file
    rows = []
    
    # Get filenames for header
    filenames = [os.path.basename(filename) for filename, _ in all_list_features]
    
    # Create header row
    header = ["Group ID"]
    
    # Add columns for each file (feature index, m/z, intensity)
    for filename in filenames:
        header.extend([
            f"{filename}_feature_index",
            f"{filename}_mz",
            f"{filename}_intensity"
        ])
    
    # Add MS/MS matching information columns
    header.extend(["MSMS_Matches", "MSMS_Details"])
    
    # Add rows for each aligned group
    for group_id, features in aligned_features.items():
        row = [f"Group_{group_id}"]
        
        # Check if features is a list (community detection) or a dictionary (clique detection)
        if isinstance(features, list):
            # Community detection format: list of dictionaries
            # Group features by dataset
            features_by_dataset = defaultdict(list)
            for feature in features:
                features_by_dataset[feature['dataset_id']].append(feature)
            
            # Add data for each file
            for dataset_id, (filename, _) in enumerate(all_list_features):
                if dataset_id in features_by_dataset:
                    # Use the first feature for this dataset (could be multiple if they're in the same community/clique)
                    feature = features_by_dataset[dataset_id][0]
                    feature_id = feature['feature_id']
                    
                    # Get original feature data
                    _, dataset_features = all_list_features[dataset_id]
                    original_feature = dataset_features[feature_id]
                    
                    # Add feature index, m/z, and intensity
                    row.extend([
                        feature_id,
                        original_feature.get('mz', ''),
                        original_feature.get('intensity', '')
                    ])
                else:
                    # No feature for this dataset in this group
                    row.extend(['', '', ''])
        else:
            # Clique detection format: dictionary mapping dataset_id to feature_id
            # Add data for each file
            for dataset_id, (filename, _) in enumerate(all_list_features):
                if dataset_id in features:
                    feature_id = features[dataset_id]
                    
                    # Get original feature data
                    _, dataset_features = all_list_features[dataset_id]
                    original_feature = dataset_features[feature_id]
                    
                    # Add feature index, m/z, and intensity
                    row.extend([
                        feature_id,
                        original_feature.get('mz', ''),
                        original_feature.get('intensity', '')
                    ])
                else:
                    # No feature for this dataset in this group
                    row.extend(['', '', ''])
        
        # Add MS/MS matching information
        msms_matches = get_msms_matching_info(features, graph)
        
        if msms_matches:
            # Format MS/MS matches summary
            msms_count = len(msms_matches)
            
            # Create detailed match information
            match_details = []
            for node1, node2, cosine_sim, shared_peaks in msms_matches:
                # Extract dataset and feature info from node IDs
                dataset1, feature1 = node1.split('_', 1)
                dataset2, feature2 = node2.split('_', 1)
                
                # Get dataset names
                dataset1_name = os.path.basename(all_list_features[int(dataset1)][0])
                dataset2_name = os.path.basename(all_list_features[int(dataset2)][0])
                
                match_detail = f"{dataset1_name}({feature1})-{dataset2_name}({feature2}):cos={cosine_sim:.3f},peaks={shared_peaks}"
                match_details.append(match_detail)
            
            row.extend([msms_count, "; ".join(match_details)])
        else:
            row.extend([0, "No MS/MS matches"])
        
        rows.append(row)
    
    # Create DataFrame and write to TSV
    df = pd.DataFrame(rows, columns=header)
    df.to_csv(output_file, sep='\t', index=False)
    
    print(f"Wrote {len(rows)} aligned feature groups to {output_file}")

def calculate_average_mz(aligned_features, feature_mzs):
    """
    Calculate the average m/z value for each aligned feature group.
    
    Parameters:
    -----------
    aligned_features : dict
        Dictionary mapping group IDs to either:
        - lists of feature dictionaries (from community detection)
        - dictionaries mapping dataset_id to feature_id (from clique detection)
    feature_mzs : dict
        Dictionary mapping (group_id, dataset_id, feature_id) to m/z values
        
    Returns:
    --------
    avg_mzs : dict
        Dictionary mapping group IDs to average m/z values
    """
    avg_mzs = {}
    
    for group_id, features in aligned_features.items():
        mz_values = []
        
        # Check if features is a list (community detection) or a dictionary (clique detection)
        if isinstance(features, list):
            # Community detection format: list of dictionaries
            for feature in features:
                dataset_id = feature['dataset_id']
                feature_id = feature['feature_id']
                
                # Get m/z value
                mz = feature_mzs.get((group_id, dataset_id, feature_id), 0)
                if mz > 0:
                    mz_values.append(mz)
        else:
            # Clique detection format: dictionary mapping dataset_id to feature_id
            for dataset_id, feature_id in features.items():
                # Get m/z value
                mz = feature_mzs.get((group_id, dataset_id, feature_id), 0)
                if mz > 0:
                    mz_values.append(mz)
        
        # Calculate average m/z
        if mz_values:
            avg_mzs[group_id] = sum(mz_values) / len(mz_values)
        else:
            avg_mzs[group_id] = 0
    
    return avg_mzs

def filter_aligned_features(aligned_features, min_datasets=2):
    """
    Filter aligned features to only include groups with features from at least min_datasets datasets.
    
    Parameters:
    -----------
    aligned_features : dict
        Dictionary mapping group IDs to lists of features
    min_datasets : int
        Minimum number of datasets that must be represented in each group
        
    Returns:
    --------
    filtered_features : dict
        Filtered dictionary of aligned features
    """
    filtered_features = {}
    
    for group_id, features in aligned_features.items():
        # Count unique datasets
        datasets = set()
        for feature in features:
            # Check if feature is a dictionary or a tuple/list
            if isinstance(feature, dict):
                if 'dataset_id' in feature:
                    datasets.add(feature['dataset_id'])
            elif isinstance(feature, (tuple, list)) and len(feature) > 0:
                # Assume first element is dataset_id
                datasets.add(feature[0])
            elif hasattr(feature, 'dataset_id'):
                # Handle object with dataset_id attribute
                datasets.add(feature.dataset_id)
            elif isinstance(feature, int):
                # If feature is just an integer, assume it's the dataset_id itself
                datasets.add(feature)
        
        # Only keep groups with features from at least min_datasets datasets
        if len(datasets) >= min_datasets:
            filtered_features[group_id] = features
    
    print(f"Filtered from {len(aligned_features)} to {len(filtered_features)} groups with at least {min_datasets} datasets")
    
    return filtered_features

def merge_similar_groups(aligned_features, feature_mzs, mz_tolerance=0.01):
    """
    Merge similar groups based on m/z similarity.
    
    Parameters:
    -----------
    aligned_features : dict
        Dictionary mapping group IDs to lists of features
    feature_mzs : dict
        Dictionary mapping (group_id, dataset_id, feature_id) to m/z values
    mz_tolerance : float
        Tolerance for m/z values to consider groups similar
        
    Returns:
    --------
    merged_features : dict
        Dictionary with merged aligned features
    """
    # Calculate average m/z for each group
    avg_mzs = calculate_average_mz(aligned_features, feature_mzs)
    
    # Sort groups by average m/z
    sorted_groups = sorted(avg_mzs.items(), key=lambda x: x[1])
    
    # Initialize merged groups
    merged_features = {}
    merged_group_id = 0
    
    # Process groups in order of increasing m/z
    i = 0
    while i < len(sorted_groups):
        group_id, avg_mz = sorted_groups[i]
        
        # Start a new merged group
        current_group = aligned_features[group_id]
        current_datasets = set(feature['dataset_id'] for feature in current_group)
        
        # Look ahead for similar groups
        j = i + 1
        while j < len(sorted_groups):
            next_group_id, next_avg_mz = sorted_groups[j]
            
            # Check if m/z difference is within tolerance
            if next_avg_mz - avg_mz > mz_tolerance:
                break
            
            # Check if groups have overlapping datasets
            next_group = aligned_features[next_group_id]
            next_datasets = set(feature['dataset_id'] for feature in next_group)
            
            if not current_datasets.intersection(next_datasets):
                # No overlap, can merge
                current_group.extend(next_group)
                current_datasets.update(next_datasets)
                
                # Remove this group from consideration
                sorted_groups.pop(j)
            else:
                j += 1
        
        # Add merged group
        merged_features[merged_group_id] = current_group
        merged_group_id += 1
        
        i += 1
    
    print(f"Merged {len(aligned_features)} groups into {len(merged_features)} groups")
    
    return merged_features 