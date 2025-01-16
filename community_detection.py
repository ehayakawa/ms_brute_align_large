import pandas as pd
import numpy as np
import community.community_louvain as community_louvain

def detect_communities(G):
    # Initial community detection using Louvain method
    partition = community_louvain.best_partition(G, weight='weight')
    
    # Post-process and validate communities
    validated_partition = validate_communities(G, partition)
    
    return validated_partition

def validate_communities(G, partition):
    """
    Validate and clean up communities based on defined criteria
    """
    # Group nodes by community
    communities = {}
    for node, comm_id in partition.items():
        if comm_id not in communities:
            communities[comm_id] = []
        communities[comm_id].append(node)
    
    # Validate each community
    valid_communities = {}
    next_valid_id = 0
    
    for comm_id, nodes in communities.items():
        if is_valid_community(G, nodes):
            # Resolve any remaining multiple connections
            cleaned_nodes = resolve_multiple_connections(G, nodes)
            if cleaned_nodes:
                valid_communities[next_valid_id] = cleaned_nodes
                next_valid_id += 1
    
    # Create new partition with only valid communities
    new_partition = {}
    for comm_id, nodes in valid_communities.items():
        for node in nodes:
            new_partition[node] = comm_id
            
    return new_partition

def is_valid_community(G, nodes):
    """
    Check if community meets validation criteria
    """
    if len(nodes) < 3:
        return False
    
    mz_values = []
    rt_values = []
    datasets = set()
    
    for node in nodes:
        # Parse the node ID string (format: "file_idx_feat_idx")
        file_index = node.split('_')[0]
        datasets.add(file_index)
        
        # Get feature data from graph node attributes
        feature_data = G.nodes[node]
        if 'mz' in feature_data:
            mz_values.append(feature_data['mz'])
        if 'rt' in feature_data:
            rt_values.append(feature_data['rt'])
    
    # Check one feature per dataset
    if len(nodes) != len(datasets):
        return False
    
    # Check m/z consistency
    if mz_values:
        mz_variance = np.var(mz_values)
        if mz_variance > 0.02:  # configurable threshold
            return False
    
    # Check RT distribution
    if rt_values:
        rt_range = max(rt_values) - min(rt_values)
        if rt_range > 1.0:  # configurable threshold
            return False
    
    return True

def resolve_multiple_connections(G, nodes):
    """
    Resolve cases where multiple features from the same dataset exist in a community
    """
    dataset_features = {}
    for node in nodes:
        file_index = node.split('_')[0]
        if file_index not in dataset_features:
            dataset_features[file_index] = []
        dataset_features[file_index].append(node)
    
    resolved_nodes = []
    for dataset, features in dataset_features.items():
        if len(features) == 1:
            resolved_nodes.extend(features)
        else:
            best_feature = select_best_feature(G, features, nodes)
            if best_feature:
                resolved_nodes.append(best_feature)
    
    return resolved_nodes if resolved_nodes else None

def select_best_feature(G, candidates, community_nodes):
    """
    Select the best feature from multiple candidates using community context
    """
    # Calculate community medians (excluding candidates)
    other_nodes = [n for n in community_nodes if n not in candidates]
    if not other_nodes:
        return candidates[0]
    
    community_mz = np.median([G.nodes[n]['precursor_mz'] for n in other_nodes])
    
    # Select feature closest to community median m/z
    best_feature = None
    min_diff = float('inf')
    
    for node in candidates:
        mz_diff = abs(G.nodes[node]['precursor_mz'] - community_mz)
        if mz_diff < min_diff:
            min_diff = mz_diff
            best_feature = node
    
    return best_feature

def generate_community_tables(G, partition, all_list_features):
    """
    Generate tables of aligned features from community detection results
    """
    aligned_features = {}
    aligned_mz = {}
    aligned_intensity = {}

    for node, community in partition.items():
        # Parse node ID string (format: "file_idx_feat_idx")
        file_index, feature_index = map(int, node.split('_'))
        feature = all_list_features[file_index][1][feature_index]
        
        if community not in aligned_features:
            aligned_features[community] = {}
            aligned_mz[community] = []
            aligned_intensity[community] = []
        
        aligned_features[community][file_index] = feature_index
        aligned_mz[community].append(feature['precursor_mz'])
        
        if 'signal_intensity' in feature:
            aligned_intensity[community].append(feature['signal_intensity'])

    return aligned_features, aligned_mz, aligned_intensity
