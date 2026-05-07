#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Module for generating detailed reports from community detection results.

This module analyzes the output from community detection and creates comprehensive
reports including statistics, quality metrics, and detailed feature information
for each detected community.

Main functions/classes:
    - generate_community_report: Creates detailed report from alignment results
    - calculate_community_statistics: Computes metrics for each community
    - assess_alignment_quality: Evaluates the quality of alignments
    - create_summary_tables: Generates summary statistics tables

Inputs:
    - TSV files from community detection alignment
    - Original feature data for detailed analysis
    - Report parameters and thresholds

Outputs:
    - Detailed HTML/text reports with community analysis
    - Quality metrics and statistics
    - Visualizations of community characteristics

Important arguments:
    - input_file: Path to community alignment TSV file
    - output_dir: Directory for report output
    - min_quality: Minimum quality threshold for reporting
"""

import os
import argparse
import pandas as pd
import networkx as nx
from collections import defaultdict

def load_graph_and_partition(output_dir):
    """
    Load the graph and partition from the output directory.
    
    Parameters:
    -----------
    output_dir : str
        Directory containing the output files
        
    Returns:
    --------
    G : networkx.Graph
        The graph
    partition : dict
        Dictionary mapping node IDs to community IDs
    """
    # Try to load the graph from a pickle file if it exists
    import pickle
    graph_file = os.path.join(output_dir, "graph.pkl")
    if os.path.exists(graph_file):
        print(f"Loading graph from {graph_file}")
        with open(graph_file, 'rb') as f:
            G = pickle.load(f)
    else:
        print(f"Graph file not found: {graph_file}")
        return None, None
    
    # Try to load the partition from a pickle file if it exists
    partition_file = os.path.join(output_dir, "partition.pkl")
    if os.path.exists(partition_file):
        print(f"Loading partition from {partition_file}")
        with open(partition_file, 'rb') as f:
            partition = pickle.load(f)
    else:
        # Try to load from the aligned features file
        aligned_file = os.path.join(output_dir, "aligned_features_community.tsv")
        if os.path.exists(aligned_file):
            print(f"Loading partition from {aligned_file}")
            df = pd.read_csv(aligned_file, sep='\t')
            partition = {}
            for _, row in df.iterrows():
                group_id = row['group_id']
                for col in df.columns:
                    if 'node_id' in col and not pd.isna(row[col]):
                        node_id = int(row[col])
                        partition[node_id] = group_id
        else:
            print(f"Partition file not found: {partition_file}")
            print(f"Aligned features file not found: {aligned_file}")
            return G, None
    
    return G, partition

def generate_community_report(G, partition, output_file):
    """
    Generate a detailed report of communities.
    
    Parameters:
    -----------
    G : networkx.Graph
        The graph
    partition : dict
        Dictionary mapping node IDs to community IDs
    output_file : str
        Path to save the report
    """
    if G is None or partition is None:
        print("Cannot generate report: Graph or partition is missing")
        return
    
    # Group nodes by community
    communities = defaultdict(list)
    for node, comm_id in partition.items():
        if node in G.nodes():
            communities[comm_id].append(node)
    
    # Sort communities by size
    sorted_communities = sorted(communities.items(), key=lambda x: len(x[1]), reverse=True)
    
    # Create report
    with open(output_file, 'w') as f:
        f.write("Community Report\n")
        f.write("===============\n\n")
        
        for comm_id, nodes in sorted_communities:
            f.write(f"Community {comm_id} ({len(nodes)} nodes)\n")
            f.write("-" * 50 + "\n")
            
            # Sort nodes by dataset_id and then by intensity (descending)
            sorted_nodes = sorted(nodes, key=lambda n: (G.nodes[n].get('dataset_id', 0), -G.nodes[n].get('intensity', 0)))
            
            for node in sorted_nodes:
                node_data = G.nodes[node]
                filename = node_data.get('filename', 'Unknown')
                if isinstance(filename, str) and os.path.basename(filename):
                    filename = os.path.basename(filename)
                
                rt = node_data.get('rt', node_data.get('retention_time', 'N/A'))
                mz = node_data.get('mz', node_data.get('precursor_mz', 'N/A'))
                intensity = node_data.get('intensity', 'N/A')
                dataset_id = node_data.get('dataset_id', 'N/A')
                
                f.write(f"  Node {node}: {filename}, RT={rt:.2f}, m/z={mz:.4f}, intensity={intensity:.1f}, dataset={dataset_id}\n")
            
            f.write("\n")
    
    print(f"Community report saved to {output_file}")

def generate_community_csv(G, partition, output_file):
    """
    Generate a CSV report of communities.
    
    Parameters:
    -----------
    G : networkx.Graph
        The graph
    partition : dict
        Dictionary mapping node IDs to community IDs
    output_file : str
        Path to save the CSV report
    """
    if G is None or partition is None:
        print("Cannot generate CSV: Graph or partition is missing")
        return
    
    # Group nodes by community
    communities = defaultdict(list)
    for node, comm_id in partition.items():
        if node in G.nodes():
            communities[comm_id].append(node)
    
    # Prepare data for CSV
    data = []
    for comm_id, nodes in communities.items():
        for node in nodes:
            node_data = G.nodes[node]
            filename = node_data.get('filename', 'Unknown')
            if isinstance(filename, str) and os.path.basename(filename):
                filename = os.path.basename(filename)
            
            rt = node_data.get('rt', node_data.get('retention_time', 'N/A'))
            mz = node_data.get('mz', node_data.get('precursor_mz', 'N/A'))
            intensity = node_data.get('intensity', 'N/A')
            dataset_id = node_data.get('dataset_id', 'N/A')
            
            data.append({
                'community_id': comm_id,
                'node_id': node,
                'filename': filename,
                'rt': rt,
                'mz': mz,
                'intensity': intensity,
                'dataset_id': dataset_id
            })
    
    # Create DataFrame and save to CSV
    df = pd.DataFrame(data)
    df.to_csv(output_file, index=False)
    print(f"Community CSV report saved to {output_file}")

def generate_simple_community_report(G, partition, output_file):
    """
    Generate a simplified community report in the format:
    community, node1 (filename, rt, mz), node2 (filename, rt, mz), ...
    
    Parameters:
    -----------
    G : networkx.Graph
        The graph
    partition : dict
        Dictionary mapping node IDs to community IDs
    output_file : str
        Path to save the report
    """
    if G is None or partition is None:
        print("Cannot generate simple report: Graph or partition is missing")
        return
    
    # Group nodes by community
    communities = defaultdict(list)
    for node, comm_id in partition.items():
        if node in G.nodes():
            communities[comm_id].append(node)
    
    # Sort communities by size
    sorted_communities = sorted(communities.items(), key=lambda x: len(x[1]), reverse=True)
    
    # Create report
    with open(output_file, 'w') as f:
        f.write("community, node1 (filename, rt, mz), node2 (filename, rt, mz), ...\n")
        
        for comm_id, nodes in sorted_communities:
            line = f"Community {comm_id}"
            
            # Sort nodes by dataset_id
            sorted_nodes = sorted(nodes, key=lambda n: G.nodes[n].get('dataset_id', 0))
            
            for node in sorted_nodes:
                node_data = G.nodes[node]
                filename = node_data.get('filename', 'Unknown')
                if isinstance(filename, str) and os.path.basename(filename):
                    filename = os.path.basename(filename)
                
                rt = node_data.get('rt', node_data.get('retention_time', 'N/A'))
                mz = node_data.get('mz', node_data.get('precursor_mz', 'N/A'))
                
                # Format rt and mz with appropriate precision
                if isinstance(rt, (int, float)):
                    rt_str = f"{rt:.2f}"
                else:
                    rt_str = str(rt)
                
                if isinstance(mz, (int, float)):
                    mz_str = f"{mz:.4f}"
                else:
                    mz_str = str(mz)
                
                line += f", Node {node} ({filename}, {rt_str}, {mz_str})"
            
            f.write(line + "\n")
    
    print(f"Simple community report saved to {output_file}")

def generate_markdown_community_report(G, partition, output_file):
    """
    Generate a markdown-formatted report of communities with tables.
    
    Parameters:
    -----------
    G : networkx.Graph
        The graph
    partition : dict
        Dictionary mapping node IDs to community IDs
    output_file : str
        Path to save the markdown report
    """
    if G is None or partition is None:
        print("Cannot generate markdown report: Graph or partition is missing")
        return
    
    # Group nodes by community
    communities = defaultdict(list)
    for node, comm_id in partition.items():
        if node in G.nodes():
            communities[comm_id].append(node)
    
    # Sort communities by size
    sorted_communities = sorted(communities.items(), key=lambda x: len(x[1]), reverse=True)
    
    # Create report
    with open(output_file, 'w') as f:
        f.write("# Community Report\n\n")
        f.write("## Summary\n\n")
        f.write("This report shows communities with features grouped by dataset. ")
        f.write("Communities with multiple features from the same dataset are highlighted.\n\n")
        f.write("Total communities: " + str(len(sorted_communities)) + "\n\n")
        
        # Count communities with multiple features from the same dataset
        multi_feature_communities = 0
        for comm_id, nodes in sorted_communities:
            # Group nodes by dataset
            nodes_by_dataset = defaultdict(list)
            for node in nodes:
                dataset_id = G.nodes[node].get('dataset_id', 'Unknown')
                nodes_by_dataset[dataset_id].append(node)
            
            # Check if any dataset has multiple features
            has_multiple = any(len(nodes) > 1 for nodes in nodes_by_dataset.values())
            if has_multiple:
                multi_feature_communities += 1
        
        f.write(f"Communities with multiple features from the same dataset: {multi_feature_communities}\n\n")
        f.write("---\n\n")
        
        for comm_id, nodes in sorted_communities:
            f.write(f"## Community {comm_id} ({len(nodes)} nodes)\n\n")
            
            # Group nodes by dataset
            nodes_by_dataset = defaultdict(list)
            for node in nodes:
                dataset_id = G.nodes[node].get('dataset_id', 'Unknown')
                nodes_by_dataset[dataset_id].append(node)
            
            # Check if any dataset has multiple features
            has_multiple = any(len(nodes) > 1 for nodes in nodes_by_dataset.values())
            
            # Add summary for this community
            f.write("### Summary\n\n")
            f.write(f"Datasets represented: {len(nodes_by_dataset)}\n\n")
            
            if has_multiple:
                f.write("⚠️ **This community contains multiple features from the same dataset(s):**\n\n")
                for dataset_id, dataset_nodes in nodes_by_dataset.items():
                    if len(dataset_nodes) > 1:
                        f.write(f"- Dataset {dataset_id}: {len(dataset_nodes)} features\n")
                f.write("\n")
            
            # Calculate average m/z and RT
            mz_values = []
            rt_values = []
            for node in nodes:
                node_data = G.nodes[node]
                mz = node_data.get('mz', node_data.get('precursor_mz', None))
                rt = node_data.get('rt', node_data.get('retention_time', None))
                if isinstance(mz, (int, float)):
                    mz_values.append(mz)
                if isinstance(rt, (int, float)):
                    rt_values.append(rt)
            
            if mz_values:
                avg_mz = sum(mz_values) / len(mz_values)
                min_mz = min(mz_values)
                max_mz = max(mz_values)
                f.write(f"Average m/z: {avg_mz:.4f} (range: {min_mz:.4f} - {max_mz:.4f})\n\n")
            
            if rt_values:
                avg_rt = sum(rt_values) / len(rt_values)
                min_rt = min(rt_values)
                max_rt = max(rt_values)
                f.write(f"Average RT: {avg_rt:.2f} (range: {min_rt:.2f} - {max_rt:.2f})\n\n")
            
            # Create table header
            f.write("### Features\n\n")
            f.write("| Dataset | Node ID | m/z | RT | Intensity | Filename | Multiple |\n")
            f.write("|---------|---------|-----|----|-----------|---------|---------|\n")
            
            # Sort datasets
            sorted_datasets = sorted(nodes_by_dataset.keys())
            
            for dataset_id in sorted_datasets:
                # Sort nodes by intensity (descending)
                sorted_nodes = sorted(nodes_by_dataset[dataset_id], 
                                     key=lambda n: -G.nodes[n].get('intensity', 0))
                
                # Check if this dataset has multiple features
                has_multiple_in_dataset = len(sorted_nodes) > 1
                
                for i, node in enumerate(sorted_nodes):
                    node_data = G.nodes[node]
                    filename = node_data.get('filename', 'Unknown')
                    if isinstance(filename, str) and os.path.basename(filename):
                        filename = os.path.basename(filename)
                    
                    rt = node_data.get('rt', node_data.get('retention_time', 'N/A'))
                    mz = node_data.get('mz', node_data.get('precursor_mz', 'N/A'))
                    intensity = node_data.get('intensity', 'N/A')
                    
                    # Format values with appropriate precision
                    if isinstance(rt, (int, float)):
                        rt_str = f"{rt:.2f}"
                    else:
                        rt_str = str(rt)
                    
                    if isinstance(mz, (int, float)):
                        mz_str = f"{mz:.4f}"
                    else:
                        mz_str = str(mz)
                    
                    if isinstance(intensity, (int, float)):
                        intensity_str = f"{intensity:.1f}"
                    else:
                        intensity_str = str(intensity)
                    
                    # Add a marker for multiple features
                    multiple_marker = "⚠️" if has_multiple_in_dataset else ""
                    
                    f.write(f"| {dataset_id} | {node} | {mz_str} | {rt_str} | {intensity_str} | {filename} | {multiple_marker} |\n")
            
            f.write("\n---\n\n")
    
    print(f"Markdown community report saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Generate community report')
    parser.add_argument('--output-dir', type=str, default='output', help='Directory containing output files')
    parser.add_argument('--report-file', type=str, default=None, help='Path to save the report (default: output_dir/community_report.txt)')
    parser.add_argument('--csv-file', type=str, default=None, help='Path to save the CSV report (default: output_dir/community_report.csv)')
    parser.add_argument('--simple-report-file', type=str, default=None, help='Path to save the simple report (default: output_dir/community_report_simple.txt)')
    parser.add_argument('--markdown-file', type=str, default=None, help='Path to save the markdown report (default: output_dir/community_report.md)')
    args = parser.parse_args()
    
    # Set default output files if not specified
    if args.report_file is None:
        args.report_file = os.path.join(args.output_dir, "community_report.txt")
    
    if args.csv_file is None:
        args.csv_file = os.path.join(args.output_dir, "community_report.csv")
    
    if args.simple_report_file is None:
        args.simple_report_file = os.path.join(args.output_dir, "community_report_simple.txt")
    
    if args.markdown_file is None:
        args.markdown_file = os.path.join(args.output_dir, "community_report.md")
    
    # Load graph and partition
    G, partition = load_graph_and_partition(args.output_dir)
    
    # Generate reports
    generate_community_report(G, partition, args.report_file)
    generate_community_csv(G, partition, args.csv_file)
    generate_simple_community_report(G, partition, args.simple_report_file)
    generate_markdown_community_report(G, partition, args.markdown_file)

if __name__ == "__main__":
    main() 