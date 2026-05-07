"""
Main module for mass spectrometry feature alignment using graph-based methods.

This module orchestrates the entire alignment pipeline for LC-MS data, identifying
the same chemical compounds across multiple datasets despite retention time drift
and m/z variations.

Main functions/classes:
    - CommunityDetector: Detects communities using Louvain algorithm for soft grouping
    - CliqueDetector: Finds maximal cliques for strict feature grouping
    - parse_arguments: Handles command-line arguments
    - main: Orchestrates the entire alignment pipeline

Inputs:
    - Excel files containing mass features (via --input-dir)
    - Optional parameters: mz_tolerance, rt_tolerance, min_datasets

Outputs:
    - TSV files with aligned features (community and clique methods)
    - PNG visualizations (graphs and heatmaps) if --visualize is specified
    - Summary statistics

Important arguments:
    --input-dir: Directory containing Excel files with mass features
    --mz-tolerance: m/z tolerance in Da (default: 0.01)
    --rt-tolerance: RT tolerance in minutes (default: 0.5)
    --min-datasets: Minimum datasets for valid group (default: 2)
    --visualize: Generate visualization plots
"""
import os
import argparse
import csv
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from collections import defaultdict
import time
import numpy as np
import logging
from pathlib import Path
from read_files import read_features, read_excel, collect_files
from graph_construction import GraphBuilder
from community_detection import detect_communities, group_features_by_community, detect_cliques, group_features_by_clique
from clique_detection import find_cliques, generate_clique_tables
from mass_feature_aligner import write_aligned_features_tsv, filter_aligned_features, calculate_average_mz, merge_similar_groups
from visualize_graph import plot_initial_graph, plot_community_graph, plot_clique_graph, visualize_subgraph, create_intensity_heatmap

class CommunityDetector:
    """
    Class for detecting communities in a graph and grouping features by community.
    """
    def __init__(self):
        self.G = None
        
    def detect_communities(self, G, hard_separation=False):
        """
        Detect communities in the graph using the Louvain method.
        
        Parameters:
        -----------
        G : networkx.Graph
            Graph to detect communities in
        hard_separation : bool
            If True, use a higher resolution and post-process communities to ensure hard separation
        """
        from community_detection import detect_communities
        self.G = G  # Store the graph
        print("Detecting communities...")
        partition = detect_communities(G, hard_separation=hard_separation)
        
        # Count communities
        communities = {}
        for node, comm_id in partition.items():
            if comm_id not in communities:
                communities[comm_id] = []
            communities[comm_id].append(node)
        
        print(f"Found {len(communities)} communities")
        
        # Print top 5 communities by size
        top_communities = sorted(communities.items(), key=lambda x: len(x[1]), reverse=True)[:5]
        for comm_id, nodes in top_communities:
            print(f"Community {comm_id}: {len(nodes)} nodes")
        
        return partition
    
    def group_features_by_community(self, partition):
        """
        Group features by community.
        """
        from community_detection import group_features_by_community
        return group_features_by_community(self.G, partition)
    
    def get_top_communities(self, partition, top_n=10):
        """
        Get the top N communities by size.
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

class CliqueDetector:
    """
    Class for detecting cliques in a graph and grouping features by clique.
    """
    def __init__(self):
        self.G = None
        
    def find_cliques(self, G):
        """
        Find cliques in the graph with optimizations to prevent excessive computation.
        """
        self.G = G  # Store the graph
        print("Finding cliques with optimizations...")
        
        # Set a maximum size for cliques to prevent excessive computation
        max_clique_size = 10
        
        # Use a more efficient approach for finding cliques
        cliques = []
        
        # First, find connected components to process them separately
        components = list(nx.connected_components(G))
        print(f"Graph has {len(components)} connected components")
        
        # Process each component
        for i, component in enumerate(components):
            if i < 5:  # Only print details for the first 5 components
                print(f"Processing component {i} with {len(component)} nodes")
            
            # Skip very large components
            if len(component) > 100:
                print(f"Skipping component {i} with {len(component)} nodes (too large)")
                continue
                
            # Extract the subgraph for this component
            subgraph = G.subgraph(component)
            
            try:
                # Find cliques in this component with a size limit
                component_cliques = []
                for clique in nx.find_cliques(subgraph):
                    if len(clique) >= 3 and len(clique) <= max_clique_size:
                        component_cliques.append(clique)
                    
                    # Limit the number of cliques per component
                    if len(component_cliques) >= 1000:
                        print(f"Reached limit of 1000 cliques for component {i}")
                        break
                
                cliques.extend(component_cliques)
                
                # Limit the total number of cliques
                if len(cliques) >= 5000:
                    print("Reached limit of 5000 total cliques")
                    break
                    
            except Exception as e:
                print(f"Error finding cliques in component {i}: {e}")
        
        print(f"Found {len(cliques)} cliques with at least 3 nodes")
        
        # Print top 5 cliques by size
        top_cliques = sorted(cliques, key=len, reverse=True)[:5]
        for i, clique in enumerate(top_cliques):
            print(f"Clique {i}: {len(clique)} nodes")
        
        return cliques, self.G
    
    def group_features_by_clique(self, cliques):
        """
        Group features by clique.
        """
        from clique_detection import group_features_by_clique
        return group_features_by_clique(cliques, self.G)
    
    def get_top_cliques(self, cliques, top_n=10):
        """
        Get the top N cliques by size.
        """
        # Sort cliques by size and return the top N
        top_cliques = sorted(cliques, key=len, reverse=True)[:top_n]
        return top_cliques

def parse_arguments():
    """
    Parse command line arguments for input and output directories
    """
    parser = argparse.ArgumentParser(description='Mass Feature Alignment')
    parser.add_argument('--input-dir', '-i', default='input_excel',
                        help='Directory containing input Excel files (default: input_excel)')
    parser.add_argument('--output-dir', '-o', default='output',
                        help='Directory for output files (default: output)')
    parser.add_argument('--mz-tolerance', type=float, default=0.01,
                        help='m/z tolerance for feature matching (default: 0.01)')
    parser.add_argument('--rt-tolerance', type=float, default=1.0,
                        help='Retention time tolerance in minutes (default: 1.0)')
    parser.add_argument('--visualize', '-v', action='store_true',
                        help='Enable graph visualization (may be resource-intensive for large datasets)')
    parser.add_argument('--max-vis-nodes', type=int, default=1000,
                        help='Maximum number of nodes to include in visualizations (default: 1000)')
    parser.add_argument('--max-vis-edges', type=int, default=5000,
                        help='Maximum number of edges to include in visualizations (default: 5000)')
    return parser.parse_args()

def test_read_excel(directory):
    """
    Test function to read Excel files from a directory
    """
    excel_files = collect_files(directory, file_extension=".xlsx")
    print(f"Found {len(excel_files)} Excel files")
    
    for excel_file in excel_files:
        try:
            list_features = read_excel(excel_file)
            print(f"Read {len(list_features)} features from {os.path.basename(excel_file)}")
        except Exception as e:
            print(f"Error reading {excel_file}: {e}")
    
    print("\n" + "="*50 + "\n")

def write_summary(all_list_features, summary_file):
    """
    Write summary statistics for processed files to a markdown file.
    
    Parameters:
    -----------
    all_list_features : list
        List of tuples (filename, features) where features is a list of dictionaries
    summary_file : str or Path
        Path to the markdown file to write
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Writing summary to {summary_file}...")
    
    # Calculate summary statistics
    summary_data = []
    total_features = 0
    
    for filename, features in all_list_features:
        # Calculate statistics
        num_features = len(features)
        total_features += num_features
        
        # Calculate average m/z, RT, and intensity
        mz_values = [feature.get('mz', 0) for feature in features if 'mz' in feature]
        rt_values = [feature.get('rt', 0) for feature in features if 'rt' in feature]
        intensity_values = [feature.get('intensity', 0) for feature in features if 'intensity' in feature]
        
        avg_mz = sum(mz_values) / len(mz_values) if mz_values else 0
        avg_rt = sum(rt_values) / len(rt_values) if rt_values else 0
        avg_intensity = sum(intensity_values) / len(intensity_values) if intensity_values else 0
        
        # Add to summary data
        summary_data.append({
            'filename': Path(filename).name,
            'num_features': num_features,
            'avg_mz': avg_mz,
            'avg_rt': avg_rt,
            'avg_intensity': avg_intensity
        })
    
    # Write summary to markdown file
    with open(summary_file, 'w') as f:
        f.write("# Mass Feature Alignment Summary\n\n")
        f.write(f"Processed {len(all_list_features)} files with a total of {total_features} features.\n\n")
        
        # Write table header
        f.write("| Filename | Features | Avg m/z | Avg RT (min) | Avg Intensity |\n")
        f.write("|----------|----------|---------|--------------|---------------|\n")
        
        # Write table rows
        for data in summary_data:
            f.write(f"| {data['filename']} | {data['num_features']} | {data['avg_mz']:.4f} | {data['avg_rt']:.2f} | {data['avg_intensity']:.2e} |\n")
    
    logger.info(f"Summary written to {summary_file}")
    logger.info(f"Processed {len(all_list_features)} files with a total of {total_features} features")

def main():
    """
    Main function for running the mass feature alignment process.
    """
    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler('mass_alignment.log')
        ]
    )
    logger = logging.getLogger(__name__)
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Mass Feature Alignment')
    parser.add_argument('--input-dir', type=str, required=True, help='Directory containing Excel files with mass features')
    parser.add_argument('--output-dir', type=str, default='output', help='Directory to save output files')
    parser.add_argument('--mz-tolerance', type=float, default=0.01, help='m/z tolerance for feature matching (in Da)')
    parser.add_argument('--rt-tolerance', type=float, default=0.5, help='RT tolerance for feature matching (in minutes)')
    parser.add_argument('--min-datasets', type=int, default=2, help='Minimum number of datasets for a valid feature group')
    parser.add_argument('--visualize', action='store_true', help='Generate visualizations')
    parser.add_argument('--hard-separation', action='store_true', help='Enable hard separation of communities for better visualization')
    parser.add_argument('--max-vis-nodes', type=int, default=1000, help='Maximum number of nodes to display in visualizations')
    parser.add_argument('--max-vis-edges', type=int, default=5000, help='Maximum number of edges to display in visualizations')
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Start timing
    start_time = time.time()
    
    # Step 1: Read Excel files and extract features
    logger.info(f"Reading Excel files from {args.input_dir}...")
    excel_files = collect_files(args.input_dir, file_extension=".xlsx")
    
    if not excel_files:
        logger.error(f"No Excel files found in {args.input_dir}")
        return
    
    logger.info(f"Found {len(excel_files)} Excel files")
    
    # Read features from each file
    all_list_features = []
    for excel_file in excel_files:
        try:
            list_features = read_excel(excel_file)
            all_list_features.append((excel_file, list_features))
            logger.info(f"Read {len(list_features)} features from {Path(excel_file).name}")
        except Exception as e:
            logger.error(f"Error reading {excel_file}: {e}")
    
    # Write summary
    summary_file = output_dir / "summary.md"
    write_summary(all_list_features, summary_file)
    
    # Step 2: Build graph from features
    graph_builder = GraphBuilder(
        mz_tolerance=args.mz_tolerance, 
        rt_tolerance=args.rt_tolerance,
        cosine_threshold=0.5,
        min_shared_peaks=3
    )
    G = graph_builder.build_graph(all_list_features)
    
    # Clean multiple connections to keep only the most likely edge between datasets
    logger.info("Cleaning multiple connections...")
    G = graph_builder.clean_multiple_connections()
    logger.info(f"Graph after cleaning: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    
    # Step 3: Detect communities
    community_detector = CommunityDetector()
    partition = community_detector.detect_communities(G, hard_separation=args.hard_separation)
    
    # Save graph and partition for later use
    import pickle
    with open(output_dir / "graph.pkl", 'wb') as f:
        pickle.dump(G, f)
    with open(output_dir / "partition.pkl", 'wb') as f:
        pickle.dump(partition, f)
    logger.info(f"Saved graph and partition to {args.output_dir}")
    
    # Step 4: Group features by community
    aligned_features_community = community_detector.group_features_by_community(partition)
    
    # Filter aligned features
    aligned_features_community = filter_aligned_features(aligned_features_community, min_datasets=args.min_datasets)
    
    # Step 5: Detect cliques
    clique_detector = CliqueDetector()
    cliques, G_cliques = clique_detector.find_cliques(G)
    
    # Step 6: Group features by clique
    aligned_features_clique = clique_detector.group_features_by_clique(cliques)
    
    # Filter aligned features
    aligned_features_clique = filter_aligned_features(aligned_features_clique, min_datasets=args.min_datasets)
    
    # Step 7: Write aligned features to TSV files
    # Calculate average m/z values for each group
    feature_mzs_community = calculate_average_mz(aligned_features_community, {})
    output_file_community = output_dir / "aligned_features_community.tsv"
    write_aligned_features_tsv(aligned_features_community, feature_mzs_community, all_list_features, output_file_community, G)
    
    feature_mzs_clique = calculate_average_mz(aligned_features_clique, {})
    output_file_clique = output_dir / "aligned_features_clique.tsv"
    write_aligned_features_tsv(aligned_features_clique, feature_mzs_clique, all_list_features, output_file_clique, G)
    
    # Step 8: Visualize results
    if args.visualize:
        logger.info("Generating visualizations...")
        
        # Plot initial graph
        pos = plot_initial_graph(G, args.output_dir)
        
        # Plot community graph
        plot_community_graph(G, partition, args.output_dir, pos, args.max_vis_nodes, args.max_vis_edges, hard_separation=args.hard_separation)
        
        # Plot clique graph
        plot_clique_graph(G, cliques, args.output_dir, pos, args.max_vis_nodes, args.max_vis_edges, hard_separation=args.hard_separation)
        
        # Create intensity heatmaps
        create_intensity_heatmap(output_file_community, args.output_dir, max_groups=50)
        create_intensity_heatmap(output_file_clique, args.output_dir, max_groups=50)
    
    # Print timing information
    elapsed_time = time.time() - start_time
    logger.info(f"Mass feature alignment completed in {elapsed_time:.2f} seconds")

if __name__ == "__main__":
    main()
