"""
Module for visualizing mass feature alignment graphs and results.

This module provides visualization functions for the feature similarity graphs,
community detection results, and clique analysis. It creates various plots
including network graphs, heatmaps, and summary visualizations.

Main functions/classes:
    - plot_initial_graph: Visualizes the initial feature similarity graph
    - plot_community_graph: Shows graph with community coloring
    - plot_clique_graph: Displays cliques with distinct visual representation
    - visualize_subgraph: Creates detailed views of specific graph regions
    - create_intensity_heatmap: Generates heatmaps of feature intensities

Inputs:
    - NetworkX graphs with feature nodes and similarity edges
    - Community/clique membership information
    - Feature intensity data across datasets
    - Visualization parameters (colors, layouts, sizes)

Outputs:
    - PNG files with graph visualizations
    - Heatmaps showing feature intensities across datasets
    - Summary plots with alignment statistics

Important arguments:
    - G: NetworkX graph to visualize
    - output_dir: Directory for saving visualization files
    - partition: Community or clique membership dictionary
    - layout: Graph layout algorithm (spring, circular, etc.)
"""
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os
import random
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from collections import defaultdict

# Global variable to store the initial layout
initial_layout = None

def get_spring_layout(G, k=None, iterations=50, seed=42):
    """
    Get a spring layout for the graph.
    
    Parameters:
    -----------
    G : networkx.Graph
        Graph to get layout for
    k : float
        Optimal distance between nodes
    iterations : int
        Number of iterations for the spring layout algorithm
    seed : int
        Random seed for reproducibility
        
    Returns:
    --------
    pos : dict
        Dictionary mapping node IDs to positions
    """
    if k is None:
        k = 1.0 / np.sqrt(G.number_of_nodes())
    return nx.spring_layout(G, k=k, iterations=iterations, seed=seed)

def get_community_layout(G, partition):
    """
    Get a layout that groups nodes by community.
    """
    # First get a layout based on the initial positions
    pos = get_spring_layout(G)
    
    # Group by community
    communities = {}
    for node, comm_id in partition.items():
        if comm_id not in communities:
            communities[comm_id] = []
        communities[comm_id].append(node)
    
    # Calculate community centers
    comm_centers = {}
    for comm_id, nodes in communities.items():
        if nodes:
            center_x = sum(pos[node][0] for node in nodes) / len(nodes)
            center_y = sum(pos[node][1] for node in nodes) / len(nodes)
            comm_centers[comm_id] = (center_x, center_y)
    
    # Adjust positions to cluster by community - make more compact
    pos_adjusted = {}
    for node, position in pos.items():
        if node in partition:
            comm_id = partition[node]
            if comm_id in comm_centers:
                # Move 80% toward community center (increased from 60%)
                center_x, center_y = comm_centers[comm_id]
                pos_adjusted[node] = (
                    position[0] * 0.2 + center_x * 0.8,  # More attraction to center
                    position[1] * 0.2 + center_y * 0.8   # More attraction to center
                )
            else:
                pos_adjusted[node] = position
        else:
            pos_adjusted[node] = position
    
    # Spread out the community centers to avoid overlap
    # Calculate the average distance between nodes in the original layout
    all_positions = list(pos.values())
    avg_distance = 0
    count = 0
    for i in range(len(all_positions)):
        for j in range(i+1, len(all_positions)):
            dist = np.sqrt((all_positions[i][0] - all_positions[j][0])**2 + 
                          (all_positions[i][1] - all_positions[j][1])**2)
            avg_distance += dist
            count += 1
    
    if count > 0:
        avg_distance /= count
        
        # Spread out community centers
        comm_centers_list = list(comm_centers.items())
        for i in range(len(comm_centers_list)):
            comm_id1, (x1, y1) = comm_centers_list[i]
            for j in range(i+1, len(comm_centers_list)):
                comm_id2, (x2, y2) = comm_centers_list[j]
                
                # Calculate distance between centers
                dist = np.sqrt((x1 - x2)**2 + (y1 - y2)**2)
                
                # If centers are too close, move them apart
                if dist < avg_distance * 3:  # Ensure communities are well separated
                    # Direction vector
                    dx = x2 - x1
                    dy = y2 - y1
                    
                    # Normalize
                    if dist > 0:
                        dx /= dist
                        dy /= dist
                    else:
                        # If centers are at the same position, move in random direction
                        angle = random.uniform(0, 2 * np.pi)
                        dx = np.cos(angle)
                        dy = np.sin(angle)
                    
                    # Move centers apart
                    move_dist = (avg_distance * 3 - dist) / 2
                    
                    # Update positions of all nodes in these communities
                    for node in communities[comm_id1]:
                        if node in pos_adjusted:
                            pos_adjusted[node] = (
                                pos_adjusted[node][0] - dx * move_dist,
                                pos_adjusted[node][1] - dy * move_dist
                            )
                    
                    for node in communities[comm_id2]:
                        if node in pos_adjusted:
                            pos_adjusted[node] = (
                                pos_adjusted[node][0] + dx * move_dist,
                                pos_adjusted[node][1] + dy * move_dist
                            )
    
    return pos_adjusted

def get_clique_layout(G, cliques):
    """
    Get a layout that groups nodes by clique.
    """
    # First get a layout based on the initial positions
    pos = get_spring_layout(G)
    
    # Create a mapping from node to clique
    node_to_clique = {}
    for i, clique in enumerate(cliques):
        for node in clique:
            if node not in node_to_clique:  # A node might be in multiple cliques, use first one
                node_to_clique[node] = i
    
    # Group by clique
    clique_groups = {}
    for i, clique in enumerate(cliques):
        clique_groups[i] = clique
    
    # Calculate clique centers
    clique_centers = {}
    for clique_id, nodes in clique_groups.items():
        if nodes:
            center_x = sum(pos[node][0] for node in nodes) / len(nodes)
            center_y = sum(pos[node][1] for node in nodes) / len(nodes)
            clique_centers[clique_id] = (center_x, center_y)
    
    # Adjust positions to cluster by clique - make more compact
    pos_adjusted = {}
    for node, position in pos.items():
        if node in node_to_clique:
            clique_id = node_to_clique[node]
            if clique_id in clique_centers:
                # Move 80% toward clique center (increased from 60%)
                center_x, center_y = clique_centers[clique_id]
                pos_adjusted[node] = (
                    position[0] * 0.2 + center_x * 0.8,  # More attraction to center
                    position[1] * 0.2 + center_y * 0.8   # More attraction to center
                )
            else:
                pos_adjusted[node] = position
        else:
            pos_adjusted[node] = position
    
    # Spread out the clique centers to avoid overlap
    # Calculate the average distance between nodes in the original layout
    all_positions = list(pos.values())
    avg_distance = 0
    count = 0
    for i in range(len(all_positions)):
        for j in range(i+1, len(all_positions)):
            dist = np.sqrt((all_positions[i][0] - all_positions[j][0])**2 + 
                          (all_positions[i][1] - all_positions[j][1])**2)
            avg_distance += dist
            count += 1
    
    if count > 0:
        avg_distance /= count
        
        # Spread out clique centers
        clique_centers_list = list(clique_centers.items())
        for i in range(len(clique_centers_list)):
            clique_id1, (x1, y1) = clique_centers_list[i]
            for j in range(i+1, len(clique_centers_list)):
                clique_id2, (x2, y2) = clique_centers_list[j]
                
                # Calculate distance between centers
                dist = np.sqrt((x1 - x2)**2 + (y1 - y2)**2)
                
                # If centers are too close, move them apart
                if dist < avg_distance * 3:  # Ensure cliques are well separated
                    # Direction vector
                    dx = x2 - x1
                    dy = y2 - y1
                    
                    # Normalize
                    if dist > 0:
                        dx /= dist
                        dy /= dist
                    else:
                        # If centers are at the same position, move in random direction
                        angle = random.uniform(0, 2 * np.pi)
                        dx = np.cos(angle)
                        dy = np.sin(angle)
                    
                    # Move centers apart
                    move_dist = (avg_distance * 3 - dist) / 2
                    
                    # Update positions of all nodes in these cliques
                    for node in clique_groups[clique_id1]:
                        if node in pos_adjusted:
                            pos_adjusted[node] = (
                                pos_adjusted[node][0] - dx * move_dist,
                                pos_adjusted[node][1] - dy * move_dist
                            )
                    
                    for node in clique_groups[clique_id2]:
                        if node in pos_adjusted:
                            pos_adjusted[node] = (
                                pos_adjusted[node][0] + dx * move_dist,
                                pos_adjusted[node][1] + dy * move_dist
                            )
    
    return pos_adjusted

def reset_layout():
    """
    Reset the global initial layout.
    This can be called if you want to start fresh with a new layout.
    """
    global initial_layout
    initial_layout = None

def plot_initial_graph(G, output_dir, pos=None, max_nodes=1000, max_edges=5000):
    """
    Plot the initial graph with nodes colored by dataset.
    
    Parameters:
    -----------
    G : networkx.Graph
        Graph to plot
    output_dir : str
        Directory to save the plot
    pos : dict, optional
        Dictionary of node positions
    max_nodes : int
        Maximum number of nodes to display
    max_edges : int
        Maximum number of edges to display
    """
    # Sample the graph if it's too large
    G_sampled = G
    if G.number_of_nodes() > max_nodes or G.number_of_edges() > max_edges:
        print(f"Graph is too large to visualize fully. Sampling to {max_nodes} nodes and {max_edges} edges.")
        G_sampled = sample_graph(G, max_nodes, max_edges)
    
    # Get layout
    if pos is None:
        pos = get_spring_layout(G_sampled)
    else:
        # Filter pos to only include nodes in G_sampled
        pos = {node: pos[node] for node in G_sampled.nodes() if node in pos}
    
    # Create figure
    plt.figure(figsize=(16, 14))
    
    # Draw edges with increased visibility
    nx.draw_networkx_edges(G_sampled, pos, width=1.0, alpha=0.7, edge_color='gray')
    
    # Draw nodes colored by dataset
    node_colors = []
    for node in G_sampled.nodes():
        dataset_id = G_sampled.nodes[node].get('dataset_id', 0)
        node_colors.append(dataset_id)
    
    nx.draw_networkx_nodes(G_sampled, pos, node_size=80, node_color=node_colors, cmap=plt.cm.rainbow, alpha=0.8)
    
    # Add labels for some nodes (e.g., high-degree nodes)
    if G_sampled.number_of_nodes() <= 100:  # Only label if not too many nodes
        nx.draw_networkx_labels(G_sampled, pos, font_size=8)
    else:
        # Label only high-degree nodes
        high_degree_nodes = sorted(G_sampled.degree, key=lambda x: x[1], reverse=True)[:20]
        labels = {node: str(node) for node, _ in high_degree_nodes}
        nx.draw_networkx_labels(G_sampled, pos, labels=labels, font_size=8)
    
    # Add a legend for dataset colors
    unique_datasets = sorted(set(node_colors))
    legend_elements = []
    for dataset_id in unique_datasets:
        color = plt.cm.rainbow(dataset_id / max(1, max(unique_datasets)))
        legend_elements.append(plt.Line2D([0], [0], marker='o', color='w', 
                               markerfacecolor=color, markersize=10, 
                               label=f'Dataset {dataset_id}'))
    
    plt.legend(handles=legend_elements, loc='upper right', title='Datasets')
    
    # Set title and remove axis
    plt.title(f"Initial Graph ({G_sampled.number_of_nodes()} nodes, {G_sampled.number_of_edges()} edges)")
    plt.axis('off')
    
    # Save figure
    output_path = os.path.join(output_dir, "initial_graph.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Initial graph visualization saved to {output_path}")

def plot_community_graph(G, partition, output_dir, pos=None, max_nodes=1000, max_edges=5000, hard_separation=False):
    """
    Plot the graph with community detection results.
    
    Steps:
    1. Sample the graph if it's too large
    2. Filter communities to ensure at most one node per dataset per community
    3. Create a new graph with ONLY intra-community edges (remove ALL inter-community edges)
    4. Apply a layout that properly separates communities
    5. Draw nodes colored by community
    6. Draw edges
    7. Draw transparent circles around communities
    
    Parameters:
    -----------
    G : networkx.Graph
        Graph to plot
    partition : dict
        Dictionary mapping node IDs to community IDs
    output_dir : str
        Directory to save the plot
    pos : dict
        Dictionary mapping node IDs to positions (ignored if hard_separation is True)
    max_nodes : int
        Maximum number of nodes to plot
    max_edges : int
        Maximum number of edges to plot
    hard_separation : bool
        If True, use a layout that strongly separates communities
    """
    if hard_separation:
        print("Plotting community graph with hard separation...")
    else:
        print("Plotting community graph...")
    
    # Filter partition to ensure at most one node per dataset per community
    filtered_partition = {}
    communities = {}
    
    # Group nodes by community
    for node, comm_id in partition.items():
        if comm_id not in communities:
            communities[comm_id] = []
        communities[comm_id].append(node)
    
    # For each community, keep only one node per dataset
    for comm_id, nodes in communities.items():
        # Group nodes by dataset
        by_dataset = {}
        for node in nodes:
            dataset_id = G.nodes[node].get('dataset_id')
            if dataset_id not in by_dataset:
                by_dataset[dataset_id] = []
            by_dataset[dataset_id].append(node)
        
        # Keep only one node per dataset (highest intensity)
        for dataset_id, dataset_nodes in by_dataset.items():
            if dataset_nodes:
                # Sort by intensity (descending)
                best_node = max(dataset_nodes, key=lambda n: G.nodes[n].get('intensity', 0))
                filtered_partition[best_node] = comm_id
    
    print(f"Filtered partition from {len(partition)} to {len(filtered_partition)} nodes")
    
    # Use the filtered partition for visualization
    partition = filtered_partition
    
    # Sample the graph if it's too large
    if G.number_of_nodes() > max_nodes or G.number_of_edges() > max_edges:
        print(f"Graph is too large ({G.number_of_nodes()} nodes, {G.number_of_edges()} edges), sampling...")
        G_sample = sample_graph(G, max_nodes, max_edges)
        print(f"Sampled graph has {G_sample.number_of_nodes()} nodes and {G_sample.number_of_edges()} edges")
        G = G_sample
        
        # Filter partition to only include nodes in the sampled graph
        partition = {node: comm_id for node, comm_id in partition.items() if node in G}
    
    # Create a new graph with only intra-community edges
    G_community = nx.Graph()
    
    # Add nodes from the partition
    for node, comm_id in partition.items():
        if node in G.nodes():
            G_community.add_node(node, **G.nodes[node])
    
    # Add edges ONLY within the same community
    for u, v in G.edges():
        if u in partition and v in partition and partition[u] == partition[v]:
            G_community.add_edge(u, v)
    
    print(f"All inter-community edges have been removed.")
    print(f"Community graph: {G_community.number_of_nodes()} nodes, {G_community.number_of_edges()} edges")
    
    # Group nodes by community
    communities = defaultdict(list)
    for node, comm_id in partition.items():
        if node in G_community.nodes():
            communities[comm_id].append(node)
    
    # Sort communities by size
    sorted_communities = sorted(communities.items(), key=lambda x: len(x[1]), reverse=True)
    
    # Create layout - always use a normal spring layout with community-aware initialization
    # First, position communities far apart from each other
    pos_init = {}
    
    # Place communities in a grid with significant spacing for initial positions
    grid_size = int(np.ceil(np.sqrt(len(sorted_communities))))
    spacing = 20.0  # Large spacing between initial community positions
    
    for i, (comm_id, nodes) in enumerate(sorted_communities):
        # Calculate grid position
        row = i // grid_size
        col = i % grid_size
        center_x = col * spacing
        center_y = row * spacing
        
        # Place all nodes in this community at the same initial position with small jitter
        for node in nodes:
            pos_init[node] = np.array([center_x, center_y]) + np.random.uniform(-0.1, 0.1, size=2)
    
    # Now run spring layout with these initial positions
    # Use a higher k value to keep communities more separated
    k_value = 5.0 if hard_separation else 2.0  # Higher k means more separation
    iterations = 200 if hard_separation else 150  # More iterations for better convergence
    
    pos = nx.spring_layout(
        G_community,
        pos=pos_init,
        k=k_value,
        iterations=iterations,
        seed=42
    )
    
    # Ensure all nodes in G_community have positions
    missing_nodes = [node for node in G_community.nodes() if node not in pos]
    if missing_nodes:
        print(f"Warning: {len(missing_nodes)} nodes have no position. Adding random positions.")
        for node in missing_nodes:
            pos[node] = np.random.uniform(-10, 10, size=2)
    
    # Create figure
    plt.figure(figsize=(16, 14))
    
    # Draw edges with increased visibility
    nx.draw_networkx_edges(G_community, pos, width=1.0, alpha=0.7, edge_color='gray')
    
    # Draw nodes colored by dataset
    node_colors = []
    for node in G_community.nodes():
        dataset_id = G_community.nodes[node].get('dataset_id', 0)
        node_colors.append(dataset_id)
    
    nx.draw_networkx_nodes(G_community, pos, node_size=80, node_color=node_colors, cmap=plt.cm.rainbow, alpha=0.8)
    
    # Add labels for some nodes (e.g., high-degree nodes)
    if G_community.number_of_nodes() <= 100:  # Only label if not too many nodes
        nx.draw_networkx_labels(G_community, pos, font_size=8)
    else:
        # Label only high-degree nodes
        high_degree_nodes = sorted(G_community.degree, key=lambda x: x[1], reverse=True)[:20]
        labels = {node: str(node) for node, _ in high_degree_nodes}
        nx.draw_networkx_labels(G_community, pos, labels=labels, font_size=8)
    
    # Add a legend for dataset colors
    unique_datasets = sorted(set(node_colors))
    legend_elements = []
    for dataset_id in unique_datasets:
        color = plt.cm.rainbow(dataset_id / max(1, max(unique_datasets)))
        legend_elements.append(plt.Line2D([0], [0], marker='o', color='w', 
                               markerfacecolor=color, markersize=10, 
                               label=f'Dataset {dataset_id}'))
    
    plt.legend(handles=legend_elements, loc='upper right', title='Datasets')
    
    # Set title and remove axis
    plt.title(f"Community Graph {'(with Hard Separation)' if hard_separation else '(with Normal Layout)'}")
    plt.axis('off')
    
    # Save figure
    output_path = os.path.join(output_dir, f"community_graph{'_hard_separation' if hard_separation else ''}.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Community graph visualization saved to {output_path}")

def plot_clique_graph(G, cliques, output_dir, pos=None, max_nodes=1000, max_edges=5000, hard_separation=False):
    """
    Plot the graph with clique detection results.
    
    Steps:
    1. Sample the graph if it's too large
    2. Filter cliques to ensure at most one node per dataset per clique
    3. Create a new graph with ONLY intra-clique edges (remove ALL inter-clique edges)
    4. Apply a layout that properly separates cliques
    5. Draw nodes colored by dataset
    6. Draw edges
    7. Save the visualization
    
    Parameters:
    -----------
    G : networkx.Graph
        Graph to plot
    cliques : list
        List of cliques (each clique is a list of node IDs)
    output_dir : str
        Directory to save the plot
    pos : dict
        Dictionary mapping node IDs to positions (ignored if hard_separation is True)
    max_nodes : int
        Maximum number of nodes to plot
    max_edges : int
        Maximum number of edges to plot
    hard_separation : bool
        If True, use a layout that strongly separates cliques
    """
    if hard_separation:
        print("Plotting clique graph with hard separation...")
    else:
        print("Plotting clique graph...")
    
    # Filter cliques to ensure at most one node per dataset per clique
    filtered_cliques = []
    for clique in cliques:
        # Group nodes by dataset
        by_dataset = {}
        for node in clique:
            dataset_id = G.nodes[node].get('dataset_id')
            if dataset_id not in by_dataset:
                by_dataset[dataset_id] = []
            by_dataset[dataset_id].append(node)
        
        # Keep only one node per dataset (highest intensity)
        filtered_clique = []
        for dataset_id, dataset_nodes in by_dataset.items():
            if dataset_nodes:
                # Sort by intensity (descending)
                best_node = max(dataset_nodes, key=lambda n: G.nodes[n].get('intensity', 0))
                filtered_clique.append(best_node)
        
        # Only keep cliques with at least 3 nodes
        if len(filtered_clique) >= 3:
            filtered_cliques.append(filtered_clique)
    
    print(f"Filtered cliques from {len(cliques)} to {len(filtered_cliques)}")
    
    # Use the filtered cliques for visualization
    cliques = filtered_cliques
    
    # Sample the graph if it's too large
    if G.number_of_nodes() > max_nodes or G.number_of_edges() > max_edges:
        print(f"Graph is too large ({G.number_of_nodes()} nodes, {G.number_of_edges()} edges), sampling...")
        G_sample = sample_graph(G, max_nodes, max_edges)
        print(f"Sampled graph has {G_sample.number_of_nodes()} nodes and {G_sample.number_of_edges()} edges")
        G = G_sample
        
        # Filter cliques to only include nodes in the sampled graph
        filtered_cliques = []
        for clique in cliques:
            filtered_clique = [node for node in clique if node in G]
            if len(filtered_clique) >= 3:  # Only keep cliques with at least 3 nodes
                filtered_cliques.append(filtered_clique)
        cliques = filtered_cliques
    
    # Create a mapping from node to clique ID for coloring
    node_to_clique = {}
    for i, clique in enumerate(cliques):
        for node in clique:
            if node not in node_to_clique:
                node_to_clique[node] = i
    
    # Create a new graph with ONLY intra-clique edges
    G_clique = nx.Graph()
    
    # Add only nodes that are in cliques
    for node in G.nodes():
        if node in node_to_clique:
            G_clique.add_node(node, **G.nodes[node])
    
    # Add ONLY intra-clique edges (edges between nodes in the same clique)
    for u, v, data in G.edges(data=True):
        if (u in node_to_clique and v in node_to_clique and 
            node_to_clique[u] == node_to_clique[v]):
            G_clique.add_edge(u, v, **data)
    
    print(f"Clique graph has {G_clique.number_of_nodes()} nodes and {G_clique.number_of_edges()} edges")
    print(f"All inter-clique edges have been removed")
    
    # Group nodes by clique
    clique_nodes = {}
    for node, clique_id in node_to_clique.items():
        if node in G_clique:
            if clique_id not in clique_nodes:
                clique_nodes[clique_id] = []
            clique_nodes[clique_id].append(node)
    
    # Sort cliques by size
    sorted_cliques = sorted(clique_nodes.items(), key=lambda x: len(x[1]), reverse=True)
    
    # Create layout - always use a normal spring layout with clique-aware initialization
    # First, position cliques far apart from each other
    pos_init = {}
    
    # Place cliques in a grid with significant spacing for initial positions
    grid_size = int(np.ceil(np.sqrt(len(sorted_cliques))))
    spacing = 20.0  # Large spacing between initial clique positions
    
    for i, (clique_id, nodes) in enumerate(sorted_cliques):
        # Calculate grid position
        row = i // grid_size
        col = i % grid_size
        center_x = col * spacing
        center_y = row * spacing
        
        # Place all nodes in this clique at the same initial position with small jitter
        for node in nodes:
            pos_init[node] = np.array([center_x, center_y]) + np.random.uniform(-0.1, 0.1, size=2)
    
    # Now run spring layout with these initial positions
    # Use a higher k value to keep cliques more separated
    k_value = 5.0 if hard_separation else 2.0  # Higher k means more separation
    iterations = 200 if hard_separation else 150  # More iterations for better convergence
    
    pos = nx.spring_layout(
        G_clique,
        pos=pos_init,
        k=k_value,
        iterations=iterations,
        seed=42
    )
    
    # Ensure all nodes in G_clique have positions
    missing_nodes = [node for node in G_clique.nodes() if node not in pos]
    if missing_nodes:
        print(f"Warning: {len(missing_nodes)} nodes have no position in clique graph. Adding random positions.")
        for node in missing_nodes:
            pos[node] = np.random.uniform(-10, 10, size=2)
    
    # Create figure
    plt.figure(figsize=(16, 14))
    
    # Draw edges with increased visibility
    nx.draw_networkx_edges(G_clique, pos, width=1.0, alpha=0.7, edge_color='gray')
    
    # Draw nodes colored by dataset
    node_colors = []
    for node in G_clique.nodes():
        dataset_id = G_clique.nodes[node].get('dataset_id', 0)
        node_colors.append(dataset_id)
    
    nx.draw_networkx_nodes(G_clique, pos, node_size=80, node_color=node_colors, cmap=plt.cm.rainbow, alpha=0.8)
    
    # Add labels for some nodes (e.g., high-degree nodes)
    if G_clique.number_of_nodes() <= 100:  # Only label if not too many nodes
        nx.draw_networkx_labels(G_clique, pos, font_size=8)
    else:
        # Label only high-degree nodes
        high_degree_nodes = sorted(G_clique.degree, key=lambda x: x[1], reverse=True)[:20]
        labels = {node: str(node) for node, _ in high_degree_nodes}
        nx.draw_networkx_labels(G_clique, pos, labels=labels, font_size=8)
    
    # Add a legend for dataset colors
    unique_datasets = sorted(set(node_colors))
    legend_elements = []
    for dataset_id in unique_datasets:
        color = plt.cm.rainbow(dataset_id / max(1, max(unique_datasets)))
        legend_elements.append(plt.Line2D([0], [0], marker='o', color='w', 
                               markerfacecolor=color, markersize=10, 
                               label=f'Dataset {dataset_id}'))
    
    plt.legend(handles=legend_elements, loc='upper right', title='Datasets')
    
    # Set title and remove axis
    plt.title(f"Clique Graph {'(with Hard Separation)' if hard_separation else '(with Normal Layout)'}")
    plt.axis('off')
    
    # Save figure
    output_path = os.path.join(output_dir, f"clique_graph{'_hard_separation' if hard_separation else ''}.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Clique graph visualization saved to {output_path}")

def visualize_subgraph(G, nodes, output_dir, filename, title=None):
    """
    Visualize a subgraph of the original graph.
    
    Parameters:
    -----------
    G : networkx.Graph
        Original graph
    nodes : list
        List of node IDs to include in the subgraph
    output_dir : str
        Directory to save the plot
    filename : str
        Filename for the plot
    title : str
        Title for the plot
    """
    # Create subgraph
    subgraph = G.subgraph(nodes)
    
    # Get layout
    pos = get_spring_layout(subgraph)
    
    # Create figure
    plt.figure(figsize=(10, 8))
    
    # Draw nodes
    nx.draw_networkx_nodes(subgraph, pos, node_size=100, node_color='skyblue', alpha=0.8)
    
    # Draw edges
    nx.draw_networkx_edges(subgraph, pos, width=1.0, alpha=0.5)
    
    # Draw labels
    nx.draw_networkx_labels(subgraph, pos, font_size=8)
    
    # Set title and remove axis
    if title:
        plt.title(title)
    else:
        plt.title(f"Subgraph ({subgraph.number_of_nodes()} nodes, {subgraph.number_of_edges()} edges)")
    plt.axis('off')
    
    # Save figure
    output_path = os.path.join(output_dir, filename)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Subgraph visualization saved to {output_path}")

def create_intensity_heatmap(aligned_features_file, output_dir, max_groups=50):
    """
    Create a heatmap of feature intensities across different files.
    
    Parameters:
    -----------
    aligned_features_file : str
        Path to the TSV file with aligned features
    output_dir : str
        Directory to save the heatmap
    max_groups : int
        Maximum number of groups to include in the heatmap
    """
    print(f"Creating intensity heatmap from {aligned_features_file}...")
    
    # Read aligned features file
    df = pd.read_csv(aligned_features_file, sep='\t')
    
    # Check if intensity columns exist
    intensity_cols = [col for col in df.columns if 'intensity' in col.lower()]
    if not intensity_cols:
        print("No intensity columns found in the aligned features file")
        return
    
    print(f"Found {len(intensity_cols)} intensity columns")
    
    # Limit to max_groups
    if len(df) > max_groups:
        print(f"Limiting to top {max_groups} groups by number of features")
        # Count non-NaN values in intensity columns for each row
        df['feature_count'] = df[intensity_cols].notna().sum(axis=1)
        df = df.sort_values('feature_count', ascending=False).head(max_groups)
        df = df.drop('feature_count', axis=1)
    
    # Create a heatmap of intensities
    intensity_data = df[intensity_cols].copy()
    
    # Convert to float
    for col in intensity_data.columns:
        intensity_data[col] = pd.to_numeric(intensity_data[col], errors='coerce')
    
    # Create heatmap with raw values
    plt.figure(figsize=(12, 10))
    sns.heatmap(intensity_data, cmap='viridis', linewidths=0.5, linecolor='gray')
    plt.title(f"Feature Intensity Heatmap (Raw Values)")
    plt.xlabel("Files")
    plt.ylabel("Feature Groups")
    
    # Save figure
    output_path = os.path.join(output_dir, "intensity_heatmap_raw.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Raw intensity heatmap saved to {output_path}")
    
    # Create heatmap with log-transformed values
    log_intensity_data = np.log1p(intensity_data)  # log(1+x) to handle zeros
    
    plt.figure(figsize=(12, 10))
    sns.heatmap(log_intensity_data, cmap='viridis', linewidths=0.5, linecolor='gray')
    plt.title(f"Feature Intensity Heatmap (Log-transformed)")
    plt.xlabel("Files")
    plt.ylabel("Feature Groups")
    
    # Save figure
    output_path = os.path.join(output_dir, "intensity_heatmap_log.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Log-transformed intensity heatmap saved to {output_path}")

def sample_graph(G, max_nodes=1000, max_edges=5000):
    """
    Sample a graph to reduce its size for visualization.
    
    Parameters:
    -----------
    G : networkx.Graph
        Graph to sample
    max_nodes : int
        Maximum number of nodes in the sampled graph
    max_edges : int
        Maximum number of edges in the sampled graph
        
    Returns:
    --------
    G_sampled : networkx.Graph
        Sampled graph
    """
    G_sampled = nx.Graph()
    
    # Sample nodes
    nodes = list(G.nodes())
    if len(nodes) > max_nodes:
        nodes = random.sample(nodes, max_nodes)
    
    # Add sampled nodes to the graph
    for node in nodes:
        G_sampled.add_node(node, **G.nodes[node])
    
    # Add edges between sampled nodes
    for u, v, data in G.edges(data=True):
        if u in G_sampled and v in G_sampled:
            G_sampled.add_edge(u, v, **data)
    
    # If still too many edges, sample edges
    if G_sampled.number_of_edges() > max_edges:
        edges = list(G_sampled.edges(data=True))
        edges = random.sample(edges, max_edges)
        
        # Create a new graph with sampled edges
        G_sampled_edges = nx.Graph()
        for u, v, data in edges:
            if u not in G_sampled_edges:
                G_sampled_edges.add_node(u, **G.nodes[u])
            if v not in G_sampled_edges:
                G_sampled_edges.add_node(v, **G.nodes[v])
            G_sampled_edges.add_edge(u, v, **data)
        
        G_sampled = G_sampled_edges
    
    print(f"Sampled graph has {G_sampled.number_of_nodes()} nodes and {G_sampled.number_of_edges()} edges")
    return G_sampled 