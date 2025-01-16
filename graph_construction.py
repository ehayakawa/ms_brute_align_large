from typing import List, Dict
import networkx as nx
import numpy as np
from scipy.spatial import KDTree

class GraphBuilder:
    def __init__(self, mz_tolerance: float = 0.01, rt_tolerance: float = 1.0):
        self.mz_tolerance = mz_tolerance
        self.rt_tolerance = rt_tolerance
        self.graph = nx.Graph()

    def calculate_similarity_score(self, feature1: Dict, feature2: Dict) -> float:
        """
        Calculate similarity score between two features
        Returns score between 0 and 1
        """
        mz_diff = abs(feature1['mz'] - feature2['mz'])
        rt_diff = abs(feature1['rt'] - feature2['rt'])
        
        # Normalize differences by tolerances
        mz_score = 1 - (mz_diff / self.mz_tolerance)
        rt_score = 1 - (rt_diff / self.rt_tolerance)
        
        # Weight m/z more heavily than RT (0.7 vs 0.3)
        return 0.7 * mz_score + 0.3 * rt_score if mz_score > 0 and rt_score > 0 else 0

    def build_graph(self, features: List[Dict]) -> nx.Graph:
        """
        Construct graph using KDTree for efficient neighbor finding
        """
        # Create KDTree for efficient searching
        feature_array = np.array([[f['mz'], f['rt']] for f in features])
        kdtree = KDTree(feature_array)
        
        # Query pairs within tolerance
        search_radius = np.sqrt(self.mz_tolerance**2 + self.rt_tolerance**2)
        pairs = kdtree.query_pairs(r=search_radius)
        
        # Add edges with similarity scores
        for i, j in pairs:
            feature1, feature2 = features[i], features[j]
            
            # Skip if features are from same dataset
            if feature1['dataset_id'] == feature2['dataset_id']:
                continue
                
            # Check if within individual tolerances
            if (abs(feature1['mz'] - feature2['mz']) <= self.mz_tolerance and 
                abs(feature1['rt'] - feature2['rt']) <= self.rt_tolerance):
                
                score = self.calculate_similarity_score(feature1, feature2)
                if score > 0:
                    self.graph.add_edge(feature1['id'], feature2['id'], 
                                      weight=score,
                                      mz_diff=abs(feature1['mz'] - feature2['mz']),
                                      rt_diff=abs(feature1['rt'] - feature2['rt']))
        
        return self.graph

    def clean_multiple_connections(self) -> nx.Graph:
        """
        Pre-process graph to resolve obvious cases of multiple connections
        """
        for node in list(self.graph.nodes()):
            neighbors = list(self.graph.neighbors(node))
            by_dataset = {}
            
            # Group neighbors by dataset
            for neighbor in neighbors:
                dataset_id = neighbor.split('_')[0]  # Assuming ID format: "dataset_featurenum"
                by_dataset.setdefault(dataset_id, []).append(neighbor)
            
            # Resolve multiple connections within each dataset
            for dataset, features in by_dataset.items():
                if len(features) > 1:
                    # Keep only the highest scoring edge
                    scores = [(f, self.graph[node][f]['weight']) for f in features]
                    best_feature = max(scores, key=lambda x: x[1])[0]
                    
                    # Remove other edges
                    for feature in features:
                        if feature != best_feature:
                            self.graph.remove_edge(node, feature)
        
        return self.graph

    def get_graph_stats(self) -> Dict:
        """
        Return basic statistics about the graph
        """
        return {
            'num_nodes': self.graph.number_of_nodes(),
            'num_edges': self.graph.number_of_edges(),
            'avg_degree': sum(dict(self.graph.degree()).values()) / self.graph.number_of_nodes(),
            'connected_components': nx.number_connected_components(self.graph)
        }
