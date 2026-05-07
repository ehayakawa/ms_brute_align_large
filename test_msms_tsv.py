#!/usr/bin/env python3
"""
Test script to demonstrate the enhanced TSV output with MS/MS matching information.

This script creates a mock example showing how the enhanced write_aligned_features_tsv
function will display MS/MS matching information in the output TSV files.
"""

import networkx as nx
import os
from collections import defaultdict

# Mock data structures to demonstrate the functionality
def create_test_data():
    """Create mock test data to demonstrate MS/MS enhancement"""
    
    # Mock aligned features (community detection format)
    aligned_features = {
        1: [
            {'dataset_id': 0, 'feature_id': 1622},
            {'dataset_id': 1, 'feature_id': 3317}, 
            {'dataset_id': 2, 'feature_id': 1261}
        ],
        2: [
            {'dataset_id': 0, 'feature_id': 1535},
            {'dataset_id': 1, 'feature_id': 3318},
            {'dataset_id': 3, 'feature_id': 892}
        ]
    }
    
    # Mock file list
    all_list_features = [
        ("dataset1.xlsx", [{'mz': 312.3265, 'intensity': 15000}, None, {'mz': 312.3265, 'intensity': 12000}]),
        ("dataset2.xlsx", [None, None, None, {'mz': 312.3271, 'intensity': 18000}]),
        ("dataset3.xlsx", [{'mz': 312.3264, 'intensity': 9000}, None, {'mz': 312.3258, 'intensity': 14000}]),
        ("dataset4.xlsx", [{'mz': 312.3267, 'intensity': 11000}])
    ]
    
    # Create mock graph with MS/MS edges
    G = nx.Graph()
    
    # Add nodes
    G.add_node("0_1622", dataset_id=0, feature_id=1622)
    G.add_node("1_3317", dataset_id=1, feature_id=3317)
    G.add_node("2_1261", dataset_id=2, feature_id=1261)
    G.add_node("0_1535", dataset_id=0, feature_id=1535)
    G.add_node("1_3318", dataset_id=1, feature_id=3318)
    G.add_node("3_892", dataset_id=3, feature_id=892)
    
    # Add MS/MS edges (these represent features matched by MS/MS similarity)
    G.add_edge("0_1622", "1_3317", 
              weight=0.85, edge_type='msms', 
              cosine_similarity=0.85, shared_peaks=12)
    G.add_edge("1_3317", "2_1261", 
              weight=0.72, edge_type='msms',
              cosine_similarity=0.72, shared_peaks=8)
    
    # Add m/z/RT edges (traditional similarity)
    G.add_edge("0_1535", "1_3318", 
              weight=0.92, edge_type='mz_rt')
    G.add_edge("1_3318", "3_892", 
              weight=0.88, edge_type='mz_rt')
    
    return aligned_features, all_list_features, G

def demonstrate_msms_matching_info():
    """Demonstrate the get_msms_matching_info function"""
    
    print("=== MS/MS Matching Information Enhancement Demo ===\\n")
    
    aligned_features, all_list_features, G = create_test_data()
    
    # Import the function (note: this would normally be from mass_feature_aligner)
    from mass_feature_aligner import get_msms_matching_info
    
    print("Test Data Overview:")
    print(f"- {len(aligned_features)} aligned feature groups")
    print(f"- {len(all_list_features)} datasets")
    print(f"- Graph has {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")
    
    # Check MS/MS edges in graph
    msms_edges = [(u, v, d) for u, v, d in G.edges(data=True) if d.get('edge_type') == 'msms']
    print(f"- {len(msms_edges)} MS/MS similarity edges")
    
    print("\\nMS/MS Edge Details:")
    for u, v, data in msms_edges:
        print(f"  {u} ↔ {v}: cosine={data['cosine_similarity']:.3f}, shared_peaks={data['shared_peaks']}")
    
    print("\\nAnalyzing MS/MS matches for each group:")
    
    for group_id, features in aligned_features.items():
        print(f"\\nGroup {group_id}:")
        print(f"  Features: {[(f['dataset_id'], f['feature_id']) for f in features]}")
        
        # Get MS/MS matching information
        msms_matches = get_msms_matching_info(features, G)
        
        if msms_matches:
            print(f"  MS/MS matches found: {len(msms_matches)}")
            for node1, node2, cosine_sim, shared_peaks in msms_matches:
                dataset1, feature1 = node1.split('_', 1)
                dataset2, feature2 = node2.split('_', 1)
                dataset1_name = os.path.basename(all_list_features[int(dataset1)][0])
                dataset2_name = os.path.basename(all_list_features[int(dataset2)][0])
                print(f"    {dataset1_name}(feature_{feature1}) ↔ {dataset2_name}(feature_{feature2})")
                print(f"    → Cosine similarity: {cosine_sim:.3f}, Shared peaks: {shared_peaks}")
        else:
            print("  No MS/MS matches (features aligned by m/z/RT similarity only)")

def show_enhanced_tsv_format():
    """Show what the enhanced TSV output will look like"""
    
    print("\\n\\n=== Enhanced TSV Output Format ===\\n")
    
    print("New columns added to TSV output:")
    print("- MSMS_Matches: Number of MS/MS similarity matches within the group")
    print("- MSMS_Details: Detailed information about each MS/MS match")
    
    print("\\nExample enhanced TSV output:")
    print("-" * 100)
    
    # Mock example of enhanced output
    header = ["Group ID", "dataset1.xlsx_feature_index", "dataset1.xlsx_mz", "dataset1.xlsx_intensity",
              "dataset2.xlsx_feature_index", "dataset2.xlsx_mz", "dataset2.xlsx_intensity", 
              "dataset3.xlsx_feature_index", "dataset3.xlsx_mz", "dataset3.xlsx_intensity",
              "MSMS_Matches", "MSMS_Details"]
    
    print("\\t".join(header))
    print("-" * 100)
    
    # Example rows
    rows = [
        ["Group_1", "1622", "312.3265", "15000", "3317", "312.3271", "18000", "1261", "312.3264", "9000",
         "2", "dataset1.xlsx(1622)-dataset2.xlsx(3317):cos=0.850,peaks=12; dataset2.xlsx(3317)-dataset3.xlsx(1261):cos=0.720,peaks=8"],
        ["Group_2", "1535", "312.3265", "12000", "3318", "312.3271", "14000", "", "", "",
         "0", "No MS/MS matches"]
    ]
    
    for row in rows:
        print("\\t".join(map(str, row)))

if __name__ == "__main__":
    try:
        demonstrate_msms_matching_info()
        show_enhanced_tsv_format()
        print("\\n✅ MS/MS enhancement demo completed successfully!")
        
    except ImportError as e:
        print(f"⚠️  Cannot import enhanced module (expected if dependencies not installed): {e}")
        print("\\nShowing conceptual output format:")
        show_enhanced_tsv_format()
        
    except Exception as e:
        print(f"❌ Error during demo: {e}")
        print("\\nShowing conceptual output format:")
        show_enhanced_tsv_format()