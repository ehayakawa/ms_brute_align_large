from read_files import read_features, read_excel, collect_files
from graph_construction import GraphBuilder
from community_detection import detect_communities, generate_community_tables
from clique_detection import find_cliques, generate_clique_tables

import csv

def test_read_excel(directory):

    """
    Test function to verify Excel file reading
    - Reads each Excel file in the directory
    - Prints sample data for verification
    """
    _, _, excel_files = collect_files(directory)
    
    if not excel_files:
        print("No Excel files found in the specified directory.")
        return

    for excel_file in excel_files:
        print(f"Reading file: {excel_file}")
        list_features = read_excel(excel_file)
        print(f"Number of features read: {len(list_features)}")
        # Print sample data for verification
        if list_features:
            print("Sample feature:")
            print(list_features[0])
            print("\nSample fragment spectrum:")
            print(list_features[0]['fragment_spectrum'][:5])  # Print first 5 peaks
        
        print("\n" + "="*50 + "\n")

def write_aligned_features_tsv(aligned_features, all_list_features, output_file):
    """
    Write aligned features to TSV file
    - Creates a table with features aligned across different files
    - Each row represents a group of aligned features
    """
    
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        
        # Create header: Group + Feature/m/z columns for each file
        header = ['Group']  # 'Group' can be either 'Community' or 'Clique'
        for file_index, _ in enumerate(all_list_features):
            header.extend([f'File_{file_index}_Feature', f'File_{file_index}_m/z'])
        writer.writerow(header)
        
        # Write data rows
        for group, file_features in aligned_features.items():
            row = [group]
            for file_index, _ in enumerate(all_list_features):
                if file_index in file_features:
                    feature_index = file_features[file_index]
                    feature = all_list_features[file_index][1][feature_index]
                    row.extend([feature_index, f"{feature['precursor_mz']:.4f}"])
                else:
                    row.extend(['', ''])# Empty cells for missing features
            writer.writerow(row)

    print(f"Aligned features have been written to {output_file}")

def main():
    """
    Main workflow of the mass feature alignment process
    """

    # Test Excel reading
    #excel_directory = 'for hayakawa'  # Update this to the correct path of your Excel files
    #test_read_excel(excel_directory)

    # SECTION 1: Data Loading

    # Read Excel files from input directory
    excel_directory = 'input_excel'  # Using the same directory as in test_read_excel
    _, _, excel_files = collect_files(excel_directory)
    # Load all features from Excel files
    all_list_features = []
    for excel_file in excel_files:
        list_features = read_excel(excel_file)
        if list_features:
            all_list_features.append((excel_file, list_features))
    
 
   
    # SECTION 2: Generate Summary Statistics
    # Calculate summary statistics for each file
 
    summary = []
    for filename, features in all_list_features:
        summary.append({
            'filename': filename,
            'num_features': len(features),
            'avg_mz': sum(f['precursor_mz'] for f in features) / len(features) if features else 0,
            'avg_rt': sum(f['retention_time'] for f in features) / len(features) if features else 0,
            'avg_intensity': sum(f['signal_intensity'] for f in features) / len(features) if features else 0
        })

    # Write summary to markdown file
    with open('all_list_features_summary.md', 'w') as f:
        f.write("# Summary of all_list_features\n\n")
        f.write("| Filename | Number of Features | Average m/z | Average RT (min) | Average Intensity |\n")
        f.write("|----------|---------------------|-------------|------------------|-------------------|\n")
        for item in summary:
            f.write(f"| {item['filename']} | {item['num_features']} | {item['avg_mz']:.2f} | {item['avg_rt']:.2f} | {item['avg_intensity']:.2f} |\n")

    print("Summary has been written to 'all_list_features_summary.md'")



    print(f"Number of Excel files processed: {len(all_list_features)}")
    print(f"Total number of features across all files: {sum(len(features) for _, features in all_list_features)}")


    # SECTION 3: Graph Construction
    # Initialize graph builder with parameters
    graph_builder = GraphBuilder(mz_tolerance=0.01, rt_tolerance=1.0)

    # Convert features to dictionary format
    all_mass_features = []
    for file_idx, (filename, features) in enumerate(all_list_features):
        for feat_idx, feature in enumerate(features):
            mass_feature = {
                'id': f"{file_idx}_{feat_idx}",
                'mz': feature['precursor_mz'],
                'rt': feature['retention_time'],
                'intensity': feature.get('signal_intensity', 0.0),
                'dataset_id': str(file_idx),
                'ms2_spectrum': feature.get('fragment_spectrum', None)
            }
            all_mass_features.append(mass_feature)

    # Build similarity graph from features
    G = graph_builder.build_graph(all_mass_features)

    # SECTION 4: Community Detection
    # Detect communities and generate alignment tables
    # Community detection
    partition = detect_communities(G)
    aligned_features_community, aligned_mz_community, aligned_intensity_community = generate_community_tables(G, partition, all_list_features)
    
    # Write community detection results
    output_file_community = 'aligned_features_community.tsv'
    write_aligned_features_tsv(aligned_features_community, all_list_features, output_file_community)

    # SECTION 5: Clique Detection
    # Find cliques and generate alignment tables

    # Clique detection
    cliques = find_cliques(G)
    aligned_features_clique, aligned_mz_clique, aligned_intensity_clique = generate_clique_tables(G, cliques, all_list_features)
    
    # Write aligned features from clique detection to TSV
    output_file_clique = 'aligned_features_clique.tsv'
    write_aligned_features_tsv(aligned_features_clique, all_list_features, output_file_clique)

    print("Aligned features using community detection:")
    print(aligned_features_community)
    print("Aligned m/z values using community detection:")
    print(aligned_mz_community)

    print("Aligned features using clique detection:")
    print(aligned_features_clique)
    print("Aligned m/z values using clique detection:")
    print(aligned_mz_clique)

if __name__ == "__main__":
    main()
