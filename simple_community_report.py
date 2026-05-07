#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Module for generating simplified community detection reports.

This module creates concise, easy-to-read reports from community detection
results, focusing on the most important metrics and findings. It provides
a streamlined alternative to the full community_report.py.

Main functions/classes:
    - generate_simple_report: Creates simplified alignment report
    - summarize_communities: Provides quick summary of communities
    - plot_basic_stats: Generates basic statistical plots

Inputs:
    - TSV files from community detection alignment
    - Basic configuration parameters

Outputs:
    - Simple text/CSV reports with key metrics
    - Basic plots showing alignment results
    - Summary statistics file

Important arguments:
    - input_file: Path to community alignment TSV file
    - output_file: Path for simplified report output
"""

import os
import argparse
import pandas as pd

def generate_report_from_tsv(aligned_file, output_file):
    """
    Generate a simple community report from the aligned features TSV file.
    
    Parameters:
    -----------
    aligned_file : str
        Path to the aligned features TSV file
    output_file : str
        Path to save the report
    """
    if not os.path.exists(aligned_file):
        print(f"Aligned features file not found: {aligned_file}")
        return
    
    # Read the aligned features file
    df = pd.read_csv(aligned_file, sep='\t')
    
    # Create report
    with open(output_file, 'w') as f:
        f.write("community, node1 (filename, rt, mz), node2 (filename, rt, mz), ...\n")
        
        for _, row in df.iterrows():
            group_id = row['group_id']
            line = f"Community {group_id}"
            
            # Find all node columns
            node_cols = [col for col in df.columns if 'node_id' in col]
            rt_cols = [col for col in df.columns if 'rt' in col and 'node' not in col]
            mz_cols = [col for col in df.columns if 'mz' in col and 'node' not in col and 'avg' not in col]
            filename_cols = [col for col in df.columns if 'filename' in col]
            
            # Add each node to the line
            for i, node_col in enumerate(node_cols):
                if pd.isna(row[node_col]):
                    continue
                
                node_id = int(row[node_col])
                
                # Get corresponding rt, mz, and filename if available
                rt = "N/A"
                mz = "N/A"
                filename = "Unknown"
                
                if i < len(rt_cols) and not pd.isna(row[rt_cols[i]]):
                    rt = f"{row[rt_cols[i]]:.2f}"
                
                if i < len(mz_cols) and not pd.isna(row[mz_cols[i]]):
                    mz = f"{row[mz_cols[i]]:.4f}"
                
                if i < len(filename_cols) and not pd.isna(row[filename_cols[i]]):
                    filename = os.path.basename(row[filename_cols[i]])
                
                line += f", Node {node_id} ({filename}, {rt}, {mz})"
            
            f.write(line + "\n")
    
    print(f"Simple community report saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Generate simple community report from aligned features')
    parser.add_argument('--output-dir', type=str, default='output', help='Directory containing output files')
    parser.add_argument('--aligned-file', type=str, default=None, help='Path to the aligned features TSV file')
    parser.add_argument('--output-file', type=str, default=None, help='Path to save the report')
    args = parser.parse_args()
    
    # Set default files if not specified
    if args.aligned_file is None:
        args.aligned_file = os.path.join(args.output_dir, "aligned_features_community.tsv")
    
    if args.output_file is None:
        args.output_file = os.path.join(args.output_dir, "community_report_simple.txt")
    
    # Generate report
    generate_report_from_tsv(args.aligned_file, args.output_file)

if __name__ == "__main__":
    main() 