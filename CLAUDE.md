# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is a graph-based mass spectrometry feature alignment system that identifies same chemical compounds across multiple LC-MS datasets. It addresses retention time drift and m/z variations using community detection or clique finding algorithms.

## Development Commands

```bash
# Install dependencies using Pipenv
pipenv install
pipenv shell

# Alternative: install with pip
pip install -r requirements.txt

# Run the main alignment process
python main.py --input-dir "path/to/excel/files" --visualize

# Common parameters:
# --mz-tolerance: m/z tolerance in Da (default: 0.01)
# --rt-tolerance: RT tolerance in minutes (default: 0.5)
# --min-datasets: Minimum datasets for valid group (default: 2)
```

## Architecture

The system follows a pipeline architecture:

1. **Data Input** (`read_files.py`): Reads Excel files with mass spectrometry data
2. **Graph Construction** (`graph_construction.py`): 
   - Creates nodes for each mass feature
   - Builds edges based on m/z/RT similarity using KD-trees
   - Cleans multiple connections between dataset pairs
3. **Feature Grouping**:
   - **Community Detection** (`community_detection.py`): Louvain algorithm for broader groups
   - **Clique Detection** (`clique_detection.py`): Maximal cliques for strict groups
4. **Output Generation** (`mass_feature_aligner.py`): Creates TSV files with aligned features
5. **Visualization** (`visualize_graph.py`): Generates graphs and heatmaps

## Key Design Decisions

- **Graph-based approach**: Features as nodes, similarities as edges
- **Dual alignment methods**: Community (soft) vs Clique (hard) detection
- **Pre-processing**: Resolves obvious multiple connections before grouping
- **KD-tree optimization**: Efficient similarity searching for large datasets

## Module Guidelines

Per `.cursorrules`:
- Stick to existing file structure
- Only add new files when necessary
- Maintain modular architecture with clear separation of concerns

## Input/Output

**Input**: Excel files with columns:
- Peak ID, Scan, RT [min], Precursor m/z, Height

**Output**: 
- TSV files with aligned features (community and clique)
- PNG visualizations (graphs and heatmaps)
- Summary statistics