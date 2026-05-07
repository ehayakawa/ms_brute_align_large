# Large-Scale Mass Spectrometry Feature Alignment

This project provides a graph-based approach for aligning mass spectrometry features across multiple datasets. It uses community detection and clique finding algorithms to identify related features and group them together.

## Features

- Read mass spectrometry features from Excel files
- Build a graph representation of features based on m/z and retention time similarity
- **MS/MS spectral matching** using cosine similarity for enhanced accuracy
- Detect communities using the Louvain algorithm
- Find cliques in the feature graph
- Align features based on community or clique membership
- Generate visualizations of the feature graph and aligned features
- Create intensity heatmaps for aligned feature groups

## Installation

1. Clone this repository:
```bash
git clone https://github.com/ehayakawa/ms_brute_align_large.git
cd ms_brute_align_large
```

2. Install the required dependencies:
```bash
pip install -r requirements.txt
```

## Usage

### Basic Usage

```bash
python main.py --input-dir "path/to/excel/files" --visualize
```

### Command Line Arguments

- `--input-dir`: Directory containing Excel files with mass features (required)
- `--output-dir`: Directory to save output files (default: "output")
- `--mz-tolerance`: m/z tolerance for feature matching in Da (default: 0.01)
- `--rt-tolerance`: RT tolerance for feature matching in minutes (default: 0.5)
- `--min-datasets`: Minimum number of datasets for a valid feature group (default: 2)
- `--visualize`: Generate visualizations (flag)

## Input Format

The input Excel files should contain mass spectrometry features with the following columns:
- `Peak ID`: Unique identifier for the feature
- `Scan`: Scan number
- `RT [min]`: Retention time in minutes
- `Precursor m/z`: Mass-to-charge ratio
- `Height`: Signal intensity
- `MSMS spectrum` (optional): Fragment spectrum in format `"m/z1 intensity1;m/z2 intensity2;..."`

## Output Files

The program generates the following output files:

- `summary.md`: Summary statistics for each input file
- `aligned_features_community.tsv`: Features aligned using community detection
- `aligned_features_clique.tsv`: Features aligned using clique detection
- `graph.pkl`: Serialized NetworkX graph object
- `partition.pkl`: Serialized community partition data
- `initial_graph.png`: Visualization of the initial feature graph (if `--visualize`)
- `community_graph.png`: Visualization of the graph with communities (if `--visualize`)
- `clique_graph.png`: Visualization of the graph with cliques (if `--visualize`)
- `intensity_heatmap_*.png`: Heatmaps of feature intensities (if `--visualize`)

The TSV output includes MS/MS matching information when available (MSMS_Matches count and MSMS_Details with cosine scores).

## Project Structure

### Core Pipeline
- `main.py`: Main script for running the alignment process
- `read_files.py`: Functions for reading Excel files and extracting features
- `graph_construction.py`: Graph building from mass spectrometry features
- `spectral_similarity.py`: MS/MS cosine similarity calculations
- `community_detection.py`: Community detection using Louvain algorithm
- `clique_detection.py`: Maximal clique finding for strict grouping
- `mass_feature_aligner.py`: Functions for aligning features and writing output
- `visualize_graph.py`: Visualization functions for graphs and heatmaps

### Standalone Tools
- `community_report.py`: CLI tool for detailed community analysis
- `simple_community_report.py`: CLI tool for simplified report generation
- `test_msms_tsv.py`: Test utilities for MS/MS TSV validation

## Example

```bash
python main.py --input-dir "input_data" --mz-tolerance 0.01 --rt-tolerance 0.5 --visualize
```

This will process all Excel files in the "input_data" directory, build a feature graph with the specified tolerances, detect communities and cliques, align features, and generate visualizations.

## License

This project is licensed under the MIT License - see the LICENSE file for details.


## Background
Mass feature is an entity of compound found in LC-MS. Mass feature always have m/z value of ionized compound as well as retention time. Mass feature might have fragment spectrum ( also known as MS2 spectrum).
In this project, there are multiple lists of mass features derived from distinct LC-MS data. One LC-MS data results in one list of mass feature. The samples analyzed on each location are essentially same, meaning chemical composition of those samples are similar. LC-MS data were acquired from different locations, so retention time may shift.


## What I want to do:

I want to build a program/system to align such lists of mass features. "Align" means that same chemical entities found in each LC-MS data which is list of mass feature are distinguished as identical chemical entity based on m/z  and retention time.  Observed m/z value of compound  may shift slightly. For example, even compound X whose theoretical m/z is 100.000 is found in two list of mass features,  m/z of mass features can be slightly different for example m/z = 100.01 or m/z = 99.98.  Therefore we must accept slight m/z error by setting tolerance.  About retention time, because LC-MS could be perfomed at different locations, so same chemical entity can be observed in different retention time.  We must accept retention time difference about 1 min, but the tolerance must be changed if necessary.

After alignment , a data table whose each row corresponds chemical entities observed over LC-MS analyses and and each column means signal intensity of the corresponding mass features in each samples will be produce.


## Workflow

A detailed workflow chart is available in [workflow_chart.md](workflow_chart.md).

The alignment process follows these main steps:

1. **Input Handling**: Read files and extract mass features (m/z, RT, spectra)
2. **Graph Construction**: Create nodes for features, connect similar features with edges
3. **Feature Grouping**: Apply community or clique detection algorithms
4. **Output Generation**: Create alignment tables with grouped features

## Project Overview

### read input files

### Create data structure

### Compare mass features in a pair of datasets

Mass features between lists A and B are compared using two key parameters:
1. m/z value (mass-to-charge ratio)
2. retention time (RT)

For example:
- Consider mass feature mf-A1 from list A with m/z = 100.01 and RT = 5.2 min
- It is compared against all features in list B using defined tolerances:
  - m/z tolerance: ±0.01
  - RT tolerance: ±1.0 min

Matching process:
1. If a mass feature in list B (e.g., mf-B1) falls within both tolerance windows:
   - m/z: 100.01 ± 0.01
   - RT: 5.2 ± 1.0 min
   Then mf-A1 and mf-B1 are considered the same chemical entity

2. Multiple matches can occur. For example:
   - mf-B1 (m/z: 100.00, RT: 5.1 min) matches mf-A1
   - mf-B4 (m/z: 100.02, RT: 5.3 min) also matches mf-A1

Graph representation:
- Each mass feature becomes a node in the graph
- Matching features are connected by edges
- Example: mf-A1 forms edges with both mf-B1 and mf-B4



### Graph representation:
- Each mass feature becomes a node in the graph
- Matching features are connected by edges
- Example network for mf-A1:
  - mf-A1 ↔ mf-B1
  - mf-A1 ↔ mf-B4
  - mf-A1 ↔ mf-C2
  - mf-A1 ↔ mf-D3
  - mf-A1 ↔ mf-D7

This extended graph structure captures all potential matches between features across multiple lists. The resulting network can become complex, requiring sophisticated clique or community detection methods to resolve the final groupings of related features.

### mass feature alignment

### Overview
This approach uses community detection algorithms with a hybrid strategy for handling multiple connections to identify groups of mass features that likely represent the same chemical entity across multiple LC-MS datasets. It's particularly suitable for:
- Large datasets (100-1000 features × 100 datasets)
- Data with RT shifts between datasets
- Cases with multiple potential matches
- Noisy or unstable measurements

### Algorithm Workflow
1. **Graph Construction**
   - Nodes: Mass features from all datasets
   - Edges: Connections based on m/z and RT similarities
   - Edge weights: Combined similarity score

2. **Pre-processing (First Phase Resolution)**
   - Remove obvious mismatches:
     - Features outside m/z tolerance
     - Features with extreme RT differences
   - Resolve clear cases of multiple connections:
     - Keep highest scoring match when difference is significant
     - Retain ambiguous cases for community analysis
   - Benefits:
     - Reduces graph complexity
     - Handles obvious cases early
     - Preserves important ambiguous connections

3. **Community Detection**
   - Algorithm: Louvain method
   - Run on partially cleaned graph
   - Maintains some ambiguous connections for context
   - Benefits:
     - Efficient for large networks
     - Handles varying community sizes
     - Natural grouping based on connection patterns

4. **Post-processing (Second Phase Resolution)**
   - Within each community:
     - Resolve remaining multiple connections
     - Use community context for better decisions
     - Apply final validation rules
   - Validation criteria:
     - One feature per dataset maximum
     - Consistent m/z values
     - Reasonable RT distribution
     - Minimum community size

### MS/MS Spectral Matching

When MS/MS data is available, the system uses a **two-case matching approach**:

1. **Case 1: No MS/MS data** - Edge weight based on m/z and RT proximity
2. **Case 2: Both features have MS/MS** - Edge weight based on cosine similarity

MS/MS edges are **prioritized** over m/z/RT edges when resolving multiple connections, ensuring structurally similar compounds are matched over proximity-based matches.

**Cosine similarity parameters:**
- `cosine_threshold`: Minimum similarity for valid MS/MS match (default: 0.5)
- `min_shared_peaks`: Minimum shared peaks required (default: 3)

### Key Parameters
- m/z tolerance: ±0.01 (adjustable)
- RT tolerance: ±1.0 min (adjustable)
- Minimum community size: 3 features
- Maximum m/z variance within community
- Cosine threshold: 0.5 (for MS/MS matching)
- Minimum shared peaks: 3 (for MS/MS matching)

### Example Resolution Process

### Advantages
1. **Scalability**
   - Efficient processing of large datasets
   - Near-linear time complexity
   - Memory-efficient implementation

2. **Robustness**
   - Tolerant to missing connections
   - Handles RT drift naturally
   - Adapts to varying data quality

3. **Flexibility**
   - Adjustable parameters
   - Can incorporate additional validation rules
   - Handles different dataset sizes






## Focus:

Since number of mass features are huge (100-1000) and number of list of mass feature is potentially large (like 100 lists), to make system that handle alignment efficiently is quite important. Also there is possibility that one mass feature in one list shows similarity to multiple features in other list. In such case, you have to judge which mass features should be aligned. This nature of complex data may gives difficulty o make alignment.

Rather than classical way of alignment, I want some clever method. For instance using graph analysis method.
Please   propose strategy to deal with alighment with such large dataset.





### Handling Multiple Feature Connections

#### Problem Description
In mass feature alignment, we often encounter situations where one feature from dataset A connects to multiple features from dataset B. For example:

- Feature mfA1 (m/z: 100.01, RT: 5.2) connects to:
  - mfB1 (m/z: 100.00, RT: 5.1)
  - mfB2 (m/z: 100.02, RT: 5.3)

This creates ambiguity that must be resolved before or during community detection.

#### Resolution Strategy

1. **Pre-processing Step**
   - Before running community detection
   - Examine all connections between datasets
   - For each feature with multiple connections to the same dataset:
     1. Calculate similarity scores
     2. Keep only the best-scoring connection
     3. Remove other connections

2. **Scoring System**
   - Primary factors:
     - m/z difference (higher weight)
     - RT difference (lower weight)
   - Additional factors (optional):
     - MS2 spectrum similarity
     - Peak intensity
     - Peak shape correlation

3. **Selection Process**
   - Calculate combined score for each connection
   - Select highest-scoring match
   - Remove other connections
   - Document the selection in output

#### Integration with Community Detection
1. Build initial graph with all possible connections
2. Apply multiple connection resolution
3. Run community detection on cleaned graph
4. Generate final feature groups

#### Benefits
- Cleaner input for community detection
- More reliable feature grouping
- Reduced ambiguity in results
- Better representation of real chemical entities

#### Example Resolution

## Feature Alignment Methods

This project implements two different approaches for mass feature alignment:

### 1. Clique Detection Method
- Finds groups where all features are directly similar to each other
- **Requirements**: All features in a group must be similar to all other features in that group
- **Example**:  ```
  Given features:
  mfA1 (m/z: 100.01, RT: 5.2)
  mfB1 (m/z: 100.00, RT: 5.1)
  mfC1 (m/z: 100.02, RT: 5.3)

  If:
  mfA1 -- mfB1 (similar)
  mfB1 -- mfC1 (similar)
  mfA1 -- mfC1 (not similar)

  Clique Output:
  Group 1: [mfA1, mfB1]
  Group 2: [mfB1, mfC1]  ```
- **Pros**:
  - Higher confidence alignments
  - Better for reference compounds
  - More precise grouping
- **Cons**:
  - May miss valid alignments
  - Sensitive to noise
  - Produces smaller groups

### 2. Community Detection Method
- Finds groups where features can be connected through intermediate features
- **Requirements**: Features in a group need only be connected through other features
- **Example**:  ```
  Using same features as above:
  Community Output:
  Group 1: [mfA1, mfB1, mfC1]  # Connected through mfB1  ```
- **Pros**:
  - More comprehensive alignments
  - Better handles RT drift
  - More robust to noise
- **Cons**:
  - May include less confident alignments
  - Requires parameter tuning
  - Can create larger, less precise groups

### Usage Recommendations
- Use **Clique Detection** when:
  - High confidence is required
  - Working with standard compounds
  - Precise alignment is critical
  
- Use **Community Detection** when:
  - Dealing with unknown compounds
  - Handling data with RT drift
  - More comprehensive alignment is needed
