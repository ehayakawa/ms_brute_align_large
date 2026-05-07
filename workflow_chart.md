# Mass Feature Alignment Workflow

```mermaid
flowchart TD
    %% Main workflow
    Start([Start]) --> InputFiles[Read MS Data Files]
    InputFiles --> FeatureExtraction[Extract Mass Features]
    FeatureExtraction --> GraphConstruction[Construct Graph]
    
    %% Graph construction details
    subgraph "Graph Construction"
        CreateNodes[Create Nodes for Features]
        BuildKDTrees[Build KD-Trees]
        CalculateSimilarity[Calculate Feature Similarities]
        AddEdges[Add Edges based on Similarity]
        ResolveMultiple[Resolve Multiple Connections]
        
        CreateNodes --> BuildKDTrees
        BuildKDTrees --> CalculateSimilarity
        CalculateSimilarity --> AddEdges
        AddEdges --> ResolveMultiple
    end
    
    GraphConstruction --> AlignmentMethods{Choose Alignment Method}
    
    %% Community detection path
    AlignmentMethods -->|Community Detection| CommunityDetection
    
    subgraph "Community Detection"
        ApplyLouvain[Apply Louvain Method]
        GroupFeaturesCommunity[Group Features by Community]
        ValidateCommunities[Validate Communities]
        
        ApplyLouvain --> GroupFeaturesCommunity
        GroupFeaturesCommunity --> ValidateCommunities
    end
    
    CommunityDetection --> GenerateCommunityTables[Generate Community Tables]
    GenerateCommunityTables --> WriteCommunityTSV[Write Community TSV]
    WriteCommunityTSV --> PrintCommunitySummary[Print Community Summary]
    
    %% Clique detection path
    AlignmentMethods -->|Clique Detection| CliqueDetection
    
    subgraph "Clique Detection"
        FindMaximalCliques[Find Maximal Cliques]
        GroupFeaturesClique[Group Features by Clique]
        ValidateCliques[Validate Cliques]
        
        FindMaximalCliques --> GroupFeaturesClique
        GroupFeaturesClique --> ValidateCliques
    end
    
    CliqueDetection --> GenerateCliqueTables[Generate Clique Tables]
    GenerateCliqueTables --> WriteCliqueTSV[Write Clique TSV]
    WriteCliqueTSV --> PrintCliqueSummary[Print Clique Summary]
    
    %% Table generation details
    subgraph "Generate Tables"
        AlignFeatures[Align Features]
        CalculateAvgMZ[Calculate Average m/z]
        CalculateAvgRT[Calculate Average RT]
        CalculateIntensity[Calculate Intensities]
        
        AlignFeatures --> CalculateAvgMZ
        CalculateAvgMZ --> CalculateAvgRT
        CalculateAvgRT --> CalculateIntensity
    end
    
    %% End of workflow
    PrintCommunitySummary --> End([End])
    PrintCliqueSummary --> End
    
    %% Styling
    classDef process fill:#f9f9f9,stroke:#333,stroke-width:1px;
    classDef decision fill:#e1f5fe,stroke:#333,stroke-width:1px;
    classDef io fill:#e8f5e9,stroke:#333,stroke-width:1px;
    classDef subgraph fill:#f5f5f5,stroke:#333,stroke-width:1px;
    
    class Start,End io;
    class InputFiles,FeatureExtraction,WriteCommunityTSV,WriteCliqueTSV io;
    class AlignmentMethods decision;
    class GraphConstruction,CommunityDetection,CliqueDetection,GenerateTables process;
```

## Workflow Steps Explanation

1. **Input Processing**
   - Read MS data files from input directory
   - Extract mass features (m/z, RT, intensity, MS2 spectra)
   - Organize features by dataset origin

2. **Graph Construction**
   - Create nodes for each mass feature
   - Build KD-Trees for efficient similarity searching
   - Calculate similarities between features based on m/z and RT
   - Add edges between similar features
   - Resolve multiple connections between datasets

3. **Feature Grouping**
   - Choose between community detection or clique detection
   - Apply selected algorithm to group similar features
   - Validate resulting groups based on criteria

4. **Output Generation**
   - Generate tables with aligned features
   - Calculate average properties for each group
   - Write results to TSV files
   - Print summary statistics 