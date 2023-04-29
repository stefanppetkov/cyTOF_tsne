# Analysis of cyTOF data acquired from CD4+ T cells by tSNE

A sample workflow which can be applied to both standard flow cytomertic and CyTOF data.
The analysis includes the following steps:

1. Data import and preprocessing
2. Transformation
3. Producing diagnostic plots
  - MDS
  - NRS
4. Cell population identification with FlowSOM and ConsensusClusterPlus
5. Visual representation with tSNE/UMAP
6. Cluster merging and annotation
7. Differential analysis
  - Differential cell population abundance
  - Differential analysis of overall marker expression
