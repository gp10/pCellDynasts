# Multimodality tests of H2BGFP dilution data
Set of tools to analyze multimodality of individual-cell H2BGFP intensity distributions.

Proliferating cells dilute their H2BGFP content two-fold upon every cell divison. Multimodality should be indicative of distinct subpopulations of cells dividing at different rates, while unimodality conforms to the single-progenitor (SP) model paradigm.

## Getting started 

Install packages: *multimode* (most calculations), *modes* (for bimodality coefficient), *foreach*, and *doMC* (for parallelised analysis) (*iterators* and *parallel* may be required too).

Each of the R scripts (rsx) reads the different datasets as variables for analysis. The variables are created from original files (typically txt) by the 'read_files.rsx' script. To start working source the file (i.e. 'source("x.rsx")')

### Main scripts/functions
* read_files.rsx
  * Create a set of variables by reading the files
* all_graph.rsx
  * Draw a small histogram of each different dataset, by animal/field of view
* modetest.rsx
  * a set of convenience functions to simplify running multimode in batch, including *matrixSearch*

### Analysis
To analyse all animals, and report the p-value of the analysis:

result <- matrixSearch(everything)

To analyse all fields of view, reporting which have a p lower than 0.05 (handy for testing synthetic datasets)

result <- matrixSearch(fov,p=0.05)

To run modality tests on a particular subset of data (e.g. data from esophagus, time=0d, mouse 6c):

result <- matrixSearch(everything$eso00mouse6c)

### Requirements
R 3.5.2		(Rstudio 1.0.136 was used)