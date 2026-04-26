**Data Analysis Workflow in R**<br>
This repository contains an R-based workflow for peptide and protein data analysis. It provides a structured pipeline from raw peptide tables to quantified protein outputs.

**Setup and Usage**
1. Configure input data In load_data.R, specify the file path to the peptide table and experiment-specific parameters
2. Define filtering criteria In filter_overview.R: set the filters appropriate for the experiment
3. Run the analysis Execute the main script to:
    * Load all required functions and data
    * Run the full analysis workflow
    * Save the results to the output folder (directory can be defined in load_data.R)

**Analysis Workflow**<br>
The analysis consists of four main parts:
1. Data preparation and overview: prepares the dataset, generates summary statistics on detected peptides and proteins
2. UPS calibration curve: constructs a calibration curve from standard scores
3. Scoring of U2OS proteins: computes scores for protein evaluation
4. Absolute protein quantification: calculates protein amounts in fmol and pg per cell

**Repository Structure**
* main : Controls the workflow execution:
    * Loads required data
    * Sources and runs analysis functions
    * Includes an additional section for sample benchmarking
* functions : Contains required packages and general helper functions
* functions_analysis:  Implements functions for the four main analysis steps
* load_data:  Handles loading of the peptide table and other required inputs
* data_overview : Applies filtering and prepares the dataset for downstream analysis
* combined_speclib:  Combines refined U2OS and predicted UPS spectral libraries
* additional_figures: Contains supplementary files for specific applications
