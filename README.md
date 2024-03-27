# AlphaFold-Multimer Local Interaction Score (AFM-LIS)
This repository contains a Jupyter notebook designed to calculate the local interaction score (LIS) from ColabFold-derived outputs, including JSON and PDB files. The LIS provides insights into the interaction strength and quality within protein complexes, specifically focusing on the predicted aligned error (PAE).

## About
This code is based on research described in the paper **Enhanced Protein-Protein Interaction Discovery via AlphaFold-Multimer** ([Kim et al., 2024](https://www.biorxiv.org/content/10.1101/2024.02.19.580970v1)). In our study, we found that although the structural accuracy of predicted protein structures by AlphaFold may be less accurate, the Local Interaction Score (LIS) and the Local Interaction Area (LIA) calculated by this strategy offer a valuable means of predicting positive protein-protein interactions (PPI). By considering both LIS and LIA, our method enables researchers to identify potential interactions, even in cases where the structural predictions may have limitations.

![Figure 2A](https://github.com/flyark/AFM-LIS/raw/main/Figure%202A.png)

> **Illustration of the algorithm for identifying the Local Interaction Area (LIA) and calculating the Local Interaction Score (LIS).** Following AlphaFold-Multimer predictions, Predicted Aligned Error (PAE) values, indicative of individual amino acid interaction likelihood, are extracted from ColabFold outputs. A PAE cutoff is applied, ignoring amino acid contacts with higher PAE values while retaining those below the threshold as part of the LIA. PAE values within the LIA are inverted to a scale of 0 to 1, where higher values denote stronger predicted interactions. The average of these inverted PAE values at the interaction interfaces yields the LIS.

## Interpreting Results for PPI Predictions

To predict positive PPI, assess LIS and LIA against these thresholds, derived from analysis of fly and human protein datasets:

**A positive PPI** is suggested if **either** of the following conditions is met:
- **Best LIS** ≥ 0.203 AND **Best LIA** ≥ 3432, or
- **Average LIS** ≥ 0.073 AND **Average LIA** ≥ 1610.

These cutoff values have been established through analysis of fly and human protein reference sets, detailed in [Tang et al. (2023)](https://www.nature.com/articles/s41467-023-37876-0) for the fly dataset and [Braun et al. (2009)](https://www.nature.com/articles/nmeth.1281) for the human dataset. For a detailed explanation of how these values were derived, please refer to [our paper](https://www.biorxiv.org/content/10.1101/2024.02.19.580970v1). Given the specificity of these values to particular reference sets, it may be necessary to adjust the cutoffs when analyzing different types of PPIs. To ensure the biological relevance of predicted interactions, **experimental validation is highly recommended.**









<img width="864" alt="image" src="https://github.com/flyark/AFM-LIS/assets/26104411/f4a2a5f5-4a8b-46f8-b1c9-e0feb5ba5557">

> **Examples of positive PPIs with flexible regions (decent LIS & low ipTM).**


## Requirements
- Python (3.6.7)
- NumPy (1.19.4)
- Pandas (1.1.4)
- Biopython (1.78)
- Multiprocessing support (fork)

Tested versions are indicated.

## Installation

For ColabFold-derived output (with json), use **alphafold_interaction_scores_github_20240209.ipynb**.

For original AlphaFold-derived output (with pkl), use **alphafold_interaction_scores_for alphafold-multimer data_v0.X.ipynb**

To use this script, follow these steps:

1. Navigate to the repository's main page.
2. Locate the Jupyter Notebook file (`.ipynb`) containing the script.
3. Download the Jupyter Notebook to your local machine.

Alternatively, if the repository is cloned or forked, ensure you have Jupyter Notebook installed on your system to run the `.ipynb` file. If you don't have Jupyter Notebook installed, it can be installed via pip:

## Features
- Extraction and processing of data from `.pdb` and `.json` files.
- Calculation of local interaction score (LIS) and area (LIA) between protein complexes.
- Support for multiprocessing to enhance processing speed.
- Generation of comprehensive data frames summarizing the interaction scores and other relevant metrics.
- Output files are saved in both CSV and Excel formats for easy analysis.

## Usage

To calculate the Local Interaction Scores using the Jupyter Notebook, follow these steps:

1. Ensure your ColabFold output (JSON and PDB files) is placed in a designated directory on your local machine. Please note that the code is originally designed to work with file names containing two protein names separated by `___`. If a different separator is used in your file names, you must either modify them to use `___` or adjust the code within the notebook to accommodate the separator you have used.
2. Open the Jupyter Notebook you downloaded or cloned from the repository.
3. In the notebook, update the `base_path` and `saving_base_path` variables to point to your input directory (where your JSON and PDB files are located) and your output directory (where you want the results to be saved), respectively.
4. Execute the cells in the Jupyter Notebook sequentially. This can be done by selecting each cell and clicking the "Run" button in the Jupyter interface or by pressing `Shift + Enter` to run the selected cell and move to the next one.
5. The notebook will process all `.pdb` files in the input directory, calculate the local interaction scores, and save the results in the specified output directory in both CSV and Excel formats for further analysis.

Make sure you have Jupyter Notebook running in an environment that has all the necessary dependencies installed. If you encounter any dependency-related issues, refer to the `Requirements` section to ensure all required packages are installed.

## Declaration of generative AI usage
This project utilized OpenAI's ChatGPT to assist in generating Python code, documentation, or other textual content.

## Reference

Kim AR, Hu Y, Comjean A, Rodiger J, Mohr SE, Perrimon N. "Enhanced Protein-Protein Interaction Discovery via AlphaFold-Multimer"
bioRxiv (2024) [paper link](https://www.biorxiv.org/content/10.1101/2024.02.19.580970v1)

FlyPredictome
https://www.flyrnai.org/tools/fly_predictome/web/
