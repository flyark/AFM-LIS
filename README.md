# AlphaFold-Multimer Local Interaction Score (AFM-LIS)

This repository contains a jupyter notebook designed to calculate the local interaction score (LIS) from ColabFold-derived outputs, including JSON and PDB files. The LIS provides insights into the interaction strength and quality within protein complexes, specifically focusing on the predicted aligned error (PAE) and other related metrics.

## Features

- Extraction and processing of data from `.pdb` and `.json` files.
- Calculation of PAE, pLDDT, and interaction scores between protein complexes.
- Support for multiprocessing to enhance processing speed.
- Generation of comprehensive dataframes summarizing the interaction scores and other relevant metrics.
- Output files are saved in both CSV and Excel formats for easy analysis.

## Requirements

- Python 3.x
- NumPy
- Pandas
- Biopython
- Multiprocessing support

## Installation

To use this script, follow these steps:

1. Navigate to the repository's main page.
2. Locate the Jupyter Notebook file (`.ipynb`) containing the script.
3. Download the Jupyter Notebook to your local machine.

Alternatively, if the repository is cloned or forked, ensure you have Jupyter Notebook installed on your system to run the `.ipynb` file. If you don't have Jupyter Notebook installed, it can be installed via pip:


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
bioRxiv (2024) doi: XXXXX

FlyPredictome
https://www.flyrnai.org/tools/fly_predictome/web/
