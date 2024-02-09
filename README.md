# Local Interaction Score Calculation for Protein Complexes

This repository contains a jupyter notebook-based python script designed to calculate the local interaction score (LIS) from ColabFold-derived outputs, including JSON and PDB files. The LIS provides insights into the interaction strength and quality within protein complexes, specifically focusing on the predicted aligned error (PAE) and other related metrics.

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

1. Ensure your ColabFold output (JSON and PDB files) is placed in a designated directory on your local machine.
2. Open the Jupyter Notebook you downloaded or cloned from the repository.
3. In the notebook, update the `base_path` and `saving_base_path` variables to point to your input directory (where your JSON and PDB files are located) and your output directory (where you want the results to be saved), respectively.
4. Execute the cells in the Jupyter Notebook sequentially. This can be done by selecting each cell and clicking the "Run" button in the Jupyter interface or by pressing `Shift + Enter` to run the selected cell and move to the next one.
5. The notebook will process all `.pdb` files in the input directory, calculate the local interaction scores, and save the results in the specified output directory in both CSV and Excel formats for further analysis.

Make sure you have Jupyter Notebook running in an environment that has all the necessary dependencies installed. If you encounter any dependency-related issues, refer to the `Requirements` section to ensure all required packages are installed.

## License

Distributed under the MIT License. See `LICENSE` for more information.

## Acknowledgments

- The jupyter notebook is designed to work with ColabFold's output.
