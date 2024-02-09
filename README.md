# Local Interaction Score Calculation for Protein Complexes

This repository contains a Python script designed to calculate the local interaction score (LIS) from ColabFold-derived outputs, including JSON and PDB files. The LIS provides insights into the interaction strength and quality within protein complexes, specifically focusing on the predicted aligned error (PAE) and other related metrics.

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
