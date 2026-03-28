# AFM-LIS: Local Interaction Score for Structure Predictions

This repository provides tools to calculate the Local Interaction Score (LIS) from structure prediction outputs for enhanced protein-protein interaction (PPI) discovery. The LIS framework accurately identifies transient and flexible interactions (e.g., involving SLiMs or IDRs) that global confidence metrics like ipTM often miss.

**Supports:** AlphaFold3, ColabFold/AlphaFold-Multimer, Boltz, Chai-1, OpenFold3

The framework includes two validated approaches:
- **iLIS (integrated Local Interaction Score)**: Our newest, robust **single metric** that combines interface confidence with a direct physical contact filter. This is the recommended score for most analyses. ([Kim et al., 2025](https://www.biorxiv.org/content/10.1101/2025.10.10.681672))
- **LIS + LIA (Local Interaction Score + Area)**: Our original dual-metric system that evaluates interactions by assessing both the _score_ (LIS) and the _area_ (LIA) of the predicted interface. ([Kim et al., 2024](https://www.biorxiv.org/content/10.1101/2024.02.19.580970v1))


## How iLIS Works: The Core Concept

The iLIS framework identifies the high-confidence local interface, filters it for direct physical contacts, and combines these scores into a single, robust metric.


**Schematic of the iLIS calculation and benchmark data**.

<img width="948" height="395" alt="image" src="https://github.com/user-attachments/assets/832cbf18-46af-4e33-83b4-6cdfcca1df8b" />

> **Figure 1A-D** from [Kim et al., 2025](https://www.biorxiv.org/content/10.1101/2025.10.10.681672). iLIS outperforms existing metrics for global structure evaluation particularly for flexible protein complexes (low pLDDT groups).



## The Evolution from LIS to iLIS
### 1. The Original Framework: LIS + LIA ###

Standard AlphaFold metrics (like ipTM) often fail on flexible interactions because non-interacting disordered regions can "dilute" the global score. Our original **Local Interaction Score (LIS)** was developed to solve this by focusing only on the confidently predicted interface, defined as all residue pairs with a Predicted Aligned Error (PAE) **≤ 12 Å**.

However, `LIS` _by itself_ had a potential failure mode: it could assign a high score to a few residue pairs that were confidently predicted but physically distant and not in actual contact. Our original (2024) framework addressed this by using `LIS` in combination with a second metric, the **Local Interaction Area (LIA)**, which measures the area in the confident interface (PAE ≤ 12 Å). In the failure case, the `LIA` score would be tiny, successfully filtering out the false positive. This dual-metric `LIS + LIA system` remains a valid approach.

<img width="864" alt="image" src="https://github.com/flyark/AFM-LIS/assets/26104411/f4a2a5f5-4a8b-46f8-b1c9-e0feb5ba5557">

> **Examples of positive PPIs with flexible regions (decent LIS & LIA & low ipTM).** All PPI examples were _experimentally_ proven previously.


### 2. The New Metric: iLIS (integrated LIS) ###

The new **iLIS** metric (2025) streamlines this process into a single, more robust score. It was designed to solve the same problem by directly integrating a physical contact filter.

iLIS is calculated as the geometric mean of two components:
- **LIS (Local Interaction Score)**: The same as the original metric; measures the average confidence of the broad LIA (PAE ≤ 12 Å).
- **cLIS (contact-filtered LIS)**: A stricter metric that measures the average confidence only for residue pairs that are in direct physical contact (PAE ≤ 12 Å **AND** Cβ–Cβ distance ≤ 8 Å).

The final score is:

`iLIS = sqrt(LIS * cLIS)`

This use of a geometric mean is critical. If no direct physical contacts are confidently predicted (`cLIS = 0`), the final `iLIS` score is forced to zero, robustly eliminating this class of false positives. It now allows single metrics to distinguish positive and negative predictions.

<img width="826" height="430" alt="image" src="https://github.com/user-attachments/assets/b04bbb9f-ad4f-49c5-b903-a9f9122a2f66" />

> **Examples of positive PPIs with flexible regions (decent iLIS & low ipTM).** All PPI examples were _experimentally_ proven previously.


## Interpreting Results & Thresholds

We provide two validated frameworks for analysis. **iLIS is the recommended single-metric approach**.

### 1. Recommended: iLIS Framework (Kim et al., 2025) ###

A high-confidence interaction is suggested if:

- **iLIS ≥ 0.223**

This threshold was established using large-scale Y2H reference sets used in yeast (Yu et al., 2008), fly (Tang et al., 2023), and human (Braun et al., 2009). Please see the details in the supplementary text in [Kim et al., 2025](https://www.biorxiv.org/content/10.1101/2025.10.10.681672).

### 2. Original: LIS + LIA Framework (Kim et al., 2024) ###

A positive PPI is suggested if either of the following conditions is met:

- **Best LIS** ≥ 0.203 **AND** **Best LIA** ≥ 3432, or

- **Average LIS** ≥ 0.073 **AND** **Average LIA** ≥ 1610.

These cutoffs were derived from specific fly (Tang et al., 2023) and human (Braun et al., 2009) reference sets. As always, optimal thresholds may vary by dataset, and experimental validation is highly recommended.


## Usage

### Command-Line Tool: `lis.py` (Recommended for batch analysis)

A Python script for calculating LIS metrics from structure prediction outputs.

```bash
# Basic usage — auto-detects platform, processes all models
python lis.py /path/to/prediction_folder/

# Process a zip file
python lis.py alphafold3_output.zip

# Parallel processing with 4 CPUs
python lis.py /path/to/predictions/ -w 4

# Custom output location
python lis.py /path/to/predictions/ -d /path/to/output/

# Show error details for failed predictions
python lis.py /path/to/predictions/ -v

# Reprocess all (ignore existing results)
python lis.py /path/to/predictions/ --no-skip-existing
```

**Features:**
- Auto-detects prediction platform from filenames
- Handles `.gz` and `.xz` compressed files transparently
- Incremental CSV output — safe to interrupt and resume
- Progress bar with elapsed time and ETA
- Parallel processing with `--workers`
- Output sorted by name and rank

**Supported platforms (auto-detected):** AlphaFold3, ColabFold, Boltz-1/2, Chai-1, OpenFold3.

### Web Tool: LIVIA

For visual analysis with no installation required:
- **[LIVIA](https://flyark.github.io/LIVIA/)** (**L**ocal **I**nteraction **VI**sualization and **A**nalysis) — drag & drop prediction files, view PAE/LIS/cLIS maps, generate ChimeraX/PyMOL scripts
- Supports all platforms listed above
- All analysis runs locally in your browser — no data leaves your device
- [Step-by-step tutorial](https://flyark.github.io/LIVIA/tutorials.html) · [GitHub](https://github.com/flyark/LIVIA)


## Output CSV Columns

| Column | Description |
|--------|-------------|
| name | Prediction name |
| rank | Rank number |
| model | Model number |
| iLIS | integrated LIS = √(LIS × cLIS) |
| iLIA | integrated LIA = √(LIA × cLIA) |
| iLISA | iLIS × iLIA |
| ipSAE | interaction prediction Score from Aligned Errors ([Dunbrack, 2025](https://www.biorxiv.org/content/10.1101/2025.02.10.637595)) |
| LIS | Local Interaction Score (PAE ≤ 12Å) |
| cLIS | contact-filtered LIS (PAE ≤ 12Å & Cβ ≤ 8Å) |
| LIA | Local Interaction Area (count of confident residue pairs) |
| cLIA | contact-filtered LIA |
| ipTM | interface predicted TM-score ([Evans et al., 2021](https://doi.org/10.1101/2021.10.04.463034)) |
| pLDDT_i/j | per-chain average pLDDT |
| pTM | predicted TM-score |
| LIR_i/j | count of Local Interaction Residues per chain |
| cLIR_i/j | count of contact-filtered LIR per chain |
| LIpLDDT_i/j | average pLDDT of LIR residues per chain |
| cLIpLDDT_i/j | average pLDDT of cLIR residues per chain |
| len_i/j | chain lengths |
| LIR_indices_i/j | LIR residue positions (compact range format) |
| cLIR_indices_i/j | cLIR residue positions (compact range format) |
| structure_file | source structure filename |


## Requirements

For `lis.py`:
- Python ≥ 3.8
- NumPy
- SciPy


## How to Cite

If you use this work, please cite the relevant paper:
- **For the iLIS metric (Recommended)**: Kim et al (2025). A Structure-Guided Kinase-Transcription Factor Interactome Atlas Reveals Docking Landscapes of the Kinome. _bioRxiv_ ([link](https://www.biorxiv.org/content/10.1101/2025.10.10.681672))
- **For the original LIS + LIA framework**: Kim et al (2024). Enhanced Protein-Protein Interaction Discovery via AlphaFold-Multimer. _bioRxiv_ ([link](https://www.biorxiv.org/content/10.1101/2024.02.19.580970v1))


## FlyPredictome

Large-scale AlphaFold-Multimer PPI predictions in *Drosophila*, scored with the LIS framework:
- **[www.flyrnai.org/tools/fly_predictome/web/](https://www.flyrnai.org/tools/fly_predictome/web/)**
- Now covers **>1.7 million** PPI predictions.


## Declaration of Generative AI Usage

This project utilized OpenAI's ChatGPT, Google's Gemini, and Anthropic's Claude to assist in generating Python code, documentation, or other textual content.
