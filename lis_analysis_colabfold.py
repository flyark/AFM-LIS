import os
import sys
import json
import argparse
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from scipy.spatial.distance import pdist, squareform

###############################################################################
# Helper: Parse the PDB once and return the structure
###############################################################################
def get_structure(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("cf", pdb_file)
    return structure

###############################################################################
# 1) Parse PDB to identify chain lengths (using the structure)
###############################################################################
def parse_pdb_chains_from_structure(structure):
    """
    Given a parsed structure, returns:
      - chain_order: list of chain IDs (e.g. ['A','B','C',...])
      - chain_lengths: list of residue counts for each chain.
    (No filtering for water residues.)
    """
    chain_order = []
    chain_lengths = []
    for chain in structure.get_chains():
        chain_order.append(chain.id)
        # Count every residue (do not skip water)
        count = sum(1 for residue in chain)
        chain_lengths.append(count)
    return chain_order, chain_lengths

###############################################################################
# 2) Read ColabFold JSON
###############################################################################
def read_colabfold_json(json_file):
    """
    Reads a typical ColabFold JSON that has keys:
      - 'pae'   => NxN array
      - 'ptm'   => float
      - 'iptm'  => float
      - 'plddt' => array of length N
    """
    with open(json_file, 'r') as f:
        data = json.load(f)

    pae = np.array(data.get('pae', []), dtype=float)
    pae = np.nan_to_num(pae)
    ptm = data.get('ptm', 0.0)
    iptm = data.get('iptm', 0.0)
    plddt = np.array(data.get('plddt', []), dtype=float)

    return {
        "pae": pae,
        "ptm": ptm,
        "iptm": iptm,
        "plddt": plddt,
    }

###############################################################################
# 3) Transform PAE => LIS
###############################################################################
def transform_pae_matrix(pae_matrix, pae_cutoff):
    """
    For each PAE value < cutoff, compute LIS = 1 - (PAE/cutoff); else 0.
    """
    transformed = np.zeros_like(pae_matrix)
    mask = (pae_matrix < pae_cutoff)
    transformed[mask] = 1.0 - (pae_matrix[mask] / pae_cutoff)
    return transformed

###############################################################################
# 4) Calculate contact map from structure
###############################################################################
def calculate_contact_map_from_structure(structure, distance_threshold=8.0):
    coords = []
    for chain in structure.get_chains():
        for residue in chain:
            if residue.has_id('CB'):
                coords.append(residue['CB'].coord)
            elif residue.has_id('CA'):
                coords.append(residue['CA'].coord)
    
    if len(coords) == 0:
        return np.zeros((0, 0), dtype=int)
    
    coords = np.array(coords)
    distmat = squareform(pdist(coords))
    contact_map = (distmat < distance_threshold).astype(int)
    return contact_map

###############################################################################
# 5) Calculate mean LIS (or cLIS) for each chain pair sub-block
###############################################################################
def calculate_mean_lis(transformed_map, subunit_number):
    """
    For each submatrix (chain_i x chain_j) in transformed_map,
    average the values > 0.
    """
    cum_lengths = np.cumsum(subunit_number)
    n_chains = len(subunit_number)
    mean_matrix = np.zeros((n_chains, n_chains))
    starts = np.concatenate(([0], cum_lengths[:-1]))
    for i in range(n_chains):
        for j in range(n_chains):
            si, ei = starts[i], cum_lengths[i]
            sj, ej = starts[j], cum_lengths[j]
            block = transformed_map[si:ei, sj:ej]
            positive = block[block > 0]
            mean_matrix[i, j] = positive.mean() if len(positive) > 0 else 0.0
    return mean_matrix

###############################################################################
# 6) Main Analysis Function: Build and return DataFrame
###############################################################################
def colabfold_lis_analysis_to_df(pdb_files, json_files, output_folder, pae_cutoff=12.0, distance_cutoff=8.0, result_save="True"):
    """
    For each (pdb, json) pair, compute metrics for every chain pair:
      - LIA: count of positions with transformed PAE > 0.
      - LIR: count of unique residue indices (rows + columns) in the submatrix.
      - cLIA and cLIR: same as above but for contact-based (cLIA) map.
      - LIS: mean of non-zero values in submatrix of transformed PAE.
      - cLIS: mean of submatrix from contact-filtered map.
      - iLIS = sqrt(LIS * cLIS)
      - iptm, ptm, and average plddt from JSON.
      - Also extracts 'rank' and 'model' from the pdb file name.
      - File name is extracted as the part before "_unrelaxed_".
      - For each chain pair, returns LIR_indices_A, LIR_indices_B, cLIR_indices_A, and cLIR_indices_B.
    """
    all_rows = []
    for pdb_file, json_file in zip(pdb_files, json_files):
        folder_name = os.path.basename(os.path.dirname(pdb_file))
        base_name = os.path.basename(pdb_file)

        # Extract file_name as the part before "_unrelaxed_"
        file_name = base_name.split("_unrelaxed_")[0] if "_unrelaxed_" in base_name else base_name

        # Extract rank
        rank = None
        if "_rank_" in base_name:
            try:
                rank = int(base_name.split("_rank_")[1].split("_")[0])
            except ValueError:
                rank = None

        # Extract model number
        model_val = None
        if "model_" in base_name:
            try:
                model_val = base_name.split("model_")[1].split("_")[0]
            except IndexError:
                model_val = None
        
        # Parse the PDB structure once
        structure = get_structure(pdb_file)
        chain_order, chain_lengths = parse_pdb_chains_from_structure(structure)
        total_res = sum(chain_lengths)
        
        # Read JSON for PAE, ptm, iptm, and plddt
        data = read_colabfold_json(json_file)
        pae_matrix = data["pae"]
        ptm = data["ptm"]
        iptm = data["iptm"]
        plddt = data["plddt"]
        confidence = 0.8 * iptm + 0.2 * ptm
        plddt_avg = float(plddt.mean()) if plddt.size > 0 else 0.0

        if pae_matrix.shape[0] != total_res:
            print(f"[Warning] PAE shape {pae_matrix.shape} != total residues {total_res} for {pdb_file}")
        
        # Transform PAE to get the LIS map (LIA map)
        transformed_pae = transform_pae_matrix(pae_matrix, pae_cutoff)
        lia_map = (transformed_pae > 0).astype(int)
        
        # Calculate contact map (from the structure) and then cLIA map
        contact_map = calculate_contact_map_from_structure(structure, distance_cutoff)
        clia_map = ((transformed_pae > 0) & (contact_map == 1)).astype(int)
        
        # Compute mean values for each chain pair submatrix
        mean_lis = calculate_mean_lis(transformed_pae, chain_lengths)
        mean_clis = calculate_mean_lis(np.where(clia_map > 0, transformed_pae, 0), chain_lengths)
        iLIS_matrix = np.sqrt(mean_lis * mean_clis)
        
        n_sub = len(chain_lengths)
        cum_lengths = np.cumsum(chain_lengths)
        starts = np.concatenate(([0], cum_lengths[:-1]))
        
        for i in range(n_sub):
            for j in range(n_sub):
                si, ei = starts[i], cum_lengths[i]
                sj, ej = starts[j], cum_lengths[j]
                sub_lia = lia_map[si:ei, sj:ej]
                sub_clia = clia_map[si:ei, sj:ej]
                
                # Count-based metrics
                LIA_val = np.count_nonzero(sub_lia)
                LIR_val = len(np.unique(np.where(sub_lia > 0)[0])) + len(np.unique(np.where(sub_lia > 0)[1]))
                cLIA_val = np.count_nonzero(sub_clia)
                cLIR_val = len(np.unique(np.where(sub_clia > 0)[0])) + len(np.unique(np.where(sub_clia > 0)[1]))
                
                # Local residue indices (1-based)
                LIR_indices_A = np.unique(np.where(sub_lia > 0)[0].flatten() + 1).tolist()
                LIR_indices_B = np.unique(np.where(sub_lia > 0)[1].flatten() + 1).tolist()
                cLIR_indices_A = np.unique(np.where(sub_clia > 0)[0].flatten() + 1).tolist()
                cLIR_indices_B = np.unique(np.where(sub_clia > 0)[1].flatten() + 1).tolist()

                LIS_val = mean_lis[i, j]
                cLIS_val = mean_clis[i, j]
                iLIS_val = np.sqrt(LIS_val * cLIS_val)
                
                row_dict = {
                    'folder_name': folder_name,
                    'file_name': file_name,
                    'chain_1': i + 1,
                    'chain_2': j + 1,
                    'rank': rank,
                    'model': model_val,
                    'iLIS': iLIS_val,
                    'LIS': LIS_val,
                    'cLIS': cLIS_val,
                    'iptm': iptm,
                    'confidence': confidence,
                    'LIA': LIA_val,
                    'LIR': LIR_val,
                    'cLIA': cLIA_val,
                    'cLIR': cLIR_val,
                    'LIR_indices_A': LIR_indices_A,
                    'LIR_indices_B': LIR_indices_B,
                    'cLIR_indices_A': cLIR_indices_A,
                    'cLIR_indices_B': cLIR_indices_B,
                    'ptm': ptm,
                    'plddt': plddt_avg
                }
                all_rows.append(row_dict)
    df_merged = pd.DataFrame(all_rows)

    # Optional grouping to average symmetric pairs (e.g. (1,2) and (2,1))
    def union_of_lists(series):
        combined = set()
        for lst in series.dropna():
            combined.update(lst)
        return sorted(combined)
    
    # Add a helper column 'chain_pair' (sorted tuple)
    df_merged['chain_pair'] = df_merged.apply(lambda row: tuple(sorted((row['chain_1'], row['chain_2']))), axis=1)
    numeric_cols = ['iLIS', 'LIS', 'cLIS', 'LIA', 'cLIA', 'LIR', 'cLIR', 'iptm', 'confidence', 'ptm', 'plddt']
    list_cols = ['LIR_indices_A', 'LIR_indices_B', 'cLIR_indices_A', 'cLIR_indices_B']
    agg_dict = {col: 'mean' for col in numeric_cols}
    agg_dict.update({col: union_of_lists for col in list_cols})
    for c in ['file_name', 'rank', 'model', 'folder_name']:
        agg_dict[c] = 'first'
    group_cols = ['chain_pair', 'file_name', 'rank', 'model', 'folder_name']
    
    # Perform the grouping operation
    df_grouped = df_merged.groupby(group_cols, as_index=False).agg(agg_dict)

    # Round values as needed
    df_grouped['LIS'] = df_grouped['LIS'].round(3)
    df_grouped['cLIS'] = df_grouped['cLIS'].round(3)
    df_grouped['iLIS'] = df_grouped['iLIS'].round(3)
    df_grouped['plddt'] = df_grouped['plddt'].round(2)
    df_grouped['chain_1'] = df_grouped['chain_pair'].apply(lambda x: x[0])
    df_grouped['chain_2'] = df_grouped['chain_pair'].apply(lambda x: x[1])
    
    final_cols = [
        'file_name', 'chain_1', 'chain_2', 'rank', 'model',
        'iLIS', 'LIS', 'cLIS', 'LIA', 'cLIA', 'LIR', 'cLIR',
        'LIR_indices_A', 'LIR_indices_B', 'cLIR_indices_A', 'cLIR_indices_B',
        'iptm', 'confidence', 'ptm', 'plddt', 'folder_name'
    ]
    df_final = df_grouped[final_cols].sort_values(['chain_1', 'chain_2'])

    if result_save == "True" and pdb_files:
        # Build output filename using the provided output_folder (without duplicating the folder path)
        output_filename = f"{file_name}_rank_00{rank}_model_{model_val}_lis_analysis.csv"
        output_path = os.path.join(output_folder, output_filename)
        df_final.to_csv(output_path, index=False)
        print("Results saved to:", output_path)

###############################################################################
# Main function: parse command-line arguments and run analysis
###############################################################################
def main():
    parser = argparse.ArgumentParser(
        description="ColabFold LIS Analysis: Compute metrics from a PDB file and its corresponding ColabFold JSON file."
    )
    parser.add_argument("pdb_file", help="Path to the PDB file")
    parser.add_argument("json_file", help="Path to the JSON file")
    parser.add_argument("-o", "--output_folder", default="lis_output", help="Output folder (default: lis_output)")
    parser.add_argument("-p", "--pae_cutoff", type=float, default=12.0, help="PAE cutoff (default: 12)")
    parser.add_argument("-d", "--distance_cutoff", type=float, default=8.0, help="Distance cutoff (default: 8)")
    args = parser.parse_args()

    # Ensure output folder exists
    if not os.path.isdir(args.output_folder):
        os.makedirs(args.output_folder)

    colabfold_lis_analysis_to_df(
        [args.pdb_file],
        [args.json_file],
        output_folder=args.output_folder,
        pae_cutoff=args.pae_cutoff,
        distance_cutoff=args.distance_cutoff,
        result_save="True"
    )

if __name__ == "__main__":
    main()
