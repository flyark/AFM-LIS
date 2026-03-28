#!/usr/bin/env python3
"""
seq_extract.py — Extract unique protein sequences from ColabFold PDB outputs
=============================================================================
Scans a folder of ColabFold predictions, extracts chain sequences from PDB
files, and outputs a deduplicated FASTA file.

Smart selection: only reads the minimum number of PDB files needed to cover
all unique proteins (based on filenames like PROT1___PROT2_unrelaxed_rank_*.pdb).

Usage:
    python seq_extract.py /path/to/predictions/ -o sequences.fasta
"""

import argparse
import gzip
import lzma
import os
import re
import sys
from collections import OrderedDict

# Standard 3-letter to 1-letter amino acid mapping
AA3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'SEC': 'U', 'PYL': 'O',
}


def read_file(path):
    """Read a file, auto-decompress .gz/.xz."""
    with open(path, 'rb') as f:
        data = f.read()
    if path.endswith('.gz'):
        data = gzip.decompress(data)
    elif path.endswith('.xz'):
        data = lzma.decompress(data)
    return data.decode('utf-8')


def extract_sequences_from_pdb(pdb_text):
    """Extract per-chain amino acid sequences from PDB ATOM records (CA only)."""
    chains = OrderedDict()
    seen = set()
    for line in pdb_text.split('\n'):
        if not line.startswith('ATOM'):
            continue
        atom_name = line[12:16].strip()
        if atom_name != 'CA':
            continue
        chain = line[21:22].strip() or 'A'
        resnum = int(line[22:26].strip())
        resname = line[17:20].strip()
        key = (chain, resnum)
        if key in seen:
            continue
        seen.add(key)
        aa = AA3TO1.get(resname, 'X')
        if chain not in chains:
            chains[chain] = []
        chains[chain].append((resnum, aa))

    # Sort by residue number and join
    result = {}
    for chain, residues in chains.items():
        residues.sort(key=lambda x: x[0])
        result[chain] = ''.join(r[1] for r in residues)
    return result


def parse_protein_names(filename):
    """Extract protein names from ColabFold filename.
    e.g., ABL1_HUMAN___KLF9_HUMAN_unrelaxed_rank_001_... -> ('ABL1_HUMAN', 'KLF9_HUMAN')
    """
    base = os.path.basename(filename)
    # Remove compression extensions
    for ext in ('.gz', '.xz'):
        if base.endswith(ext):
            base = base[:-len(ext)]

    # Match ColabFold pattern: PREFIX_unrelaxed_rank_*
    m = re.match(r'^(.+)_unrelaxed_rank_\d+', base)
    if not m:
        return None, None

    prefix = m.group(1)

    # Split by ___ (triple underscore) separator
    if '___' in prefix:
        parts = prefix.split('___', 1)
        return parts[0], parts[1]

    return None, None


def find_pdb_files(folder):
    """Find all PDB files in folder, return list of (path, protA, protB)."""
    pdbs = []
    for fname in os.listdir(folder):
        if not (fname.endswith('.pdb') or fname.endswith('.pdb.gz') or fname.endswith('.pdb.xz')):
            continue
        if fname.startswith('._') or '__MACOSX' in fname:
            continue

        prot_a, prot_b = parse_protein_names(fname)
        if prot_a and prot_b:
            pdbs.append((os.path.join(folder, fname), prot_a, prot_b))

    return pdbs


def smart_extract(folder, separator='___'):
    """Extract unique sequences with minimum PDB reads."""
    pdbs = find_pdb_files(folder)
    if not pdbs:
        print('[seq_extract] ERROR: No ColabFold PDB files found', file=sys.stderr)
        sys.exit(1)

    # Collect all unique protein names
    all_proteins = set()
    for _, prot_a, prot_b in pdbs:
        all_proteins.add(prot_a)
        all_proteins.add(prot_b)

    print(f'[seq_extract] {len(all_proteins)} unique proteins, {len(pdbs)} PDB files')

    # Smart selection: pick PDBs that cover the most missing proteins
    sequences = {}  # protein_name -> sequence
    remaining = set(all_proteins)
    files_read = 0

    # Sort PDBs to prioritize rank_001 (best model)
    pdbs.sort(key=lambda x: x[0])
    rank1_pdbs = [p for p in pdbs if '_rank_001_' in p[0] or '_rank_1_' in p[0]]
    other_pdbs = [p for p in pdbs if p not in rank1_pdbs]

    for pdb_path, prot_a, prot_b in rank1_pdbs + other_pdbs:
        if not remaining:
            break

        # Skip if both proteins already have sequences
        need_a = prot_a in remaining
        need_b = prot_b in remaining
        if not need_a and not need_b:
            continue

        # Read PDB and extract sequences
        try:
            pdb_text = read_file(pdb_path)
            chain_seqs = extract_sequences_from_pdb(pdb_text)
            files_read += 1
        except Exception as e:
            print(f'[seq_extract] WARNING: {os.path.basename(pdb_path)}: {e}', file=sys.stderr)
            continue

        chains = sorted(chain_seqs.keys())
        if len(chains) >= 2:
            if need_a and chains[0] in chain_seqs:
                sequences[prot_a] = chain_seqs[chains[0]]
                remaining.discard(prot_a)
            if need_b and chains[1] in chain_seqs:
                sequences[prot_b] = chain_seqs[chains[1]]
                remaining.discard(prot_b)
        elif len(chains) == 1:
            # Single chain — assign to whichever is needed
            if need_a:
                sequences[prot_a] = chain_seqs[chains[0]]
                remaining.discard(prot_a)

        if files_read % 10 == 0 or not remaining:
            print(f'\r[seq_extract] Read {files_read} PDBs, '
                  f'{len(sequences)}/{len(all_proteins)} proteins covered',
                  end='', flush=True)

    print(f'\r[seq_extract] Read {files_read} PDBs, '
          f'{len(sequences)}/{len(all_proteins)} proteins covered')

    if remaining:
        print(f'[seq_extract] WARNING: {len(remaining)} proteins not found: '
              f'{", ".join(sorted(remaining)[:5])}{"..." if len(remaining) > 5 else ""}',
              file=sys.stderr)

    return sequences


def main():
    parser = argparse.ArgumentParser(
        description='Extract unique protein sequences from ColabFold PDB outputs',
    )
    parser.add_argument('path', help='Folder with ColabFold prediction PDB files')
    parser.add_argument('--output', '-o', default=None,
                        help='Output FASTA file (default: <folder>_sequences.fasta)')
    parser.add_argument('--separator', default='___',
                        help='Separator between protein names in filename (default: ___)')
    args = parser.parse_args()

    if not os.path.isdir(args.path):
        print(f'[seq_extract] ERROR: {args.path} is not a directory', file=sys.stderr)
        sys.exit(1)

    # Determine output path early for skip check
    output = args.output
    if output is None:
        basename = os.path.basename(args.path.rstrip('/'))
        output = os.path.join(args.path, f'{basename}_sequences.fasta')

    # Skip if already exists and covers all proteins
    if os.path.exists(output):
        pdbs = find_pdb_files(args.path)
        all_proteins = set()
        for _, a, b in pdbs:
            all_proteins.add(a)
            all_proteins.add(b)
        # Count sequences in existing FASTA
        existing_count = 0
        with open(output, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    existing_count += 1
        if existing_count >= len(all_proteins):
            print(f'[seq_extract] Up to date: {output} ({existing_count} sequences, {len(all_proteins)} proteins)')
            return
        else:
            print(f'[seq_extract] Updating: {existing_count} sequences but {len(all_proteins)} proteins found')

    sequences = smart_extract(args.path, args.separator)

    # Verify: unique sequences match expected count
    unique_seqs = set(sequences.values())
    n_names = len(sequences)
    n_unique = len(unique_seqs)
    if n_names != n_unique:
        # Find which proteins share the same sequence
        seq_to_names = {}
        for name, seq in sequences.items():
            seq_to_names.setdefault(seq, []).append(name)
        duplicates = {tuple(names): len(seq) for seq, names in seq_to_names.items() if len(names) > 1}
        print(f'[seq_extract] NOTE: {n_names} protein names but {n_unique} unique sequences')
        for names, seqlen in sorted(duplicates.items()):
            print(f'  Same sequence ({seqlen} aa): {", ".join(names)}')
    else:
        print(f'[seq_extract] Verified: {n_names} names = {n_unique} unique sequences')

    # Write FASTA sorted by protein name
    with open(output, 'w') as f:
        for name in sorted(sequences.keys()):
            seq = sequences[name]
            f.write(f'>{name}\n')
            # Wrap at 80 characters
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')

    print(f'[seq_extract] Saved {n_names} sequences to {output}')


if __name__ == '__main__':
    main()
