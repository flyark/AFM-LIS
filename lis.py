#!/usr/bin/env python3
"""
LIS -- Local Interaction Score Analysis (CLI)
==============================================
Calculate LIS/cLIS/iLIS metrics from structure prediction outputs.

Supports: AlphaFold3, ColabFold, Boltz, Chai-1, OpenFold3, Generic
Input:    folder or zip file with prediction outputs
Output:   CSV with one row per model per chain pair

Architecture: scan folder, process each structure+PAE pair independently,
append results to CSV one model at a time. No grouping, no averaging.

Usage:
    python lis.py <path> [options]
      path              folder or zip file with prediction outputs
      --output / -o     output CSV filename
      --output-dir / -d output directory (default: input folder)
      --pae-cutoff      PAE cutoff (default: 12)
      --cb-cutoff       Cb distance cutoff in Angstroms (default: 8)
      --workers / -w    parallel workers (default: 1)
      --no-skip-existing  reprocess everything

Dependencies: numpy, scipy
"""

import argparse
import gzip
import io
import json
import lzma
import math
import os
import re
import shutil
import sys
import tarfile
import tempfile
import zipfile
from collections import OrderedDict

import numpy as np
from scipy.spatial.distance import pdist, squareform

# ============================================================================
# Constants
# ============================================================================

ION_NAMES = {'ZN', 'CA', 'MG', 'MN', 'FE', 'CU', 'NA', 'K', 'CL', 'NI', 'CO', 'CD'}

CSV_HEADER = (
    'name,rank,model,chain_i,chain_j,iLIS,iLIA,iLISA,ipSAE,actifpTM,LIS,cLIS,LIA,cLIA,'
    'ipTM,pLDDT_i,pLDDT_j,pTM,LIR_i,LIR_j,cLIR_i,cLIR_j,'
    'LIpLDDT_i,LIpLDDT_j,cLIpLDDT_i,cLIpLDDT_j,'
    'len_i,len_j,LIR_indices_i,LIR_indices_j,cLIR_indices_i,cLIR_indices_j,'
    'structure_file'
)


# ============================================================================
# File I/O -- read files from folders, zips, with .gz/.xz decompression
# ============================================================================

def _strip_compression_ext(name):
    """Strip .gz / .xz extension from a filename."""
    if name.endswith('.gz'):
        return name[:-3]
    if name.endswith('.xz'):
        return name[:-3]
    return name


def _decode_content(name, data):
    """Decompress if needed and decode text files to str, leave binary as bytes."""
    basename = os.path.basename(name)
    if basename.endswith('.gz'):
        data = gzip.decompress(data)
        basename = basename[:-3]
    if basename.endswith('.xz'):
        data = lzma.decompress(data)
        basename = basename[:-3]
    if any(basename.endswith(ext) for ext in ('.json', '.pdb', '.cif', '.txt', '.csv')):
        try:
            return data.decode('utf-8')
        except UnicodeDecodeError:
            return data
    return data


def scan_files(path):
    """Scan a folder or zip for filenames without reading file contents.

    Returns (filenames, read_fn) where:
        filenames: list of relative file paths (with .gz/.xz extensions stripped)
        read_fn(name): reads and returns the content of a single file on demand
    """
    path = str(path)

    if _is_tar_zstd(path):
        file_map = {}
        _scan_tar_zstd(path, file_map)
        filenames = [k for k in file_map.keys() if k != '__tmpdir__']
        read_fn = _make_tar_zstd_reader(path, file_map)
    elif zipfile.is_zipfile(path):
        file_map = {}
        _scan_zip(path, file_map, prefix='')
        filenames = list(file_map.keys())
        read_fn = _make_zip_reader(file_map)
    elif os.path.isdir(path):
        file_map = {}
        _scan_dir(path, file_map)
        filenames = list(file_map.keys())
        read_fn = _make_dir_reader(file_map)
    else:
        print(f"[LIS] ERROR: {path} is not a valid folder, zip, or tar.zstd file", file=sys.stderr)
        sys.exit(1)

    return filenames, read_fn, file_map


def _scan_zip(zip_path, file_map, prefix='', nesting_chain=None):
    """Recursively scan a zip file (handles nested zips)."""
    if nesting_chain is None:
        nesting_chain = []
    with zipfile.ZipFile(zip_path) as zf:
        for name in zf.namelist():
            if name.endswith('/'):
                continue
            full_name = prefix + name
            basename = os.path.basename(name)
            if basename.endswith('.zip'):
                data = zf.read(name)
                with tempfile.NamedTemporaryFile(suffix='.zip', delete=False) as tmp:
                    tmp.write(data)
                    tmp_path = tmp.name
                try:
                    if nesting_chain:
                        new_chain = nesting_chain + [(name,)]
                    else:
                        new_chain = [(zip_path, name)]
                    _scan_zip(tmp_path, file_map,
                              prefix=full_name.replace('.zip', '/'),
                              nesting_chain=new_chain)
                finally:
                    os.unlink(tmp_path)
                continue
            clean_name = _strip_compression_ext(full_name)
            if nesting_chain:
                file_map[clean_name] = nesting_chain + [(name,)]
            else:
                file_map[clean_name] = (zip_path, name)


def _is_tar_zstd(path):
    """Check if path is a tar.zstd or tar.zst file."""
    return os.path.isfile(path) and any(
        path.endswith(ext) for ext in ('.tar.zstd', '.tar.zst', '.tar.zstandard')
    )


def _extract_tar_zstd_to_tmpdir(tar_path):
    """Extract a zstd-compressed tar to a temp directory. Returns the temp dir path."""
    import subprocess
    tmpdir = tempfile.mkdtemp(prefix='lis_tar_')
    print(f"[LIS] Extracting {os.path.basename(tar_path)} to temp dir...", file=sys.stderr)
    proc2 = subprocess.Popen(['zstd', '-d', '-c', tar_path], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    proc3 = subprocess.Popen(['tar', 'xf', '-', '-C', tmpdir], stdin=proc2.stdout, stderr=subprocess.DEVNULL)
    proc2.stdout.close()
    proc3.communicate()
    print(f"[LIS] Extraction complete.", file=sys.stderr)
    return tmpdir


def _scan_tar_zstd(tar_path, file_map):
    """Extract tar.zstd to temp dir and scan as a regular directory."""
    tmpdir = _extract_tar_zstd_to_tmpdir(tar_path)
    # Store tmpdir path so reader can find it
    file_map['__tmpdir__'] = tmpdir
    _scan_dir(tmpdir, file_map)


def _make_tar_zstd_reader(tar_path, file_map):
    """Create a reader using the extracted temp directory."""
    return _make_dir_reader(file_map)


def _scan_dir(dirpath, file_map):
    """Recursively scan a directory (handles nested zips)."""
    for root, _dirs, fnames in os.walk(dirpath):
        for fname in fnames:
            fpath = os.path.join(root, fname)
            relname = os.path.relpath(fpath, dirpath)
            if fname.endswith('.zip') and zipfile.is_zipfile(fpath):
                _scan_zip(fpath, file_map, prefix=relname.replace('.zip', '/'))
                continue
            clean_name = _strip_compression_ext(relname)
            file_map[clean_name] = fpath


def _make_dir_reader(file_map):
    """Create a reader function for directory-based file_map."""
    def read_fn(name):
        fpath = file_map.get(name)
        if fpath is None:
            return None
        with open(fpath, 'rb') as f:
            data = f.read()
        return _decode_content(fpath, data)
    return read_fn


def _make_zip_reader(file_map):
    """Create a reader function for zip-based file_map."""
    def read_fn(name):
        entry = file_map.get(name)
        if entry is None:
            return None
        if isinstance(entry, tuple) and len(entry) == 2:
            zip_path, raw_name = entry
            with zipfile.ZipFile(zip_path) as zf:
                data = zf.read(raw_name)
            return _decode_content(raw_name, data)
        # Nested zip case
        chain = entry
        outer_zip_path, first_nested_entry = chain[0]
        remaining = chain[1:]
        with zipfile.ZipFile(outer_zip_path) as zf:
            data = zf.read(first_nested_entry)
        tmp_paths = []
        try:
            for step in remaining:
                entry_name = step[0]
                with tempfile.NamedTemporaryFile(suffix='.zip', delete=False) as tmp:
                    tmp.write(data)
                    tmp_paths.append(tmp.name)
                with zipfile.ZipFile(tmp_paths[-1]) as zf:
                    data = zf.read(entry_name)
        finally:
            for tp in tmp_paths:
                try:
                    os.unlink(tp)
                except OSError:
                    pass
        final_name = remaining[-1][0] if remaining else first_nested_entry
        return _decode_content(final_name, data)
    return read_fn


# ============================================================================
# Platform Detection
# ============================================================================

def detect_platform(filenames, read_fn):
    """Auto-detect the prediction platform from filenames."""
    basenames = [os.path.basename(f) for f in filenames]

    # Tamarind AF3 / OpenFold3: result_sample_N_model.pdb + result_sample_N_confidences.json
    has_tamarind_model = any(re.match(r'^result_sample_\d+_model\.pdb$', b) for b in basenames)
    has_tamarind_conf = any(re.match(r'^result_sample_\d+_confidences\.json$', b) for b in basenames)
    if has_tamarind_model and has_tamarind_conf:
        for fname in filenames:
            if os.path.basename(fname) == 'experiment_config.json':
                content = read_fn(fname)
                if isinstance(content, str):
                    try:
                        ec = json.loads(content)
                        if ec.get('inference_ckpt_path', '') and 'of3' in ec['inference_ckpt_path']:
                            return 'openfold3'
                    except (json.JSONDecodeError, KeyError):
                        pass
        return 'alphafold3'

    # Standard AF3: *_model_N.cif + *_summary_confidences_N.json + *_full_data_N.json
    has_af3_model = any(re.search(r'_model_\d+\.cif$', b) for b in basenames)
    has_af3_summary = any(re.search(r'_summary_confidences_\d+\.json$', b) for b in basenames)
    has_af3_full = any(re.search(r'_full_data_\d+\.json$', b) for b in basenames)
    if has_af3_model and has_af3_summary and has_af3_full:
        return 'alphafold3'

    # AF3 Server output: seed-N_sample-M/model.cif + seed-N_sample-M/confidences.json
    has_server_model = any(re.search(r'seed-\d+_sample-\d+/model\.cif$', f) for f in filenames)
    has_server_conf = any(re.search(r'seed-\d+_sample-\d+/confidences\.json$', f) for f in filenames)
    if has_server_model and has_server_conf:
        return 'alphafold3'

    # ColabFold
    has_cf_model = any(re.search(r'_unrelaxed_rank_\d+', b) and (b.endswith('.pdb') or b.endswith('.cif'))
                       for b in basenames)
    has_cf_scores = any(re.search(r'_scores_rank_\d+', b) and b.endswith('.json') for b in basenames)
    if has_cf_model and has_cf_scores:
        return 'colabfold'

    # Boltz
    has_boltz_conf = any(b.startswith('confidence_') and b.endswith('.json') for b in basenames)
    has_boltz_struct = any(b.endswith('.cif') or b.endswith('.pdb') for b in basenames)
    if has_boltz_conf and has_boltz_struct:
        return 'boltz'

    # Chai-1
    has_chai_model = any(re.match(r'^pred\.(rank_\d+|model_idx_\d+)\.cif$', b) for b in basenames)
    has_chai_scores = any(re.match(r'^scores\.(rank_\d+|model_idx_\d+)\.(json|npz)$', b) for b in basenames)
    if has_chai_model or (has_chai_scores and any(b.startswith('pred.') for b in basenames)):
        return 'chai'

    return 'generic'


# ============================================================================
# Model Discovery -- one flat list of (name, rank, model, struct, pae, scores, fmt)
# ============================================================================

def _get_toplevel_name(filepath):
    """Extract top-level directory name from a relative path as the prediction name.
    e.g. 'result-p53-mdm2_chai-1/pred.model_idx_1.cif' → 'result-p53-mdm2_chai-1'
    e.g. 'pred.model_idx_1.cif' → None (no parent directory)
    """
    parts = filepath.split('/')
    if len(parts) >= 2 and parts[0] and not parts[0].startswith('__'):
        return parts[0]
    return None


def find_models(filenames, platform, read_fn):
    """Scan filenames and yield model tuples.

    Yields: (name, rank, model, struct_path, pae_path, scores_path, fmt)
    Each tuple represents one independent structure+PAE pair to analyze.
    """
    # Filter out macOS resource fork files
    filenames = [f for f in filenames if '__MACOSX' not in f and not os.path.basename(f).startswith('._')]
    basenames_map = {os.path.basename(f): f for f in filenames}

    if platform == 'colabfold':
        yield from _find_colabfold(filenames, basenames_map)
    elif platform in ('alphafold3', 'openfold3'):
        yield from _find_af3(filenames, basenames_map)
    elif platform == 'boltz':
        yield from _find_boltz(filenames, basenames_map)
    elif platform == 'chai':
        yield from _find_chai(filenames, basenames_map)
    else:
        yield from _find_generic(filenames, basenames_map)


def _find_colabfold(filenames, basenames_map):
    """ColabFold: *_unrelaxed_rank_N_model_M.pdb + *_scores_rank_N_model_M.json"""
    for name in filenames:
        base = os.path.basename(name)
        m = re.match(r'^(.+)_unrelaxed_rank_(\d+)(.*)\.(pdb|cif)$', base)
        if not m:
            continue
        prefix = m.group(1)
        rank = str(int(m.group(2)))  # strip leading zeros: 001 -> 1
        rest = m.group(3)
        fmt = m.group(4)

        # Use the full PDB filename as the model identifier
        model_str = base

        # Find matching scores file (must match same prefix)
        scores_path = None
        padded = rank.zfill(3)
        for fn in filenames:
            fb = os.path.basename(fn)
            if fb.startswith(prefix) and '_scores_rank_' in fb and f'rank_{padded}' in fb and fb.endswith('.json'):
                scores_path = fn
                break
        if not scores_path:
            for fn in filenames:
                fb = os.path.basename(fn)
                if fb.startswith(prefix) and 'scores' in fb and f'rank_{rank}' in fb and fb.endswith('.json'):
                    scores_path = fn
                    break
        if not scores_path:
            for fn in filenames:
                fb = os.path.basename(fn)
                if fb.startswith(prefix) and 'scores' in fb and fb.endswith('.json'):
                    scores_path = fn
                    break

        # ColabFold: PAE is inside the scores JSON under 'pae' key
        yield (prefix, rank, model_str, name, scores_path, scores_path, fmt)


def _find_af3(filenames, basenames_map):
    """AlphaFold3 standard: *_model_N.cif + *_full_data_N.json + *_summary_confidences_N.json
    Tamarind AF3: result_sample_N_model.pdb + result_sample_N_confidences.json
    """
    # Standard AF3
    for name in filenames:
        base = os.path.basename(name)
        m = re.match(r'^(.+)_model_(\d+)\.cif$', base)
        if not m:
            continue
        prefix = m.group(1)
        idx = m.group(2)
        full_data = basenames_map.get(f'{prefix}_full_data_{idx}.json')
        summary = basenames_map.get(f'{prefix}_summary_confidences_{idx}.json')
        if not full_data or not summary:
            continue
        yield (prefix, idx, os.path.basename(name), name, full_data, summary, 'cif')

    # Tamarind AF3
    for name in filenames:
        base = os.path.basename(name)
        m = re.match(r'^result_sample_(\d+)_model\.pdb$', base)
        if not m:
            continue
        idx = m.group(1)
        conf_path = basenames_map.get(f'result_sample_{idx}_confidences.json')
        agg_path = basenames_map.get(f'result_sample_{idx}_confidences_aggregated.json')
        pred_name = _get_toplevel_name(name) or 'prediction'
        yield (pred_name, idx, base, name, conf_path, agg_path or conf_path, 'pdb')

    # AF3 Server output: prediction_name/seed-N_sample-M/model.cif + confidences.json
    filenames_set = set(filenames)
    for name in filenames:
        base = os.path.basename(name)
        if base != 'model.cif':
            continue
        dirpath = os.path.dirname(name)
        sample_dir = os.path.basename(dirpath)
        m = re.match(r'seed-(\d+)_sample-(\d+)', sample_dir)
        if not m:
            continue
        seed_idx = m.group(1)
        sample_idx = m.group(2)
        # Prediction name from grandparent directory
        pred_dir = os.path.dirname(dirpath)
        pred_name = os.path.basename(pred_dir)
        if not pred_name or pred_name in ('AF3_outputs', ''):
            pred_name = 'af3_prediction'
        # Find confidences.json in same directory via direct path lookup
        conf_path = os.path.join(dirpath, 'confidences.json')
        if conf_path not in filenames_set:
            conf_path = None
        summary_path = os.path.join(dirpath, 'summary_confidences.json')
        if summary_path not in filenames_set:
            summary_path = None
        rank = f'{seed_idx}_{sample_idx}'
        if conf_path:
            yield (pred_name, rank, base, name, conf_path, summary_path or conf_path, 'cif')


def _find_boltz(filenames, basenames_map):
    """Boltz: *_model_N.pdb/.cif + confidence_*_model_N.json + pae_*_model_N.npz

    Groups files by parent directory to handle multiple predictions.
    """
    # Group files by parent directory
    from collections import defaultdict
    dir_groups = defaultdict(list)
    for f in filenames:
        parent = os.path.dirname(f)
        dir_groups[parent].append(f)

    # If all files in one directory, treat as single prediction
    if len(dir_groups) <= 1:
        dir_groups = {'': filenames}

    for dirpath, dir_files in sorted(dir_groups.items()):
        struct_files = sorted([f for f in dir_files if f.endswith('.cif') or f.endswith('.pdb')],
                              key=lambda f: os.path.basename(f))
        conf_files = [f for f in dir_files if os.path.basename(f).startswith('confidence') and f.endswith('.json')]
        pae_files = [f for f in dir_files if os.path.basename(f).startswith('pae') and f.endswith('.npz')]

        if not struct_files:
            continue

        # Extract prediction name from directory or top-level folder
        pred_name = os.path.basename(dirpath) if dirpath else ''
        # Use top-level directory name if current name is generic
        if not pred_name or pred_name in ('Boltz1_outputs', 'boltz_outputs', 'predictions', 'result', 'results'):
            toplevel = _get_toplevel_name(struct_files[0]) if struct_files else None
            pred_name = toplevel or 'boltz'

        for i, sf in enumerate(struct_files):
            fmt = 'cif' if sf.endswith('.cif') else 'pdb'
            # Extract model index from filename
            sb = os.path.basename(sf)
            model_idx = str(i)
            m = re.search(r'model_(\d+)', sb)
            if m:
                model_idx = m.group(1)

            # Match confidence and pae by model index
            conf_path = None
            for cf in conf_files:
                if f'model_{model_idx}' in os.path.basename(cf):
                    conf_path = cf
                    break
            if not conf_path and i < len(conf_files):
                conf_path = conf_files[i]

            pae_path = None
            for pf in pae_files:
                if f'model_{model_idx}' in os.path.basename(pf):
                    pae_path = pf
                    break
            if not pae_path and i < len(pae_files):
                pae_path = pae_files[i]

            yield (pred_name, model_idx, sb, sf, pae_path or conf_path, conf_path, fmt)


def _find_chai(filenames, basenames_map):
    """Chai-1: pred.rank_N.cif + scores.rank_N.json + pae.rank_N.npy/.npz"""
    for name in filenames:
        base = os.path.basename(name)
        m = re.match(r'^pred\.(rank_(\d+)|model_idx_(\d+))\.cif$', base)
        if not m:
            continue
        rank_or_idx = m.group(2) or m.group(3)
        is_rank = m.group(2) is not None

        score_path = (
            basenames_map.get(f'scores.rank_{rank_or_idx}.json')
            or basenames_map.get(f'scores.model_idx_{rank_or_idx}.json')
        )
        if not score_path:
            for fn in filenames:
                if 'scores' in os.path.basename(fn) and fn.endswith('.json'):
                    score_path = fn
                    break

        pae_path = (
            basenames_map.get(f'pae.rank_{rank_or_idx}.npy')
            or basenames_map.get(f'pae.rank_{rank_or_idx}.npz')
            or basenames_map.get(f'pae.model_idx_{rank_or_idx}.npy')
            or basenames_map.get(f'pae.model_idx_{rank_or_idx}.npz')
        )

        pred_name = _get_toplevel_name(name) or 'prediction'
        yield (pred_name, rank_or_idx, base, name, pae_path or score_path, score_path, 'cif')


def _find_generic(filenames, basenames_map):
    """Generic: any .pdb/.cif + matching .json or .npz"""
    struct_files = sorted(
        [f for f in filenames if f.endswith('.cif') or f.endswith('.pdb')],
        key=lambda f: os.path.basename(f)
    )
    json_files = [f for f in filenames if f.endswith('.json')]
    npz_files = [f for f in filenames if f.endswith('.npz')]

    for i, sf in enumerate(struct_files):
        struct_base = os.path.basename(sf)
        fmt = 'cif' if sf.endswith('.cif') else 'pdb'
        name_prefix = re.sub(r'\.(cif|pdb)$', '', struct_base)

        pae_source = None
        idx_match = re.search(r'(\d+)', name_prefix)
        idx = idx_match.group(1) if idx_match else None

        for jf in json_files:
            if name_prefix in os.path.basename(jf):
                pae_source = jf
                break
        if not pae_source and idx:
            for jf in json_files:
                fb = os.path.basename(jf)
                if idx in fb and 'aggregated' not in fb:
                    pae_source = jf
                    break
        if not pae_source:
            pae_source = json_files[i] if i < len(json_files) else (json_files[0] if json_files else None)
        if not pae_source and npz_files:
            pae_source = npz_files[i] if i < len(npz_files) else npz_files[0]

        yield ('generic', str(i), struct_base, sf, pae_source, pae_source, fmt)


# ============================================================================
# PAE Extraction
# ============================================================================

def extract_pae(pae_source, read_fn):
    """Extract PAE matrix from a file. Returns 2D numpy array or None."""
    if not pae_source:
        return None
    content = read_fn(pae_source)
    if content is None:
        return None

    basename = os.path.basename(pae_source)

    # .npy file (Chai-1)
    if basename.endswith('.npy'):
        if isinstance(content, str):
            content = content.encode('utf-8')
        arr = np.load(io.BytesIO(content))
        if arr.ndim == 3:
            arr = arr[0]
        return arr.astype(np.float32)

    # .npz file (Boltz, Chai-1)
    if basename.endswith('.npz'):
        if isinstance(content, str):
            content = content.encode('utf-8')
        npz = np.load(io.BytesIO(content))
        if 'pae' in npz:
            arr = npz['pae']
        else:
            arr = None
            for key in npz.files:
                v = npz[key]
                if v.ndim >= 2:
                    arr = v
                    break
        if arr is None:
            return None
        if arr.ndim == 3:
            arr = arr[0]
        return arr.astype(np.float32)

    # JSON file
    if not isinstance(content, str):
        return None
    try:
        data = json.loads(content)
    except json.JSONDecodeError:
        return None

    # AF3 full_data: {pae: [[...]]}
    if 'pae' in data and isinstance(data['pae'], list):
        pae = data['pae']
        if isinstance(pae[0], list):
            return np.array(pae, dtype=np.float32)
        n = int(round(math.sqrt(len(pae))))
        if n * n == len(pae):
            return np.array(pae, dtype=np.float32).reshape(n, n)

    # AlphaFold DB format: [{predicted_aligned_error: [[...]]}]
    if isinstance(data, list) and data and 'predicted_aligned_error' in (data[0] if isinstance(data[0], dict) else {}):
        return np.array(data[0]['predicted_aligned_error'], dtype=np.float32)

    # Direct predicted_aligned_error field
    if 'predicted_aligned_error' in data:
        pae = data['predicted_aligned_error']
        if isinstance(pae, list) and isinstance(pae[0], list):
            return np.array(pae, dtype=np.float32)

    # Boltz confidence JSON pae_matrix
    for key in ('pae_matrix', 'predicted_aligned_error_matrix'):
        if key in data:
            return np.array(data[key], dtype=np.float32)

    # Tamarind AF3 / OpenFold3: pde key
    if 'pde' in data and isinstance(data['pde'], list):
        pde = data['pde']
        if isinstance(pde[0], list):
            return np.array(pde, dtype=np.float32)
        n = int(round(math.sqrt(len(pde))))
        if n * n == len(pde):
            return np.array(pde, dtype=np.float32).reshape(n, n)

    return None


# ============================================================================
# Confidence Score Extraction
# ============================================================================

def _unwrap(v):
    """Unwrap single-element arrays (Chai-1 wraps scalars in arrays)."""
    if isinstance(v, (list, np.ndarray)) and len(v) == 1:
        return float(v[0])
    return v


def extract_confidence_scores(confidence_path, read_fn):
    """Extract pTM, ipTM, chain_pair_iptm, pLDDT, etc. from confidence JSON."""
    if not confidence_path:
        return {}
    content = read_fn(confidence_path)
    if not content or not isinstance(content, str):
        return {}
    try:
        data = json.loads(content)
    except json.JSONDecodeError:
        return {}

    scores = {}

    # AF3 summary_confidences / ColabFold / Chai-1
    if 'ptm' in data:
        scores['pTM'] = _unwrap(data['ptm'])
    if 'iptm' in data:
        scores['ipTM'] = _unwrap(data['iptm'])
    if 'chain_pair_iptm' in data:
        scores['chainPairIptm'] = data['chain_pair_iptm']
    if 'atom_plddts' in data:
        plddts = data['atom_plddts']
        scores['pLDDT'] = sum(plddts) / len(plddts) if plddts else 0

    # ColabFold plddt array
    if 'plddt' in data:
        plddts = data['plddt']
        if isinstance(plddts, list) and plddts:
            scores.setdefault('pLDDT', sum(plddts) / len(plddts))

    # Chai-1 per_chain_pair_iptm
    if 'per_chain_pair_iptm' in data and 'chainPairIptm' not in scores:
        raw = data['per_chain_pair_iptm']
        if (isinstance(raw, list) and raw
                and isinstance(raw[0], list) and raw[0]
                and isinstance(raw[0][0], list)):
            scores['chainPairIptm'] = raw[0]
        else:
            scores['chainPairIptm'] = raw
    if 'aggregate_score' in data:
        scores['aggregateScore'] = _unwrap(data['aggregate_score'])

    # Tamarind AF3 aggregated format
    if 'avg_plddt' in data and 'pLDDT' not in scores:
        scores['pLDDT'] = data['avg_plddt']
    if 'sample_ranking_score' in data:
        scores['rankingScore'] = data['sample_ranking_score']
    if 'iptm_by_asym_id_pair' in data and 'chainPairIptm' not in scores:
        raw = data['iptm_by_asym_id_pair']
        all_ids = set()
        for key in raw:
            m2 = re.match(r'\((\d+),\s*(\d+)\)', key)
            if m2:
                all_ids.add(m2.group(1))
                all_ids.add(m2.group(2))
        ids = sorted(all_ids, key=int)
        if ids:
            matrix = []
            for i in ids:
                row = []
                for j in ids:
                    row.append(raw.get(f'({i}, {j})', raw.get(f'({j}, {i})', 0)))
                matrix.append(row)
            scores['chainPairIptm'] = matrix

    # Boltz confidence
    if 'confidence_score' in data:
        scores['confidence'] = _unwrap(data['confidence_score'])
    if 'ptm_score' in data:
        scores['pTM'] = _unwrap(data['ptm_score'])
    if 'iptm_score' in data:
        scores['ipTM'] = _unwrap(data['iptm_score'])
    if 'pair_chains_iptm' in data:
        raw = data['pair_chains_iptm']
        if isinstance(raw, list):
            scores['chainPairIptm'] = raw
        elif isinstance(raw, dict):
            keys = sorted(raw.keys(), key=lambda x: int(x))
            matrix = []
            for row_key in keys:
                row = raw[row_key]
                col_keys = sorted(row.keys(), key=lambda x: int(x))
                matrix.append([row[ck] for ck in col_keys])
            scores['chainPairIptm'] = matrix
    if 'complex_plddt' in data:
        scores.setdefault('pLDDT', data['complex_plddt'])

    return scores


# ============================================================================
# Structure Parsing -- PDB and mmCIF
# ============================================================================

def parse_pdb_coords(pdb_text):
    """Extract one Cb (Ca for GLY, P for nucleic) coordinate per residue from PDB text."""
    residues = OrderedDict()
    for line in pdb_text.split('\n'):
        if not line.startswith('ATOM') and not line.startswith('HETATM'):
            continue
        if len(line) < 54:
            continue
        atom_name = line[12:16].strip()
        comp_id = line[17:20].strip()
        chain = line[21:22].strip() or 'A'
        try:
            resnum = int(line[22:26].strip())
        except ValueError:
            continue
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        key = f'{chain}:{resnum}'

        if atom_name == 'CB':
            residues[key] = {'chain': chain, 'resnum': resnum, 'x': x, 'y': y, 'z': z, 'has_p': False}
        elif atom_name == 'CA' and comp_id == 'GLY' and key not in residues:
            residues[key] = {'chain': chain, 'resnum': resnum, 'x': x, 'y': y, 'z': z, 'has_p': False}
        elif atom_name == 'P' and key not in residues:
            residues[key] = {'chain': chain, 'resnum': resnum, 'x': x, 'y': y, 'z': z, 'has_p': True}
        elif line.startswith('HETATM') and comp_id in ION_NAMES and f'{chain}:1' not in residues:
            residues[f'{chain}:1'] = {'chain': chain, 'resnum': 1, 'x': x, 'y': y, 'z': z, 'has_p': False}

    return list(residues.values())


def parse_cif_coords(cif_text):
    """Extract one Cb (Ca for GLY, P for nucleic) coordinate per residue from mmCIF text."""
    residues = OrderedDict()
    in_atom_site = False
    col_names = []

    for line in cif_text.split('\n'):
        if line.startswith('_atom_site.'):
            in_atom_site = True
            col_names.append(line.strip().split('.')[1])
            continue
        if in_atom_site and not line.startswith('_atom_site.') and not line.startswith('#') and line.strip():
            if line.startswith('loop_') or line.startswith('_'):
                in_atom_site = False
                continue
            parts = line.strip().split()
            if len(parts) < len(col_names):
                continue

            def get_col(name, _parts=parts, _cols=col_names):
                idx = _cols.index(name) if name in _cols else -1
                return _parts[idx] if idx >= 0 else ''

            group_pdb = get_col('group_PDB')
            if group_pdb not in ('ATOM', 'HETATM'):
                continue

            atom_name = get_col('label_atom_id')
            comp_id = get_col('label_comp_id')
            chain = get_col('label_asym_id')
            res_seq = get_col('label_seq_id')

            try:
                x = float(get_col('Cartn_x'))
                y = float(get_col('Cartn_y'))
                z = float(get_col('Cartn_z'))
            except ValueError:
                continue

            if group_pdb == 'HETATM' and comp_id in ION_NAMES:
                ion_key = f'{chain}:1'
                if ion_key not in residues:
                    residues[ion_key] = {'chain': chain, 'resnum': 1, 'x': x, 'y': y, 'z': z, 'has_p': False}
                continue

            try:
                resnum = int(res_seq)
            except ValueError:
                continue
            key = f'{chain}:{resnum}'

            if atom_name == 'CB':
                residues[key] = {'chain': chain, 'resnum': resnum, 'x': x, 'y': y, 'z': z, 'has_p': False}
            elif atom_name == 'CA' and comp_id == 'GLY' and key not in residues:
                residues[key] = {'chain': chain, 'resnum': resnum, 'x': x, 'y': y, 'z': z, 'has_p': False}
            elif atom_name == 'P' and key not in residues:
                residues[key] = {'chain': chain, 'resnum': resnum, 'x': x, 'y': y, 'z': z, 'has_p': True}

    return list(residues.values())


def parse_structure_coords(text, fmt):
    """Parse structure coordinates from PDB or CIF format."""
    return parse_pdb_coords(text) if fmt == 'pdb' else parse_cif_coords(text)


# ============================================================================
# Chain Boundary Extraction
# ============================================================================

def get_chains_from_pdb(pdb_text):
    """Extract chain names, sizes, and types from PDB ATOM records."""
    chain_order = []
    chain_counts = OrderedDict()
    seen_residues = set()

    for line in pdb_text.split('\n'):
        if not line.startswith('ATOM'):
            continue
        atom_name = line[12:16].strip()
        if atom_name not in ('CA', 'P'):
            continue
        chain = line[21:22].strip() or 'A'
        resnum = line[22:26].strip()
        rkey = f'{chain}:{resnum}'
        if rkey in seen_residues:
            continue
        seen_residues.add(rkey)
        if chain not in chain_counts:
            chain_order.append(chain)
            chain_counts[chain] = 0
        chain_counts[chain] += 1

    return {
        'names': chain_order,
        'sizes': [chain_counts[c] for c in chain_order],
        'types': ['protein'] * len(chain_order),
    }


def get_chains_from_cif(cif_text):
    """Extract chain names, sizes, and types from CIF atom_site records."""
    chain_order = []
    chain_counts = OrderedDict()
    chain_types = {}
    seen_residues = set()
    in_atom_site = False
    col_names = []

    for line in cif_text.split('\n'):
        if line.startswith('_atom_site.'):
            in_atom_site = True
            col_names.append(line.strip().split('.')[1])
            continue
        if in_atom_site and not line.startswith('_atom_site.') and not line.startswith('#') and line.strip():
            if line.startswith('loop_') or line.startswith('_'):
                in_atom_site = False
                continue
            parts = line.strip().split()
            if len(parts) < len(col_names):
                continue

            def get_col(name, _parts=parts, _cols=col_names):
                idx = _cols.index(name) if name in _cols else -1
                return _parts[idx] if idx >= 0 else ''

            group_pdb = get_col('group_PDB')
            atom_name = get_col('label_atom_id')
            comp_id = get_col('label_comp_id')
            chain = get_col('label_asym_id')
            res_seq = get_col('label_seq_id')

            counted = False
            if group_pdb == 'ATOM' and atom_name in ('CA', 'P'):
                rkey = f'{chain}:{res_seq}'
                if rkey not in seen_residues:
                    seen_residues.add(rkey)
                    counted = True
                    if atom_name == 'P':
                        chain_types[chain] = 'dna' if comp_id.startswith('D') else 'rna'
                    elif chain not in chain_types:
                        chain_types[chain] = 'protein'
            elif group_pdb == 'HETATM' and comp_id in ION_NAMES:
                rkey = f'{chain}:ion'
                if rkey not in seen_residues:
                    seen_residues.add(rkey)
                    counted = True
                    chain_types[chain] = 'ion'

            if counted:
                if chain not in chain_counts:
                    chain_order.append(chain)
                    chain_counts[chain] = 0
                chain_counts[chain] += 1

    return {
        'names': chain_order,
        'sizes': [chain_counts[c] for c in chain_order],
        'types': [chain_types.get(c, 'protein') for c in chain_order],
    }


def get_chains_from_structure(text, fmt):
    """Get chain info from structure file."""
    return get_chains_from_pdb(text) if fmt == 'pdb' else get_chains_from_cif(text)


# ============================================================================
# B-factor / pLDDT Parsing
# ============================================================================

def parse_bfactors_per_residue(text, fmt):
    """Parse B-factors (pLDDT) per residue from structure file.

    Returns dict of {chain:resnum: bfactor} using CA atoms.
    """
    bfactors = {}

    if fmt == 'pdb':
        for line in text.split('\n'):
            if not line.startswith('ATOM') or len(line) < 66:
                continue
            if line[12:16].strip() != 'CA':
                continue
            chain = line[21:22].strip() or 'A'
            try:
                rn = int(line[22:26].strip())
                bf = float(line[60:66].strip())
            except ValueError:
                continue
            bfactors[f'{chain}:{rn}'] = bf
    else:
        in_atom_site = False
        col_names = []
        for line in text.split('\n'):
            if line.startswith('_atom_site.'):
                in_atom_site = True
                col_names.append(line.strip().split('.')[1])
                continue
            if in_atom_site and not line.startswith('_atom_site.') and not line.startswith('#') and line.strip():
                if line.startswith('loop_') or line.startswith('_'):
                    in_atom_site = False
                    continue
                parts = line.strip().split()
                if len(parts) < len(col_names):
                    continue

                def get_col(name, _parts=parts, _cols=col_names):
                    idx = _cols.index(name) if name in _cols else -1
                    return _parts[idx] if idx >= 0 else ''

                if get_col('group_PDB') != 'ATOM' or get_col('label_atom_id') != 'CA':
                    continue
                chain = get_col('label_asym_id')
                try:
                    rn = int(get_col('label_seq_id'))
                    bf = float(get_col('B_iso_or_equiv'))
                except ValueError:
                    continue
                bfactors[f'{chain}:{rn}'] = bf

    return bfactors


def compute_chain_plddt(struct_text, fmt):
    """Compute per-chain average pLDDT from B-factors."""
    bfs = parse_bfactors_per_residue(struct_text, fmt)
    chain_vals = {}
    for key, val in bfs.items():
        chain = key.split(':')[0]
        chain_vals.setdefault(chain, []).append(val)
    return {c: sum(v) / len(v) for c, v in chain_vals.items() if v}


# ============================================================================
# Contact Map
# ============================================================================

def compute_contact_map(coords, threshold=8):
    """Compute NxN binary contact map from Cb coordinates.

    Uses Cb-Cb distance <= threshold (with 4A adjustment for phosphorus atoms).
    """
    n = len(coords)
    if n == 0:
        return np.zeros((0, 0), dtype=np.uint8), 0

    xyz = np.array([[c['x'], c['y'], c['z']] for c in coords])
    has_p = np.array([c['has_p'] for c in coords])

    distances = squareform(pdist(xyz))

    p_adjustment = np.zeros_like(distances)
    p_mask = has_p[:, np.newaxis] | has_p[np.newaxis, :]
    p_adjustment[p_mask] = -4.0
    adjusted = distances + p_adjustment

    contact = (adjusted < threshold).astype(np.uint8)
    return contact, n


# ============================================================================
# PAE Transform
# ============================================================================

def transform_pae_matrix(pae, pae_cutoff=12):
    """Symmetrize and transform PAE to confidence scores.

    Symmetrization: (PAE[i,j] + PAE[j,i]) / 2
    Transform: confidence = 1 - PAE/cutoff if PAE < cutoff, else 0
    """
    sym = (pae + pae.T) / 2.0
    transformed = np.zeros_like(sym)
    mask = sym < pae_cutoff
    transformed[mask] = 1.0 - sym[mask] / pae_cutoff
    return transformed


# ============================================================================
# ipSAE Calculation (Dunbrack et al. 2025)
# ============================================================================

def calc_ipsae(pae, si, ei, sj, ej, pae_cutoff):
    """Calculate ipSAE (Dunbrack d0res method) for a chain pair.

    d0 from per-residue count of inter-chain residues with PAE < cutoff.
    No distance filter. PAE cutoff only.
    Returns max of two asymmetric scores.
    """
    len_i = ei - si
    len_j = ej - sj
    if len_i == 0 or len_j == 0:
        return 0.0

    def _ipsae_one_direction(block):
        """block[r, s]: PAE from residues r->s. Returns max per-residue score."""
        mask = block < pae_cutoff
        good_counts = mask.sum(axis=1)  # per residue
        has_good = good_counts > 0
        if not has_good.any():
            return 0.0
        d0 = np.maximum(1.0, 1.24 * (np.maximum(good_counts[has_good], 27) - 15) ** (1.0 / 3.0) - 1.8)
        d0sq = d0 * d0
        # For each residue with good contacts, compute mean TM-score
        block_good = block[has_good]
        mask_good = mask[has_good]
        tm_vals = np.where(mask_good, 1.0 / (1.0 + (block_good ** 2) / d0sq[:, None]), 0.0)
        scores = tm_vals.sum(axis=1) / good_counts[has_good]
        return float(scores.max())

    block_ij = pae[si:ei, sj:ej].astype(np.float64)
    block_ji = pae[sj:ej, si:ei].astype(np.float64)
    return max(_ipsae_one_direction(block_ij), _ipsae_one_direction(block_ji))


def calc_ipsae_d0chn(pae, si, ei, sj, ej, pae_cutoff):
    """Calculate ipSAE with d0 from chain pair length (Dunbrack d0chn variant).

    d0 = f(len_chain1 + len_chain2). No distance filter (matches Dunbrack's d0chn).
    For each residue i in chain1, score = mean of TM(PAE[i,j]) for all j in chain2
    where PAE[i,j] < pae_cutoff. ipSAE = max over all residues.
    Returns max of two asymmetric scores (A->B, B->A).
    """
    len_i = ei - si
    len_j = ej - sj
    if len_i == 0 or len_j == 0:
        return 0.0

    d0 = max(1.0, 1.24 * (max(len_i + len_j, 27) - 15) ** (1.0 / 3.0) - 1.8)
    d0sq = d0 * d0

    # ipSAE_d0chn(I->J): for each residue in chain I, score against chain J
    max_score_ij = 0.0
    for ri in range(si, ei):
        score_sum = 0.0
        count = 0
        for rj in range(sj, ej):
            v = pae[ri, rj] if ri < pae.shape[0] and rj < pae.shape[1] else 31.0
            if v < pae_cutoff:
                score_sum += 1.0 / (1.0 + (v * v) / d0sq)
                count += 1
        if count > 0:
            score = score_sum / count
            if score > max_score_ij:
                max_score_ij = score

    # ipSAE_d0chn(J->I)
    max_score_ji = 0.0
    for rj in range(sj, ej):
        score_sum = 0.0
        count = 0
        for ri in range(si, ei):
            v = pae[rj, ri] if rj < pae.shape[0] and ri < pae.shape[1] else 31.0
            if v < pae_cutoff:
                score_sum += 1.0 / (1.0 + (v * v) / d0sq)
                count += 1
        if count > 0:
            score = score_sum / count
            if score > max_score_ji:
                max_score_ji = score

    return max(max_score_ij, max_score_ji)


# ============================================================================
# actifpTM Calculation (Varga & Ovchinnikov 2025)
# Approximation from PAE matrix + CB contact map (without distogram logits)
# ============================================================================

def calc_actifptm(pae, contact, n_use, si, ei, sj, ej):
    """Calculate approximate actifpTM (actual interface pTM) for a chain pair.

    Uses TM-score transform on PAE values, weighted by binary CB contacts.
    d0 is computed from full complex length (matching official implementation).
    Returns max per-residue weighted TM score across all interface residues.
    """
    n_total = pae.shape[0]
    clipped = max(n_total, 19)
    d0 = 1.24 * (clipped - 15) ** (1.0 / 3.0) - 1.8
    d0sq = d0 * d0

    def _actifptm_one_direction(r_start, r_end, s_start, s_end):
        ri_end = min(r_end, n_use)
        si_end = min(s_end, n_use)
        if ri_end <= r_start or si_end <= s_start:
            return 0.0
        contact_block = contact[r_start:ri_end, s_start:si_end]
        weight_sums = contact_block.sum(axis=1)
        has_contact = weight_sums > 0
        if not has_contact.any():
            return 0.0
        pae_block = pae[r_start:ri_end, s_start:si_end].astype(np.float64)
        tm_vals = np.where(contact_block, 1.0 / (1.0 + (pae_block ** 2) / d0sq), 0.0)
        scores = tm_vals[has_contact].sum(axis=1) / weight_sums[has_contact]
        return float(scores.max())

    return max(_actifptm_one_direction(si, ei, sj, ej),
               _actifptm_one_direction(sj, ej, si, ei))


# ============================================================================
# Single-Model Analysis
# ============================================================================

def _avg_bfactor(res_set, chain, bfactors):
    """Average B-factor for a set of residue numbers in a chain."""
    vals = []
    for r in res_set:
        v = bfactors.get(f'{chain}:{r}')
        if v is not None:
            vals.append(v)
    return sum(vals) / len(vals) if vals else float('nan')


def analyze_single_model(struct_text, pae_matrix, scores, fmt, platform,
                         pae_path, read_fn, pae_cutoff=12, cb_cutoff=8):
    """Compute LIS metrics for one model.

    Returns list of dicts (one per symmetrized chain pair), each containing
    all fields needed for one CSV row.
    """
    # Get chain info from structure
    chain_info = get_chains_from_structure(struct_text, fmt)

    # For AF3, try to get chain IDs from full_data JSON (token_chain_ids)
    chain_ids = None
    if platform in ('alphafold3', 'openfold3') and pae_path:
        pae_content = read_fn(pae_path)
        if isinstance(pae_content, str):
            try:
                fd = json.loads(pae_content)
                if 'token_chain_ids' in fd:
                    chain_ids = fd['token_chain_ids']
            except json.JSONDecodeError:
                pass

    # Build chain boundaries
    if chain_ids:
        chain_map = OrderedDict()
        for c in chain_ids:
            chain_map[c] = chain_map.get(c, 0) + 1
        chain_names = list(chain_map.keys())
        sizes = list(chain_map.values())
        n_total = len(chain_ids)
        struct_types = {}
        if chain_info.get('types'):
            for cname, t in zip(chain_info['names'], chain_info['types']):
                struct_types[cname] = t
    else:
        chain_names = chain_info['names']
        sizes = chain_info['sizes']
        n_total = sum(sizes)

    # Ensure PAE size matches
    pae = pae_matrix
    pae_size = pae.shape[0]
    if pae_size != n_total:
        n_total = pae_size

    # Cumulative sums for chain boundaries
    cum_sum = np.cumsum(sizes)
    starts = np.concatenate(([0], cum_sum[:-1]))

    # Build transformed confidence map (symmetrized PAE)
    transformed = transform_pae_matrix(pae[:n_total, :n_total], pae_cutoff)

    # Contact map + distance matrix
    coords = parse_structure_coords(struct_text, fmt)
    contact, n_coords = compute_contact_map(coords, cb_cutoff)
    n_use = min(n_total, n_coords)

    # Distance matrix for ipSAE (15Å cutoff)
    dist_matrix = None
    if len(coords) > 0:
        xyz = np.array([[c['x'], c['y'], c['z']] for c in coords])
        from scipy.spatial.distance import cdist
        dist_matrix = cdist(xyz, xyz)

    # ipTM
    iptm_matrix = scores.get('chainPairIptm')
    global_iptm = scores.get('ipTM', 0)

    # B-factors for LIpLDDT
    bfs = parse_bfactors_per_residue(struct_text, fmt)
    chain_plddt = compute_chain_plddt(struct_text, fmt)

    # Per chain pair metrics (asymmetric: both (i,j) and (j,i))
    nc = len(chain_names)
    pairs = {}

    for i in range(nc):
        for j in range(nc):
            if i == j:
                continue
            si, ei = int(starts[i]), int(min(cum_sum[i], n_total))
            sj, ej = int(starts[j]), int(min(cum_sum[j], n_total))

            # Vectorized LIS/cLIS computation
            t_block = transformed[si:ei, sj:ej]
            t_pos = t_block > 0

            lis_sum = float(t_block[t_pos].sum())
            lis_count_avg = int(t_pos.sum())

            # LIR: residues with at least one positive transformed value
            lir_i_mask = t_pos.any(axis=1)
            lir_j_mask = t_pos.any(axis=0)
            lir_i = set(np.where(lir_i_mask)[0] + 1)
            lir_j = set(np.where(lir_j_mask)[0] + 1)

            # Contact-weighted (cLIS)
            c_ei = min(ei, n_use)
            c_ej = min(ej, n_use)
            c_si = si
            c_sj = sj
            if c_ei > c_si and c_ej > c_sj:
                contact_block = contact[c_si:c_ei, c_sj:c_ej].astype(bool)
                t_contact = t_block[:c_ei-c_si, :c_ej-c_sj]
                ct_pos = t_pos[:c_ei-c_si, :c_ej-c_sj] & contact_block
                clis_sum = float(t_contact[ct_pos].sum())
                clis_count_avg = int(ct_pos.sum())
                clir_i = set(np.where(ct_pos.any(axis=1))[0] + 1)
                clir_j = set(np.where(ct_pos.any(axis=0))[0] + 1)
            else:
                clis_sum = 0.0
                clis_count_avg = 0
                clir_i = set()
                clir_j = set()

            # LIA counts (asymmetric PAE < cutoff)
            pae_ij = pae[si:ei, sj:ej]
            pae_ji = pae[sj:ej, si:ei]
            lis_count_ab = int((pae_ij < pae_cutoff).sum())
            lis_count_ba = int((pae_ji < pae_cutoff).sum())

            if c_ei > c_si and c_ej > c_sj:
                pae_ij_c = pae_ij[:c_ei-c_si, :c_ej-c_sj]
                pae_ji_c = pae_ji[:c_ej-c_sj, :c_ei-c_si]
                clis_count_ab = int(((pae_ij_c < pae_cutoff) & contact_block).sum())
                clis_count_ba = int(((pae_ji_c < pae_cutoff) & contact_block.T).sum())
            else:
                clis_count_ab = 0
                clis_count_ba = 0

            lis_val = lis_sum / lis_count_avg if lis_count_avg > 0 else 0.0
            clis_val = clis_sum / clis_count_avg if clis_count_avg > 0 else 0.0
            ilis_val = math.sqrt(lis_val * clis_val)

            iptm_val = global_iptm
            if iptm_matrix:
                try:
                    if isinstance(iptm_matrix, dict):
                        # OpenFold3 format: {'(A, B)': 0.235}
                        key = f'({chain_names[i]}, {chain_names[j]})'
                        if key in iptm_matrix:
                            iptm_val = float(iptm_matrix[key])
                    elif isinstance(iptm_matrix, (list, np.ndarray)):
                        if i < len(iptm_matrix) and j < len(iptm_matrix[i]):
                            v = iptm_matrix[i][j]
                            if v is not None:
                                iptm_val = float(v)
                except (KeyError, IndexError, TypeError):
                    pass

            try:
                ipsae = calc_ipsae(pae, si, min(ei, pae.shape[0]),
                                    sj, min(ej, pae.shape[1]), 10)
            except Exception:
                ipsae = 0.0

            try:
                actifptm = calc_actifptm(pae, contact, n_use, si, min(ei, pae.shape[0]),
                                          sj, min(ej, pae.shape[1]))
            except Exception:
                actifptm = 0.0

            lia_count = lis_count_ab + lis_count_ba
            clia_count = clis_count_ab + clis_count_ba
            ilia_val = math.sqrt(lia_count * clia_count)
            ilisa_val = ilis_val * ilia_val

            liplddt_i = _avg_bfactor(lir_i, chain_names[i], bfs)
            liplddt_j = _avg_bfactor(lir_j, chain_names[j], bfs)
            cliplddt_i = _avg_bfactor(clir_i, chain_names[i], bfs)
            cliplddt_j = _avg_bfactor(clir_j, chain_names[j], bfs)

            key = f'{chain_names[i]},{chain_names[j]}'
            pairs[key] = {
                'ci': chain_names[i],
                'cj': chain_names[j],
                'LIS': lis_val,
                'cLIS': clis_val,
                'iLIS': ilis_val,
                'iLIA': ilia_val,
                'iLISA': ilisa_val,
                'ipTM': iptm_val,
                'ipSAE': ipsae,
                'actifpTM': actifptm,
                'LIA': lia_count,
                'cLIA': clia_count,
                'LIpLDDT_i': liplddt_i,
                'LIpLDDT_j': liplddt_j,
                'cLIpLDDT_i': cliplddt_i,
                'cLIpLDDT_j': cliplddt_j,
                'lirI': lir_i,
                'lirJ': lir_j,
                'clirI': clir_i,
                'clirJ': clir_j,
                'lenI': sizes[i],
                'lenJ': sizes[j],
            }

    # Symmetrize: combine (i,j) and (j,i) into one entry
    symmetric = {}
    seen = set()
    for key, val in pairs.items():
        ci, cj = key.split(',')
        canon = ','.join(sorted([ci, cj]))
        if canon in seen:
            continue
        seen.add(canon)

        rev_key = f'{cj},{ci}'
        if rev_key in pairs:
            rv = pairs[rev_key]
            s = {
                'ci': ci, 'cj': cj,
                'LIS': (val['LIS'] + rv['LIS']) / 2,
                'cLIS': (val['cLIS'] + rv['cLIS']) / 2,
                'ipTM': (val['ipTM'] + rv['ipTM']) / 2,
                'ipSAE': max(val['ipSAE'], rv['ipSAE']),
                'actifpTM': max(val['actifpTM'], rv['actifpTM']),
                'LIA': val['LIA'] + rv['LIA'],
                'cLIA': val['cLIA'] + rv['cLIA'],
                'lenI': val['lenI'], 'lenJ': val['lenJ'],
                'lirI': val['lirI'], 'lirJ': rv['lirI'],
                'clirI': val['clirI'], 'clirJ': rv['clirI'],
                'LIpLDDT_i': val['LIpLDDT_i'], 'LIpLDDT_j': rv['LIpLDDT_i'],
                'cLIpLDDT_i': val['cLIpLDDT_i'], 'cLIpLDDT_j': rv['cLIpLDDT_i'],
            }
            s['iLIS'] = math.sqrt(s['LIS'] * s['cLIS'])
            s['iLIA'] = math.sqrt(s['LIA'] * s['cLIA'])
            s['iLISA'] = s['iLIS'] * s['iLIA']
            symmetric[f'{ci}-{cj}'] = s
        else:
            val_copy = dict(val)
            val_copy['iLIS'] = math.sqrt(val['LIS'] * val['cLIS'])
            val_copy['iLIA'] = math.sqrt(val['LIA'] * val['cLIA'])
            val_copy['iLISA'] = val_copy['iLIS'] * val_copy['iLIA']
            symmetric[f'{ci}-{cj}'] = val_copy

    # Build result list with per-chain pLDDT from structure
    results = []
    for key, v in sorted(symmetric.items()):
        plddt_i = chain_plddt.get(v['ci'])
        plddt_j = chain_plddt.get(v['cj'])
        if plddt_i is None and 'pLDDT' in scores:
            plddt_i = scores['pLDDT']
        if plddt_j is None and 'pLDDT' in scores:
            plddt_j = scores['pLDDT']

        results.append({
            **v,
            'pLDDT_i': plddt_i,
            'pLDDT_j': plddt_j,
            'pTM': scores.get('pTM'),
        })

    return results


# ============================================================================
# CSV Output
# ============================================================================

def format_indices(res_set):
    """Format a set of residue indices as a compact range string for CSV.
    e.g., {1,2,3,4,10,11,13} -> '"1-4,10-11,13"'
    Quoted because commas are used inside.
    """
    if not res_set:
        return ''
    sorted_pos = sorted(res_set)
    ranges = []
    start = end = sorted_pos[0]
    for p in sorted_pos[1:]:
        if p == end + 1:
            end = p
        else:
            ranges.append(f'{start}-{end}' if start != end else str(start))
            start = end = p
    ranges.append(f'{start}-{end}' if start != end else str(start))
    return '"[' + ','.join(ranges) + ']"'


def _extract_model_num(struct_filename):
    """Extract model number (integer) from structure filename."""
    m = re.search(r'model_(\d+)', struct_filename)
    if m:
        return m.group(1)
    m = re.search(r'sample_(\d+)', struct_filename)
    if m:
        return m.group(1)
    m = re.search(r'model_idx_(\d+)', struct_filename)
    if m:
        return m.group(1)
    return ''


def format_row(name, rank, struct_file, pair):
    """Format one CSV row from a pair dict."""
    def fmt_plddt(v):
        if v is None or (isinstance(v, float) and math.isnan(v)):
            return ''
        return f'{v:.1f}'

    model_num = _extract_model_num(struct_file)

    row = [
        name, rank, model_num, pair['ci'], pair['cj'],
        f"{pair['iLIS']:.4f}", f"{pair['iLIA']:.1f}", f"{pair['iLISA']:.1f}",
        f"{pair['ipSAE']:.4f}", f"{pair.get('actifpTM', 0):.4f}",
        f"{pair['LIS']:.4f}", f"{pair['cLIS']:.4f}",
        f"{pair['LIA']:.1f}", f"{pair['cLIA']:.1f}",
        f"{pair['ipTM']:.4f}",
        fmt_plddt(pair.get('pLDDT_i')),
        fmt_plddt(pair.get('pLDDT_j')),
        f"{pair['pTM']:.3f}" if pair.get('pTM') is not None else '',
        str(len(pair['lirI'])), str(len(pair['lirJ'])),
        str(len(pair['clirI'])), str(len(pair['clirJ'])),
        fmt_plddt(pair.get('LIpLDDT_i')),
        fmt_plddt(pair.get('LIpLDDT_j')),
        fmt_plddt(pair.get('cLIpLDDT_i')),
        fmt_plddt(pair.get('cLIpLDDT_j')),
        str(pair['lenI']), str(pair['lenJ']),
        format_indices(pair['lirI']),
        format_indices(pair['lirJ']),
        format_indices(pair['clirI']),
        format_indices(pair['clirJ']),
        struct_file,
    ]
    return ','.join(row)


# ============================================================================
# Worker Functions (top-level for multiprocessing pickling)
# ============================================================================

def _do_process(model_tuple, read_fn, detected, pae_cutoff, cb_cutoff, verbose=False):
    """Process a single model. Returns (name, rank, rows, error_msg) tuple."""
    name, rank, model_label, struct_path, pae_path, scores_path, fmt = model_tuple

    struct_text = read_fn(struct_path)
    if not struct_text or not isinstance(struct_text, str):
        return name, rank, None, f'structure file unreadable: {struct_path}'

    pae = extract_pae(pae_path, read_fn)
    if pae is None:
        return name, rank, None, f'PAE not found or unreadable: {pae_path}'
    pae = np.nan_to_num(pae, nan=31.0)

    scores = extract_confidence_scores(scores_path, read_fn)
    if scores_path != pae_path and pae_path:
        full_scores = extract_confidence_scores(pae_path, read_fn)
        for k, v in full_scores.items():
            if k not in scores:
                scores[k] = v

    try:
        pairs = analyze_single_model(
            struct_text, pae, scores, fmt, detected,
            pae_path, read_fn, pae_cutoff, cb_cutoff)
    except Exception as e:
        return name, rank, None, f'analysis error: {e}'

    rows = [format_row(name, rank, model_label, pair) for pair in pairs]
    return name, rank, rows, None


def _process_one_sequential(model_tuple, read_fn, detected, pae_cutoff, cb_cutoff, verbose=False):
    """Sequential wrapper."""
    return _do_process(model_tuple, read_fn, detected, pae_cutoff, cb_cutoff, verbose)


def _mp_worker(args):
    """Multiprocessing worker — creates its own read_fn from file_map."""
    model_tuple, detected, pae_cutoff, cb_cutoff, file_map, verbose = args
    read_fn = _make_dir_reader(file_map)
    return _do_process(model_tuple, read_fn, detected, pae_cutoff, cb_cutoff, verbose)


def _sort_csv(filepath):
    """Sort a CSV file by name (col 0), then rank (col 1)."""
    try:
        import csv as _csv
        with open(filepath, 'r') as f:
            reader = _csv.reader(f)
            header = next(reader)
            rows = list(reader)
        if not rows:
            return
        rows.sort(key=lambda r: (r[0], r[1].zfill(5) if len(r) > 1 and r[1] else ''))
        with open(filepath, 'w', newline='') as f:
            writer = _csv.writer(f)
            writer.writerow(header)
            writer.writerows(rows)
        pass  # sorted silently
    except Exception:
        pass


# ============================================================================
# Main Pipeline
# ============================================================================

PLATFORM_LABELS = {
    'alphafold3': 'AlphaFold3',
    'colabfold': 'ColabFold',
    'boltz': 'Boltz',
    'chai': 'Chai-1',
    'openfold3': 'OpenFold3',
    'generic': 'Generic',
}


def run(path, output=None, output_dir=None, pae_cutoff=12, cb_cutoff=8,
        platform=None, skip_existing=True, workers=1, verbose=False):
    """Run the LIS analysis pipeline.

    Processes each model independently, one at a time, appending CSV rows
    incrementally. No grouping, no averaging.
    """
    path = str(path)
    print(f'[LIS] Scanning files from: {path}')
    filenames, read_fn, file_map = scan_files(path)
    print(f'[LIS] Found {len(filenames)} files')

    if not filenames:
        print('[LIS] ERROR: No files found', file=sys.stderr)
        sys.exit(1)

    # Detect platform
    if platform:
        detected = platform
    else:
        detected = detect_platform(filenames, read_fn)
    print(f'[LIS] Platform: {PLATFORM_LABELS.get(detected, detected)}')

    # Find all models
    models = list(find_models(filenames, detected, read_fn))
    # Filter macOS resource forks
    models = [m for m in models if not m[0].startswith('._')]

    if not models:
        print('[LIS] ERROR: No prediction models found', file=sys.stderr)
        sys.exit(1)

    print(f'[LIS] {len(models)} model(s) found')

    # Determine output path
    if output is None:
        basename = os.path.basename(path.rstrip('/'))
        basename = re.sub(r'\.(zip|tar\.gz|tar\.zstd|tar\.zst|tgz)$', '', basename)
        output = f'{basename}_lis_analysis.csv'
    if output_dir is None:
        output_dir = path if os.path.isdir(path) else os.path.dirname(os.path.abspath(path))
    os.makedirs(output_dir, exist_ok=True)
    output = os.path.join(output_dir, os.path.basename(output))

    # Load existing (name,rank) combos for skip-existing
    existing_keys = set()
    if skip_existing and os.path.exists(output):
        with open(output, 'r') as f:
            for line in f:
                parts = line.split(',')
                if len(parts) >= 2 and parts[0] != 'name':
                    existing_keys.add((parts[0], parts[1]))
        if existing_keys:
            print(f'[LIS] Found existing CSV with {len(existing_keys)} (name,rank) combo(s) -- will skip these')

    # Filter to models that need processing
    models_todo = [m for m in models if (m[0], m[1]) not in existing_keys]
    total_skipped = len(models) - len(models_todo)

    if not models_todo:
        print(f'[LIS] All {total_skipped} model(s) already processed.')
        # Still sort existing CSV
        _sort_csv(output)
        return output

    print(f'[LIS] {len(models_todo)} to process' +
          (f', {total_skipped} skipped' if total_skipped > 0 else ''))

    # Write header if file is new
    if not os.path.exists(output) or not existing_keys:
        with open(output, 'w') as f:
            f.write(CSV_HEADER + '\n')

    import time as _time
    total_done = 0
    total_failed = 0
    total = len(models_todo)
    t_start = _time.time()

    def _progress(name, rank, ok, err_msg=None):
        nonlocal total_done, total_failed
        if ok:
            total_done += 1
        else:
            total_failed += 1
        n = total_done + total_failed
        elapsed = _time.time() - t_start
        per_item = elapsed / n if n > 0 else 0
        eta = per_item * (total - n)
        elapsed_str = f'{int(elapsed // 60)}m{int(elapsed % 60):02d}s' if elapsed >= 60 else f'{int(elapsed)}s'
        eta_str = f'{int(eta // 60)}m{int(eta % 60):02d}s' if eta >= 60 else f'{int(eta)}s'
        pct = n * 100 // total
        bar_len = 30
        filled = bar_len * n // total
        bar = '█' * filled + '░' * (bar_len - filled)
        status = 'OK' if ok else 'FAIL'
        print(f'\r[LIS] {bar} {pct}% ({n}/{total}) {elapsed_str} elapsed, ETA {eta_str} | {name} {status}      ', end='', flush=True)
        if not ok and err_msg and verbose:
            print(f'\n      >> {err_msg}', flush=True)

    # Check if input is a folder (needed for multiprocessing — zip read_fn can't be pickled)
    is_folder = os.path.isdir(path)

    if workers > 1 and is_folder:
        from multiprocessing import Pool

        # file_map is picklable (dict of str → str for folders)
        task_args = [
            (m, detected, pae_cutoff, cb_cutoff, file_map, verbose)
            for m in models_todo
        ]

        with Pool(processes=workers) as pool:
            for name, rank, result_rows, err_msg in pool.imap_unordered(_mp_worker, task_args):
                if result_rows:
                    with open(output, 'a') as f:
                        f.write('\n'.join(result_rows) + '\n')
                _progress(name, rank, result_rows is not None, err_msg)
    else:
        if workers > 1 and not is_folder:
            print('[LIS] Note: parallel mode only works with folders (not zips). Using sequential.')
        for gi, model_tuple in enumerate(models_todo):
            name, rank = model_tuple[0], model_tuple[1]
            name, rank, result_rows, err_msg = _process_one_sequential(
                model_tuple, read_fn, detected, pae_cutoff, cb_cutoff, verbose)
            if result_rows:
                with open(output, 'a') as f:
                    f.write('\n'.join(result_rows) + '\n')
            _progress(name, rank, result_rows is not None, err_msg)

    print()  # newline after progress bar
    _sort_csv(output)

    elapsed = _time.time() - t_start
    elapsed_str = f'{int(elapsed // 60)}m{int(elapsed % 60):02d}s' if elapsed >= 60 else f'{elapsed:.1f}s'
    summary = f'{total_done} done'
    if total_skipped > 0:
        summary += f', {total_skipped} skipped'
    if total_failed > 0:
        summary += f', {total_failed} failed'
    print(f'[LIS] {summary} in {elapsed_str} — {output}\n')

    # Clean up temp directory from tar.zstd extraction
    tmpdir = file_map.get('__tmpdir__')
    if tmpdir and os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir, ignore_errors=True)
        print(f'[LIS] Cleaned up temp directory')

    return output


# ============================================================================
# CLI
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='LIS -- Local Interaction Score Analysis (CLI)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Supported platforms (auto-detected):
  AlphaFold3   *_model_N.cif + *_full_data_N.json + *_summary_confidences_N.json
  ColabFold    *_unrelaxed_rank_N*.pdb + *_scores_rank_N*.json
  Boltz        *.pdb/.cif + confidence_*.json + pae_*.npz
  Chai-1       pred.rank_N.cif + scores.model_idx_N.json + pae.*.npy/.npz
  OpenFold3    result_sample_N_model.pdb + result_sample_N_confidences*.json
  Generic      any .pdb/.cif + PAE .json

Examples:
  python lis.py alphafold3_output.zip
  python lis.py colabfold_results/ -o results.csv
  python lis.py prediction.zip --pae-cutoff 10 --platform boltz
  python lis.py /path/to/predictions/ -d /path/to/output/
  python lis.py /path/to/predictions/ --no-skip-existing
  python lis.py /path/to/predictions/ -w 4
        """,
    )
    parser.add_argument('path', help='Path to folder or zip file with prediction outputs')
    parser.add_argument('--output', '-o', default=None,
                        help='Output CSV filename (default: <name>_lis_analysis.csv)')
    parser.add_argument('--output-dir', '-d', default=None,
                        help='Output directory (default: input folder)')
    parser.add_argument('--pae-cutoff', type=float, default=12,
                        help='PAE cutoff for confidence transform (default: 12)')
    parser.add_argument('--cb-cutoff', type=float, default=8,
                        help='Cb distance cutoff in Angstroms (default: 8)')
    parser.add_argument('--platform', default=None,
                        choices=['alphafold3', 'colabfold', 'boltz', 'chai', 'openfold3', 'generic'],
                        help='Force platform detection (default: auto-detect)')
    parser.add_argument('--workers', '-w', type=int, default=1,
                        help='Number of parallel workers (default: 1)')
    parser.add_argument('--no-skip-existing', action='store_true',
                        help='Reprocess all predictions even if already in the output CSV')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Show error details for failed predictions')
    args = parser.parse_args()

    if not os.path.exists(args.path):
        print(f'[LIS] ERROR: Path not found: {args.path}', file=sys.stderr)
        sys.exit(1)

    run(args.path, output=args.output, output_dir=args.output_dir,
        pae_cutoff=args.pae_cutoff, cb_cutoff=args.cb_cutoff,
        platform=args.platform,
        skip_existing=not args.no_skip_existing,
        workers=args.workers, verbose=args.verbose)


if __name__ == '__main__':
    main()
