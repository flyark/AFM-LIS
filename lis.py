#!/usr/bin/env python3
"""
LIS -- Local Interaction Score Analysis (CLI)
==============================================
Calculate LIS/cLIS/iLIS metrics from structure prediction outputs.

Supports: AlphaFold3, ColabFold, AlphaPulldown (AF-Multimer), Boltz, Chai-1, OpenFold3, Generic
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
    'ipTM,pLDDT_i,pLDDT_j,pLDDT,pTM,LIR_i,LIR_j,cLIR_i,cLIR_j,'
    'LIpLDDT_i,LIpLDDT_j,LIpLDDT,cLIpLDDT_i,cLIpLDDT_j,cLIpLDDT,'
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

    # Protenix-v2: *_sample_N.cif + *_summary_confidence_sample_N.json + *_full_data_sample_N.json
    has_prot_model = any(re.search(r'_sample_\d+\.cif$', b) for b in basenames)
    has_prot_summary = any(re.search(r'_summary_confidence_sample_\d+\.json$', b) for b in basenames)
    has_prot_full = any(re.search(r'_full_data_sample_\d+\.json$', b) for b in basenames)
    if has_prot_model and has_prot_summary and has_prot_full:
        return 'alphafold3'  # use the AF3 discovery path; _find_af3 handles the Protenix layout

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

    # AlphaPulldown / AlphaFold-Multimer: pae_model_N_*.json + confidence_model_N_*.json
    # (+ ranked_N.pdb / unrelaxed_model_N_*.pdb / ranking_debug.json)
    # Must come before Boltz because both have confidence_*.json files.
    has_ap_pae = any(re.match(r'^pae_model_\d+_.+\.json$', b) for b in basenames)
    has_ap_conf = any(re.match(r'^confidence_model_\d+_.+\.json$', b) for b in basenames)
    if has_ap_pae and has_ap_conf:
        return 'alphapulldown'

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

    # ESMFold2 (biohub API): *_model.pdb + *_pae.json
    has_esm_model = any(re.search(r'_model\.pdb$', b) for b in basenames)
    has_esm_pae = any(re.search(r'_pae\.json$', b) for b in basenames)
    if has_esm_model and has_esm_pae:
        return 'esmfold2'

    # Folder-per-pair pattern (filename-agnostic, format-flexible). Any folder that
    # contains at least one structure file (.pdb / .cif / .pdb.gz / .cif.gz / .pdb.xz /
    # .cif.xz) AND at least one PAE-style array file (.npy / .npz) is treated as a
    # single-model prediction. An optional .json carries the scores. The folder
    # basename becomes the prediction name. Mixed redundancy (e.g. both .pdb AND
    # .pdb.gz in the same folder) is fine — _find_esmfold2_native picks one by
    # priority (uncompressed > compressed, .npy > .npz).
    from collections import defaultdict as _dd
    _by_dir = _dd(list)
    for _f in filenames:
        _by_dir[os.path.dirname(_f)].append(os.path.basename(_f))
    _STRUCT_EXTS = ('.pdb', '.cif', '.pdb.gz', '.cif.gz', '.pdb.xz', '.cif.xz')
    _PAE_EXTS = ('.npy', '.npz')
    for _d, _names in _by_dir.items():
        _has_struct = any(n.endswith(_STRUCT_EXTS) for n in _names)
        _has_pae    = any(n.endswith(_PAE_EXTS) for n in _names)
        if _has_struct and _has_pae:
            return 'esmfold2_native'

    # Local AF3 (official google-deepmind/alphafold3 output): name-prefixed, no _N index, e.g.
    #   job/job_model.cif + job_confidences.json (PAE) + job_summary_confidences.json + job_data.json (input),
    #   plus per-sample job/seed-S_sample-M/job_seed-S_sample-M_confidences.json subdirs. The indexed
    #   ('_model_N.cif') and bare ('seed-N_sample-M/model.cif') detectors miss it. Placed LAST, so it only
    #   upgrades a folder that every other platform rejected — it can never re-classify one that already works.
    has_local_conf = any(b.endswith('_confidences.json') and not b.endswith('_summary_confidences.json') for b in basenames)
    has_local_summary = any(b.endswith('_summary_confidences.json') for b in basenames)
    has_local_model = any(b.endswith('_model.cif') for b in basenames)
    if has_local_conf and has_local_summary and has_local_model:
        return 'alphafold3'

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
        got = False
        for t in _find_af3(filenames, basenames_map):
            got = True
            yield t
        if not got:                      # official local AlphaFold 3 output (name-prefixed, no _N index)
            yield from _find_af3_local(filenames, basenames_map)
    elif platform == 'alphapulldown':
        yield from _find_alphapulldown(filenames, basenames_map, read_fn)
    elif platform == 'esmfold2':
        yield from _find_esmfold2(filenames, basenames_map)
    elif platform == 'esmfold2_native':
        yield from _find_esmfold2_native(filenames, basenames_map)
    elif platform == 'boltz':
        yield from _find_boltz(filenames, basenames_map)
    elif platform == 'chai':
        yield from _find_chai(filenames, basenames_map)
    else:
        yield from _find_generic(filenames, basenames_map, read_fn)


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

        # Find matching scores file. Require the prefix to be followed by
        # '_scores_' so we don't grab a longer prediction's file when one
        # name is a prefix of another (e.g. AAA___BBB vs AAA___BBB___CCC).
        scores_path = None
        padded = rank.zfill(3)
        scores_prefix = prefix + '_scores_'
        for fn in filenames:
            fb = os.path.basename(fn)
            if fb.startswith(scores_prefix) and f'rank_{padded}' in fb and fb.endswith('.json'):
                scores_path = fn
                break
        if not scores_path:
            for fn in filenames:
                fb = os.path.basename(fn)
                if fb.startswith(scores_prefix) and f'rank_{rank}' in fb and fb.endswith('.json'):
                    scores_path = fn
                    break
        if not scores_path:
            for fn in filenames:
                fb = os.path.basename(fn)
                if fb.startswith(scores_prefix) and fb.endswith('.json'):
                    scores_path = fn
                    break

        # ColabFold: PAE is inside the scores JSON under 'pae' key
        yield (prefix, rank, model_str, name, scores_path, scores_path, fmt)


def _find_af3(filenames, basenames_map):
    """AlphaFold3 standard: *_model_N.cif + *_full_data_N.json + *_summary_confidences_N.json
    Tamarind AF3: result_sample_N_model.pdb + result_sample_N_confidences.json
    Protenix-v2: *_sample_N.cif + *_full_data_sample_N.json + *_summary_confidence_sample_N.json
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

    # Protenix-v2 (AF3-like layout, slightly different basenames)
    for name in filenames:
        base = os.path.basename(name)
        m = re.match(r'^(.+)_sample_(\d+)\.cif$', base)
        if not m:
            continue
        prefix = m.group(1)
        idx = m.group(2)
        full_data = basenames_map.get(f'{prefix}_full_data_sample_{idx}.json')
        summary = basenames_map.get(f'{prefix}_summary_confidence_sample_{idx}.json')
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


def _find_af3_local(filenames, basenames_map):
    """Official local AlphaFold 3 output (github.com/google-deepmind/alphafold3).

    Files are name-prefixed with no _N index, which the indexed ('_model_N.cif') and
    the bare ('seed-N_sample-M/model.cif') AF3 finders both miss:
        job/job_model.cif                 job/job_confidences.json          (PAE: 'pae' key)
        job/job_summary_confidences.json  job/job_data.json                 (input; no PAE)
        job/seed-S_sample-M/job_seed-S_sample-M_model.cif  (+ _confidences.json, _summary_confidences.json)

    Prefers the per-sample ensemble (one rank each); the top-level model is a copy of
    the best sample, so it is used only when no seed-S_sample-M subdirs are present.
    """
    filenames_set = set(filenames)
    seed_models, top_models = [], []
    for name in filenames:
        base = os.path.basename(name)
        if not base.endswith('_model.cif'):
            continue
        prefix = base[:-len('_model.cif')]                       # 'job' or 'job_seed-S_sample-M'
        dirpath = os.path.dirname(name)
        conf = os.path.join(dirpath, f'{prefix}_confidences.json')
        if conf not in filenames_set:                            # a PAE sibling is required
            continue
        summ = os.path.join(dirpath, f'{prefix}_summary_confidences.json')
        summ = summ if summ in filenames_set else conf
        m = re.search(r'seed-(\d+)_sample-(\d+)', prefix) or re.search(r'seed-(\d+)_sample-(\d+)', dirpath)
        pred = re.sub(r'_seed-\d+_sample-\d+$', '', prefix) or _get_toplevel_name(name) or 'af3_prediction'
        if m:
            seed_models.append((pred, f'{m.group(1)}_{m.group(2)}', base, name, conf, summ))
        else:
            top_models.append((pred, '0', base, name, conf, summ))
    for pred, rank, base, name, conf, summ in (seed_models or top_models):
        yield (pred, rank, base, name, conf, summ, 'cif')


def _find_alphapulldown(filenames, basenames_map, read_fn):
    """AlphaPulldown / AlphaFold-Multimer (via AlphaPulldown).

    AlphaFold/AlphaPulldown writes one bundle per protein pair:
      - pae_<model_id>.json         AlphaFold pae_json output:
                                    [{"predicted_aligned_error": [[..]],
                                      "max_predicted_aligned_error": float}]
      - confidence_<model_id>.json  AlphaFold confidence_json output:
                                    {"residueNumber": [..],
                                     "confidenceScore": [..],  ← per-residue pLDDT
                                     "confidenceCategory": [..]}
      - unrelaxed_<model_id>.pdb    per-model coordinates (B-factors = pLDDT)
      - relaxed_<model_id>.pdb      optional, AMBER-relaxed
      - ranked_N.pdb                rank-sorted copies (N=0 is best)
      - ranking_debug.json          {"iptm+ptm": {model_id: score},
                                     "order": [model_id, ...]}

    Note: per-model scalar pTM and ipTM are only stored inside
    result_<model_id>.pkl (which we do not parse). The CSV's ipTM/pTM
    columns will therefore be empty for AlphaPulldown — iLIS/cLIS/ipSAE
    from PAE are the primary screening metrics here.

    <model_id> is typically `model_M_multimer_v3_pred_K` (M = 1..5 model
    variant, K = prediction index). Each prediction pair usually lives in
    its own subdirectory, so we group files by directory and process each
    group independently.
    """
    from collections import defaultdict
    by_dir = defaultdict(list)
    for f in filenames:
        by_dir[os.path.dirname(f)].append(f)

    for d, files in by_dir.items():
        local = {os.path.basename(f): f for f in files}

        pae_by_id = {}
        for b, f in local.items():
            m = re.match(r'^pae_(model_\d+.*?)\.json$', b)
            if m:
                pae_by_id[m.group(1)] = f

        conf_by_id = {}
        for b, f in local.items():
            m = re.match(r'^confidence_(model_\d+.*?)\.json$', b)
            if m:
                conf_by_id[m.group(1)] = f

        # Prefer relaxed over unrelaxed when both exist for the same model_id
        struct_by_id = {}
        for b, f in local.items():
            m = re.match(r'^unrelaxed_(model_\d+.*?)\.(pdb|cif)$', b)
            if m:
                struct_by_id.setdefault(m.group(1), f)
        for b, f in local.items():
            m = re.match(r'^relaxed_(model_\d+.*?)\.(pdb|cif)$', b)
            if m:
                struct_by_id[m.group(1)] = f

        ranked_by_rank = {}
        for b, f in local.items():
            m = re.match(r'^ranked_(\d+)\.(pdb|cif)$', b)
            if m:
                ranked_by_rank[int(m.group(1))] = f

        if not pae_by_id:
            continue
        if not struct_by_id and not ranked_by_rank:
            continue

        # Read ranking_debug.json → ordered list of model_ids (best first)
        rank_order = []
        if 'ranking_debug.json' in local:
            try:
                content = read_fn(local['ranking_debug.json'])
                if isinstance(content, str):
                    dbg = json.loads(content)
                    raw = dbg.get('order', [])
                    if isinstance(raw, list):
                        rank_order = [str(x) for x in raw]
            except (json.JSONDecodeError, ValueError, KeyError, TypeError):
                pass

        pred_name = os.path.basename(d) or 'prediction'

        if rank_order and struct_by_id:
            # Best path: rank index → model_id → per-model struct + PAE
            for rank_idx, model_id in enumerate(rank_order):
                if model_id not in pae_by_id or model_id not in struct_by_id:
                    continue
                struct_path = struct_by_id[model_id]
                pae_path = pae_by_id[model_id]
                conf_path = conf_by_id.get(model_id, pae_path)
                fmt = 'cif' if struct_path.lower().endswith('.cif') else 'pdb'
                yield (pred_name, str(rank_idx), os.path.basename(struct_path),
                       struct_path, pae_path, conf_path, fmt)
        elif struct_by_id:
            # No ranking info — iterate by sorted model_id with 0-indexed rank
            sorted_ids = sorted(set(struct_by_id.keys()) & set(pae_by_id.keys()))
            for rank_idx, model_id in enumerate(sorted_ids):
                struct_path = struct_by_id[model_id]
                pae_path = pae_by_id[model_id]
                conf_path = conf_by_id.get(model_id, pae_path)
                fmt = 'cif' if struct_path.lower().endswith('.cif') else 'pdb'
                yield (pred_name, str(rank_idx), os.path.basename(struct_path),
                       struct_path, pae_path, conf_path, fmt)
        else:
            # Only ranked_N.pdb available — need ranking_debug to match PAE
            if not rank_order:
                continue
            for rank_idx, model_id in enumerate(rank_order):
                if rank_idx not in ranked_by_rank or model_id not in pae_by_id:
                    continue
                struct_path = ranked_by_rank[rank_idx]
                pae_path = pae_by_id[model_id]
                conf_path = conf_by_id.get(model_id, pae_path)
                fmt = 'cif' if struct_path.lower().endswith('.cif') else 'pdb'
                # Use model_id label so the CSV's `model` column carries the
                # underlying AF-Multimer model identifier (not the rank).
                yield (pred_name, str(rank_idx), model_id,
                       struct_path, pae_path, conf_path, fmt)


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


def _find_esmfold2(filenames, basenames_map):
    """ESMFold2 (biohub API): single model — *_model.pdb + *_pae.json"""
    for name in filenames:
        base = os.path.basename(name)
        m = re.match(r'^(.+)_model\.pdb$', base)
        if not m:
            continue
        prefix = m.group(1)
        pae_path = basenames_map.get(f'{prefix}_pae.json')
        if not pae_path:
            continue
        # PAE JSON also carries ptm/iptm/plddt → use it as the confidence source too
        yield (prefix, '0', os.path.basename(name), name, pae_path, pae_path, 'pdb')


def _find_esmfold2_native(filenames, basenames_map):
    """Folder-per-pair (and optionally per-sample-within-pair).

    Naming is fully agnostic: a folder is a prediction set whenever it contains
    one or more structure files and one or more PAE files. Multi-sample matching
    uses the longest-common-prefix trick — whatever prefix is shared by all
    struct stems is treated as common, and the remaining suffix is the sample ID.
    The same is done for PAE stems and JSON stems. IDs that appear in both
    struct- and PAE-sets are paired. This handles cases like:

      complex_s1.pdb.gz + pae_s1.npz + metrics_s1.json       (suffix _s1)
      model_1.pdb + pae_1.npz + scores_1.json                (suffix _1)
      protein_pair_0.pdb + pae_0.npy + scores_0.json         (mixed prefixes)
      Fsn___SkpB_0.pdb + Fsn___SkpB_0.npz + Fsn___SkpB_0.json
      model1.pdb + pae1.npz                                  (no underscore)

    Accepted structure formats: .pdb, .cif, .pdb.gz, .cif.gz, .pdb.xz, .cif.xz
    Accepted PAE formats:       .npy, .npz   (.npy preferred — direct array)
    For each struct/PAE/json candidate, uncompressed wins over compressed when
    both exist for the same sample ID.
    """
    import re
    from collections import defaultdict
    by_dir = defaultdict(list)
    for f in filenames:
        by_dir[os.path.dirname(f)].append(f)

    STRUCT_EXT_RX = re.compile(r'\.(pdb|cif)(\.gz|\.xz)?$', re.IGNORECASE)
    PAE_EXT_RX = re.compile(r'\.(npy|npz)$', re.IGNORECASE)
    EXT_STRIP_RX = re.compile(r'\.(pdb|cif|npy|npz|json)(\.gz|\.xz)?$', re.IGNORECASE)

    def stem(path):
        return EXT_STRIP_RX.sub('', os.path.basename(path))

    def longest_common_prefix(strs):
        if not strs:
            return ''
        s1, s2 = min(strs), max(strs)
        i = 0
        while i < len(s1) and i < len(s2) and s1[i] == s2[i]:
            i += 1
        return s1[:i]

    def ids_by_lcp(files):
        """Return {sample_id: file} where sample_id = stem - common_prefix.
        With 1 file, sample_id is '' (treated as a single-model marker)."""
        stems = [stem(f) for f in files]
        if len(files) == 1:
            return {'': files[0]}
        pref = longest_common_prefix(stems)
        out = {}
        for f, st in zip(files, stems):
            sid = st[len(pref):].lstrip('_-.')  # trim common leading separators
            out[sid] = f
        return out

    def struct_priority_key(p):
        b = os.path.basename(p).lower()
        return 0 if not (b.endswith('.gz') or b.endswith('.xz')) else 1

    def pae_priority_key(p):
        return 0 if p.lower().endswith('.npy') else 1

    for d, files in by_dir.items():
        struct_files = [f for f in files if STRUCT_EXT_RX.search(os.path.basename(f))]
        pae_files = [f for f in files if PAE_EXT_RX.search(os.path.basename(f))]
        json_files = [f for f in files if os.path.basename(f).lower().endswith('.json')]
        if not struct_files or not pae_files:
            continue

        pred_name = os.path.basename(d) or 'prediction'

        # Pre-sort by priority so duplicates resolve cleanly
        struct_files = sorted(struct_files, key=struct_priority_key)
        pae_files = sorted(pae_files, key=pae_priority_key)

        # If only one struct or one PAE, treat as single-model
        if len(struct_files) == 1 or len(pae_files) == 1:
            sp = struct_files[0]
            pp = pae_files[0]
            jp = json_files[0] if json_files else pp
            fmt = 'cif' if '.cif' in os.path.basename(sp).lower() else 'pdb'
            yield (pred_name, '0', os.path.basename(sp), sp, pp, jp, fmt)
            continue

        # Multi-sample: derive sample IDs by stripping the longest common prefix
        # within each file-type group, then pair up matching IDs.
        struct_ids = ids_by_lcp(struct_files)
        pae_ids = ids_by_lcp(pae_files)
        json_ids = ids_by_lcp(json_files) if json_files else {}

        common = set(struct_ids) & set(pae_ids)
        if not common:
            # Could not split by ID — fall back to single pair
            sp = struct_files[0]
            pp = pae_files[0]
            jp = json_files[0] if json_files else pp
            fmt = 'cif' if '.cif' in os.path.basename(sp).lower() else 'pdb'
            yield (pred_name, '0', os.path.basename(sp), sp, pp, jp, fmt)
            continue

        for sid in sorted(common):
            sp = struct_ids[sid]
            pp = pae_ids[sid]
            jp = json_ids.get(sid) or (json_files[0] if json_files else pp)
            fmt = 'cif' if '.cif' in os.path.basename(sp).lower() else 'pdb'
            rank = sid or '0'
            yield (pred_name, rank, os.path.basename(sp), sp, pp, jp, fmt)


def _find_generic(filenames, basenames_map, read_fn=None):
    """Generic: any .pdb/.cif + a JSON/NPZ carrying PAE.

    Keeps the original name -> index -> positional pairing (so folders that already work
    are untouched), with two additions: obvious non-PAE JSONs (AF3 '*_data.json' input,
    ranking/summary/config files) are dropped from the candidate pool so they cannot
    shadow the real PAE file, and if the chosen file still has no readable PAE the
    remaining candidates are scanned for one that does. Both only ever rescue a case that
    used to fail; a previously-correct pick is left exactly as before.
    """
    struct_files = sorted(
        [f for f in filenames if f.endswith('.cif') or f.endswith('.pdb')],
        key=lambda f: os.path.basename(f)
    )
    all_json = [f for f in filenames if f.endswith('.json')]
    npz_files = [f for f in filenames if f.endswith('.npz')]

    def is_non_pae(f):
        b = os.path.basename(f).lower()
        return any(k in b for k in (
            '_data.json', 'ranking_scores', 'ranking_debug', 'config',
            '_summary_confidences.json', 'summary_confidence', 'metrics.json'))
    # Drop obvious non-PAE JSONs, but keep original order and never end up with nothing.
    json_files = [f for f in all_json if not is_non_pae(f)] or all_json

    for i, sf in enumerate(struct_files):
        struct_base = os.path.basename(sf)
        fmt = 'cif' if sf.endswith('.cif') else 'pdb'
        name_prefix = re.sub(r'\.(cif|pdb)$', '', struct_base)
        idx_match = re.search(r'(\d+)', name_prefix)
        idx = idx_match.group(1) if idx_match else None

        pae_source = None
        for jf in json_files:                       # (1) shares the structure's base name
            if name_prefix in os.path.basename(jf):
                pae_source = jf
                break
        if not pae_source and idx:                  # (2) shares the numeric index
            for jf in json_files:
                fb = os.path.basename(jf)
                if idx in fb and 'aggregated' not in fb:
                    pae_source = jf
                    break
        if not pae_source:                          # (3) positional
            pae_source = json_files[i] if i < len(json_files) else (json_files[0] if json_files else None)
        if not pae_source and npz_files:
            pae_source = npz_files[i] if i < len(npz_files) else npz_files[0]

        # Generosity: only if the chosen file has no readable PAE, look for one that does.
        if read_fn is not None and pae_source is not None and extract_pae(pae_source, read_fn) is None:
            for cf in json_files + npz_files:
                if cf == pae_source:
                    continue
                try:
                    if extract_pae(cf, read_fn) is not None:
                        pae_source = cf
                        break
                except Exception:
                    continue

        yield ('generic', str(i), struct_base, sf, pae_source, pae_source, fmt)


# ============================================================================
# Manifest — the explicit "declare your data" escape hatch
# ============================================================================

def load_manifest(filenames, read_fn, manifest_arg=None):
    """Load an optional lis.json manifest for layouts auto-detection doesn't recognise.

    Looked up from --manifest, else a lis.json / lis_manifest.json among the inputs.
    Returns the parsed dict (with '_dir' = the manifest's directory, for resolving relative
    paths) or None. Every field is optional:
        platform   force the finder (same as --platform)
        structure  glob for structure files            pae      glob for the PAE file
        pae_key    JSON key holding the PAE matrix      summary  glob for the scores file
        models     [ {name, structure, pae, summary} ] explicit per-prediction mapping
    """
    content, src = None, None
    if manifest_arg:
        try:
            content, src = open(manifest_arg).read(), manifest_arg
        except OSError:
            return None
    else:
        for f in filenames:
            if os.path.basename(f).lower() in ('lis.json', 'lis_manifest.json'):
                content, src = read_fn(f), f
                break
    if not isinstance(content, str):
        return None
    try:
        m = json.loads(content)
    except (json.JSONDecodeError, TypeError):
        return None
    if not isinstance(m, dict):
        return None
    m['_dir'] = os.path.dirname(src)
    return m


def manifest_has_layout(manifest):
    """True if the manifest says how to find models (not merely which platform)."""
    return bool(manifest) and bool(manifest.get('models') or manifest.get('pae'))


def _find_from_manifest(manifest, filenames):
    """Yield model tuples from a manifest. Two modes: an explicit 'models' list, or globs
    (structure/pae/summary) paired per structure. Carries pae_key as an 8th tuple element
    when declared, so extract_pae can read a non-standard PAE key."""
    import fnmatch
    base = manifest.get('_dir', '')
    pae_key = manifest.get('pae_key')

    def resolve(pattern):
        if not pattern:
            return []
        for cand in ([f'{base}/{pattern}'] if base else []) + [pattern]:
            if cand in filenames:
                return [cand]
        bn = os.path.basename(pattern)                       # glob on basename, in file order
        return [f for f in filenames if fnmatch.fnmatch(os.path.basename(f), bn)]

    def emit(name, rank, s, p, sm):
        fmt = 'cif' if s.endswith('.cif') else 'pdb'
        t = (name, str(rank), os.path.basename(s), s, p, sm or p, fmt)
        return t + (pae_key,) if pae_key else t

    if manifest.get('models'):
        for i, e in enumerate(manifest['models']):
            s = resolve(e.get('structure'))
            p = resolve(e.get('pae') or e.get('structure'))
            if not s or not p:
                continue
            sm = resolve(e.get('summary'))
            yield emit(e.get('name', f'model_{i}'), e.get('rank', i), s[0], p[0], sm[0] if sm else None)
        return

    structs = sorted(resolve(manifest.get('structure')) or
                     [f for f in filenames if f.endswith(('.cif', '.pdb'))], key=os.path.basename)
    pae_all = resolve(manifest.get('pae'))
    summ_all = resolve(manifest.get('summary'))

    def sibling(files, sdir, sname, i):
        if not files:
            return None
        same = [f for f in files if os.path.dirname(f) == sdir]
        for pool in (same, files):
            for f in pool:
                if sname in os.path.basename(f):
                    return f
        return (same or files)[i] if i < len(same or files) else (same or files)[0]

    for i, s in enumerate(structs):
        sdir = os.path.dirname(s)
        sname = re.sub(r'\.(cif|pdb)$', '', os.path.basename(s))
        yield emit(sname, i, s, sibling(pae_all, sdir, sname, i), sibling(summ_all, sdir, sname, i))


# ============================================================================
# PAE Extraction
# ============================================================================

def extract_pae(pae_source, read_fn, pae_key=None, allow_pickle=False):
    """Extract PAE matrix from a file. Returns 2D numpy array or None.

    pae_key, if given (from a lis.json manifest), is tried first so a non-standard JSON
    key can be read; otherwise the usual keys (pae / predicted_aligned_error / ...) apply.

    allow_pickle enables reading PAE from a .pkl/.pickle (raw AlphaFold2/-Multimer
    result_*.pkl). Off by default because unpickling executes arbitrary code; it is only
    reachable via an explicit manifest 'pae' entry plus the --allow-pickle flag.
    """
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

    # .pkl / .pickle (raw AlphaFold2 / AlphaFold-Multimer result_*.pkl, which stores the
    # PAE under 'predicted_aligned_error'). Unpickling runs arbitrary code, so this path is
    # gated behind allow_pickle (--allow-pickle) and only reached via an explicit manifest
    # 'pae' entry — auto-detection never selects a pickle.
    if basename.endswith(('.pkl', '.pickle')):
        if not allow_pickle:
            return None  # caller emits a --allow-pickle hint
        if isinstance(content, str):
            content = content.encode('latin-1')
        import pickle
        try:
            obj = pickle.loads(content)
        except Exception:
            return None
        # Usually the AF2 result dict; occasionally just the bare array.
        candidates = []
        if isinstance(obj, dict):
            if pae_key and pae_key in obj:
                candidates.append(obj[pae_key])
            for k in ('predicted_aligned_error', 'pae', 'pae_matrix',
                      'predicted_aligned_error_matrix'):
                if k in obj:
                    candidates.append(obj[k])
        else:
            candidates.append(obj)
        for v in candidates:
            try:
                arr = np.asarray(v, dtype=np.float32)
            except (ValueError, TypeError):
                continue
            if arr.ndim == 3:
                arr = arr[0]
            if arr.ndim == 2 and arr.shape[0] == arr.shape[1] and arr.shape[0] > 0:
                return arr
        return None

    # JSON file
    if not isinstance(content, str):
        return None
    try:
        data = json.loads(content)
    except json.JSONDecodeError:
        return None

    # Manifest-declared PAE key wins (lets a non-standard key be read).
    if pae_key and isinstance(data, dict) and pae_key in data:
        v = data[pae_key]
        if isinstance(v, list) and v:
            if isinstance(v[0], list):
                return np.array(v, dtype=np.float32)
            n = int(round(math.sqrt(len(v))))
            if n * n == len(v):
                return np.array(v, dtype=np.float32).reshape(n, n)

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

    # Protenix-v2: token_pair_pae (also exposes token_pair_pde + contact_probs)
    if 'token_pair_pae' in data and isinstance(data['token_pair_pae'], list):
        pae = data['token_pair_pae']
        if isinstance(pae[0], list):
            return np.array(pae, dtype=np.float32)

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
    # ESMFold2 native uses 'interface_ptm' as the global ipTM key
    if 'interface_ptm' in data and 'ipTM' not in scores:
        scores['ipTM'] = _unwrap(data['interface_ptm'])
    if 'chain_pair_iptm' in data:
        scores['chainPairIptm'] = data['chain_pair_iptm']
    # ESMFold2 native uses 'pair_chains_iptm' for the per-chain-pair matrix
    if 'pair_chains_iptm' in data and 'chainPairIptm' not in scores:
        scores['chainPairIptm'] = data['pair_chains_iptm']
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

    # AlphaFold confidence_json format (AlphaPulldown's confidence_<model_id>.json):
    # {"residueNumber": [..], "confidenceScore": [..], "confidenceCategory": [..]}
    # The per-residue confidenceScore is pLDDT — average it for a global pLDDT.
    if 'confidenceScore' in data and isinstance(data['confidenceScore'], list):
        plddts = data['confidenceScore']
        if plddts:
            scores.setdefault('pLDDT', sum(plddts) / len(plddts))

    # AF-Multimer ranking_debug.json composite (0.8·iPTM + 0.2·pTM); kept as
    # rankingScore — not aliased to ipTM to avoid mixing semantics.
    if 'ranking_confidence' in data:
        scores['rankingScore'] = _unwrap(data['ranking_confidence'])

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
    """Extract one Cb (Ca for GLY, P for nucleic) coordinate per polymer residue,
    plus one coordinate per non-polymer HETATM atom, from PDB text. See
    parse_cif_coords for why per-atom HETATM coords keep the PAE alignment."""
    residues = OrderedDict()
    het_idx = 0
    for line in pdb_text.split('\n'):
        if not line.startswith('ATOM') and not line.startswith('HETATM'):
            continue
        if len(line) < 54:
            continue
        atom_name = line[12:16].strip()
        comp_id = line[17:20].strip()
        chain = line[21:22].strip() or 'A'
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])

        if line.startswith('HETATM'):
            # Ion = one monatomic token; glycan/ligand = one token per atom.
            if comp_id in ION_NAMES:
                if f'{chain}:1' not in residues:
                    residues[f'{chain}:1'] = {'chain': chain, 'resnum': 1, 'x': x, 'y': y, 'z': z, 'has_p': False}
            else:
                residues[f'het:{chain}:{het_idx}'] = {'chain': chain, 'resnum': 1, 'x': x, 'y': y, 'z': z, 'has_p': 'P' in atom_name}
                het_idx += 1
            continue

        try:
            resnum = int(line[22:26].strip())
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


def parse_cif_coords(cif_text):
    """Extract one Cb (Ca for GLY, P for nucleic) coordinate per polymer residue,
    plus one coordinate per non-polymer HETATM atom, from mmCIF text.

    AlphaFold3 tokenizes every glycan/ligand atom separately, so the PAE matrix
    has one token per HETATM atom. Emitting one coordinate per atom keeps the
    contact map aligned with the PAE tokens; otherwise every chain after a glycan
    loses its contact-based metrics (cLIS -> iLIS, actifpTM)."""
    residues = OrderedDict()
    het_idx = 0
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
            # Boltz CIFs wrap the trailing pdbx_PDB_model_num field onto the next
            # line for some atoms (data line has 18 fields + a lone "1" on the
            # next line). Don't enforce strict column-count match — just gate
            # on the first token being ATOM/HETATM. get_col() returns '' for
            # any column past parts.length, which is harmless: we never read
            # the wrapped model_num column.
            if not parts or parts[0] not in ('ATOM', 'HETATM'):
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

            if group_pdb == 'HETATM':
                # Ion = one monatomic token; glycan/ligand = one token per atom.
                if comp_id in ION_NAMES:
                    ion_key = f'{chain}:1'
                    if ion_key not in residues:
                        residues[ion_key] = {'chain': chain, 'resnum': 1, 'x': x, 'y': y, 'z': z, 'has_p': False}
                else:
                    residues[f'het:{chain}:{het_idx}'] = {'chain': chain, 'resnum': 1, 'x': x, 'y': y, 'z': z, 'has_p': 'P' in atom_name}
                    het_idx += 1
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
            # Boltz CIFs wrap the trailing pdbx_PDB_model_num field onto the next
            # line for some atoms (data line has 18 fields + a lone "1" on the
            # next line). Don't enforce strict column-count match — just gate
            # on the first token being ATOM/HETATM. get_col() returns '' for
            # any column past parts.length, which is harmless: we never read
            # the wrapped model_num column.
            if not parts or parts[0] not in ('ATOM', 'HETATM'):
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
            elif group_pdb == 'HETATM':
                # Glycan / ligand (non-ion HETATM): one token per atom in AF3.
                chain_types[chain] = 'glycan'
                counted = True

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

    # Auto-scale 0–1 pLDDT to 0–100 (e.g. ESMFold2 native PDB stores pLDDT on 0–1).
    # Use max value as a heuristic: AlphaFold/ColabFold use 0–100, ESMFold uses 0–1.
    if bfactors:
        mx = max(bfactors.values())
        if 0 < mx <= 1.0:
            bfactors = {k: v * 100.0 for k, v in bfactors.items()}

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
    """Transform PAE to confidence scores. Asymmetric, per-direction.

    Per-direction LIS = 1 - mean(PAE<cutoff)/cutoff is recovered by leaving
    PAE[i,j] and PAE[j,i] independent. The (i,j) and (j,i) chain-pair entries
    are then averaged downstream in analyze_single_model's symmetrize step.
    """
    pae = np.asarray(pae, dtype=np.float64)
    transformed = np.zeros_like(pae)
    mask = pae < pae_cutoff
    transformed[mask] = 1.0 - pae[mask] / pae_cutoff
    return transformed


def calc_pae_chain_pair_iptm(pae, starts, ends):
    """PAE-based per-chain-pair iPTM matrix (used when the platform doesn't ship one).

    Matches AlphaFold3's chain_pair_iptm definition:
      For chain pair (a, b):
        - Off-diagonal: d0 = f(N_a + N_b); for each residue i in chain_a ∪ chain_b
          compute row average of TM-score over the OPPOSITE chain only; iPTM = max row.
        - Diagonal (a == b): d0 = f(N_a); compute pTM for that chain.
      d0 = 1.24·(max(N, 19) − 15)^(1/3) − 1.8.

    The exact AF formula uses E[1/(1+(PAE/d0)²)] over the saved PAE-bin distribution;
    we only have E[PAE], so f(E[PAE]) underestimates the true value by Jensen's
    inequality. Empirical bias vs AF3 ground truth (PRC2, MIS12C): mean ‑0.01..‑0.03,
    max |Δ| ≈ 0.06 — ranking-preserving and useful for screening.
    """
    pae = np.asarray(pae, dtype=np.float64)
    nc = len(starts)
    out = np.zeros((nc, nc), dtype=np.float64)
    for a in range(nc):
        sa, ea = int(starts[a]), int(ends[a])
        N_a = ea - sa
        if N_a <= 0:
            continue
        for b in range(nc):
            sb, eb = int(starts[b]), int(ends[b])
            N_b = eb - sb
            if N_b <= 0:
                continue
            if a == b:
                N_for_d0 = max(N_a, 19)
            else:
                N_for_d0 = max(N_a + N_b, 19)
            d0 = 1.24 * (N_for_d0 - 15) ** (1.0/3.0) - 1.8
            if d0 <= 0:
                continue
            d0sq = d0 * d0
            if a == b:
                block = pae[sa:ea, sa:ea]
                tm = 1.0 / (1.0 + (block * block) / d0sq)
                row_avgs = tm.mean(axis=1)
            else:
                block_ab = pae[sa:ea, sb:eb]
                tm_ab = 1.0 / (1.0 + (block_ab * block_ab) / d0sq)
                row_a = tm_ab.mean(axis=1)
                block_ba = pae[sb:eb, sa:ea]
                tm_ba = 1.0 / (1.0 + (block_ba * block_ba) / d0sq)
                row_b = tm_ba.mean(axis=1)
                row_avgs = np.concatenate([row_a, row_b])
            out[a, b] = float(row_avgs.max()) if row_avgs.size else 0.0
    return out


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

    # Build transformed confidence map (asymmetric PAE, per-direction)
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

    # PAE-based per-chain-pair iPTM fallback (Boltz2 / AlphaPulldown / Generic etc.
    # ship only a global scalar). For 2-chain complexes per-pair == global so we
    # skip; for 3+ chains the global value is misleading on every row.
    # Underbias ~0.01–0.07 vs ground truth (Jensen's), but ranking-preserving.
    if iptm_matrix is None and len(chain_names) >= 3:
        chain_ends = list(cum_sum)
        chain_starts = list(starts)
        iptm_matrix = calc_pae_chain_pair_iptm(pae[:n_total, :n_total], chain_starts, chain_ends)

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

            # Vectorized LIS/cLIS computation (transformed is asymmetric per-direction)
            t_block = transformed[si:ei, sj:ej]
            t_pos = t_block > 0

            lis_sum = float(t_block[t_pos].sum())
            lis_count_avg = int(t_pos.sum())

            # LIR: cell-based union — a residue is in LIR if either direction is confident
            # on any of its cells (paired across the block's row/col).
            t_block_rev = transformed[sj:ej, si:ei]
            either_pos = t_pos | (t_block_rev.T > 0)
            lir_i = set(np.where(either_pos.any(axis=1))[0] + 1)
            lir_j = set(np.where(either_pos.any(axis=0))[0] + 1)

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
                # cLIR: cell-based union + contact
                either_contact = either_pos[:c_ei-c_si, :c_ej-c_sj] & contact_block
                clir_i = set(np.where(either_contact.any(axis=1))[0] + 1)
                clir_j = set(np.where(either_contact.any(axis=0))[0] + 1)
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
                # val['LIA'] already covers both AB and BA directions, so no extra sum
                'LIA': val['LIA'],
                'cLIA': val['cLIA'],
                'lenI': val['lenI'], 'lenJ': val['lenJ'],
                # chain-preserved cell-based union — every LIR_A residue has a paired LIR_B
                'lirI': val['lirI'] | rv['lirJ'],
                'lirJ': val['lirJ'] | rv['lirI'],
                'clirI': val['clirI'] | rv['clirJ'],
                'clirJ': val['clirJ'] | rv['clirI'],
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


def _wmean_plddt(vi, vj, wi, wj):
    """Weight-averaged per-pair pLDDT from two per-chain values, weighted by residue count
    (chain length for overall pLDDT, interface-residue count for LIpLDDT/cLIpLDDT). This equals
    the mean over the pooled residues. Returns None if either value is missing or weights are 0."""
    if vi is None or vj is None:
        return None
    if (isinstance(vi, float) and math.isnan(vi)) or (isinstance(vj, float) and math.isnan(vj)):
        return None
    w = (wi or 0) + (wj or 0)
    if w <= 0:
        return None
    return (wi * vi + wj * vj) / w


def format_row(name, rank, struct_file, pair):
    """Format one CSV row from a pair dict."""
    def fmt_plddt(v):
        if v is None or (isinstance(v, float) and math.isnan(v)):
            return ''
        return f'{v:.1f}'

    model_num = _extract_model_num(struct_file)

    # Per-pair pLDDT: length-weighted mean of the two chains' pLDDT; the interface variants are
    # weighted by interface-residue count (LIR/cLIR) rather than chain length.
    plddt_pair = _wmean_plddt(pair.get('pLDDT_i'), pair.get('pLDDT_j'), pair['lenI'], pair['lenJ'])
    liplddt_pair = _wmean_plddt(pair.get('LIpLDDT_i'), pair.get('LIpLDDT_j'), len(pair['lirI']), len(pair['lirJ']))
    cliplddt_pair = _wmean_plddt(pair.get('cLIpLDDT_i'), pair.get('cLIpLDDT_j'), len(pair['clirI']), len(pair['clirJ']))

    row = [
        name, rank, model_num, pair['ci'], pair['cj'],
        f"{pair['iLIS']:.4f}", f"{pair['iLIA']:.1f}", f"{pair['iLISA']:.1f}",
        f"{pair['ipSAE']:.4f}", f"{pair.get('actifpTM', 0):.4f}",
        f"{pair['LIS']:.4f}", f"{pair['cLIS']:.4f}",
        f"{pair['LIA']:.1f}", f"{pair['cLIA']:.1f}",
        f"{pair['ipTM']:.4f}",
        fmt_plddt(pair.get('pLDDT_i')),
        fmt_plddt(pair.get('pLDDT_j')),
        fmt_plddt(plddt_pair),
        f"{pair['pTM']:.3f}" if pair.get('pTM') is not None else '',
        str(len(pair['lirI'])), str(len(pair['lirJ'])),
        str(len(pair['clirI'])), str(len(pair['clirJ'])),
        fmt_plddt(pair.get('LIpLDDT_i')),
        fmt_plddt(pair.get('LIpLDDT_j')),
        fmt_plddt(liplddt_pair),
        fmt_plddt(pair.get('cLIpLDDT_i')),
        fmt_plddt(pair.get('cLIpLDDT_j')),
        fmt_plddt(cliplddt_pair),
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

def _extract_alphapulldown_scalars(struct_path, model_label, read_fn):
    """Read exact per-model iPTM and pTM from a sibling ranking_debug.json.

    AlphaPulldown writes (for multimer mode):
        {"iptm+ptm": {model_id: composite},   # 0.8·iPTM + 0.2·pTM
         "order":    [model_id, ...],
         "iptm":     {model_id: iptm}}        # raw per-model iPTM

    AF-Multimer's ranking score is composite = 0.8·iPTM + 0.2·pTM, so
        pTM = 5·composite − 4·iPTM  (algebraically exact, not an approximation).

    Returns a dict with optional keys: ipTM, pTM, rankingScore. Empty on miss.
    """
    ranking_path = os.path.join(os.path.dirname(struct_path), 'ranking_debug.json')
    content = read_fn(ranking_path)
    if not content or not isinstance(content, str):
        return {}
    try:
        rd = json.loads(content)
    except json.JSONDecodeError:
        return {}

    # Pull the AF-Multimer model_id out of any label form we yield
    # (e.g. 'unrelaxed_model_3_multimer_v3_pred_0.pdb' → 'model_3_multimer_v3_pred_0').
    m = re.search(r'model_\d+(?:_[A-Za-z0-9]+)*', model_label)
    if not m:
        return {}
    model_id = m.group(0)

    out = {}
    iptm_dict = rd.get('iptm', {}) or {}
    composite_dict = rd.get('iptm+ptm', {}) or {}

    if model_id in iptm_dict:
        iptm = float(iptm_dict[model_id])
        out['ipTM'] = iptm
        if model_id in composite_dict:
            composite = float(composite_dict[model_id])
            out['rankingScore'] = composite
            ptm = 5.0 * composite - 4.0 * iptm
            # Tiny floating-point excursions outside [0, 1] are OK to clamp.
            if -0.005 <= ptm <= 1.005:
                out['pTM'] = max(0.0, min(1.0, ptm))
    return out


def _do_process(model_tuple, read_fn, detected, pae_cutoff, cb_cutoff, verbose=False, allow_pickle=False):
    """Process a single model. Returns (name, rank, rows, error_msg) tuple.

    Any unexpected exception (corrupt gz, truncated file, unreadable JSON, etc.)
    is captured and returned as err_msg so one bad prediction doesn't kill the
    whole batch.
    """
    name, rank, model_label, struct_path, pae_path, scores_path, fmt = model_tuple[:7]
    pae_key = model_tuple[7] if len(model_tuple) > 7 else None   # optional, from a lis.json manifest

    try:
        struct_text = read_fn(struct_path)
        if not struct_text or not isinstance(struct_text, str):
            return name, rank, None, f'structure file unreadable: {struct_path}'

        pae = extract_pae(pae_path, read_fn, pae_key, allow_pickle=allow_pickle)
        if pae is None:
            if pae_path and os.path.basename(pae_path).endswith(('.pkl', '.pickle')) and not allow_pickle:
                return name, rank, None, (f'PAE is a pickle ({os.path.basename(pae_path)}); re-run with '
                                          f'--allow-pickle to load it (only for files you trust).')
            return name, rank, None, (f'PAE not found or unreadable: {pae_path}. If this layout is '
                                      f'unrecognised, declare it in a lis.json manifest (pae / pae_key) or pass --manifest.')
        pae = np.nan_to_num(pae, nan=31.0)

        scores = extract_confidence_scores(scores_path, read_fn)
        if scores_path != pae_path and pae_path:
            full_scores = extract_confidence_scores(pae_path, read_fn)
            for k, v in full_scores.items():
                if k not in scores:
                    scores[k] = v

        # AlphaPulldown: pull exact per-model iPTM and recover pTM from
        # ranking_debug.json (composite = 0.8·iPTM + 0.2·pTM → pTM = 5·comp − 4·iPTM).
        if detected == 'alphapulldown':
            ap_scalars = _extract_alphapulldown_scalars(struct_path, model_label, read_fn)
            for k, v in ap_scalars.items():
                if k not in scores or scores.get(k) is None:
                    scores[k] = v

        try:
            pairs = analyze_single_model(
                struct_text, pae, scores, fmt, detected,
                pae_path, read_fn, pae_cutoff, cb_cutoff)
        except Exception as e:
            return name, rank, None, f'analysis error: {e}'

        rows = [format_row(name, rank, model_label, pair) for pair in pairs]
        return name, rank, rows, None
    except EOFError as e:
        return name, rank, None, f'corrupt gz (EOFError): {pae_path}: {e}'
    except (OSError, zipfile.BadZipFile) as e:
        return name, rank, None, f'I/O error reading {pae_path}: {type(e).__name__}: {e}'
    except Exception as e:
        return name, rank, None, f'unexpected error on {pae_path}: {type(e).__name__}: {e}'


def _process_one_sequential(model_tuple, read_fn, detected, pae_cutoff, cb_cutoff, verbose=False, allow_pickle=False):
    """Sequential wrapper."""
    return _do_process(model_tuple, read_fn, detected, pae_cutoff, cb_cutoff, verbose, allow_pickle)


def _mp_worker(args):
    """Multiprocessing worker — creates its own read_fn from file_map."""
    model_tuple, detected, pae_cutoff, cb_cutoff, file_map, verbose, allow_pickle = args
    read_fn = _make_dir_reader(file_map)
    return _do_process(model_tuple, read_fn, detected, pae_cutoff, cb_cutoff, verbose, allow_pickle)


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
    'alphapulldown': 'AlphaPulldown (AF-Multimer)',
    'colabfold': 'ColabFold',
    'boltz': 'Boltz',
    'chai': 'Chai-1',
    'openfold3': 'OpenFold3',
    'generic': 'Generic',
}


def run(path, output=None, output_dir=None, pae_cutoff=12, cb_cutoff=8,
        platform=None, skip_existing=True, workers=1, verbose=False, manifest_path=None,
        allow_pickle=False):
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

    # An optional lis.json manifest is the explicit override for layouts detection cannot
    # parse. It may declare the models outright (or via globs), or just force the platform.
    manifest = load_manifest(filenames, read_fn, manifest_path)
    if manifest and manifest_has_layout(manifest):
        detected = manifest.get('platform') or 'manifest'
        print(f'[LIS] Platform: {PLATFORM_LABELS.get(detected, detected)} (from lis.json manifest)')
        models = list(_find_from_manifest(manifest, filenames))
    else:
        if manifest and manifest.get('platform') and not platform:
            platform = manifest['platform']
        if platform:
            detected = platform
        else:
            detected = detect_platform(filenames, read_fn)
        print(f'[LIS] Platform: {PLATFORM_LABELS.get(detected, detected)}')
        models = list(find_models(filenames, detected, read_fn))
    # Filter macOS resource forks
    models = [m for m in models if not m[0].startswith('._')]

    if not models:
        print('[LIS] ERROR: No prediction models found', file=sys.stderr)
        if not (manifest and manifest_has_layout(manifest)):
            print('[LIS]   Unrecognised layout? Add a lis.json manifest (or pass --manifest)\n'
                  '[LIS]   declaring "structure"/"pae"/"pae_key", or an explicit "models" list.', file=sys.stderr)
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
        if not ok and err_msg:
            # Always surface failures so the bad-file path lands in the log,
            # not just under --verbose. Goes to stderr so it's visible even
            # when stdout is captured by a progress collector.
            print(f'\n[LIS] FAIL {name} rank={rank}: {err_msg}', file=sys.stderr, flush=True)

    # Check if input is a folder (needed for multiprocessing — zip read_fn can't be pickled)
    is_folder = os.path.isdir(path)

    if workers > 1 and is_folder:
        from multiprocessing import Pool

        # file_map is picklable (dict of str → str for folders)
        task_args = [
            (m, detected, pae_cutoff, cb_cutoff, file_map, verbose, allow_pickle)
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
                model_tuple, read_fn, detected, pae_cutoff, cb_cutoff, verbose, allow_pickle)
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
  AlphaFold3     *_model_N.cif + *_full_data_N.json + *_summary_confidences_N.json
  ColabFold      *_unrelaxed_rank_N*.pdb + *_scores_rank_N*.json
  AlphaPulldown  ranked_N.pdb / unrelaxed_model_N_*.pdb + pae_model_N_*.json
                 + confidence_model_N_*.json (+ ranking_debug.json)
  Boltz          *.pdb/.cif + confidence_*.json + pae_*.npz
  Chai-1         pred.rank_N.cif + scores.model_idx_N.json + pae.*.npy/.npz
  OpenFold3      result_sample_N_model.pdb + result_sample_N_confidences*.json
  Generic        any .pdb/.cif + PAE .json

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
                        choices=['alphafold3', 'alphapulldown', 'colabfold', 'boltz', 'chai', 'openfold3', 'generic'],
                        help='Force platform detection (default: auto-detect)')
    parser.add_argument('--manifest', default=None,
                        help='Path to a lis.json manifest declaring the data (structure/pae/pae_key/summary '
                             'globs, or an explicit models list) — an override for unrecognised layouts')
    parser.add_argument('--workers', '-w', type=int, default=1,
                        help='Number of parallel workers (default: 1)')
    parser.add_argument('--no-skip-existing', action='store_true',
                        help='Reprocess all predictions even if already in the output CSV')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Show error details for failed predictions')
    parser.add_argument('--allow-pickle', action='store_true',
                        help='Permit reading PAE from a .pkl/.pickle file declared in a manifest '
                             '(e.g. a raw AlphaFold2/-Multimer result_*.pkl). Unpickling executes '
                             'arbitrary code — only use on files you trust.')
    args = parser.parse_args()

    if not os.path.exists(args.path):
        print(f'[LIS] ERROR: Path not found: {args.path}', file=sys.stderr)
        sys.exit(1)

    run(args.path, output=args.output, output_dir=args.output_dir,
        pae_cutoff=args.pae_cutoff, cb_cutoff=args.cb_cutoff,
        platform=args.platform, manifest_path=args.manifest,
        skip_existing=not args.no_skip_existing,
        workers=args.workers, verbose=args.verbose, allow_pickle=args.allow_pickle)


if __name__ == '__main__':
    main()
