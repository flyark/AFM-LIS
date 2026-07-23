#!/usr/bin/env python3
"""Regression tests for the official *local* AlphaFold 3 output (issue #14) and the
generic PAE-source generosity that came with it.

Run:  python tests/test_local_af3.py      (exit 0 = all pass, 1 = a check failed)
"""
import os, sys, csv, json, tempfile, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
import lis

failures = []
def check(cond, msg):
    print(('  ok   ' if cond else '  FAIL ') + msg)
    if not cond:
        failures.append(msg)

# ---- 1. Detection + finder are filename-only (no PAE read needed) --------------------
LOCAL = [
    'job/job_model.cif', 'job/job_confidences.json', 'job/job_summary_confidences.json',
    'job/job_data.json', 'job/job_ranking_scores.csv',
    'job/seed-1_sample-0/job_seed-1_sample-0_model.cif',
    'job/seed-1_sample-0/job_seed-1_sample-0_confidences.json',
    'job/seed-1_sample-0/job_seed-1_sample-0_summary_confidences.json',
    'job/seed-1_sample-1/job_seed-1_sample-1_model.cif',
    'job/seed-1_sample-1/job_seed-1_sample-1_confidences.json',
    'job/seed-1_sample-1/job_seed-1_sample-1_summary_confidences.json',
]
none = lambda p: None
check(lis.detect_platform(LOCAL, none) == 'alphafold3', 'local AF3 layout detects as alphafold3')
models = list(lis.find_models(LOCAL, 'alphafold3', none))
check(len(models) == 2, 'per-sample ensemble yields one model per seed_sample (got %d)' % len(models))
check(all(m[4].endswith('_confidences.json') and not m[4].endswith('_data.json') for m in models),
      "PAE source is *_confidences.json, never the *_data.json input")

# ---- 2. Regression: existing formats must still detect as before ----------------------
check(lis.detect_platform(['x_model_0.cif', 'x_full_data_0.json', 'x_summary_confidences_0.json'], none)
      == 'alphafold3', 'indexed (standard) AF3 still detects as alphafold3')
check(lis.detect_platform(['a_unrelaxed_rank_001.pdb', 'a_scores_rank_001.json'], none)
      == 'colabfold', 'ColabFold still detects as colabfold')
check(lis.detect_platform(['complex.pdb', 'x.json'], none) == 'generic', 'plain cif+json still generic')

# ---- 3. Generic generosity: a decoy *_data.json must not shadow the real PAE json -----
gj = [{'predicted_aligned_error': [[0.5, 5.0], [5.0, 0.5]]}]  # PAE-bearing
files = {'p.cif': 'x', 'p_data.json': json.dumps({'name': 'p'}), 'p_scores.json': json.dumps(gj)}
def rd(path):
    return files.get(os.path.basename(path))
gmodels = list(lis.find_models(list(files), 'generic', rd))
check(gmodels and os.path.basename(gmodels[0][4]) == 'p_scores.json',
      'generic picks the PAE json (p_scores.json), not the *_data.json input')

# ---- 4. End-to-end: a real minimal local-AF3 dimer produces a valid iLIS --------------
tmp = tempfile.mkdtemp()
subprocess.run([sys.executable, os.path.join(HERE, 'make_local_af3_fixture.py'), tmp],
               check=True, capture_output=True)
job = os.path.join(tmp, 'dimer_A_B')
subprocess.run([sys.executable, os.path.join(ROOT, 'lis.py'), job], check=True, capture_output=True)
csv_path = os.path.join(job, 'dimer_A_B_lis_analysis.csv')
rows = list(csv.DictReader(open(csv_path))) if os.path.exists(csv_path) else []
check(bool(rows) and all(float(r['iLIS']) > 0 for r in rows),
      'end-to-end local-AF3 run produces a positive iLIS (%s)' % [r.get('iLIS') for r in rows])

# ---- 5. Manifest: the explicit "declare your data" escape hatch ----------------------
MFN = ['d/x.cif', 'd/x_ae.json', 'd/x_data.json']
glob_mani = {'structure': '*.cif', 'pae': '*_ae.json', 'pae_key': 'my_pae', '_dir': 'd'}
mm = list(lis._find_from_manifest(glob_mani, MFN))
check(len(mm) == 1 and os.path.basename(mm[0][4]) == 'x_ae.json' and mm[0][-1] == 'my_pae',
      'manifest glob mode pairs the declared PAE file and carries pae_key')
expl_mani = {'models': [{'name': 'p', 'structure': 'd/x.cif', 'pae': 'd/x_ae.json'}], '_dir': 'd'}
mm2 = list(lis._find_from_manifest(expl_mani, MFN))
check(bool(mm2) and mm2[0][0] == 'p' and mm2[0][3] == 'd/x.cif' and mm2[0][4] == 'd/x_ae.json',
      'manifest explicit "models" maps structure/pae directly')
check(lis.manifest_has_layout(glob_mani) and not lis.manifest_has_layout({'platform': 'boltz'}),
      'manifest_has_layout distinguishes a layout from a platform-only hint')
pk_files = {'x_ae.json': json.dumps({'my_pae': [[0.5, 5.0], [5.0, 0.5]]})}
pae = lis.extract_pae('x_ae.json', lambda p: pk_files.get(os.path.basename(p)), 'my_pae')
check(pae is not None and tuple(pae.shape) == (2, 2), 'extract_pae reads a manifest-declared pae_key')

# ---- 6. .pkl PAE (raw AF2/-Multimer result_*.pkl): gated behind --allow-pickle ---------
import pickle as _pickle
pdir = tempfile.mkdtemp()
# minimal 2-chain PDB, close interface so iLIS > 0 (chain A CB y=1.5, chain B CB y=8.5 -> 7 A)
_lines, _aid = [], 0
for _ch in ('A', 'B'):
    for _res in range(1, 4):
        for _atom in ('N', 'CA', 'C', 'O', 'CB'):
            _aid += 1
            _y = (0.0 if _ch == 'A' else 5.0) + {'N': -0.5, 'CA': 0.0, 'C': 0.5, 'O': 1.0,
                                                 'CB': 1.5 if _ch == 'A' else 3.5}[_atom]
            _lines.append('ATOM  %5d  %-3s ALA %s%4d    %8.3f%8.3f%8.3f  1.00 90.00           %s'
                          % (_aid, _atom, _ch, _res, _res * 4.0, _y, 0.0, _atom[0]))
open(os.path.join(pdir, 'complex.pdb'), 'w').write('\n'.join(_lines) + '\nEND\n')
# AF2-style result pickle: PAE under predicted_aligned_error (like the real result dict)
_tok = [c for c in 'AB' for _ in range(3)]
_pae = [[0.5 if i == j else (2.0 if _tok[i] == _tok[j] else 5.0) for j in range(6)] for i in range(6)]
_pkl = os.path.join(pdir, 'result_model_1_multimer_v3_pred_0.pkl')
with open(_pkl, 'wb') as _f:
    _pickle.dump({'predicted_aligned_error': _pae, 'ptm': 0.8, 'iptm': 0.85}, _f)
open(os.path.join(pdir, 'lis.json'), 'w').write(
    json.dumps({'structure': '*.pdb', 'pae': 'result_*.pkl', 'pae_key': 'predicted_aligned_error'}))

# unit: extract_pae refuses a .pkl unless allow_pickle=True
_rd = lambda p: (open(p, 'rb').read() if p and p.endswith(('.pkl', '.pickle')) else open(p).read()) if p else None
check(lis.extract_pae(_pkl, _rd, 'predicted_aligned_error', allow_pickle=False) is None,
      'extract_pae refuses a .pkl unless allow_pickle=True')
_arr = lis.extract_pae(_pkl, _rd, 'predicted_aligned_error', allow_pickle=True)
check(_arr is not None and tuple(_arr.shape) == (6, 6),
      'extract_pae reads predicted_aligned_error from a .pkl with allow_pickle=True')

# end-to-end: without the flag it fails with a --allow-pickle hint; with it, positive iLIS
_r1 = subprocess.run([sys.executable, os.path.join(ROOT, 'lis.py'), pdir, '--no-skip-existing'],
                     capture_output=True, text=True)
check('--allow-pickle' in (_r1.stderr + _r1.stdout),
      'without --allow-pickle the .pkl fails with a hint to --allow-pickle')
subprocess.run([sys.executable, os.path.join(ROOT, 'lis.py'), pdir, '--allow-pickle', '--no-skip-existing'],
               check=True, capture_output=True)
_csv = os.path.join(pdir, os.path.basename(pdir) + '_lis_analysis.csv')
_rows = list(csv.DictReader(open(_csv))) if os.path.exists(_csv) else []
check(bool(_rows) and all(float(r['iLIS']) > 0 for r in _rows),
      'end-to-end .pkl PAE with --allow-pickle produces a positive iLIS (%s)' % [r.get('iLIS') for r in _rows])

print(('\nALL PASS' if not failures else '\nFAILED: %d check(s)' % len(failures)))
sys.exit(1 if failures else 0)
