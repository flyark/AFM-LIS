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

print(('\nALL PASS' if not failures else '\nFAILED: %d check(s)' % len(failures)))
sys.exit(1 if failures else 0)
