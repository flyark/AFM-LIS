#!/usr/bin/env python3
"""Generate a minimal but valid *official local AlphaFold 3* output for a 2-chain dimer,
used as a regression fixture for the local-AF3 detection fix (issue #14). Two chains of
3 residues each; chains placed so several Cbeta pairs are < 8 A (cLIR) and inter-chain
PAE is low (LIR), so lis.py produces a non-trivial iLIS.

Usage: python tests/make_local_af3_fixture.py [outdir]   (default: tests/fixtures)
"""
import os, sys, json

OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), 'fixtures')
JOB = 'dimer_A_B'
job_dir = os.path.join(OUT, JOB)
os.makedirs(job_dir, exist_ok=True)

# 2 chains x 3 residues (ALA). Chain A near y=0, chain B near y=5 -> interface Cbeta < 8 A.
CHAINS = ['A', 'B']
NRES = 3
def coords(chain, res, atom):
    x = res * 4.0
    y = 0.0 if chain == 'A' else 5.0
    dy = {'N': -0.5, 'CA': 0.0, 'C': 0.5, 'O': 1.0, 'CB': 1.5 if chain == 'A' else 3.5}[atom]
    return x, y + dy, 0.0

def make_cif():
    hdr = ['group_PDB', 'id', 'type_symbol', 'label_atom_id', 'label_alt_id', 'label_comp_id',
           'label_asym_id', 'label_entity_id', 'label_seq_id', 'pdbx_PDB_ins_code',
           'Cartn_x', 'Cartn_y', 'Cartn_z', 'occupancy', 'B_iso_or_equiv',
           'auth_asym_id', 'auth_seq_id', 'pdbx_PDB_model_num']
    lines = ['data_%s' % JOB, '#', 'loop_'] + ['_atom_site.%s' % h for h in hdr]
    aid = 0
    for ci, ch in enumerate(CHAINS, start=1):
        for res in range(1, NRES + 1):
            for atom in ['N', 'CA', 'C', 'O', 'CB']:
                aid += 1
                x, y, z = coords(ch, res, atom)
                el = atom[0]
                lines.append('ATOM %d %s %s . ALA %s %d %d ? %.3f %.3f %.3f 1.00 90.00 %s %d 1'
                             % (aid, el, atom, ch, ci, res, x, y, z, ch, res))
    lines.append('#')
    return '\n'.join(lines) + '\n'

def make_confidences():
    n = len(CHAINS) * NRES  # 6 tokens
    tok_chain = [c for c in CHAINS for _ in range(NRES)]
    tok_res = [r for _ in CHAINS for r in range(1, NRES + 1)]
    pae = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            pae[i][j] = 0.5 if i == j else (2.0 if tok_chain[i] == tok_chain[j] else 5.0)  # all <=12
    return {
        'pae': pae,
        'token_chain_ids': tok_chain,
        'token_res_ids': tok_res,
        'atom_plddts': [90.0] * (n * 5),
        'contact_probs': [[1.0 if i != j else 0.0 for j in range(n)] for i in range(n)],
    }

def make_summary():
    return {
        'iptm': 0.85, 'ptm': 0.80, 'ranking_score': 0.85,
        'chain_iptm': [0.85, 0.85], 'chain_ptm': [0.80, 0.80],
        'chain_pair_iptm': [[0.80, 0.85], [0.85, 0.80]],
        'chain_pair_pae_min': [[0.5, 5.0], [5.0, 0.5]],
        'fraction_disordered': 0.0, 'has_clash': 0.0, 'num_recycles': 10,
    }

def make_data():  # the AF3 *input* json — deliberately carries NO PAE
    return {'name': JOB, 'modelSeeds': [1], 'dialect': 'alphafold3', 'version': 2,
            'sequences': [{'protein': {'id': 'A', 'sequence': 'AAA'}},
                          {'protein': {'id': 'B', 'sequence': 'AAA'}}]}

cif, conf, summ, data = make_cif(), make_confidences(), make_summary(), make_data()
# top-level (best model) + two seed-sample subdirs
w = lambda p, s: open(p, 'w').write(s if isinstance(s, str) else json.dumps(s))
w(os.path.join(job_dir, f'{JOB}_model.cif'), cif)
w(os.path.join(job_dir, f'{JOB}_confidences.json'), conf)
w(os.path.join(job_dir, f'{JOB}_summary_confidences.json'), summ)
w(os.path.join(job_dir, f'{JOB}_data.json'), data)
w(os.path.join(job_dir, f'{JOB}_ranking_scores.csv'), 'seed,sample,ranking_score\n1,0,0.85\n1,1,0.83\n')
for s in (0, 1):
    d = os.path.join(job_dir, f'seed-1_sample-{s}')
    os.makedirs(d, exist_ok=True)
    pre = f'{JOB}_seed-1_sample-{s}'
    w(os.path.join(d, f'{pre}_model.cif'), cif)
    w(os.path.join(d, f'{pre}_confidences.json'), conf)
    w(os.path.join(d, f'{pre}_summary_confidences.json'), summ)
print('wrote fixture ->', job_dir)
