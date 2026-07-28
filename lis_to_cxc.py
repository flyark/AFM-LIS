#!/usr/bin/env python3
"""
lis_to_cxc.py — Generate ChimeraX .cxc visualization scripts directly from
lis.py output (`*_lis_analysis.csv`).

Works with any predictor lis.py supports — ColabFold (AlphaFold-Multimer),
AlphaFold 3, Boltz-2, Chai-1, OpenFold3 — and handles multi-chain folds
where AF3-class predictors produce N chains per fold.

Input :  one or more *_lis_analysis.csv files + a folder with
         the .cif / .pdb model files referenced in `structure_file`.
         Files may live in nested subfolders under --pdb-root; lookup is
         recursive.

Outputs (written to --out/, created if missing):
  <fold>__rank<N>_model<M>_interface.cxc   one per (fold, rank, model)
                                            — open in ChimeraX (double-click
                                            or `open <path>` from CLI)
  summary.tsv  (only when --summary passed)  one row per chain-pair-in-cxc;
                                              ~21 columns + any --metadata
                                              passthrough columns. Default:
                                              not emitted — the lis.py CSV
                                              is usually enough for
                                              downstream pandas work.

Expected CSV columns (lis.py current schema):
  name (str)         fold identifier, e.g. "CCNQ_HUMAN___PITX1_HUMAN"
  rank, model (str)  per-model rank index and model index
  chain_i, chain_j   single-letter chain IDs (A, B, ...)
  iLIS, LIS, cLIS    interaction scores (float, 0..1)
  ipTM               AlphaFold ipTM (float)
  LIR_indices_i/j    residue indices in compact form "[16-20,22-23,...]"
                     or numpy-array form "[ 16 17 18 ...]"
  cLIR_indices_i/j   contact-residue indices (same format)
  structure_file     basename of .cif/.pdb model; can live nested under
                     --pdb-root (recursive search)
  len_i, len_j       chain lengths (optional; goes into manifest)
  pLDDT_i, pLDDT_j   per-chain mean pLDDT (optional; goes into manifest)

  OLD-schema *_lis_all_ranks.csv files are auto-detected; iLIS is
  synthesized as sqrt(LIS * cLIS) on interface scores since OLD schema
  has no iLIS column.

Multi-chain handling (AF3 / Boltz / Chai etc.):
  A 3-chain fold has rows for A-B, A-C, B-C — three chain pairs per model.
  We emit a single cxc per (fold, rank, model) that colors LIR/cLIR for
  every chain pair simultaneously, so opening the cxc loads the structure
  once and shows every interface highlighted. Coloring is per-chain (LIVIA
  convention): chain A has one color across every pair it appears in.

iLIS cutoff (--ilis-cutoff, default 0.223 = FPR 10% on the ColabFold
calibration):
  Chain-pair rows with iLIS below the cutoff are dropped. A (fold, rank,
  model) where every pair is below the cutoff is skipped entirely — no
  cxc emitted, structure not opened in ChimeraX. Pass --ilis-cutoff 0 to
  disable filtering and emit every pair lis.py produced. NOTE: thresholds
  are calibrated on ColabFold (AF2-multimer) Y2H reference sets; AF3,
  Boltz-2, Chai-1, OpenFold3 may need different thresholds.

CLI examples:
    # Smallest invocation: one CSV, structures next to it, default cutoff
    python lis_to_cxc.py --csv X_lis_analysis.csv --pdb-root . --out cxc/

    # Batch: a directory of *_lis_analysis.csv files (e.g. ColabFold chunks)
    python lis_to_cxc.py --csv-dir lis_results/ --pdb-root structures/ --out cxc/

    # Cross-platform: loosen cutoff (the default 0.223 is ColabFold-calibrated;
    # AF3 / Chai-1 / Boltz-2 / OpenFold3 may score lower for the same biology)
    python lis_to_cxc.py --csv X.csv --pdb-root . --out cxc/ --ilis-cutoff 0.15

    # Disable the cutoff entirely (emit every pair lis.py produced)
    python lis_to_cxc.py --csv X.csv --pdb-root . --out cxc/ --ilis-cutoff 0

    # Selective: only specific folds of interest (exact match on `name`)
    python lis_to_cxc.py --csv X.csv --pdb-root . --out cxc/ \\
        --pairs "FOO_HUMAN___BAR_HUMAN,BAZ_HUMAN___QUX_HUMAN"
    python lis_to_cxc.py --csv X.csv --pdb-root . --out cxc/ \\
        --pairs-file shortlist.txt    # one fold name per line; # = comment

    # Top-N folds by best iLIS, with symbol metadata attached
    python lis_to_cxc.py --csv-dir lis_results/ --pdb-root structures/ \\
        --out cxc_top100/ --top-n 100 --metadata pair_symbols.tsv

    # Visual customization
    python lis_to_cxc.py --csv X.csv --pdb-root . --out cxc/ --single-color
    python lis_to_cxc.py --csv X.csv --pdb-root . --out cxc/ \\
        --palette purple-gold
    python lis_to_cxc.py --csv X.csv --pdb-root . --out cxc/ \\
        --colors "#aa0000,#0066cc,#2ca02c"
    python lis_to_cxc.py --csv X.csv --pdb-root . --out cxc/ --gap-fill 0

    # Fully independent LIR and cLIR colors per chain (no lighten derivation):
    python lis_to_cxc.py --csv X.csv --pdb-root . --out cxc/ \\
        --clir-colors "#aa0000,#0066cc" \\
        --lir-colors  "#ffaa00,#88ccff"

    # Custom cLIR but keep --lir-lighten derivation for LIR:
    python lis_to_cxc.py --csv X.csv --pdb-root . --out cxc/ \\
        --clir-colors "#aa0000,#0066cc"

    # Inspect without writing (verifies inputs resolve, structures found)
    python lis_to_cxc.py --csv X.csv --pdb-root . --out cxc/ --dry-run

    # Discover built-in palettes
    python lis_to_cxc.py --list-palettes

--metadata TSV format (tab-separated, first column must be `name`):
    name\\tbait_symbol\\tprey_symbol\\tbait_uniprot\\tprey_uniprot\\torganism
    FOO_HUMAN___BAR_HUMAN\\tFOO\\tBAR\\tP12345\\tQ67890\\thuman
  Remaining columns are passed through verbatim into both the cxc header
  (as a `# Metadata:` block) and the manifest (as appended columns).

Default color palette (matches the LIVIA web tool):
    LIR chain i:  #80CBC4 (light teal)   cLIR chain i: #00897B (dark teal)
    LIR chain j:  #FFAB91 (light orange) cLIR chain j: #E64A19 (dark red)

For 3+ chain folds, additional chains cycle through tab10 colors.
Run `python lis_to_cxc.py --list-palettes` for all built-in palettes.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from pathlib import Path
from typing import Iterable

# ----------------------------------------------------------------------
# Residue-index parsing
# ----------------------------------------------------------------------

# lis.py NEW schema encodes residue indices as bracket-wrapped, comma-
# separated mix of single positions and ranges, e.g. "[16-20,22-23,25-32,34]".
# Older OLD-schema CSVs encode as numpy-array-string "[ 16 17 18 ...]".

def parse_residue_indices(raw: str) -> list[int]:
    """Parse either NEW-schema range format or OLD-schema array format."""
    if not isinstance(raw, str) or not raw.strip() or raw.strip() in ("[]", "nan"):
        return []
    inner = raw.strip().strip("[]").replace("\n", " ").strip()
    if not inner:
        return []
    out: list[int] = []
    # NEW schema: tokens separated by commas; each token is either an int
    # or a range "a-b".
    if "," in inner or "-" in inner and " " not in inner:
        for tok in inner.split(","):
            tok = tok.strip()
            if not tok:
                continue
            if "-" in tok:
                a, b = tok.split("-", 1)
                try:
                    a, b = int(a), int(b)
                    out.extend(range(a, b + 1))
                except ValueError:
                    continue
            else:
                try:
                    out.append(int(tok))
                except ValueError:
                    continue
        if out:
            return sorted(set(out))
    # OLD schema fallback: whitespace-separated ints
    for tok in re.split(r"\s+", inner):
        try:
            out.append(int(tok))
        except ValueError:
            continue
    return sorted(set(out))


def fill_gaps(positions: list[int], max_gap: int = 5) -> list[int]:
    """Bridge ≤ max_gap residues so cartoon coloring looks continuous."""
    if len(positions) < 2:
        return positions
    out = [positions[0]]
    for prev, curr in zip(positions, positions[1:]):
        gap = curr - prev - 1
        if 0 < gap <= max_gap:
            out.extend(range(prev + 1, curr))
        out.append(curr)
    return out


def prune_short_segments(positions: list[int], min_len: int = 2) -> list[int]:
    """Drop any contiguous run of residues shorter than `min_len`.

    With `min_len=2` (default): isolated single residues are dropped but
    every 2+ run survives. With `min_len=3`: also drops dimer runs.
    Useful for cleaning a LIR set before exporting a partial structure
    for Foldseek-multimer — singleton residues add noise without
    contributing a recognizable structural fragment. Pass `min_len=1`
    to disable (every residue kept).
    """
    if min_len <= 1 or not positions:
        return positions
    positions = sorted(set(positions))
    runs: list[list[int]] = []
    cur: list[int] = [positions[0]]
    for p in positions[1:]:
        if p == cur[-1] + 1:
            cur.append(p)
        else:
            runs.append(cur)
            cur = [p]
    runs.append(cur)
    return [p for run in runs if len(run) >= min_len for p in run]


def positions_to_compact_ranges(positions: list[int], sep: str = ",") -> str:
    """`[1, 2, 3, 10, 11, 15]` → `'1-3,10-11,15'`. Range-compression only."""
    if not positions:
        return ""
    positions = sorted(set(positions))
    ranges: list[str] = []
    a = b = positions[0]
    for p in positions[1:]:
        if p == b + 1:
            b = p
        else:
            ranges.append(f"{a}-{b}" if a != b else f"{a}")
            a = b = p
    ranges.append(f"{a}-{b}" if a != b else f"{a}")
    return sep.join(ranges)


def positions_to_resspec(positions: list[int], chain: str) -> str:
    """`[1, 2, 3, 10, 11, 15]` → `'/A:1-3,10-11,15'` for ChimeraX."""
    compressed = positions_to_compact_ranges(positions)
    return f"/{chain}:{compressed}" if compressed else ""


# ----------------------------------------------------------------------
# Color palettes
# ----------------------------------------------------------------------
#
# Design (matches LIVIA web tool, https://flyark.github.io/LIVIA): each
# CHAIN gets one *base* color (used for cLIR). The corresponding LIR color
# is derived by lightening (`lighten(base, 0.5)`) — i.e. two shades of the
# same hue per chain. With `--single-color` (LIVIA's "Solid" mode), LIR
# and cLIR collapse to the same base color.
#
# When a fold has > len(palette) chains, palettes cycle. For 3+ chain
# folds the LIVIA web tool defaults to "tab10-gradient"; we mirror that.

# LIVIA-derived base names (used for tab10-tail extension; the toolkit's
# whole palette set is already LIVIA-flavored so the prefix would just
# be noise).
_LIVIA_BASE_NAMES = [
    "teal-coral", "blue-orange", "purple-gold", "slate-rose",
    "indigo-tangerine", "green-red",
    "cyan-navy-peach-crimson", "sky-indigo-salmon-red",
    "lime-forest-gold-maroon", "mint-darkteal-apricot-burnt",
]

# Each palette: list of base (cLIR) hex colors, one per chain in
# alphabetical order. First N chains get the first N colors.
PALETTES: dict[str, list[str]] = {
    # 2-chain LIVIA presets (the cLIR-dark of each preset pair → base
    # colors for chain A then chain B). For 3+ chains the list extends
    # with Tab10 colors so additional chains are still distinguishable.
    "teal-coral":     ["#00897B", "#E64A19"],   # default
    "blue-orange":    ["#2471A3", "#E67E22"],
    "purple-gold":    ["#5e35b1", "#f9a825"],
    "slate-rose":     ["#546e7a", "#c2185b"],
    "indigo-tangerine": ["#303f9f", "#e65100"],
    "green-red":      ["#2e7d32", "#c62828"],
    "cyan-navy-peach-crimson":     ["#1A237E", "#B71C1C", "#4DD0E1", "#FFAB91"],
    "sky-indigo-salmon-red":       ["#283593", "#C62828", "#81D4FA", "#EF9A9A"],
    "lime-forest-gold-maroon":     ["#1B5E20", "#B71C1C", "#AED581", "#FFD54F"],
    "mint-darkteal-apricot-burnt": ["#004D40", "#BF360C", "#80CBC4", "#FFCC80"],
    # Generic n-chain palettes
    "tab10": [
        "#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
        "#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf",
    ],
    "tab20": [
        "#1f77b4","#aec7e8","#ff7f0e","#ffbb78","#2ca02c","#98df8a","#d62728","#ff9896",
        "#9467bd","#c5b0d5","#8c564b","#c49c94","#e377c2","#f7b6d2","#7f7f7f","#c7c7c7",
        "#bcbd22","#dbdb8d","#17becf","#9edae5",
    ],
    "mono": ["#37474F"],   # single color for every chain — solid mono look
}
# Tab10 tail appended to every LIVIA 2-color preset so >2-chain folds
# still get unique colors per chain instead of recycling chain A's color.
_TAB10_TAIL = [c for c in PALETTES["tab10"][2:]]
for _name in _LIVIA_BASE_NAMES:
    if len(PALETTES[_name]) < 10:
        PALETTES[_name] = PALETTES[_name] + _TAB10_TAIL[: 10 - len(PALETTES[_name])]

# Backwards-compat: accept the old `livia-*` prefixed names silently.
PALETTE_ALIASES: dict[str, str] = {f"livia-{n}": n for n in _LIVIA_BASE_NAMES}


def lighten(hex_color: str, amount: float = 0.5) -> str:
    """LIVIA's `lighten`: mix `hex_color` with white by `amount` (0..1)."""
    h = hex_color.lstrip("#")
    if len(h) != 6:
        return hex_color
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    lr = round(r + (255 - r) * amount)
    lg = round(g + (255 - g) * amount)
    lb = round(b + (255 - b) * amount)
    return f"#{lr:02X}{lg:02X}{lb:02X}"


def assign_chain_colors(
    chains: list[str],
    palette: list[str],
    *,
    single_color: bool = False,
    lir_lighten: float = 0.5,
    lir_palette: list[str] | None = None,
    clir_palette: list[str] | None = None,
) -> dict[str, tuple[str, str]]:
    """Map each chain → (lir_hex, clir_hex).

    Per-CHAIN (not per-pair) — so chain A always gets the same color
    whether it appears in pair A-B or pair A-C. Each palette cycles if
    there are more chains than colors.

    Default mode: cLIR = `palette[i]`; LIR = `lighten(cLIR, lir_lighten)`.
    `single_color=True`: LIR == cLIR == base (LIVIA's "Solid" mode).

    Overrides (independent):
      `clir_palette`: per-chain cLIR colors. Replaces `palette` for the
        cLIR (dark) shade. LIR is still derived unless `lir_palette` is
        also given.
      `lir_palette`: per-chain LIR colors. Replaces the lighten()
        derivation entirely; LIR is no longer tied to the cLIR base.
        Overrides `single_color` for the LIR shade.

    Either or both can be passed. When both are given, LIR and cLIR are
    fully independent (no derivation between them).
    """
    out: dict[str, tuple[str, str]] = {}
    for i, ch in enumerate(chains):
        # cLIR base: explicit override > palette
        if clir_palette is not None:
            clir_hex = clir_palette[i % len(clir_palette)]
        else:
            clir_hex = palette[i % len(palette)]
        # LIR: explicit override > single-color collapse > lighten derivation
        if lir_palette is not None:
            lir_hex = lir_palette[i % len(lir_palette)]
        elif single_color:
            lir_hex = clir_hex
        else:
            lir_hex = lighten(clir_hex, lir_lighten)
        out[ch] = (lir_hex, clir_hex)
    return out


# ----------------------------------------------------------------------
# Core cxc emitter
# ----------------------------------------------------------------------

CXC_HEADER = """# {fold}  rank={rank}  model={model}
# Generated by lis_to_cxc.py — one cxc per (fold, rank, model), one
# section per chain-pair interface. Open in ChimeraX (double-click .cxc
# or `open <this file>`); the structure loads automatically and each
# chain-pair's LIR / cLIR residues are colored per the LIVIA palette.
#
# Provenance (absolute paths at generation time):
#   Source CSV:       {csv_abspath}
#   Structure file:   {structure_abspath}
#   Output cxc:       {cxc_abspath}
{path_note}
#
# iLIS thresholds (0.223 / 0.339 / 0.551 = FPR 10% / 5% / 1%) are calibrated
# on large-scale Y2H reference sets in yeast, fly, and human predicted using
# ColabFold (AF2-multimer). See Kim et al., 2026 for benchmark details.
# Different platforms (AF3, Boltz-2, Chai-1, OpenFold3) may require different
# thresholds — validate against known positives for platform-specific claims.
#
# Interface summary ({n_pairs} chain pair{n_pairs_plural} in this model):
{interface_summary}

# Clear any 2D labels left over from a previous cxc in the same session
# so this cxc's title/pair overlay starts from a clean slate.
# A fresh ChimeraX (e.g. 1.12) only registers the `~2dlabels` delete form
# after `2dlabels` has been used once in the session, so calling it first
# errors "~2dlabels is an unknown command" (flyark/AFM-LIS#15). Create a
# throwaway label to load the command, then `~2dlabels` drops ALL labels
# (ours are namespaced `lis2cxc_*`; a future ChimeraX with wildcard delete
# could narrow this to `2dlabels delete lis2cxc_*`).
2dlabels create lis2cxc_reg text " "
~2dlabels

close
# Path is quoted so spaces / parentheses in directory names (e.g. macOS
# `Library/CloudStorage`, Dropbox, iCloud folders, "Image Lab (Bio-RaD)")
# don't break ChimeraX's `open` tokenizer.
open "{structure_openpath}"

# Print the interface summary to ChimeraX's log so it's visible at-run-
# time (comments above are parsed-then-discarded by ChimeraX, only useful
# in a text editor). `echo` writes a line to the log panel.
echo ============================================================
echo {fold}  rank={rank}  model={model}
echo ============================================================
{interface_echo}
echo ============================================================

# Whole-structure starting state: all chains hidden, all atoms white.
# The per-chain LIR phase below shows only the residues that
# participate in surviving (post-cutoff) chain pairs. Non-LIR
# residues stay hidden — focused LIVIA-style interface view.
hide cartoon
color white

lighting soft
graphics silhouettes true
~select

# 2D label overlay on the viewport itself (so the interface info is
# visible on the structure, not just in the log). To hide: `~2dlabels`.
{interface_2dlabels}
"""

CXC_PAIR_COMMENT = """
# ---- {chain_i} ↔ {chain_j}  iLIS={ilis:.3f}  LIS={lis:.3f}  cLIS={clis:.3f}  ipTM={iptm:.2f} ----
# chain {chain_i}: LIR={lir_i_str}
# chain {chain_j}: LIR={lir_j_str}
#       cLIR (in contact) chain {chain_i}: {clir_i_str}
#       cLIR (in contact) chain {chain_j}: {clir_j_str}
"""

# ----------------------------------------------------------------------
# LIR-only partial structure extraction (for Foldseek-multimer etc.)
# ----------------------------------------------------------------------

def _filter_pdb_by_lir(structure_path: Path,
                       lir_by_chain: dict[str, set[int]],
                       out_path: Path) -> bool:
    """Filter a PDB file to keep only ATOM records whose (chain, resnum)
    is in `lir_by_chain`. TER records are inserted between chains.
    Returns True if at least one atom was written."""
    out_lines: list[str] = []
    chains_seen: set[str] = set()
    last_chain: str | None = None
    serial = 1
    with open(structure_path) as f:
        for line in f:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            chain = line[21:22].strip()
            try:
                resnum = int(line[22:26])
            except ValueError:
                continue
            if chain not in lir_by_chain or resnum not in lir_by_chain[chain]:
                continue
            if last_chain is not None and chain != last_chain:
                out_lines.append("TER\n")
            # Renumber serial in cols 7-11 to keep PDB valid after deletion.
            out_lines.append(line[:6] + f"{serial:5d}" + line[11:])
            serial += 1
            last_chain = chain
            chains_seen.add(chain)
    if not out_lines:
        return False
    with open(out_path, "w") as f:
        f.write("REMARK   1 LIR-filtered partial structure (lis_to_cxc.py --save-lir-pdb)\n")
        f.write(f"REMARK   1 Source: {structure_path.name}\n")
        f.write(f"REMARK   1 Chains: {','.join(sorted(chains_seen))}\n")
        for ch in sorted(chains_seen):
            n = len(lir_by_chain[ch])
            f.write(f"REMARK   1 Chain {ch}: {n} LIR residue(s)\n")
        f.writelines(out_lines)
        f.write("TER\nEND\n")
    return True


def _filter_cif_by_lir_to_pdb(structure_path: Path,
                              lir_by_chain: dict[str, set[int]],
                              out_path: Path) -> bool:
    """Parse an mmCIF atom_site loop, filter by LIR, write PDB ATOM records.
    Output is always PDB (single canonical format for Foldseek)."""
    with open(structure_path) as f:
        content = f.read()
    # Find the atom_site loop header.
    m = re.search(r'(?ms)^loop_\s*\n((?:_atom_site\.[^\n]+\n)+)', content)
    if not m:
        return False
    col_names = [ln.strip().split('.', 1)[1].strip()
                 for ln in m.group(1).splitlines() if ln.strip()]
    def idx(*candidates):
        for c in candidates:
            if c in col_names:
                return col_names.index(c)
        return None
    i_group = idx('group_PDB')
    i_atom  = idx('label_atom_id')
    i_alt   = idx('label_alt_id')
    i_res   = idx('label_comp_id', 'auth_comp_id')
    i_chain = idx('auth_asym_id', 'label_asym_id')
    i_resn  = idx('auth_seq_id', 'label_seq_id')
    i_x     = idx('Cartn_x')
    i_y     = idx('Cartn_y')
    i_z     = idx('Cartn_z')
    i_occ   = idx('occupancy')
    i_b     = idx('B_iso_or_equiv')
    i_elem  = idx('type_symbol')
    if None in (i_group, i_atom, i_res, i_chain, i_resn, i_x, i_y, i_z):
        return False
    # Data section starts after the column block; ends at next loop_ / data_ / save_.
    data_start = m.end()
    tail = content[data_start:]
    end_m = re.search(r'\n(?:loop_|data_|save_|#\s*$)', tail)
    if end_m:
        tail = tail[:end_m.start()]
    out_lines: list[str] = []
    chains_seen: set[str] = set()
    last_chain: str | None = None
    serial = 1
    for raw in tail.splitlines():
        s = raw.strip()
        if not s or s.startswith('#'):
            continue
        toks = s.split()
        if len(toks) < len(col_names):
            continue
        if toks[i_group] not in ('ATOM', 'HETATM'):
            continue
        chain = toks[i_chain]
        try:
            resnum = int(toks[i_resn])
        except ValueError:
            continue
        if chain not in lir_by_chain or resnum not in lir_by_chain[chain]:
            continue
        atom_name = toks[i_atom].strip('"').strip("'")
        resname = toks[i_res][:3]
        try:
            x = float(toks[i_x]); y = float(toks[i_y]); z = float(toks[i_z])
        except ValueError:
            continue
        occ = float(toks[i_occ]) if i_occ is not None else 1.0
        bfac = float(toks[i_b]) if i_b is not None else 0.0
        elem = (toks[i_elem] if i_elem is not None else atom_name[:1]).upper()[:2]
        # PDB atom-name column rule: 1-3 chars left-padded with one space.
        if len(atom_name) >= 4:
            atom_field = atom_name[:4]
        else:
            atom_field = f" {atom_name:<3s}"
        if last_chain is not None and chain != last_chain:
            out_lines.append("TER\n")
        out_lines.append(
            f"ATOM  {serial:5d} {atom_field} {resname:>3s} {chain[:1]:1s}"
            f"{resnum:4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
            f"{occ:6.2f}{bfac:6.2f}          {elem:>2s}\n"
        )
        serial += 1
        last_chain = chain
        chains_seen.add(chain)
    if not out_lines:
        return False
    with open(out_path, "w") as f:
        f.write("REMARK   1 LIR-filtered partial structure (lis_to_cxc.py --save-lir-pdb)\n")
        f.write(f"REMARK   1 Source: {structure_path.name} (converted from mmCIF)\n")
        f.write(f"REMARK   1 Chains: {','.join(sorted(chains_seen))}\n")
        for ch in sorted(chains_seen):
            n = len(lir_by_chain[ch])
            f.write(f"REMARK   1 Chain {ch}: {n} LIR residue(s)\n")
        f.writelines(out_lines)
        f.write("TER\nEND\n")
    return True


def extract_lir_partial_pdb(structure_path: Path,
                            rows: list[dict],
                            gap_fill: int,
                            out_path: Path,
                            min_segment: int = 5) -> bool:
    """Write a partial PDB with only LIR-region atoms from `rows`.

    Cleanup pipeline before filtering:
      1. parse_residue_indices    (raw LIR set per pair)
      2. union across all chain-pair rows the chain participates in
      3. fill_gaps(..., gap_fill) (bridge short gaps)
      4. prune_short_segments(..., min_segment)
         (drop isolated single residues by default — they're noise as
         Foldseek input)

    Output is always PDB; CIF inputs are converted to PDB. Useful as
    direct input for Foldseek-multimer or any structural-similarity
    search that should operate on the confident interaction interface
    only, not the full multimer.
    """
    lir_by_chain: dict[str, set[int]] = {}
    for r in rows:
        ci = (r.get("chain_i") or "A").strip()
        cj = (r.get("chain_j") or "B").strip()
        lir_by_chain.setdefault(ci, set()).update(
            parse_residue_indices(r.get("LIR_indices_i", ""))
        )
        lir_by_chain.setdefault(cj, set()).update(
            parse_residue_indices(r.get("LIR_indices_j", ""))
        )
    # Apply gap-fill + segment pruning per chain so cleanup respects
    # each chain's contiguity independently.
    cleaned: dict[str, set[int]] = {}
    for ch, positions in lir_by_chain.items():
        filled = fill_gaps(sorted(positions), gap_fill)
        pruned = prune_short_segments(filled, min_segment)
        if pruned:
            cleaned[ch] = set(pruned)
    if not cleaned:
        return False
    ext = structure_path.suffix.lower()
    if ext in (".pdb", ".ent"):
        return _filter_pdb_by_lir(structure_path, cleaned, out_path)
    if ext in (".cif", ".mmcif"):
        return _filter_cif_by_lir_to_pdb(structure_path, cleaned, out_path)
    return False


CXC_FOOTER = """
~select

# Optional alternate coloring (uncomment to use):
# color bfactor palette alphafold     # pLDDT spectrum
# color bfactor palette ChimeraX/palettes/pLDDT-AF2     # ChimeraX built-in

# Optional: save the current viewport as a publication-quality PNG with
# transparent background (uncomment, edit filename, run from the command
# line in ChimeraX). Useful for paper figures and slide decks.
# save myfigure.png transparentBackground true

# Print the same tip to ChimeraX's Log panel so it's visible at runtime
# (`#` comments above are silently ignored by the parser — only `echo`
# lines surface in the Log).
echo ------------------------------------------------------------
echo Tip: to save this view as a transparent-background PNG for figures:
echo     save myfigure.png transparentBackground true
echo ------------------------------------------------------------
"""


def emit_cxc(
    rows: list[dict],
    structure_path: Path,
    out_path: Path,
    *,
    gap_fill: int = 5,
    palette: list[str],
    single_color: bool = False,
    csv_path: Path | None = None,
    metadata: dict[str, str] | None = None,
    plddt: bool = False,
    lir_palette: list[str] | None = None,
    clir_palette: list[str] | None = None,
    lir_lighten: float = 0.5,
    min_segment: int = 5,
    relative: bool = False,
) -> bool:
    """Emit one cxc covering every chain-pair row in this (fold, rank, model).

    Filtering by --ilis-cutoff / --pairs / --top-N happens upstream in
    main(); by the time we get here, `rows` is the surviving pair set
    for this group. Returns True if the cxc was written; False only if
    `rows` is empty.

    Coloring is per-chain (LIVIA convention): each chain in the fold gets
    one base color; LIR = lighten(base, lir_lighten), cLIR = base (gradient
    mode). With `single_color=True`, LIR == cLIR == base. The palette
    cycles if there are more chains than entries.

    With `plddt=True` the chain palette is unused — the whole structure
    is rendered as cartoon and colored by pLDDT (`color bfactor palette
    alphafold`). Interface info (title, 2D labels, log echo, header
    summary) is still emitted; only the cartoon color phases are
    replaced.

    2D label IDs are namespaced with the `lis2cxc_` prefix so they
    don't collide with labels from other tools or hand-created labels.
    """
    if not rows:
        return False
    kept = rows

    # Per-chain color assignment based on all chains in this group.
    all_chains = sorted({
        (r.get("chain_i") or "?").strip() for r in kept
    } | {
        (r.get("chain_j") or "?").strip() for r in kept
    })
    chain_colors = assign_chain_colors(
        all_chains, palette,
        single_color=single_color,
        lir_lighten=lir_lighten,
        lir_palette=lir_palette,
        clir_palette=clir_palette,
    )

    # `open` path — absolute by default (reliable double-click on the machine that
    # generated it), or relative to the cxc's own folder with --relative. ChimeraX
    # resolves a relative `open` against the command file's location, so a relative
    # path makes the cxc portable across machines/OS (e.g. generate on Linux, view
    # on Windows) as long as the cxc + structure move together (flyark/AFM-LIS#17).
    if relative:
        try:
            open_path = os.path.relpath(structure_path.resolve(), out_path.resolve().parent)
        except ValueError:                       # different Windows drive: no relative path exists
            open_path = str(structure_path.resolve())
        open_path = open_path.replace(os.sep, '/')   # forward slashes so a Linux-made path also opens on Windows
        path_note = ("# The `open` line below uses a path RELATIVE to this cxc's location, so the\n"
                     "# bundle is portable: move or share the cxc together with its structure\n"
                     "# (keeping their relative layout) and it opens on any machine/OS. ChimeraX\n"
                     "# resolves the relative `open` against the cxc's own folder (--relative).")
    else:
        open_path = str(structure_path.resolve())
        path_note = ("# The `open` line below uses the structure's ABSOLUTE path — reliable for\n"
                     "# double-clicking on the machine that generated it. For a portable or\n"
                     "# cross-OS bundle (e.g. generate on Linux, view on Windows), regenerate\n"
                     "# with --relative (flyark/AFM-LIS#17).")

    # Build the at-a-glance interface summary block. One line per chain pair
    # with scores; two indented lines per pair listing LIR residue ranges
    # on each side. LIR (broader confident-interface region) reads as a
    # contiguous boundary; cLIR (direct contacts) is too granular for an
    # at-a-glance summary. Sorted by iLIS descending so strongest
    # interfaces float to the top.
    sorted_kept = sorted(kept, key=lambda r: -_float(r.get("iLIS"), 0.0))
    summary_lines: list[str] = []
    echo_lines: list[str] = []
    for r in sorted_kept:
        ci = (r.get("chain_i") or "?").strip()
        cj = (r.get("chain_j") or "?").strip()
        lir_i = parse_residue_indices(r.get("LIR_indices_i", ""))
        lir_j = parse_residue_indices(r.get("LIR_indices_j", ""))
        # Re-compress to range form regardless of input format. lis.py
        # variants emit either "[16-20,22-23]" (already compressed) or
        # "[29, 30, 32, 33]" (expanded); we always want the compact form
        # in the summary for human readability.
        lir_i_ranges = positions_to_compact_ranges(lir_i) or "(none)"
        lir_j_ranges = positions_to_compact_ranges(lir_j) or "(none)"
        header_text = (
            f"{ci} <-> {cj}  "
            f"iLIS={_float(r.get('iLIS')):.3f}  "
            f"LIS={_float(r.get('LIS')):.3f}  "
            f"cLIS={_float(r.get('cLIS')):.3f}  "
            f"ipTM={_float(r.get('ipTM')):.2f}"
        )
        lir_i_text = f"LIR {ci} ({len(lir_i)} residue{'s' if len(lir_i)!=1 else ''}): {lir_i_ranges}"
        lir_j_text = f"LIR {cj} ({len(lir_j)} residue{'s' if len(lir_j)!=1 else ''}): {lir_j_ranges}"
        # Comment form for in-file inspection.
        summary_lines.append(f"#   {header_text}")
        summary_lines.append(f"#     {lir_i_text}")
        summary_lines.append(f"#     {lir_j_text}")
        # `echo` form so the same info prints to ChimeraX's log when the
        # cxc runs (comments are silently consumed by the parser).
        echo_lines.append(f"echo   {header_text}")
        echo_lines.append(f"echo     {lir_i_text}")
        echo_lines.append(f"echo     {lir_j_text}")
    interface_summary = "\n".join(summary_lines)
    n_pairs_plural = "s" if len(kept) != 1 else ""
    interface_echo = "\n".join(
        [f"echo Interface summary ({len(kept)} chain pair{n_pairs_plural} in this model):"]
        + echo_lines
    )

    head = rows[0]

    # Per-chain unions of LIR and cLIR positions across every pair this
    # chain participates in. Built here so both the 2dlabels overlay
    # below and the per-chain emit phases further down can reuse them.
    # Cleanup applied: (1) gap_fill bridges short gaps so the cartoon
    # band reads as a continuous boundary; (2) min_segment drops short
    # runs (default 5) — without this, isolated single residues float
    # disconnected from the main structure in the 3D view.
    lir_by_chain: dict[str, set[int]] = {}
    clir_by_chain: dict[str, set[int]] = {}
    for r in kept:
        chain_i = (r.get("chain_i") or "A").strip()
        chain_j = (r.get("chain_j") or "B").strip()
        lir_by_chain.setdefault(chain_i, set()).update(
            parse_residue_indices(r.get("LIR_indices_i", ""))
        )
        lir_by_chain.setdefault(chain_j, set()).update(
            parse_residue_indices(r.get("LIR_indices_j", ""))
        )
        clir_by_chain.setdefault(chain_i, set()).update(
            parse_residue_indices(r.get("cLIR_indices_i", ""))
        )
        clir_by_chain.setdefault(chain_j, set()).update(
            parse_residue_indices(r.get("cLIR_indices_j", ""))
        )
    # Apply per-chain gap-fill + short-segment pruning to LIR.
    for ch, positions in list(lir_by_chain.items()):
        filled = fill_gaps(sorted(positions), gap_fill)
        pruned = prune_short_segments(filled, min_segment)
        lir_by_chain[ch] = set(pruned)
    # cLIR ⊆ LIR (cLIR is the subset of LIR with direct contacts). After
    # pruning LIR, drop any cLIR residue whose LIR neighborhood was just
    # pruned away — otherwise the 3D would color cLIR onto a hidden
    # cartoon. Keeps cLIR consistent with what's actually shown.
    for ch in list(clir_by_chain.keys()):
        if ch in lir_by_chain:
            clir_by_chain[ch] = clir_by_chain[ch] & lir_by_chain[ch]
        else:
            clir_by_chain[ch] = set()

    # 2D label overlay rendered on the viewport itself. Title at the top,
    # then one line per chain pair (iLIS + ipTM only — per-pair residue
    # lists are too verbose for an on-screen overlay; they stay in the
    # echo/log and in the header comments). xpos/ypos are normalized 0–1
    # with origin at the lower-left. Cap visible pair labels at
    # MAX_VISIBLE to avoid overrunning the bottom of the viewport on
    # large complexes (e.g. 22-chain folds with 21 pairs).
    MAX_VISIBLE = 10
    # Font sizes — smaller than the original to leave more room for the
    # structure underneath. Title 14 / pair 11 / LIR 10 keep enough
    # legibility while shrinking the label panel ~20%.
    SIZE_TITLE = 14
    SIZE_PAIR  = 11
    SIZE_LIR   = 10
    # ypos delta between successive lines, tuned so adjacent bgColor-white
    # rectangles touch without overlap. Size N text is roughly N * 0.0017
    # of normalized height in a typical ChimeraX viewport.
    LINE_DY = 0.022       # between pair lines (size 11)
    LIR_DY  = 0.020       # between LIR sub-lines (size 10)
    # margin 0 removes the default bgColor padding so rectangles butt up.
    LABEL_TAIL = "margin 0"
    label_lines: list[str] = []
    # Label as `iLIS_rank` when the rank was reassigned by iLIS
    # (synthetic-rank platforms: AF3 / Boltz-2 / Chai-1 / OpenFold3 —
    # detected via the `_orig_rank` flag set during reassignment).
    # ColabFold keeps the predictor's native confidence rank → plain `rank`.
    rank_label = "iLIS_rank" if head.get("_orig_rank") is not None else "rank"
    title_text = (
        f"{head.get('name','?')}  "
        f"{rank_label}={head.get('rank','?')}  "
        f"model={head.get('model','?')}"
    )
    label_lines.append(
        f'2dlabels create lis2cxc_title text "{title_text}" '
        f'xpos 0.02 ypos 0.96 size {SIZE_TITLE} color black bgColor white {LABEL_TAIL}'
    )
    # First pair sits one title-height (~0.025) below the title.
    ypos = 0.935
    # Per-pair labels are split into 4 sub-labels so each chain letter
    # can carry that chain's color (matching the 3D coloring).
    # CHAR_DX scaled for the smaller pair font (size 11).
    CHAR_DX = 0.010
    for i, r in enumerate(sorted_kept[:MAX_VISIBLE]):
        ci = (r.get("chain_i") or "?").strip()
        cj = (r.get("chain_j") or "?").strip()
        ilis = _float(r.get("iLIS"))
        iptm = _float(r.get("ipTM"))
        # In pLDDT mode the 3D view is colored by confidence, not chain,
        # so the per-chain label colors wouldn't map to anything in the
        # viewport. Keep label letters black there.
        col_i = "black" if plddt else chain_colors[ci][1]
        col_j = "black" if plddt else chain_colors[cj][1]
        x = 0.02
        label_lines.append(
            f'2dlabels create lis2cxc_pair{i}_i text "{ci}" '
            f'xpos {x:.4f} ypos {ypos:.3f} size {SIZE_PAIR} color {col_i} bgColor white {LABEL_TAIL}'
        )
        x += CHAR_DX * max(1, len(ci))
        label_lines.append(
            f'2dlabels create lis2cxc_pair{i}_sep text "-" '
            f'xpos {x:.4f} ypos {ypos:.3f} size {SIZE_PAIR} color black bgColor white {LABEL_TAIL}'
        )
        x += CHAR_DX
        label_lines.append(
            f'2dlabels create lis2cxc_pair{i}_j text "{cj}" '
            f'xpos {x:.4f} ypos {ypos:.3f} size {SIZE_PAIR} color {col_j} bgColor white {LABEL_TAIL}'
        )
        x += CHAR_DX * max(1, len(cj))
        score_text = f"  iLIS={ilis:.3f}  ipTM={iptm:.2f}"
        label_lines.append(
            f'2dlabels create lis2cxc_pair{i}_rest text "{score_text}" '
            f'xpos {x:.4f} ypos {ypos:.3f} size {SIZE_PAIR} color black bgColor white {LABEL_TAIL}'
        )
        ypos -= LINE_DY
    if len(sorted_kept) > MAX_VISIBLE:
        remaining = len(sorted_kept) - MAX_VISIBLE
        label_lines.append(
            f'2dlabels create lis2cxc_pairmore text "... +{remaining} more pair(s) - see log" '
            f'xpos 0.02 ypos {ypos:.3f} size {SIZE_LIR} color #555555 bgColor white {LABEL_TAIL}'
        )

    # Per-chain LIR residue lines in the label overlay (LIVIA shows the
    # binding residues in the side panel; we put them in the on-canvas
    # label area). LIR (broader confident-interface region, gap-filled)
    # gives an aggregated boundary view that reads cleanly on screen;
    # cLIR (contact-only) is too scattered for an at-a-glance overlay.
    # Each line uses the chain's base color so it visually maps to the
    # 3D coloring. Long lists truncate aggressively (MAX_LIR_CHARS = 40)
    # so the text doesn't run wide enough to overlap the structure.
    MAX_LIR_CHARS = 40
    for ch in all_chains:
        lir_positions = sorted(lir_by_chain.get(ch, set()))
        if not lir_positions:
            continue
        col = chain_colors[ch][1]  # base (dark) color — legible on white bg
        ranges = positions_to_compact_ranges(lir_positions)
        if len(ranges) > MAX_LIR_CHARS:
            ranges = ranges[:MAX_LIR_CHARS - 3] + "..."
        text_part = f"LIR {ch} ({len(lir_positions)}): {ranges}"
        label_lines.append(
            f'2dlabels create lis2cxc_lir_{ch} text "{text_part}" '
            f'xpos 0.02 ypos {ypos:.3f} size {SIZE_LIR} color {col} bgColor white {LABEL_TAIL}'
        )
        ypos -= LIR_DY
    interface_2dlabels = "\n".join(label_lines)
    text = CXC_HEADER.format(
        fold=head.get("name", "?"),
        rank=head.get("rank", "?"),
        model=head.get("model", "?"),
        structure_openpath=open_path,
        path_note=path_note,
        structure_abspath=str(structure_path.resolve()),
        csv_abspath=str(csv_path.resolve()) if csv_path else "(not provided)",
        cxc_abspath=str(out_path.resolve()),
        n_pairs=len(kept),
        n_pairs_plural=n_pairs_plural,
        interface_summary=interface_summary,
        interface_echo=interface_echo,
        interface_2dlabels=interface_2dlabels,
    )
    # Note if rank was recomputed (synthetic rank=model from lis.py on
    # non-ColabFold platforms — AF3 / Boltz-2 / Chai-1 / OpenFold3).
    orig_rank = head.get("_orig_rank")
    if orig_rank is not None and str(orig_rank) != str(head.get("rank")):
        text += (f"\n# Rank reassigned by lis_to_cxc.py from descending iLIS "
                 f"(lis.py CSV had rank={orig_rank}, equal to model — typical of "
                 f"AF3 / Boltz-2 / Chai-1 / OpenFold3 platforms which have no native rank).\n")

    # Record LIVIA-style gap-fill setting used for LIR rendering.
    # (Irrelevant in --plddt mode but cheap to record for the manifest.)
    if gap_fill > 0:
        text += f"\n# LIR gap-fill: bridged gaps <= {gap_fill} residues for visual continuity (--gap-fill {gap_fill}).\n"
    else:
        text += "\n# LIR gap-fill: disabled (--gap-fill 0); exact residue ranges only.\n"

    # Optional metadata block (from --metadata TSV in main).
    if metadata:
        text += "\n# Metadata (from --metadata TSV):\n"
        for k, v in metadata.items():
            if v:
                text += f"#   {k}: {v}\n"

    # Color legend in comments so user can see which chain is which —
    # only meaningful in LIVIA mode; in --plddt mode the structure is
    # colored by confidence, not chain, so the chain palette is unused.
    if not plddt:
        text += "\n# Chain color assignments (LIVIA per-chain convention):\n"
        for ch in all_chains:
            lir_c, clir_c = chain_colors[ch]
            if single_color:
                text += f"#   chain {ch}: {clir_c} (LIR = cLIR, solid mode)\n"
            else:
                text += f"#   chain {ch}: cLIR={clir_c} | LIR={lir_c} (50%-lightened)\n"
    else:
        text += "\n# Color mode: pLDDT spectrum (--plddt) — chain palette unused.\n"
        text += "#   per-residue: red/orange < 50, yellow 50-70, light blue 70-90, dark blue > 90\n"

    # Per-pair comment blocks. Actual color/show commands are emitted
    # once per chain further down, using the per-chain unions built
    # earlier (lir_by_chain / clir_by_chain). Per-chain emission order
    # — ALL LIR first, then ALL cLIR — guarantees cLIR always wins in
    # the final paint regardless of which pair declared it; matches
    # LIVIA's rendering convention.
    for r in kept:
        chain_i = (r.get("chain_i") or "A").strip()
        chain_j = (r.get("chain_j") or "B").strip()
        text += CXC_PAIR_COMMENT.format(
            chain_i=chain_i, chain_j=chain_j,
            ilis=_float(r.get("iLIS")),
            lis=_float(r.get("LIS")),
            clis=_float(r.get("cLIS")),
            iptm=_float(r.get("ipTM")),
            lir_i_str=r.get("LIR_indices_i", ""),
            lir_j_str=r.get("LIR_indices_j", ""),
            clir_i_str=r.get("cLIR_indices_i", ""),
            clir_j_str=r.get("cLIR_indices_j", ""),
        )

    # --plddt mode: show the same per-chain LIR scope as LIVIA modes
    # (matching user expectation — same residues visible across all
    # color modes), but color them by pLDDT confidence instead of the
    # chain palette. pLDDT is encoded in the bfactor field for
    # AlphaFold-derived PDBs.
    if plddt:
        text += "\n# pLDDT-spectrum rendering on per-chain LIR residues\n"
        text += "# (same scope as LIVIA modes; coloring by per-residue pLDDT).\n"
        for ch in all_chains:
            lir_positions = sorted(lir_by_chain.get(ch, set()))
            if not lir_positions:
                continue
            resspec = positions_to_resspec(lir_positions, ch)
            text += f"show {resspec} cartoon\n"
        text += "color bfactor palette alphafold\n"
        text += "\n~select\n"
        out_path.write_text(text)
        return True

    # Phase 1: per-chain LIR (show cartoon + light color).
    text += "\n# Per-chain LIR: show cartoon + light-shade fill.\n"
    for ch in all_chains:
        lir_positions = sorted(lir_by_chain.get(ch, set()))
        if not lir_positions:
            continue
        col_lir, _ = chain_colors[ch]
        resspec = positions_to_resspec(lir_positions, ch)
        text += f"show {resspec} cartoon\n"
        text += f"color {resspec} {col_lir}\n"

    # Phase 2: per-chain cLIR override (dark shade — applied last so it
    # always wins over LIR coloring, regardless of pair ordering).
    text += "\n# Per-chain cLIR: dark-shade override on contact residues.\n"
    for ch in all_chains:
        clir_positions = sorted(clir_by_chain.get(ch, set()))
        if not clir_positions:
            continue
        _, col_clir = chain_colors[ch]
        resspec = positions_to_resspec(clir_positions, ch)
        text += f"color {resspec} {col_clir}\n"

    text += CXC_FOOTER
    out_path.write_text(text)
    return True


def parse_colors_flag(value: str, flag_name: str = "--colors") -> list[str]:
    """`#80CBC4,#E64A19,#1f77b4` → ['#80CBC4','#E64A19','#1f77b4'].

    `flag_name` is used in error messages only — pass the actual CLI flag
    being parsed so the user sees a useful diagnostic.
    """
    out = []
    for tok in value.split(","):
        tok = tok.strip()
        if not tok:
            continue
        if not tok.startswith("#"):
            tok = "#" + tok
        if len(tok) not in (4, 7):  # #abc or #aabbcc
            raise SystemExit(f"!! {flag_name}: '{tok}' is not a valid hex color")
        out.append(tok.upper() if len(tok) == 7 else tok)
    if not out:
        raise SystemExit(f"!! {flag_name}: at least one color required")
    return out


def _float(v, default: float = float("nan")) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


# ----------------------------------------------------------------------
# CSV reader (NEW schema; auto-detects OLD via column-name fallbacks)
# ----------------------------------------------------------------------

def _strip_platform_suffix(name: str) -> str:
    """`fold_x_afm3 (AlphaFold3)` → `fold_x_afm3`. Strips a single trailing
    ` (...)` parenthetical, which the LIVIA web tool appends to `name`
    when it exports lis_analysis CSVs client-side from AF3 / ColabFold /
    Boltz-2 / Chai-1 / OpenFold3 zips."""
    return re.sub(r"\s*\([^)]+\)\s*$", "", str(name))


def _normalize_rank_model(r: dict) -> dict:
    """Some lis-CSV variants only have `model` (e.g. `"model_0"`) without
    a separate `rank` column (LIVIA web tool export). Synthesize `rank`
    from `model` when missing so the (fold, rank, model) grouping works.
    """
    if "rank" not in r or r.get("rank") in (None, "", "nan"):
        m = re.search(r"(\d+)", str(r.get("model", "")))
        r["rank"] = m.group(1) if m else r.get("model", "0")
    return r


def read_lis_csv(path: Path) -> list[dict]:
    """Return rows as dicts with NEW-schema keys, regardless of input schema.

    Handles three variants seen in the wild:
      1. Modern lis.py CLI output: full NEW schema with `rank` + `structure_file`.
      2. LIVIA web tool export: NEW schema MINUS `rank` and `structure_file`;
         `name` has a ` (Platform)` suffix.
      3. OLD-schema `_lis_all_ranks.csv` from older AF2-multimer batches:
         long column names, `Protein_1`/`Protein_2`.
    """
    with open(path) as f:
        rows = list(csv.DictReader(f))
    if not rows:
        return []
    cols = set(rows[0].keys())
    if "name" in cols and "chain_i" in cols:
        # Variant 1 or 2 (NEW schema). Strip platform suffix from `name`
        # and synthesize `rank` if missing.
        for r in rows:
            r["name"] = _strip_platform_suffix(r.get("name", ""))
            _normalize_rank_model(r)
        return rows
    # OLD schema heuristic: long verbose column names.
    if "Protein_1" in cols and "Local Interaction Score (Interface)" in cols:
        renamed = []
        for r in rows:
            # OLD schema doesn't carry an iLIS column — lis.py computes it
            # downstream as the geometric mean of LIS and cLIS on the
            # interface scores. Compute it here so downstream tooling
            # (e.g. --ilis-cutoff filtering) Just Works on OLD-schema CSVs
            # without the user needing to re-run lis.py.
            lis  = _float(r.get("Local Interaction Score (Interface)", ""), 0.0)
            clis = _float(r.get("Contact Local Interaction Score (Interface)", ""), 0.0)
            ilis = (lis * clis) ** 0.5 if (lis > 0 and clis > 0) else 0.0
            renamed.append({
                "name":            f"{r.get('Protein_1','?')}___{r.get('Protein_2','?')}",
                "rank":            r.get("Rank", ""),
                "model":           r.get("Rank", ""),
                "chain_i":         "A",
                "chain_j":         "B",
                "iLIS":            f"{ilis:.6f}",  # synthesized from LIS × cLIS
                "LIS":             r.get("Local Interaction Score (Interface)", ""),
                "cLIS":            r.get("Contact Local Interaction Score (Interface)", ""),
                "ipTM":            r.get("ipTM", ""),
                "LIR_indices_i":   r.get("Local Interaction Residue Indice A", ""),
                "LIR_indices_j":   r.get("Local Interaction Residue Indice B", ""),
                "cLIR_indices_i":  r.get("Contact Local Interaction Residue Indice A", ""),
                "cLIR_indices_j":  r.get("Contact Local Interaction Residue Indice B", ""),
                "structure_file":  r.get("pdb_file", ""),
                "_old_schema":     "1",
            })
        return renamed
    raise SystemExit(f"!! unrecognized schema in {path.name}; columns sample: {sorted(cols)[:8]}")


# ----------------------------------------------------------------------
# main
# ----------------------------------------------------------------------

def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    # Allow `--list-palettes` without requiring input files
    p.add_argument("--list-palettes", action="store_true",
                   help="print available named palettes and exit")
    src = p.add_mutually_exclusive_group(required=False)
    src.add_argument("--csv", type=Path, help="one *_lis_analysis.csv file")
    src.add_argument("--csv-dir", type=Path, help="folder of *_lis_analysis.csv files; processed together")
    p.add_argument("--pdb-root", type=Path,
                   help="folder containing the .cif / .pdb model files referenced by `structure_file`")
    p.add_argument("--out", type=Path,
                   help="output folder for emitted .cxc files (created if missing)")
    p.add_argument("--gap-fill", type=int, default=5,
                   help=("LIVIA-style LIR gap-fill: bridge gaps <= N residues "
                         "in the LIR cartoon coloring so the highlighted band "
                         "looks continuous over short missing stretches "
                         "(default 5). Set to 0 to disable. Applies to LIR "
                         "only; cLIR (in-contact residues) is always rendered "
                         "exactly. Value is recorded in the cxc header + "
                         "manifest for reproducibility."))
    p.add_argument("--palette", default="teal-coral",
                   help=("named palette of per-chain base colors. Default: teal-coral. "
                         "Old `livia-*` prefixed names are still accepted as aliases. "
                         "Run with --list-palettes to see all. Use --colors to supply custom hexes."))
    p.add_argument("--colors",
                   help=("comma-separated custom hex colors (e.g. '#80CBC4,#E64A19,#1f77b4') — "
                         "overrides --palette. One color per chain; cycles if there are more chains. "
                         "Acts as the cLIR (dark) base; LIR is derived via --lir-lighten unless "
                         "--lir-colors is also given. For fully independent LIR vs cLIR, use "
                         "--lir-colors + --clir-colors instead."))
    p.add_argument("--clir-colors",
                   help=("comma-separated custom hex colors for cLIR (e.g. '#aa0000,#0066cc'). "
                         "Per-chain in chain-letter order; cycles if fewer colors than chains. "
                         "Overrides --colors / --palette for the cLIR (dark) shade only. "
                         "LIR is still derived via --lir-lighten unless --lir-colors is also "
                         "given. Use together with --lir-colors for fully independent shades."))
    p.add_argument("--lir-colors",
                   help=("comma-separated custom hex colors for LIR (e.g. '#ffaa00,#88ccff'). "
                         "Per-chain in chain-letter order; cycles if fewer colors than chains. "
                         "Overrides the --lir-lighten derivation — LIR is no longer tied to "
                         "the cLIR base. Combine with --clir-colors for fully independent "
                         "LIR / cLIR coloring, or with --palette / --colors to keep palette-"
                         "derived cLIR + custom LIR."))
    p.add_argument("--single-color", action="store_true",
                   help=("'Solid' mode (LIVIA): LIR and cLIR use the same color per chain "
                         "(no light/dark distinction). Default is gradient (LIR = 50%% lightened cLIR). "
                         "Overridden by --lir-colors if both are given."))
    p.add_argument("--lir-lighten", type=float, default=0.5,
                   help=("amount to lighten the base color for LIR display (0.0..1.0; default 0.5). "
                         "Ignored in --single-color mode or when --lir-colors is given."))
    p.add_argument("--plddt", action="store_true",
                   help=("color the whole structure by per-residue pLDDT "
                         "(AlphaFold spectrum: red/orange < 50, yellow 50-70, "
                         "light blue 70-90, dark blue > 90) instead of LIVIA "
                         "LIR/cLIR per-chain coloring. Interface info (title, "
                         "2D labels, log echo, header summary) is still "
                         "emitted; pLDDT just replaces the cartoon color "
                         "phases. Useful for confidence-first views. "
                         "Mutually exclusive with --single-color / --colors / "
                         "--palette (palette settings are silently ignored)."))
    p.add_argument("--pairs",
                   help=("comma-separated list of fold names of interest "
                         "(matched against the CSV `name` column, exact match). "
                         "Only matching folds get a cxc. Useful when a CSV "
                         "contains many folds but you only want a few. "
                         "Combine with --pairs-file to extend the list."))
    p.add_argument("--pairs-file", type=Path,
                   help=("path to a text file with one fold name per line "
                         "(matched against the CSV `name` column, exact match). "
                         "Lines starting with `#` and blank lines are ignored. "
                         "Combine with --pairs to extend the list."))
    p.add_argument("--top-n", type=int, default=0,
                   help=("keep only the top-N folds by best iLIS across ranks "
                         "(0 = no top-N filter, default). Adds two passes over "
                         "the CSV(s): pass 1 builds per-fold max iLIS, pass 2 "
                         "emits cxc only for the top-N fold names. Combines "
                         "with --pairs / --pairs-file (union)."))
    p.add_argument("--best-rank-only", action="store_true",
                   help=("for each fold, keep only the (rank, model) with the "
                         "highest iLIS — 1 cxc per fold instead of one per "
                         "(fold, rank, model) group. Drastically reduces "
                         "clutter on multi-rank predictions (ColabFold 5 ranks, "
                         "AF3 5 models). For synthetic-rank platforms "
                         "(AF3/Boltz/Chai/OF3) the comparison uses the "
                         "iLIS-reordered rank, so 'rank 1' is always the most "
                         "confident model. Combines with --pairs / --pairs-file / "
                         "--top-n (applied after those filters)."))
    p.add_argument("--save-lir-pdb", nargs="?", const="", default=None, metavar="DIR",
                   help=("ALSO emit a partial PDB per cxc containing only "
                         "LIR-region atoms (gap-filled + short-segment-pruned "
                         "per-chain union). Useful as direct input to "
                         "Foldseek-multimer or any downstream search that "
                         "needs the confident interface only, not the whole "
                         "multimer. Bare `--save-lir-pdb` → written alongside "
                         "cxc files in --out. `--save-lir-pdb DIR` → written "
                         "to DIR (created if missing). Filename: <fold>__"
                         "rank<N>_model<M>_LIR.pdb. Output is always PDB "
                         "regardless of input format (PDB→PDB filter; "
                         "mmCIF→atom_site parse→PDB)."))
    p.add_argument("--save-lir-cxc", nargs="?", const="", default=None, metavar="DIR",
                   help=("ALSO emit a companion .cxc that opens the LIR-only "
                         "partial PDB directly, so the confident interface can be "
                         "viewed as a standalone structure in ChimeraX (not just "
                         "colored on the whole multimer). Implies writing that "
                         "partial PDB — see --save-lir-pdb. Same per-chain "
                         "LIR/cLIR coloring; residue numbers are preserved. Bare "
                         "`--save-lir-cxc` → written alongside cxc files in --out. "
                         "`--save-lir-cxc DIR` → written to DIR. Filename: "
                         "<fold>__rank<N>_model<M>_LIR.cxc."))
    p.add_argument("--lir-min-segment", type=int, default=5,
                   help=("drop any contiguous run of LIR residues shorter than "
                         "this length AFTER gap-filling (default 5). Applies "
                         "to BOTH the 3D cxc coloring (so isolated 1-4 residue "
                         "fragments don't float disconnected from the main "
                         "interface in the viewport) AND --save-lir-pdb output "
                         "(cleaner Foldseek-multimer queries). Pass 1 to "
                         "disable. cLIR residues whose LIR neighborhood was "
                         "pruned are also dropped (cLIR ⊆ LIR by definition)."))
    p.add_argument("--metadata", type=Path,
                   help=("TSV file with metadata to attach to each cxc + the "
                         "manifest. First column header must be `name` and "
                         "values must match the CSV `name` column. Remaining "
                         "columns are passed through verbatim — e.g. "
                         "`bait_symbol`, `prey_symbol`, `bait_uniprot`, "
                         "`organism`. Each matching cxc gets a `# Metadata:` "
                         "comment block in the header."))
    p.add_argument("--summary", nargs="?", const="summary.tsv", default=None,
                   metavar="FILE",
                   help=("opt-in TSV summary of emitted cxc files. Pass "
                         "`--summary` (bare) for the default name "
                         "`summary.tsv`, or `--summary <name>.tsv` for a "
                         "custom name. Written next to the cxc/ folder; one "
                         "row per chain-pair-in-cxc. Columns: source_csv, "
                         "name, rank (recomputed by iLIS for AF3/Boltz/Chai/"
                         "OF3 platforms), orig_rank, model, chain_i/j, "
                         "iLIS / LIS / cLIS / ipTM, len_i/j, pLDDT_i/j, "
                         "cLIR counts + compact ranges, gap_fill, "
                         "color_mode, structure_file (abs), cxc_file (abs), "
                         "plus any --metadata columns. Useful as the "
                         "input→output bridge for pandas filtering. "
                         "Default: no summary emitted (the lis.py CSV is "
                         "enough for most downstream work)."))
    p.add_argument("--no-manifest", action="store_true",
                   help=argparse.SUPPRESS)  # deprecated alias — no-op now
    p.add_argument("--ilis-cutoff", type=float, default=0.223,
                   help=("only emit chain-pair sections whose iLIS >= this value "
                         "(default 0.223 = FPR 10%% from the iLIS xlsx pipeline ladder; "
                         "stricter tiers: 0.339 = FPR 5%%, 0.551 = FPR 1%%; "
                         "pass 0.0 to disable filtering entirely and emit every pair lis.py produced). "
                         "A (fold,rank,model) with zero surviving pairs is skipped entirely. "
                         "NOTE: thresholds are calibrated on large-scale Y2H reference sets in yeast, "
                         "fly, and human predicted using ColabFold (AF2-multimer); see Kim et al., 2026 "
                         "for iLIS benchmark details. Different platforms (AF3, Boltz-2, Chai-1, "
                         "OpenFold3) may require different thresholds."))
    p.add_argument("--relative", action="store_true",
                   help=("write the structure path in the cxc's `open` line RELATIVE to "
                         "the cxc's own folder instead of absolute. ChimeraX resolves a "
                         "relative `open` against the command file's location, so this "
                         "makes the bundle portable across machines and OSes (e.g. generate "
                         "on Linux, view the cxc on Windows) as long as the cxc and its "
                         "structure move together. Default is absolute (reliable for "
                         "double-clicking on the machine that generated it). Applies to the "
                         "main cxc and the --save-lir-cxc companion."))
    p.add_argument("--dry-run", action="store_true",
                   help="report what would be emitted without writing files")
    args = p.parse_args()

    if args.list_palettes:
        print("Available palettes (one base color per chain; LIR derived by lightening unless --single-color):")
        for name, cols in PALETTES.items():
            n_unique = len(cols)
            sample = " ".join(cols[:4]) + ("  ..." if n_unique > 4 else "")
            print(f"  {name:40s}  {n_unique:2d} colors:  {sample}")
        print("\nOr supply custom hex colors directly via --colors '#hex1,#hex2,...'.")
        return

    if args.csv is None and args.csv_dir is None:
        p.error("one of --csv or --csv-dir is required (or use --list-palettes)")
    if args.pdb_root is None or args.out is None:
        p.error("--pdb-root and --out are required (unless --list-palettes)")

    # Resolve palette (cLIR base). `palette` is used when --clir-colors
    # is NOT given, and also as the legacy single-source-of-truth.
    if args.colors:
        palette = parse_colors_flag(args.colors, "--colors")
        palette_source = f"--colors ({len(palette)} custom)"
    else:
        # Accept the old `livia-*` prefix transparently.
        palette_name = PALETTE_ALIASES.get(args.palette, args.palette)
        if palette_name in PALETTES:
            palette = PALETTES[palette_name]
            palette_source = f"--palette {palette_name}"
            if palette_name != args.palette:
                palette_source += f" (alias of `{args.palette}`)"
        else:
            p.error(f"unknown palette '{args.palette}'. Run with --list-palettes to see options.")

    # Optional per-chain overrides for LIR / cLIR (fully independent).
    lir_palette = parse_colors_flag(args.lir_colors, "--lir-colors") if args.lir_colors else None
    clir_palette = parse_colors_flag(args.clir_colors, "--clir-colors") if args.clir_colors else None
    if clir_palette is not None:
        palette_source += f"  | --clir-colors ({len(clir_palette)})"
    if lir_palette is not None:
        palette_source += f"  | --lir-colors ({len(lir_palette)})"

    # Pair-of-interest filter. Union of --pairs (inline) and --pairs-file
    # (one per line, # and blanks ignored). Empty set = no filter, every
    # fold in the CSV is processed.
    pairs_of_interest: set[str] = set()
    if args.pairs:
        pairs_of_interest.update(p.strip() for p in args.pairs.split(",") if p.strip())
    if args.pairs_file:
        if not args.pairs_file.exists():
            sys.exit(f"!! --pairs-file not found: {args.pairs_file}")
        for line in args.pairs_file.read_text().splitlines():
            s = line.strip()
            if s and not s.startswith("#"):
                pairs_of_interest.add(s)

    # Metadata pass-through: TSV with header row; first column must be
    # `name` and match the CSV `name` column. Remaining columns are
    # passed through to the manifest + cxc headers verbatim.
    metadata: dict[str, dict[str, str]] = {}
    metadata_cols: list[str] = []
    if args.metadata:
        if not args.metadata.exists():
            sys.exit(f"!! --metadata file not found: {args.metadata}")
        with open(args.metadata) as f:
            md_reader = csv.DictReader(f, delimiter="\t")
            if not md_reader.fieldnames or md_reader.fieldnames[0] != "name":
                sys.exit(f"!! --metadata first column must be 'name', got "
                         f"{md_reader.fieldnames[:2] if md_reader.fieldnames else 'empty'}")
            metadata_cols = list(md_reader.fieldnames[1:])
            for r in md_reader:
                metadata[r["name"]] = {c: r.get(c, "") for c in metadata_cols}
        print(f"metadata: loaded {len(metadata)} entries with columns {metadata_cols}")

    # --top-N pre-scan: walk all input CSVs once to compute per-fold max
    # iLIS, sort, take top N fold names, union into pairs_of_interest.
    # Adds an extra read pass but lets the main loop stay simple.
    if args.top_n > 0:
        _csvs_for_scan = ([args.csv] if args.csv
                          else sorted(args.csv_dir.glob("*_lis_analysis.csv")))
        fold_max_ilis: dict[str, float] = {}
        for cp in _csvs_for_scan:
            for r in read_lis_csv(cp):
                name = r.get("name", "?")
                ilis = _float(r.get("iLIS"), 0.0)
                if ilis > fold_max_ilis.get(name, -1.0):
                    fold_max_ilis[name] = ilis
        top_n_names = [name for name, _ in
                       sorted(fold_max_ilis.items(), key=lambda kv: -kv[1])[:args.top_n]]
        pairs_of_interest.update(top_n_names)
        print(f"--top-n {args.top_n}: kept {len(top_n_names)} fold name(s) "
              f"from {len(fold_max_ilis)} total (max iLIS = "
              f"{fold_max_ilis[top_n_names[0]]:.3f}, "
              f"min kept iLIS = {fold_max_ilis[top_n_names[-1]]:.3f})")

    args.out.mkdir(parents=True, exist_ok=True)
    if args.lir_colors and args.clir_colors:
        mode = "independent (LIR + cLIR both explicit per chain)"
    elif args.lir_colors:
        mode = f"LIR explicit + cLIR from palette ({len(palette)})"
    elif args.clir_colors:
        mode = "cLIR explicit + LIR derived" + (
            f" (lighten {args.lir_lighten})" if not args.single_color else " (single-color: LIR = cLIR)"
        )
    elif args.single_color:
        mode = "single-color (LIR = cLIR)"
    else:
        mode = f"gradient (LIR = lighten(base, {args.lir_lighten}))"
    print(f"palette: {palette_source} → {len(palette)} base color(s)  |  mode: {mode}")
    if pairs_of_interest:
        print(f"pair filter: keeping only {len(pairs_of_interest)} fold name(s) of interest")
    if args.ilis_cutoff > 0:
        print(f"iLIS cutoff: emit only chain pairs with iLIS ≥ {args.ilis_cutoff}  (default FPR 10%; pass --ilis-cutoff 0 to disable)")
        print(f"  thresholds calibrated on Y2H reference sets in yeast/fly/human via ColabFold (Kim et al. 2026);")
        print(f"  AF3 / Boltz-2 / Chai-1 / OpenFold3 may require different thresholds — validate against known positives.")
    else:
        print(f"iLIS cutoff: none (every chain pair in the CSV gets a section)")

    csvs = [args.csv] if args.csv else sorted(args.csv_dir.glob("*_lis_analysis.csv"))
    if not csvs:
        sys.exit(f"!! no *_lis_analysis.csv found at {args.csv or args.csv_dir}")

    total_rows = total_groups = total_written = total_skipped_missing = 0
    total_pairs_filtered = total_groups_filtered = 0
    total_groups_not_in_pairs = 0
    all_seen_names: set[str] = set()
    manifest_rows: list[dict] = []  # one per chain-pair-in-cxc
    for csv_path in csvs:
        rows = read_lis_csv(csv_path)
        total_rows += len(rows)
        print(f"\n=== {csv_path.name}  ({len(rows)} rows) ===")

        # For non-ColabFold platforms, lis.py emits synthetic rank
        # (rank == model — no native confidence ordering). Recompute
        # rank per fold by sorting models on max-iLIS so rank=1 means
        # "best iLIS model" on AF3 / Boltz-2 / Chai-1 / OpenFold3.
        # ColabFold's native rank (from `_rank_NNN_` filename) differs
        # from model and is preserved.
        per_fold = {}
        for r in rows:
            per_fold.setdefault(r.get("name", "?"), []).append(r)
        recomputed_folds: set[str] = set()
        for fold, fold_rows in per_fold.items():
            synthetic = all(str(r.get("rank")) == str(r.get("model"))
                            for r in fold_rows)
            if not synthetic:
                continue
            # Score each (rank, model) group by its best iLIS, then
            # reassign ranks 1..N descending.
            by_key = {}
            for r in fold_rows:
                k = (str(r.get("rank")), str(r.get("model")))
                by_key.setdefault(k, []).append(r)
            ordered = sorted(by_key.keys(),
                             key=lambda k: -max(_float(r.get("iLIS"), 0.0)
                                                for r in by_key[k]))
            for new_rank, k in enumerate(ordered, start=1):
                for r in by_key[k]:
                    r["_orig_rank"] = r["rank"]  # remember CSV value
                    r["rank"] = str(new_rank)
            recomputed_folds.add(fold)
        if recomputed_folds:
            print(f"  recomputed rank by iLIS for {len(recomputed_folds)} fold(s) "
                  f"with synthetic rank=model (AF3/Boltz/Chai/OF3 style).")

        # Group rows by (name, rank, model) — one cxc per group.
        groups: dict[tuple[str, str, str], list[dict]] = {}
        for r in rows:
            key = (r.get("name", "?"), r.get("rank", "?"), r.get("model", "?"))
            groups.setdefault(key, []).append(r)
            all_seen_names.add(key[0])

        # --best-rank-only: keep only the (rank, model) per fold whose
        # max iLIS is highest. Applied AFTER rank reassignment above so
        # the comparison uses the iLIS-reordered rank for synthetic-rank
        # platforms (AF3/Boltz/Chai/OF3).
        if args.best_rank_only:
            best_per_fold: dict[str, tuple[tuple[str, str, str], float]] = {}
            for key, grp in groups.items():
                fold = key[0]
                max_ilis = max(_float(r.get("iLIS"), 0.0) for r in grp)
                if fold not in best_per_fold or max_ilis > best_per_fold[fold][1]:
                    best_per_fold[fold] = (key, max_ilis)
            kept_keys = {kv[0] for kv in best_per_fold.values()}
            dropped_n = len(groups) - len(kept_keys)
            groups = {k: groups[k] for k in kept_keys}
            print(f"  --best-rank-only: kept best (rank, model) per fold "
                  f"({len(groups)} folds, dropped {dropped_n} lower-ranked groups)")
        total_groups += len(groups)

        for (name, rank, model), grp in sorted(groups.items()):
            # Pair-of-interest filter — skip silently if a list was provided
            # and this fold isn't on it. Tracked separately from the iLIS
            # cutoff so the summary footer reports "not in pair list" vs
            # "below cutoff" as distinct skip reasons.
            if pairs_of_interest and name not in pairs_of_interest:
                total_groups_not_in_pairs += 1
                continue
            # Apply iLIS cutoff at the per-pair level: drop weak pairs but
            # keep the (fold,rank,model) group around if at least one pair
            # survives. If none do, skip the whole group — no point loading
            # a structure with no interface to highlight.
            if args.ilis_cutoff > 0:
                kept = [r for r in grp if _float(r.get("iLIS"), 0.0) >= args.ilis_cutoff]
                dropped = len(grp) - len(kept)
                if dropped:
                    total_pairs_filtered += dropped
                if not kept:
                    total_groups_filtered += 1
                    print(f"  - skip {name}  rank={rank}  model={model}  (all {len(grp)} pair{'s' if len(grp)!=1 else ''} < iLIS {args.ilis_cutoff})")
                    continue
                grp = kept

            # All rows in a group share the same structure_file. lis.py
            # records the **basename** only — the actual file may live
            # in a nested subfolder under --pdb-root, especially for
            # Boltz-2 / Chai-1 / ColabFold / OpenFold3 extracted zips
            # (e.g. `result-X/predictions/result/result_model_0.pdb`).
            # Try in order:
            #   1. Direct `pdb_root / sf`        (AF3 flat extraction)
            #   2. Recursive `pdb_root.rglob(sf)` (everything else)
            #   3. Convention-based candidates    (no structure_file in CSV)
            #   4. Fuzzy recursive glob           (last-ditch)
            sf = grp[0].get("structure_file", "")
            cif_path = args.pdb_root / sf if sf else None
            if cif_path is not None and not cif_path.exists():
                # Recursive basename search — handles every platform's
                # nested zip layout without needing to know its structure.
                recursive = list(args.pdb_root.rglob(sf)) if sf else []
                cif_path = recursive[0] if recursive else None
            if cif_path is None or not cif_path.exists():
                # Most AF3 model files: `<name>_model_<N>.cif`. Extract
                # numeric model index regardless of whether `model` is
                # "0" or "model_0".
                m_idx = re.search(r"(\d+)", str(model))
                model_n = m_idx.group(1) if m_idx else str(model)
                candidates = [
                    f"{name}_model_{model_n}.cif",
                    f"{name}_model_{model_n}.pdb",
                ]
                cif_path = None
                for cand in candidates:
                    # Loop var renamed from `p` to avoid shadowing the
                    # argparse parser bound at the top of main().
                    cand_path = args.pdb_root / cand
                    if cand_path.exists():
                        cif_path = cand_path
                        break
                if cif_path is None:
                    # Last-ditch fuzzy recursive glob.
                    matches = list(args.pdb_root.rglob(f"*{name}*model*{model_n}*.cif")) + \
                              list(args.pdb_root.rglob(f"*{name}*model*{model_n}*.pdb")) + \
                              list(args.pdb_root.rglob(f"*{name}*rank*{rank}*.pdb")) + \
                              list(args.pdb_root.rglob(f"*{name}*rank*{rank}*.pdb.gz"))
                    cif_path = matches[0] if matches else None
            if cif_path is None or not cif_path.exists():
                total_skipped_missing += 1
                print(f"  ! missing structure for {name}  rank={rank}  model={model}  (looked for: {sf!r}, recursive {sf!r}, or {name}_model_*.cif)")
                continue

            # Normalize model index for the filename: handle both raw int
            # ("0") and LIVIA-web string ("model_0") forms cleanly so we
            # don't end up with `__rank0_modelmodel_0_interface.cxc`.
            m_idx_match = re.search(r"(\d+)", str(model))
            model_clean = m_idx_match.group(1) if m_idx_match else str(model)
            out_name = f"{name}__rank{rank}_model{model_clean}_interface.cxc"
            out_path = args.out / out_name

            if args.dry_run:
                print(f"  [dry-run] would emit {out_name}  ({len(grp)} chain pairs)")
                continue

            md_row = metadata.get(name, {}) if metadata else {}
            wrote = emit_cxc(grp, cif_path, out_path,
                             gap_fill=args.gap_fill,
                             palette=palette,
                             single_color=args.single_color,
                             lir_lighten=args.lir_lighten,
                             lir_palette=lir_palette,
                             clir_palette=clir_palette,
                             min_segment=args.lir_min_segment,
                             csv_path=csv_path,
                             metadata=md_row if md_row else None,
                             plddt=args.plddt,
                             relative=args.relative)
            if wrote:
                total_written += 1
                print(f"  ✓ {out_name}  ({len(grp)} chain pair{'s' if len(grp) != 1 else ''})")
                # --save-lir-pdb: emit a partial structure with only the
                # gap-filled + pruned LIR atoms. Foldseek-multimer / etc.
                # use this as the "confident interface" portable input.
                # --save-lir-cxc implies writing the partial PDB it opens, so the
                # partial is extracted when EITHER flag is set.
                if args.save_lir_pdb is not None or args.save_lir_cxc is not None:
                    lir_dir = (Path(args.save_lir_pdb) if args.save_lir_pdb
                               else Path(args.save_lir_cxc) if args.save_lir_cxc
                               else args.out)
                    lir_dir.mkdir(parents=True, exist_ok=True)
                    lir_pdb_name = f"{name}__rank{rank}_model{model_clean}_LIR.pdb"
                    lir_pdb_path = lir_dir / lir_pdb_name
                    ok = extract_lir_partial_pdb(
                        cif_path, grp,
                        gap_fill=args.gap_fill,
                        out_path=lir_pdb_path,
                        min_segment=args.lir_min_segment,
                    )
                    if ok:
                        print(f"    + LIR-only partial PDB: {lir_pdb_name}")
                        # --save-lir-cxc: companion cxc that opens the partial PDB
                        # directly (standalone interface view). emit_cxc uses the
                        # structure only as the `open` target; all LIR/cLIR ranges
                        # come from the CSV rows and are preserved in the partial.
                        if args.save_lir_cxc is not None:
                            lir_cxc_name = f"{name}__rank{rank}_model{model_clean}_LIR.cxc"
                            lir_cxc_path = lir_dir / lir_cxc_name
                            emit_cxc(grp, lir_pdb_path, lir_cxc_path,
                                     gap_fill=args.gap_fill,
                                     palette=palette,
                                     single_color=args.single_color,
                                     lir_lighten=args.lir_lighten,
                                     lir_palette=lir_palette,
                                     clir_palette=clir_palette,
                                     min_segment=args.lir_min_segment,
                                     csv_path=csv_path,
                                     metadata=md_row if md_row else None,
                                     plddt=args.plddt,
                                     relative=args.relative)
                            print(f"    + LIR-only cxc: {lir_cxc_name}")
                    else:
                        print(f"    ! could not extract LIR PDB for {out_name} "
                              f"(unsupported format or empty LIR after pruning)")
                # Collect one manifest row per chain-pair in this group.
                for r in grp:
                    ci = (r.get("chain_i") or "?").strip()
                    cj = (r.get("chain_j") or "?").strip()
                    clir_i = parse_residue_indices(r.get("cLIR_indices_i", ""))
                    clir_j = parse_residue_indices(r.get("cLIR_indices_j", ""))
                    row = {
                        "source_csv":      csv_path.name,
                        "name":            name,
                        "rank":            rank,
                        "model":           model,
                        "chain_i":         ci,
                        "chain_j":         cj,
                        "iLIS":            f"{_float(r.get('iLIS'), 0.0):.4f}",
                        "LIS":             f"{_float(r.get('LIS'), 0.0):.4f}",
                        "cLIS":            f"{_float(r.get('cLIS'), 0.0):.4f}",
                        "ipTM":            f"{_float(r.get('ipTM'), 0.0):.4f}",
                        "len_i":           r.get("len_i", ""),
                        "len_j":           r.get("len_j", ""),
                        "pLDDT_i":         r.get("pLDDT_i", ""),
                        "pLDDT_j":         r.get("pLDDT_j", ""),
                        "cLIR_i_n":        len(clir_i),
                        "cLIR_j_n":        len(clir_j),
                        "cLIR_i_ranges":   positions_to_compact_ranges(clir_i),
                        "cLIR_j_ranges":   positions_to_compact_ranges(clir_j),
                        "orig_rank":       r.get("_orig_rank", ""),
                        "gap_fill":        args.gap_fill,
                        "color_mode":      ("plddt" if args.plddt else
                                            "single-color" if args.single_color else
                                            "gradient"),
                        "structure_file":  str(cif_path.resolve()),
                        "cxc_file":        str(out_path.resolve()),
                    }
                    # Append any metadata columns for this name.
                    for c in metadata_cols:
                        row[c] = md_row.get(c, "")
                    manifest_rows.append(row)

    # Write summary TSV alongside the cxc/ folder — only if --summary
    # was passed. (Default is no summary; the lis.py CSV is usually
    # sufficient as the input-side analysis table.)
    if args.summary and manifest_rows:
        summary_path = args.out / args.summary
        fieldnames = list(manifest_rows[0].keys())
        with open(summary_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(manifest_rows)
        print(f"\nwrote summary: {summary_path}  ({len(manifest_rows)} rows × {len(fieldnames)} cols)")

    print(f"\n=== summary ===")
    print(f"  CSV files processed:           {len(csvs)}")
    print(f"  total rows (chain-pair × rank × model): {total_rows}")
    print(f"  (fold, rank, model) groups:    {total_groups}")
    print(f"  cxc files written:             {total_written}")
    if pairs_of_interest:
        print(f"  groups skipped (not in --pairs/--pairs-file): {total_groups_not_in_pairs}")
        unmatched = sorted(pairs_of_interest - all_seen_names)
        if unmatched:
            print(f"  WARN: {len(unmatched)} fold name(s) in --pairs not found in any CSV:")
            for name in unmatched[:10]:
                print(f"    - {name}")
            if len(unmatched) > 10:
                print(f"    ... and {len(unmatched) - 10} more")
    if total_skipped_missing:
        print(f"  groups skipped (no structure): {total_skipped_missing}")
    if args.ilis_cutoff > 0:
        print(f"  pairs filtered (iLIS < {args.ilis_cutoff}):    {total_pairs_filtered}")
        print(f"  groups skipped (all pairs < cutoff): {total_groups_filtered}")


if __name__ == "__main__":
    main()
