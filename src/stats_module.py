"""
stats_module.py — Chemical shift statistics, reference database, and predictor.

Provides:
- compute_shift_stats(): aggregate shifts by (residue, atom, ss_class)
- ShiftPredictor: given PDB + DSSP, predict expected chemical shifts
- random_coil_reference: built-in RC values for deviation analysis
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional, Dict, List, Tuple


# ---------------------------------------------------------------------------
# Random Coil Reference Values (Wishart 1995 / BMRB averages)
# These are the baseline — deviations from these indicate structure.
# ---------------------------------------------------------------------------

RANDOM_COIL = {
    # residue -> {atom -> rc_shift_ppm}
    'A': {'CA': 52.3, 'CB': 19.1, 'N': 123.8, 'H': 8.25, 'C': 177.8},
    'R': {'CA': 56.1, 'CB': 30.9, 'N': 120.5, 'H': 8.27, 'C': 176.3},
    'N': {'CA': 53.1, 'CB': 38.9, 'N': 118.7, 'H': 8.38, 'C': 175.2},
    'D': {'CA': 54.2, 'CB': 41.1, 'N': 120.4, 'H': 8.37, 'C': 176.3},
    'C': {'CA': 58.2, 'CB': 28.0, 'N': 118.8, 'H': 8.32, 'C': 174.6},
    'Q': {'CA': 55.8, 'CB': 29.4, 'N': 119.8, 'H': 8.27, 'C': 175.9},
    'E': {'CA': 56.6, 'CB': 30.2, 'N': 120.2, 'H': 8.36, 'C': 176.6},
    'G': {'CA': 45.1, 'CB': None, 'N': 108.8, 'H': 8.33, 'C': 173.8},
    'H': {'CA': 55.9, 'CB': 29.4, 'N': 118.2, 'H': 8.41, 'C': 174.1},
    'I': {'CA': 61.1, 'CB': 38.8, 'N': 120.9, 'H': 8.22, 'C': 175.8},
    'L': {'CA': 55.2, 'CB': 42.4, 'N': 121.8, 'H': 8.16, 'C': 177.6},
    'K': {'CA': 56.5, 'CB': 33.1, 'N': 120.4, 'H': 8.25, 'C': 176.6},
    'M': {'CA': 55.6, 'CB': 33.1, 'N': 119.6, 'H': 8.28, 'C': 176.3},
    'F': {'CA': 57.7, 'CB': 39.6, 'N': 120.3, 'H': 8.30, 'C': 175.8},
    'P': {'CA': 63.3, 'CB': 32.1, 'N': 136.5, 'H': None, 'C': 177.3},
    'S': {'CA': 58.3, 'CB': 63.8, 'N': 115.7, 'H': 8.31, 'C': 174.6},
    'T': {'CA': 61.8, 'CB': 69.8, 'N': 113.6, 'H': 8.24, 'C': 174.7},
    'W': {'CA': 57.5, 'CB': 29.9, 'N': 121.3, 'H': 8.18, 'C': 176.1},
    'Y': {'CA': 57.9, 'CB': 38.8, 'N': 120.3, 'H': 8.18, 'C': 175.9},
    'V': {'CA': 62.2, 'CB': 32.9, 'N': 119.9, 'H': 8.19, 'C': 176.3},
}


# ---------------------------------------------------------------------------
# Statistics computation
# ---------------------------------------------------------------------------

def compute_shift_stats(df: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate chemical shift statistics grouped by (residue, atom, ss_class).

    Input df must have columns: residue, atom, shift, ss_class
    Returns df with: residue | atom | ss_class | count | mean | median | std | min | max
    """
    required = {'residue', 'atom', 'shift'}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"DataFrame missing columns: {missing}")

    if 'ss_class' not in df.columns:
        df = df.copy()
        df['ss_class'] = 'unknown'

    stats = (
        df.groupby(['residue', 'atom', 'ss_class'])['shift']
        .agg(
            count='count',
            mean='mean',
            median='median',
            std='std',
            min='min',
            max='max',
        )
        .reset_index()
    )
    return stats.sort_values(['residue', 'atom', 'ss_class']).reset_index(drop=True)


def add_random_coil_deviation(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add a 'rc_deviation' column = observed shift - random coil reference.
    Positive CA deviation → helical tendency.
    Negative CA deviation → beta-sheet tendency.
    """
    df = df.copy()
    deviations = []
    for _, row in df.iterrows():
        res  = row['residue']
        atom = row['atom']
        shift = row.get('shift') or row.get('mean')
        rc = RANDOM_COIL.get(res, {}).get(atom)
        if rc is not None and shift is not None:
            deviations.append(float(shift) - float(rc))
        else:
            deviations.append(None)
    df['rc_deviation'] = deviations
    return df


# ---------------------------------------------------------------------------
# Shift predictor
# ---------------------------------------------------------------------------

class ShiftPredictor:
    """
    Predict NMR chemical shifts for a protein given its sequence + secondary structure.

    Usage:
        predictor = ShiftPredictor(stats_df)
        predictions = predictor.predict(sequence, ss_map)
    """

    def __init__(self, stats_df: pd.DataFrame):
        """
        stats_df: output of compute_shift_stats() — the reference database.
        """
        self.stats = stats_df
        self._build_lookup()

    def _build_lookup(self):
        """Pre-index stats for fast lookup."""
        self._lookup = {}
        for _, row in self.stats.iterrows():
            key = (row['residue'], row['atom'], row['ss_class'])
            self._lookup[key] = {
                'mean':   row['mean'],
                'median': row.get('median', row['mean']),
                'std':    row.get('std', 1.0),
                'count':  row.get('count', 0),
            }

    def predict_shift(
        self,
        residue: str,
        atom: str,
        ss_class: str,
        method: str = 'mean'
    ) -> Optional[float]:
        """
        Predict shift for one (residue, atom, ss_class) triple.
        method: 'mean', 'median', or 'sample' (draws from distribution)
        Falls back to random coil if not in stats.
        """
        key = (residue, atom, ss_class)
        entry = self._lookup.get(key)

        if entry is None:
            # fallback: ignore ss_class
            for ss in ['coil', 'unknown']:
                entry = self._lookup.get((residue, atom, ss))
                if entry:
                    break

        if entry is None:
            # last resort: random coil
            rc = RANDOM_COIL.get(residue, {}).get(atom)
            return rc

        if method == 'mean':
            return entry['mean']
        elif method == 'median':
            return entry['median']
        elif method == 'sample':
            std = entry['std'] or 1.0
            return float(np.random.normal(entry['mean'], std))
        else:
            raise ValueError(f"Unknown method: {method}")

    def predict(
        self,
        sequence: str,
        ss_map: Dict[int, str],
        atoms: List[str] = None,
        method: str = 'mean'
    ) -> pd.DataFrame:
        """
        Predict shifts for every (residue, atom) in a full sequence.

        sequence: one-letter amino acid string
        ss_map: {residue_number (1-based) -> ss_class}
        atoms: list of atom types to predict (default: CA, CB, N)
        method: 'mean', 'median', or 'sample'

        Returns DataFrame: seq_id | residue | atom | ss_class | predicted_shift | rc_ref | rc_deviation
        """
        if atoms is None:
            atoms = ['CA', 'CB', 'N']

        rows = []
        for i, aa in enumerate(sequence, 1):
            ss = ss_map.get(i, 'coil')
            for atom in atoms:
                pred = self.predict_shift(aa, atom, ss, method=method)
                rc   = RANDOM_COIL.get(aa, {}).get(atom)
                dev  = (pred - rc) if (pred is not None and rc is not None) else None
                rows.append({
                    'seq_id':          i,
                    'residue':         aa,
                    'atom':            atom,
                    'ss_class':        ss,
                    'predicted_shift': pred,
                    'rc_ref':          rc,
                    'rc_deviation':    dev,
                })

        return pd.DataFrame(rows)

    def predict_uncertainty(
        self,
        sequence: str,
        ss_map: Dict[int, str],
        n_samples: int = 100,
        atoms: List[str] = None,
    ) -> pd.DataFrame:
        """
        Predict shifts with uncertainty estimates via Monte Carlo sampling.
        Returns mean ± std of sampled predictions.
        """
        if atoms is None:
            atoms = ['CA', 'CB', 'N']

        all_samples = {(i+1, aa, atom): [] for i, aa in enumerate(sequence) for atom in atoms}

        for _ in range(n_samples):
            for i, aa in enumerate(sequence, 1):
                ss = ss_map.get(i, 'coil')
                for atom in atoms:
                    val = self.predict_shift(aa, atom, ss, method='sample')
                    all_samples[(i, aa, atom)].append(val)

        rows = []
        for (seq_id, aa, atom), samples in all_samples.items():
            valid = [s for s in samples if s is not None]
            if not valid:
                continue
            rows.append({
                'seq_id':       seq_id,
                'residue':      aa,
                'atom':         atom,
                'ss_class':     ss_map.get(seq_id, 'coil'),
                'pred_mean':    np.mean(valid),
                'pred_std':     np.std(valid),
                'pred_95_lo':   np.percentile(valid, 2.5),
                'pred_95_hi':   np.percentile(valid, 97.5),
            })

        return pd.DataFrame(rows).sort_values(['seq_id', 'atom']).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Assignment scoring (Phase 5 groundwork)
# ---------------------------------------------------------------------------

def score_peak_assignment(
    peak_shift: float,
    residue: str,
    atom: str,
    ss_class: str,
    stats_df: pd.DataFrame,
) -> float:
    """
    Compute likelihood score for assigning a peak to (residue, atom, ss_class).
    Uses Gaussian probability: p(shift | residue, atom, ss_class).
    Returns log-probability (higher = better).
    """
    match = stats_df[
        (stats_df['residue'] == residue) &
        (stats_df['atom'] == atom) &
        (stats_df['ss_class'] == ss_class)
    ]
    if match.empty:
        # Fallback to RC reference with broad sigma
        rc = RANDOM_COIL.get(residue, {}).get(atom)
        if rc is None:
            return -np.inf
        sigma = 2.0
        return -0.5 * ((peak_shift - rc) / sigma) ** 2

    mean = match.iloc[0]['mean']
    std  = max(match.iloc[0]['std'], 0.3)  # floor std to avoid division issues
    log_p = -0.5 * ((peak_shift - mean) / std) ** 2 - np.log(std)
    return float(log_p)


def rank_assignments(
    peak_shift: float,
    atom: str,
    sequence: str,
    ss_map: Dict[int, str],
    stats_df: pd.DataFrame,
    top_n: int = 5,
) -> pd.DataFrame:
    """
    Given an experimental peak at peak_shift ppm for a given atom type,
    rank all possible residue assignments by likelihood.

    Returns top_n candidates with columns:
        seq_id | residue | ss_class | log_prob | rank
    """
    candidates = []
    for i, aa in enumerate(sequence, 1):
        ss = ss_map.get(i, 'coil')
        score = score_peak_assignment(peak_shift, aa, atom, ss, stats_df)
        candidates.append({
            'seq_id':   i,
            'residue':  aa,
            'ss_class': ss,
            'log_prob': score,
        })

    df = pd.DataFrame(candidates)
    df = df.sort_values('log_prob', ascending=False).reset_index(drop=True)
    df['rank'] = df.index + 1
    return df.head(top_n)
