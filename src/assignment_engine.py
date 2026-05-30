"""
assignment_engine.py — Phase 5: Automated NMR Peak Assignment Engine

Solves the peak assignment problem:
  Given a list of experimental 13C peak positions (ppm),
  find the best one-to-one mapping to (residue, atom, ss_class) slots
  in the protein sequence.

Algorithm
---------
1. Build a cost matrix: rows = experimental peaks, cols = residue-atom slots
   Cost = -log_probability(peak | residue, atom, ss_class)  [Gaussian model]

2. Apply hard constraints:
   - Sequence composition: count of each residue type limits assignments
   - Atom-type filtering: CA peaks only compete against CA slots, etc.
   - Impossible assignments get cost = +inf

3. Solve with the Hungarian algorithm (scipy.optimize.linear_sum_assignment)
   → optimal minimum-cost matching in O(n^3)

4. For multi-atom joint assignment (CA+CB pairs):
   - Build a joint cost using both dimensions simultaneously
   - Dramatically improves specificity (two coordinates >> one)

5. Output: assignment table with confidence scores per peak

Usage
-----
    from assignment_engine import AssignmentEngine
    engine = AssignmentEngine(stats_df, sequence, ss_map)
    result = engine.assign(experimental_peaks, atom='CA')
    result = engine.assign_joint(ca_peaks, cb_peaks)  # 2D CA-CB matching
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from scipy.optimize import linear_sum_assignment
from collections import Counter


# ---------------------------------------------------------------------------
# Random coil references (duplicated here for self-contained module)
# ---------------------------------------------------------------------------

RANDOM_COIL = {
    'A': {'CA': 52.3, 'CB': 19.1, 'N': 123.8},
    'C': {'CA': 58.2, 'CB': 28.0, 'N': 118.8},
    'D': {'CA': 54.2, 'CB': 41.1, 'N': 120.4},
    'E': {'CA': 56.6, 'CB': 30.2, 'N': 120.2},
    'F': {'CA': 57.7, 'CB': 39.6, 'N': 120.3},
    'G': {'CA': 45.1, 'CB': None, 'N': 108.8},
    'H': {'CA': 55.9, 'CB': 29.4, 'N': 118.2},
    'I': {'CA': 61.1, 'CB': 38.8, 'N': 120.9},
    'K': {'CA': 56.5, 'CB': 33.1, 'N': 120.4},
    'L': {'CA': 55.2, 'CB': 42.4, 'N': 121.8},
    'M': {'CA': 55.6, 'CB': 33.1, 'N': 119.6},
    'N': {'CA': 53.1, 'CB': 38.9, 'N': 118.7},
    'P': {'CA': 63.3, 'CB': 32.1, 'N': 136.5},
    'Q': {'CA': 55.8, 'CB': 29.4, 'N': 119.8},
    'R': {'CA': 56.1, 'CB': 30.9, 'N': 120.5},
    'S': {'CA': 58.3, 'CB': 63.8, 'N': 115.7},
    'T': {'CA': 61.8, 'CB': 69.8, 'N': 113.6},
    'V': {'CA': 62.2, 'CB': 32.9, 'N': 119.9},
    'W': {'CA': 57.5, 'CB': 29.9, 'N': 121.3},
    'Y': {'CA': 57.9, 'CB': 38.8, 'N': 120.3},
}

INF_COST = 1e9   # cost for forbidden assignments


# ---------------------------------------------------------------------------
# Assignment slot: one (seq_id, residue, atom, ss_class) position
# ---------------------------------------------------------------------------

class AssignmentSlot:
    __slots__ = ['seq_id', 'residue', 'atom', 'ss_class']

    def __init__(self, seq_id: int, residue: str, atom: str, ss_class: str):
        self.seq_id   = seq_id
        self.residue  = residue
        self.atom     = atom
        self.ss_class = ss_class

    def __repr__(self):
        return f"{self.residue}{self.seq_id}-{self.atom}({self.ss_class})"


# ---------------------------------------------------------------------------
# Core engine
# ---------------------------------------------------------------------------

class AssignmentEngine:
    """
    Automated NMR peak assignment engine.

    Parameters
    ----------
    stats_df : DataFrame
        Reference database from pipeline (residue | atom | ss_class | mean | std | count)
    sequence : str
        One-letter amino acid sequence
    ss_map : dict
        {seq_id (1-based) -> ss_class ('helix'/'strand'/'coil')}
    """

    def __init__(
        self,
        stats_df: pd.DataFrame,
        sequence: str,
        ss_map: Dict[int, str],
    ):
        self.stats_df = stats_df
        self.sequence = sequence
        self.ss_map   = ss_map
        self._build_lookup()
        self._build_slots()
        self._composition = Counter(sequence)

    def _build_lookup(self):
        """Pre-index stats_df for O(1) lookup."""
        self._lookup: Dict[tuple, dict] = {}
        for _, row in self.stats_df.iterrows():
            key = (row['residue'], row['atom'], row['ss_class'])
            self._lookup[key] = {
                'mean':  float(row['mean']),
                'std':   max(float(row.get('std', 1.0) or 1.0), 0.3),
                'count': int(row.get('count', 1)),
            }

    def _build_slots(self):
        """Build all (seq_id, residue, atom, ss_class) slots from sequence + ss_map."""
        self._slots_by_atom: Dict[str, List[AssignmentSlot]] = {}
        for i, aa in enumerate(self.sequence, 1):
            ss = self.ss_map.get(i, 'coil')
            for atom in ['CA', 'CB', 'N', 'H']:
                if aa == 'G' and atom == 'CB':
                    continue  # Gly has no CB
                if aa == 'P' and atom in ('N', 'H'):
                    continue  # Pro has no amide NH
                slot = AssignmentSlot(i, aa, atom, ss)
                self._slots_by_atom.setdefault(atom, []).append(slot)

    def _peak_slot_cost(self, peak_ppm: float, slot: AssignmentSlot) -> float:
        """
        Cost = -log P(peak | slot distribution).
        Lower cost = better match.
        Uses Gaussian log-probability with fallback to RC reference.
        """
        key = (slot.residue, slot.atom, slot.ss_class)
        entry = self._lookup.get(key)

        if entry is None:
            # Fallback: try coil, then unknown
            for ss_fb in ('coil', 'unknown'):
                entry = self._lookup.get((slot.residue, slot.atom, ss_fb))
                if entry:
                    break

        if entry is None:
            # Last resort: random coil reference with broad sigma
            rc = RANDOM_COIL.get(slot.residue, {}).get(slot.atom)
            if rc is None:
                return INF_COST
            sigma = 2.0
            return 0.5 * ((peak_ppm - rc) / sigma) ** 2

        mean  = entry['mean']
        sigma = entry['std']
        # Negative log-likelihood of Gaussian (ignoring constant terms)
        return 0.5 * ((peak_ppm - mean) / sigma) ** 2 + np.log(sigma)

    # ------------------------------------------------------------------ #
    # 1D assignment (single atom type)
    # ------------------------------------------------------------------ #

    def build_cost_matrix(
        self,
        experimental_peaks: List[float],
        atom: str,
        max_ppm_deviation: float = 8.0,
    ) -> Tuple[np.ndarray, List[AssignmentSlot]]:
        """
        Build cost matrix for Hungarian algorithm.

        Rows    = experimental peaks (n_peaks)
        Columns = residue-atom slots (n_slots)

        Returns (cost_matrix, slots_list)
        """
        slots = self._slots_by_atom.get(atom, [])
        if not slots:
            raise ValueError(f"No slots found for atom {atom}. Check sequence and ss_map.")

        n_peaks = len(experimental_peaks)
        n_slots = len(slots)

        cost = np.full((n_peaks, n_slots), INF_COST, dtype=float)

        for i, peak in enumerate(experimental_peaks):
            for j, slot in enumerate(slots):
                # Hard cutoff: peak too far from any reasonable shift for this atom
                rc = RANDOM_COIL.get(slot.residue, {}).get(slot.atom)
                if rc is not None and abs(peak - rc) > max_ppm_deviation:
                    continue  # leave as INF_COST
                cost[i, j] = self._peak_slot_cost(peak, slot)

        return cost, slots

    def assign(
        self,
        experimental_peaks: List[float],
        atom: str = 'CA',
        max_ppm_deviation: float = 8.0,
        min_confidence: float = 0.0,
    ) -> pd.DataFrame:
        """
        Assign experimental peaks to residue-atom slots (1D, single atom type).

        Parameters
        ----------
        experimental_peaks : list of ppm values (from experimental spectrum)
        atom : 'CA', 'CB', or 'N'
        max_ppm_deviation : peaks further than this from RC reference are forbidden
        min_confidence : filter output to rows with confidence >= this value (0–1)

        Returns
        -------
        DataFrame with columns:
            peak_id | exp_shift | seq_id | residue | atom | ss_class |
            predicted_shift | delta_ppm | cost | confidence | status
        """
        peaks = sorted(experimental_peaks)  # sort for cleaner output
        cost_matrix, slots = self.build_cost_matrix(peaks, atom, max_ppm_deviation)

        n_peaks, n_slots = cost_matrix.shape

        # Pad to square if needed (Hungarian requires square or more cols than rows)
        if n_peaks > n_slots:
            # More peaks than slots — some peaks will be unassigned (dummy columns)
            pad = np.full((n_peaks, n_peaks - n_slots), INF_COST * 0.9)
            cost_padded = np.hstack([cost_matrix, pad])
            dummy_slots = n_peaks - n_slots
        else:
            cost_padded = cost_matrix
            dummy_slots = 0

        # Run Hungarian algorithm
        row_ind, col_ind = linear_sum_assignment(cost_padded)

        rows = []
        for r, c in zip(row_ind, col_ind):
            peak = peaks[r]
            is_dummy = (dummy_slots > 0 and c >= n_slots)

            if is_dummy or cost_padded[r, c] >= INF_COST * 0.8:
                rows.append({
                    'peak_id':          r + 1,
                    'exp_shift':        round(peak, 3),
                    'seq_id':           None,
                    'residue':          None,
                    'atom':             atom,
                    'ss_class':         None,
                    'predicted_shift':  None,
                    'delta_ppm':        None,
                    'cost':             None,
                    'confidence':       0.0,
                    'status':           'unassigned',
                })
                continue

            slot = slots[c]
            raw_cost = cost_matrix[r, c]

            # Compute predicted shift for this slot
            key = (slot.residue, slot.atom, slot.ss_class)
            entry = self._lookup.get(key)
            if entry is None:
                for ss_fb in ('coil', 'unknown'):
                    entry = self._lookup.get((slot.residue, slot.atom, ss_fb))
                    if entry: break
            pred_shift = entry['mean'] if entry else RANDOM_COIL.get(slot.residue, {}).get(slot.atom)
            delta = round(peak - pred_shift, 3) if pred_shift is not None else None

            rows.append({
                'peak_id':          r + 1,
                'exp_shift':        round(peak, 3),
                'seq_id':           slot.seq_id,
                'residue':          slot.residue,
                'atom':             atom,
                'ss_class':         slot.ss_class,
                'predicted_shift':  round(pred_shift, 3) if pred_shift else None,
                'delta_ppm':        delta,
                'cost':             round(raw_cost, 4),
                'confidence':       None,   # filled below
                'status':           'assigned',
            })

        result = pd.DataFrame(rows)

        # Compute confidence scores (convert cost to 0–1 probability)
        assigned = result[result['status'] == 'assigned'].copy()
        if not assigned.empty:
            costs = assigned['cost'].values.astype(float)
            # Softmax-style: lower cost → higher confidence
            exp_neg = np.exp(-costs)
            # Normalize against all assigned peaks
            total = exp_neg.sum()
            confidences = exp_neg / total * len(exp_neg) / len(result)
            confidences = np.clip(confidences, 0, 1)
            result.loc[result['status'] == 'assigned', 'confidence'] = np.round(confidences, 3)

        result.loc[result['status'] == 'unassigned', 'confidence'] = 0.0

        if min_confidence > 0:
            result = result[result['confidence'] >= min_confidence].reset_index(drop=True)

        return result.sort_values('seq_id', na_position='last').reset_index(drop=True)

    # ------------------------------------------------------------------ #
    # 2D joint assignment (CA + CB simultaneously)
    # ------------------------------------------------------------------ #

    def assign_joint(
        self,
        ca_peaks: List[float],
        cb_peaks: List[float],
        max_ppm_deviation: float = 8.0,
        ca_weight: float = 1.0,
        cb_weight: float = 1.2,   # CB slightly up-weighted — more discriminating
    ) -> pd.DataFrame:
        """
        Joint CA-CB assignment using 2D cost matrix.

        Each experimental (CA, CB) pair is matched to a residue position.
        The cost is the sum of CA and CB costs — both must be consistent
        with the same residue and SS class.

        ca_peaks, cb_peaks : lists of ppm values, same length
                             (i-th CA peak and i-th CB peak come from the same residue)

        Returns DataFrame with same columns as assign(), plus 'cb_exp_shift' and 'cb_delta_ppm'.
        """
        if len(ca_peaks) != len(cb_peaks):
            raise ValueError(
                f"ca_peaks ({len(ca_peaks)}) and cb_peaks ({len(cb_peaks)}) must have equal length.\n"
                "Each pair should be the CA and CB shifts from the same residue."
            )

        # Slots: only residues that have both CA and CB
        ca_slots = {s.seq_id: s for s in self._slots_by_atom.get('CA', [])}
        cb_slots = {s.seq_id: s for s in self._slots_by_atom.get('CB', [])}
        joint_slots = [ca_slots[sid] for sid in sorted(ca_slots) if sid in cb_slots]

        n_peaks = len(ca_peaks)
        n_slots = len(joint_slots)

        # Build joint cost matrix
        cost = np.full((n_peaks, n_slots), INF_COST, dtype=float)

        for i, (ca, cb) in enumerate(zip(ca_peaks, cb_peaks)):
            for j, slot in enumerate(joint_slots):
                ca_cost = self._peak_slot_cost(ca, slot)
                # Temporary CB slot using same seq_id
                cb_slot = AssignmentSlot(slot.seq_id, slot.residue, 'CB', slot.ss_class)
                cb_cost = self._peak_slot_cost(cb, cb_slot)

                if ca_cost >= INF_COST or cb_cost >= INF_COST:
                    continue

                cost[i, j] = ca_weight * ca_cost + cb_weight * cb_cost

        # Pad and solve
        if n_peaks > n_slots:
            pad = np.full((n_peaks, n_peaks - n_slots), INF_COST * 0.9)
            cost_padded = np.hstack([cost, pad])
            dummy_slots = n_peaks - n_slots
        else:
            cost_padded = cost
            dummy_slots = 0

        row_ind, col_ind = linear_sum_assignment(cost_padded)

        rows = []
        for r, c in zip(row_ind, col_ind):
            ca  = ca_peaks[r]
            cb  = cb_peaks[r]
            is_dummy = (dummy_slots > 0 and c >= n_slots)

            if is_dummy or cost_padded[r, c] >= INF_COST * 0.8:
                rows.append({
                    'peak_id': r + 1,
                    'ca_exp': round(ca, 3), 'cb_exp': round(cb, 3),
                    'seq_id': None, 'residue': None,
                    'ss_class': None,
                    'ca_pred': None, 'cb_pred': None,
                    'ca_delta': None, 'cb_delta': None,
                    'joint_cost': None, 'confidence': 0.0, 'status': 'unassigned',
                })
                continue

            slot = joint_slots[c]
            ca_entry = self._lookup.get((slot.residue, 'CA', slot.ss_class)) or \
                       self._lookup.get((slot.residue, 'CA', 'coil'))
            cb_entry = self._lookup.get((slot.residue, 'CB', slot.ss_class)) or \
                       self._lookup.get((slot.residue, 'CB', 'coil'))

            ca_pred = ca_entry['mean'] if ca_entry else RANDOM_COIL.get(slot.residue, {}).get('CA')
            cb_pred = cb_entry['mean'] if cb_entry else RANDOM_COIL.get(slot.residue, {}).get('CB')

            rows.append({
                'peak_id':   r + 1,
                'ca_exp':    round(ca, 3),
                'cb_exp':    round(cb, 3),
                'seq_id':    slot.seq_id,
                'residue':   slot.residue,
                'ss_class':  slot.ss_class,
                'ca_pred':   round(ca_pred, 3) if ca_pred else None,
                'cb_pred':   round(cb_pred, 3) if cb_pred else None,
                'ca_delta':  round(ca - ca_pred, 3) if ca_pred else None,
                'cb_delta':  round(cb - cb_pred, 3) if cb_pred else None,
                'joint_cost': round(cost[r, c], 4),
                'confidence': None,
                'status':    'assigned',
            })

        result = pd.DataFrame(rows)

        # Confidence
        assigned = result[result['status'] == 'assigned']
        if not assigned.empty:
            costs = assigned['joint_cost'].values.astype(float)
            exp_neg = np.exp(-costs / costs.mean())
            confidences = np.clip(exp_neg / exp_neg.sum() * len(exp_neg) / len(result), 0, 1)
            result.loc[result['status'] == 'assigned', 'confidence'] = np.round(confidences, 3)
        result.loc[result['status'] == 'unassigned', 'confidence'] = 0.0

        return result.sort_values('seq_id', na_position='last').reset_index(drop=True)

    # ------------------------------------------------------------------ #
    # Ranked candidates (for any single peak — diagnostic use)
    # ------------------------------------------------------------------ #

    def rank_candidates(
        self,
        peak_ppm: float,
        atom: str,
        top_n: int = 10,
    ) -> pd.DataFrame:
        """
        For a single experimental peak, rank all possible assignments by cost.
        Useful for inspecting ambiguous peaks or validating the engine.

        Returns DataFrame: rank | seq_id | residue | ss_class | predicted_shift | delta_ppm | cost | confidence
        """
        slots = self._slots_by_atom.get(atom, [])
        rows = []
        for slot in slots:
            cost = self._peak_slot_cost(peak_ppm, slot)
            entry = self._lookup.get((slot.residue, slot.atom, slot.ss_class)) or \
                    self._lookup.get((slot.residue, slot.atom, 'coil'))
            pred = entry['mean'] if entry else RANDOM_COIL.get(slot.residue, {}).get(slot.atom)
            rows.append({
                'seq_id': slot.seq_id, 'residue': slot.residue,
                'ss_class': slot.ss_class,
                'predicted_shift': round(pred, 3) if pred else None,
                'delta_ppm': round(peak_ppm - pred, 3) if pred else None,
                'cost': round(cost, 4),
            })

        df = pd.DataFrame(rows).sort_values('cost').head(top_n).reset_index(drop=True)
        df['rank'] = df.index + 1

        # Confidence relative to top_n candidates
        costs = df['cost'].values.astype(float)
        exp_neg = np.exp(-costs + costs.min())
        df['confidence'] = np.round(exp_neg / exp_neg.sum(), 3)

        return df

    # ------------------------------------------------------------------ #
    # Assignment quality metrics
    # ------------------------------------------------------------------ #

    def score_assignment(self, assignment_df: pd.DataFrame) -> dict:
        """
        Compute quality metrics for a completed assignment.

        Returns dict with:
        - n_assigned: number of peaks assigned
        - n_unassigned: peaks with no match
        - mean_delta_ppm: average |predicted - experimental| for assigned peaks
        - rmsd_ppm: RMS deviation
        - mean_confidence: average confidence score
        - high_confidence_frac: fraction of peaks with confidence >= 0.5
        - composition_violations: residue types assigned more than sequence allows
        """
        assigned = assignment_df[assignment_df['status'] == 'assigned']
        unassigned = assignment_df[assignment_df['status'] == 'unassigned']

        deltas = assigned['delta_ppm'].dropna().abs()
        confidences = assigned['confidence'].dropna()

        # Check composition constraints
        if 'residue' in assigned.columns:
            assigned_counts = Counter(assigned['residue'].dropna())
            violations = {
                res: (count, self._composition[res])
                for res, count in assigned_counts.items()
                if count > self._composition.get(res, 0)
            }
        else:
            violations = {}

        return {
            'n_assigned':            len(assigned),
            'n_unassigned':          len(unassigned),
            'assignment_rate':       round(len(assigned) / max(len(assignment_df), 1), 3),
            'mean_delta_ppm':        round(deltas.mean(), 3) if len(deltas) else None,
            'rmsd_ppm':              round(np.sqrt((deltas**2).mean()), 3) if len(deltas) else None,
            'mean_confidence':       round(confidences.mean(), 3) if len(confidences) else None,
            'high_confidence_frac':  round((confidences >= 0.5).mean(), 3) if len(confidences) else None,
            'composition_violations': violations,
        }


# ---------------------------------------------------------------------------
# Convenience runner — called from pipeline.py
# ---------------------------------------------------------------------------

def run_assignment_engine(
    experimental_peaks: Dict[str, List[float]],
    stats_df: pd.DataFrame,
    sequence: str,
    ss_map: Dict[int, str],
    results_dir: Path,
    ca_cb_pairs: List[Tuple[float, float]] = None,
) -> dict:
    """
    Full Phase 5 assignment run.

    Parameters
    ----------
    experimental_peaks : {'CA': [ppm, ...], 'CB': [ppm, ...], 'N': [ppm, ...]}
        Experimental peak positions per atom type.
        Pass whatever you have — CA only is fine, CA+CB is better.
    stats_df : reference database (from pipeline)
    sequence : one-letter sequence string
    ss_map : {seq_id -> ss_class} (from DSSP)
    results_dir : where to save output CSVs
    ca_cb_pairs : list of (ca_ppm, cb_ppm) tuples for joint 2D assignment (optional)
                  If None, 1D assignments are run for each atom type separately.

    Returns
    -------
    dict with keys: '1d_results', '2d_result', 'metrics', 'engine'
    """
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    print("\n" + "=" * 60)
    print("  PHASE 5 — AUTOMATED ASSIGNMENT ENGINE")
    print("=" * 60)
    print(f"  Sequence length:    {len(sequence)} residues")
    print(f"  Reference DB rows:  {len(stats_df)}")
    for atom, peaks in experimental_peaks.items():
        print(f"  Experimental peaks ({atom}): {len(peaks)}")

    engine = AssignmentEngine(stats_df, sequence, ss_map)
    output = {'engine': engine, '1d_results': {}, 'metrics': {}}

    # ── 1D assignments ────────────────────────────────────────────────────
    for atom, peaks in experimental_peaks.items():
        if not peaks:
            continue
        print(f"\n--- 1D Assignment: {atom} ({len(peaks)} peaks) ---")
        try:
            result = engine.assign(peaks, atom=atom)
            assigned = result[result['status'] == 'assigned']
            metrics  = engine.score_assignment(result)

            print(f"  Assigned:       {metrics['n_assigned']} / {len(peaks)}")
            print(f"  RMSD:           {metrics['rmsd_ppm']} ppm")
            print(f"  Mean |Δδ|:      {metrics['mean_delta_ppm']} ppm")
            print(f"  Mean confidence:{metrics['mean_confidence']}")
            if metrics['composition_violations']:
                print(f"  ⚠ Composition violations: {metrics['composition_violations']}")

            # Show top assignments
            print(f"\n  Top assignments (sorted by seq_id):")
            display_cols = ['peak_id', 'exp_shift', 'seq_id', 'residue',
                            'ss_class', 'predicted_shift', 'delta_ppm', 'confidence', 'status']
            print(result[display_cols].head(20).to_string(index=False))

            out_path = results_dir / f"assignment_{atom}.csv"
            result.to_csv(out_path, index=False)
            print(f"\n  Saved → {out_path}")

            output['1d_results'][atom] = result
            output['metrics'][atom] = metrics

        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback; traceback.print_exc()

    # ── 2D joint CA-CB assignment ─────────────────────────────────────────
    if ca_cb_pairs:
        print(f"\n--- 2D Joint CA-CB Assignment ({len(ca_cb_pairs)} pairs) ---")
        try:
            ca_list = [p[0] for p in ca_cb_pairs]
            cb_list = [p[1] for p in ca_cb_pairs]
            result_2d = engine.assign_joint(ca_list, cb_list)

            assigned_2d = result_2d[result_2d['status'] == 'assigned']
            n_total = len(result_2d)
            print(f"  Assigned:   {len(assigned_2d)} / {n_total} peak pairs")
            print(f"  Unassigned: {n_total - len(assigned_2d)}")

            if not assigned_2d.empty:
                ca_deltas = assigned_2d['ca_delta'].dropna().abs()
                cb_deltas = assigned_2d['cb_delta'].dropna().abs()
                print(f"  CA RMSD:    {np.sqrt((ca_deltas**2).mean()):.3f} ppm")
                print(f"  CB RMSD:    {np.sqrt((cb_deltas**2).mean()):.3f} ppm")
                print(f"  Confidence: {assigned_2d['confidence'].mean():.3f} (mean)")

            display_2d = ['peak_id', 'ca_exp', 'cb_exp', 'seq_id', 'residue',
                          'ss_class', 'ca_pred', 'ca_delta', 'cb_delta', 'confidence', 'status']
            print(f"\n  Top 2D assignments:")
            print(result_2d[display_2d].head(20).to_string(index=False))

            out_path = results_dir / "assignment_2d_cacb.csv"
            result_2d.to_csv(out_path, index=False)
            print(f"\n  Saved → {out_path}")

            output['2d_result'] = result_2d

        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback; traceback.print_exc()

    # ── Summary ───────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("  ASSIGNMENT SUMMARY")
    print("=" * 60)
    for atom, m in output['metrics'].items():
        print(f"  {atom}: {m['n_assigned']} assigned, RMSD={m['rmsd_ppm']} ppm, "
              f"conf={m['mean_confidence']}")
    if 'CA' in output['metrics']:
        m = output['metrics']['CA']
        if m['rmsd_ppm'] is not None and m['rmsd_ppm'] < 1.5:
            print(f"\n  ✓ CA RMSD < 1.5 ppm — assignment quality is GOOD")
        elif m['rmsd_ppm'] is not None and m['rmsd_ppm'] < 3.0:
            print(f"\n  ⚠ CA RMSD 1.5–3.0 ppm — assignment quality is MODERATE")
            print(f"    Consider: more BMRB entries in reference DB (run --batch)")
        elif m['rmsd_ppm'] is not None:
            print(f"\n  ✗ CA RMSD > 3.0 ppm — assignment quality is POOR")
            print(f"    Likely cause: reference DB too small (only 1 entry?)")
            print(f"    Fix: python src/pipeline.py --batch  (builds 10-entry reference DB)")

    return output
