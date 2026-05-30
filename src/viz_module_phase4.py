"""
viz_module_phase4.py — Phase 4: Enhanced spectrum simulation.

New capabilities over Phase 3:
  1. Voigt lineshape (Lorentzian × Gaussian convolution) — more realistic ssNMR peaks
  2. Relative intensity scaling by carbon multiplicity (CH, CH2, CH3)
  3. Side-chain carbons (CB, CG, CD, CE) alongside backbone CA/N
  4. 2D CA-CB correlation spectrum simulation
  5. NMRPipe-compatible text export (.tab format)
  6. Overlay: experimental vs predicted spectrum comparison

All functions are standalone additions — drop into viz_module.py or import directly.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
from pathlib import Path
from typing import Optional, List, Dict, Tuple
from scipy.signal import fftconvolve
from scipy.special import voigt_profile


# ---------------------------------------------------------------------------
# Carbon multiplicity table
# Number of carbons of each type per residue (for relative intensity scaling).
# Values reflect 13C natural abundance × number of equivalent sites.
# In DEPT/CP-MAS experiments, CH > CH2 > CH3 ≈ quaternary C (near zero).
# For direct 13C, all are roughly equal; this table handles DEPT editing.
# ---------------------------------------------------------------------------

CARBON_MULTIPLICITY = {
    # residue -> {atom -> (n_carbons, carbon_type)}
    # carbon_type: 'CH', 'CH2', 'CH3', 'quat' (quaternary), 'none'
    'A': {'CA': (1, 'CH'),  'CB': (1, 'CH3'), 'C': (1, 'quat')},
    'C': {'CA': (1, 'CH'),  'CB': (1, 'CH2'), 'C': (1, 'quat')},
    'D': {'CA': (1, 'CH'),  'CB': (1, 'CH2'), 'CG': (1, 'quat'), 'C': (1, 'quat')},
    'E': {'CA': (1, 'CH'),  'CB': (1, 'CH2'), 'CG': (1, 'CH2'), 'CD': (1, 'quat'), 'C': (1, 'quat')},
    'F': {'CA': (1, 'CH'),  'CB': (1, 'CH2'), 'CG': (1, 'quat'), 'CD': (2, 'CH'), 'CE': (2, 'CH'), 'CZ': (1, 'CH'), 'C': (1, 'quat')},
    'G': {'CA': (1, 'CH2'), 'C': (1, 'quat')},
    'H': {'CA': (1, 'CH'),  'CB': (1, 'CH2'), 'CG': (1, 'quat'), 'CD': (1, 'CH'), 'CE': (1, 'CH'), 'C': (1, 'quat')},
    'I': {'CA': (1, 'CH'),  'CB': (1, 'CH'),  'CG1': (1, 'CH2'), 'CG2': (1, 'CH3'), 'CD1': (1, 'CH3'), 'C': (1, 'quat')},
    'K': {'CA': (1, 'CH'),  'CB': (1, 'CH2'), 'CG': (1, 'CH2'), 'CD': (1, 'CH2'), 'CE': (1, 'CH2'), 'C': (1, 'quat')},
    'L': {'CA': (1, 'CH'),  'CB': (1, 'CH2'), 'CG': (1, 'CH'),  'CD1': (1, 'CH3'), 'CD2': (1, 'CH3'), 'C': (1, 'quat')},
    'M': {'CA': (1, 'CH'),  'CB': (1, 'CH2'), 'CG': (1, 'CH2'), 'CE': (1, 'CH3'), 'C': (1, 'quat')},
    'N': {'CA': (1, 'CH'),  'CB': (1, 'CH2'), 'CG': (1, 'quat'), 'C': (1, 'quat')},
    'P': {'CA': (1, 'CH'),  'CB': (1, 'CH2'), 'CG': (1, 'CH2'), 'CD': (1, 'CH2'), 'C': (1, 'quat')},
    'Q': {'CA': (1, 'CH'),  'CB': (1, 'CH2'), 'CG': (1, 'CH2'), 'CD': (1, 'quat'), 'C': (1, 'quat')},
    'R': {'CA': (1, 'CH'),  'CB': (1, 'CH2'), 'CG': (1, 'CH2'), 'CD': (1, 'CH2'), 'CZ': (1, 'quat'), 'C': (1, 'quat')},
    'S': {'CA': (1, 'CH'),  'CB': (1, 'CH2'), 'C': (1, 'quat')},
    'T': {'CA': (1, 'CH'),  'CB': (1, 'CH'),  'CG2': (1, 'CH3'), 'C': (1, 'quat')},
    'V': {'CA': (1, 'CH'),  'CB': (1, 'CH'),  'CG1': (1, 'CH3'), 'CG2': (1, 'CH3'), 'C': (1, 'quat')},
    'W': {'CA': (1, 'CH'),  'CB': (1, 'CH2'), 'CG': (1, 'quat'), 'CD1': (1, 'CH'), 'CD2': (1, 'quat'), 'CE2': (1, 'quat'), 'CE3': (1, 'CH'), 'CZ2': (1, 'CH'), 'CZ3': (1, 'CH'), 'CH2': (1, 'CH'), 'C': (1, 'quat')},
    'Y': {'CA': (1, 'CH'),  'CB': (1, 'CH2'), 'CG': (1, 'quat'), 'CD': (2, 'CH'), 'CE': (2, 'CH'), 'CZ': (1, 'quat'), 'C': (1, 'quat')},
}

# DEPT intensity weights by carbon type
# DEPT-135: CH and CH3 positive, CH2 negative, quaternary C absent
# Direct 13C CP-MAS: all present, intensity ≈ proportional to n_carbons
DEPT_WEIGHTS = {
    'CH':   1.0,
    'CH2': -1.0,   # negative in DEPT-135
    'CH3':  1.0,
    'quat': 0.0,   # absent in DEPT
    'none': 0.0,
}

DIRECT_WEIGHTS = {
    'CH':   1.0,
    'CH2':  2.0,   # twice the signal (2 carbons)
    'CH3':  3.0,   # three times the signal
    'quat': 1.0,
    'none': 0.0,
}

# Random coil reference shifts for side-chain atoms (Wishart 1995 + BMRB averages)
SIDECHAIN_RC = {
    'A': {'CB': 19.1},
    'C': {'CB': 28.0},
    'D': {'CB': 41.1, 'CG': 178.4},
    'E': {'CB': 30.2, 'CG': 35.6, 'CD': 182.0},
    'F': {'CB': 39.6, 'CG': 137.3, 'CD': 131.8, 'CE': 130.9, 'CZ': 129.4},
    'G': {},
    'H': {'CB': 29.4, 'CG': 131.9, 'CD': 120.0, 'CE': 138.0},
    'I': {'CB': 38.8, 'CG1': 27.2, 'CG2': 17.4, 'CD1': 13.7},
    'K': {'CB': 33.1, 'CG': 24.7, 'CD': 29.2, 'CE': 42.0},
    'L': {'CB': 42.4, 'CG': 26.9, 'CD1': 24.8, 'CD2': 23.2},
    'M': {'CB': 33.1, 'CG': 31.5, 'CE': 17.2},
    'N': {'CB': 38.9, 'CG': 176.4},
    'P': {'CB': 32.1, 'CG': 27.3, 'CD': 50.5},
    'Q': {'CB': 29.4, 'CG': 33.7, 'CD': 178.5},
    'R': {'CB': 30.9, 'CG': 27.3, 'CD': 43.5, 'CZ': 159.5},
    'S': {'CB': 63.8},
    'T': {'CB': 69.8, 'CG2': 21.5},
    'V': {'CB': 32.9, 'CG1': 21.6, 'CG2': 21.6},
    'W': {'CB': 29.9, 'CG': 112.0, 'CD1': 126.1, 'CD2': 129.2, 'CE2': 138.0, 'CE3': 121.3, 'CZ2': 114.5, 'CZ3': 121.3, 'CH2': 123.9},
    'Y': {'CB': 38.8, 'CG': 128.4, 'CD': 132.9, 'CE': 117.5, 'CZ': 155.9},
}


# ---------------------------------------------------------------------------
# 1. Voigt lineshape
# ---------------------------------------------------------------------------

def voigt_peak(x: np.ndarray, center: float, sigma: float, gamma: float) -> np.ndarray:
    """
    Compute a Voigt profile at positions x.

    sigma: Gaussian width (inhomogeneous broadening, typical ssNMR: 0.1–0.5 ppm)
    gamma: Lorentzian half-width (natural + homogeneous broadening, typical: 0.1–0.3 ppm)

    The Voigt profile is more realistic for solid-state NMR than pure Lorentzian:
    - Lorentzian component = T2 relaxation (homogeneous)
    - Gaussian component = chemical shift distribution (inhomogeneous / disorder)
    """
    z = (x - center + 1j * gamma) / (sigma * np.sqrt(2))
    return voigt_profile(x - center, sigma, gamma)


def lorentzian_peak(x: np.ndarray, center: float, gamma: float) -> np.ndarray:
    """Pure Lorentzian (for comparison / simpler spectra)."""
    return (gamma / 2) ** 2 / ((x - center) ** 2 + (gamma / 2) ** 2)


def gaussian_peak(x: np.ndarray, center: float, sigma: float) -> np.ndarray:
    """Pure Gaussian."""
    return np.exp(-0.5 * ((x - center) / sigma) ** 2)


# ---------------------------------------------------------------------------
# 2. Core 1D spectrum simulator (Phase 4 — enhanced)
# ---------------------------------------------------------------------------

def simulate_1d_spectrum(
    peaks: List[Dict],
    ppm_range: Tuple[float, float] = (0, 200),
    n_points: int = 16000,
    sigma: float = 0.15,       # Gaussian width (ppm) — inhomogeneous broadening
    gamma: float = 0.15,       # Lorentzian half-width (ppm) — homogeneous broadening
    lineshape: str = 'voigt',  # 'voigt', 'lorentzian', or 'gaussian'
    experiment: str = 'direct',  # 'direct' (CP-MAS) or 'dept135'
    normalize: bool = True,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Simulate a 1D 13C NMR spectrum.

    Parameters
    ----------
    peaks : list of dicts, each with:
        'shift' (float, ppm)
        'residue' (str, one-letter)
        'atom' (str, e.g. 'CA')
        'intensity' (float, optional — overrides auto-scaling)
    ppm_range : (lo, hi) in ppm
    n_points : number of frequency points
    sigma : Gaussian broadening (ppm) — models disorder/inhomogeneity
    gamma : Lorentzian broadening (ppm) — models T2 relaxation
    lineshape : 'voigt' (default), 'lorentzian', 'gaussian'
    experiment : 'direct' (all carbons) or 'dept135' (CH+CH3 positive, CH2 negative)
    normalize : scale max intensity to 1.0

    Returns
    -------
    (x, spectrum) arrays
    """
    x = np.linspace(ppm_range[0], ppm_range[1], n_points)
    spectrum = np.zeros_like(x)

    weight_table = DEPT_WEIGHTS if experiment == 'dept135' else DIRECT_WEIGHTS

    for peak in peaks:
        shift = peak.get('shift')
        if shift is None or np.isnan(shift):
            continue
        if not (ppm_range[0] <= shift <= ppm_range[1]):
            continue

        # Determine intensity
        if 'intensity' in peak:
            intensity = float(peak['intensity'])
        else:
            residue = peak.get('residue', 'X')
            atom = peak.get('atom', 'CA')
            mult_info = CARBON_MULTIPLICITY.get(residue, {}).get(atom)
            if mult_info:
                n_carbons, ctype = mult_info
                intensity = n_carbons * weight_table.get(ctype, 1.0)
            else:
                intensity = 1.0

        if intensity == 0:
            continue

        # Lineshape
        if lineshape == 'voigt':
            peak_shape = voigt_peak(x, shift, sigma, gamma)
        elif lineshape == 'lorentzian':
            peak_shape = lorentzian_peak(x, shift, gamma)
        else:
            peak_shape = gaussian_peak(x, shift, sigma)

        spectrum += intensity * peak_shape

    if normalize and spectrum.max() > 0:
        spectrum = spectrum / spectrum.max()

    return x, spectrum


# ---------------------------------------------------------------------------
# 3. Build peak list from predictions DataFrame
# ---------------------------------------------------------------------------

def predictions_to_peaks(
    predictions_df: pd.DataFrame,
    atoms: List[str] = None,
    include_sidechains: bool = True,
    sequence: str = None,
) -> List[Dict]:
    """
    Convert a predictions DataFrame into a list of peak dicts for simulate_1d_spectrum().

    If include_sidechains=True and sequence is provided, adds side-chain carbons
    using SIDECHAIN_RC random-coil references adjusted by the predicted CA deviation
    as a first-order structural correction.
    """
    if atoms is None:
        atoms = ['CA', 'CB', 'C']

    peaks = []

    # Backbone peaks from predictions
    for _, row in predictions_df.iterrows():
        if row['atom'] not in atoms:
            continue
        shift = row.get('predicted_shift')
        if shift is None or (isinstance(shift, float) and np.isnan(shift)):
            continue
        peaks.append({
            'shift': float(shift),
            'residue': row['residue'],
            'atom': row['atom'],
            'seq_id': row['seq_id'],
            'ss_class': row.get('ss_class', 'coil'),
        })

    # Side-chain carbons (approximated from RC + structural correction)
    if include_sidechains and sequence:
        # Build a CA deviation map from predictions (structural correction proxy)
        ca_dev = {}
        for _, row in predictions_df[predictions_df['atom'] == 'CA'].iterrows():
            dev = row.get('rc_deviation')
            if dev is not None and not np.isnan(float(dev)):
                ca_dev[int(row['seq_id'])] = float(dev)

        for i, aa in enumerate(sequence, 1):
            sc_shifts = SIDECHAIN_RC.get(aa, {})
            dev = ca_dev.get(i, 0.0)
            # Side-chain structural correction: approximately 0.3-0.5× backbone CA deviation
            sc_correction = dev * 0.35

            for sc_atom, rc_shift in sc_shifts.items():
                if sc_atom == 'CB':
                    continue  # CB already in backbone predictions
                adjusted_shift = rc_shift + sc_correction
                peaks.append({
                    'shift': adjusted_shift,
                    'residue': aa,
                    'atom': sc_atom,
                    'seq_id': i,
                    'ss_class': 'sidechain',
                })

    return peaks


# ---------------------------------------------------------------------------
# 4. 2D CA-CB correlation spectrum
# ---------------------------------------------------------------------------

def simulate_2d_cacb(
    predictions_df: pd.DataFrame,
    ppm_range_ca: Tuple[float, float] = (40, 75),
    ppm_range_cb: Tuple[float, float] = (10, 75),
    n_points: int = 512,
    sigma_ca: float = 0.2,
    sigma_cb: float = 0.25,
    gamma_ca: float = 0.15,
    gamma_cb: float = 0.20,
    output_path: Optional[str] = None,
    title: str = "Simulated CA-CB 2D Correlation",
    color_by_ss: bool = True,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Simulate a 2D CA-CB correlation spectrum (analogous to 2D DARR/PDSD in ssNMR).

    Each residue contributes one cross-peak at (CA shift, CB shift).
    Glycine (no CB) is omitted.

    Returns
    -------
    (ca_axis, cb_axis, spectrum_2d)
    Also saves a plot if output_path is given.
    """
    ca_axis = np.linspace(ppm_range_ca[0], ppm_range_ca[1], n_points)
    cb_axis = np.linspace(ppm_range_cb[0], ppm_range_cb[1], n_points)
    spectrum = np.zeros((n_points, n_points))

    # Build (CA, CB) pairs per residue
    ca_shifts = predictions_df[predictions_df['atom'] == 'CA'].set_index('seq_id')['predicted_shift']
    cb_shifts = predictions_df[predictions_df['atom'] == 'CB'].set_index('seq_id')['predicted_shift']
    ss_map = predictions_df[predictions_df['atom'] == 'CA'].set_index('seq_id')['ss_class']
    res_map = predictions_df[predictions_df['atom'] == 'CA'].set_index('seq_id')['residue']

    cross_peaks = []
    for seq_id in ca_shifts.index:
        if seq_id not in cb_shifts.index:
            continue  # Gly or missing CB
        ca = ca_shifts[seq_id]
        cb = cb_shifts[seq_id]
        if ca is None or cb is None:
            continue
        if np.isnan(float(ca)) or np.isnan(float(cb)):
            continue

        cross_peaks.append({
            'seq_id': seq_id,
            'ca': float(ca),
            'cb': float(cb),
            'ss_class': ss_map.get(seq_id, 'coil'),
            'residue': res_map.get(seq_id, 'X'),
        })

        # Add peak to 2D matrix (Gaussian in both dimensions)
        ca_idx = np.argmin(np.abs(ca_axis - ca))
        cb_idx = np.argmin(np.abs(cb_axis - cb))

        # 2D Voigt as outer product of 1D profiles
        ca_profile = gaussian_peak(ca_axis, ca, sigma_ca)
        cb_profile = gaussian_peak(cb_axis, cb, sigma_cb)
        spectrum += np.outer(cb_profile, ca_profile)

    # Normalize
    if spectrum.max() > 0:
        spectrum = spectrum / spectrum.max()

    # Plot
    if output_path or True:  # always build fig for return
        SS_COLORS = {
            'helix': '#E74C3C', 'strand': '#3498DB',
            'coil': '#95A5A6', 'unknown': '#BDC3C7', 'sidechain': '#2ECC71',
        }

        fig, ax = plt.subplots(figsize=(9, 8))

        # Contour plot
        levels = np.linspace(0.05, 1.0, 12)
        ax.contourf(ca_axis, cb_axis, spectrum, levels=levels, cmap='Blues', alpha=0.6)
        ax.contour(ca_axis, cb_axis, spectrum, levels=levels, colors='navy', linewidths=0.4, alpha=0.5)

        # Scatter cross-peaks colored by SS class
        for cp in cross_peaks:
            color = SS_COLORS.get(cp['ss_class'], 'black')
            ax.scatter(cp['ca'], cp['cb'], c=color, s=25, zorder=5, alpha=0.9)
            ax.annotate(
                f"{cp['residue']}{cp['seq_id']}",
                (cp['ca'], cp['cb']),
                fontsize=5, ha='left', va='bottom', color='#2C3E50', alpha=0.8,
            )

        ax.set_xlabel('CA Chemical Shift (ppm)', fontsize=12)
        ax.set_ylabel('CB Chemical Shift (ppm)', fontsize=12)
        ax.set_title(title, fontsize=13, fontweight='bold')
        ax.invert_xaxis()
        ax.invert_yaxis()

        # Legend
        patches = [
            mpatches.Patch(color='#E74C3C', label='Helix'),
            mpatches.Patch(color='#3498DB', label='Strand'),
            mpatches.Patch(color='#95A5A6', label='Coil'),
        ]
        ax.legend(handles=patches, fontsize=9, loc='upper right')
        ax.grid(True, alpha=0.2)
        plt.tight_layout()

        if output_path:
            fig.savefig(output_path, dpi=150, bbox_inches='tight')
            print(f"[VIZ Phase 4] Saved 2D CA-CB spectrum → {output_path}")
        plt.close(fig)

    return ca_axis, cb_axis, spectrum, cross_peaks


# ---------------------------------------------------------------------------
# 5. Publication-quality 1D spectrum plot
# ---------------------------------------------------------------------------

def plot_1d_spectrum_phase4(
    predictions_df: pd.DataFrame,
    sequence: str = None,
    ppm_range: Tuple[float, float] = (0, 200),
    sigma: float = 0.15,
    gamma: float = 0.15,
    include_sidechains: bool = True,
    experiment: str = 'direct',
    output_path: Optional[str] = None,
    title: str = "Simulated 13C Solid-State NMR Spectrum",
    show_regions: bool = True,
    bmrb_id: str = '',
) -> plt.Figure:
    """
    Full publication-quality 1D 13C spectrum with:
    - Voigt lineshape
    - Intensity-scaled peaks
    - Side-chain carbons
    - Annotated chemical shift regions
    - Color-coded peak markers by SS class
    """
    peaks = predictions_to_peaks(
        predictions_df,
        atoms=['CA', 'CB', 'C'],
        include_sidechains=include_sidechains,
        sequence=sequence,
    )

    x, spectrum = simulate_1d_spectrum(
        peaks,
        ppm_range=ppm_range,
        sigma=sigma,
        gamma=gamma,
        lineshape='voigt',
        experiment=experiment,
    )

    fig, (ax_main, ax_ca) = plt.subplots(
        2, 1, figsize=(16, 8),
        gridspec_kw={'height_ratios': [3, 1]},
        sharex=True,
    )

    # ── Main spectrum ─────────────────────────────────────────────────────
    ax_main.plot(x, spectrum, color='#2C3E50', linewidth=1.0, label='Simulated spectrum')
    ax_main.fill_between(x, spectrum, alpha=0.1, color='#3498DB')

    # Annotated chemical shift regions (ssNMR reference ranges)
    if show_regions:
        regions = [
            (170, 185, '#FFF3CD', 'Carbonyl C'),
            (110, 145, '#E8F5E9', 'Aromatic C'),
            (60, 80,  '#E3F2FD', 'Cα + Ser/Thr Cβ'),
            (40, 60,  '#FCE4EC', 'Cα aliphatic'),
            (10, 40,  '#F3E5F5', 'Cβ/Cγ aliphatic'),
        ]
        for lo, hi, color, label in regions:
            if hi > ppm_range[0] and lo < ppm_range[1]:
                ax_main.axvspan(lo, hi, alpha=0.15, color=color, label=label)

    ax_main.set_ylabel('Intensity (a.u.)', fontsize=11)
    ax_main.set_title(title, fontsize=13, fontweight='bold')
    ax_main.set_xlim(ppm_range[1], ppm_range[0])  # NMR convention: right→left
    ax_main.grid(True, alpha=0.15)
    ax_main.legend(fontsize=8, loc='upper left')

    # ── CA-only zoom panel ────────────────────────────────────────────────
    SS_COLORS = {
        'helix': '#E74C3C', 'strand': '#3498DB',
        'coil': '#95A5A6', 'unknown': '#BDC3C7', 'sidechain': '#2ECC71',
    }

    ca_peaks = [p for p in peaks if p.get('atom') == 'CA']
    for p in ca_peaks:
        color = SS_COLORS.get(p.get('ss_class', 'coil'), 'grey')
        ax_ca.axvline(p['shift'], color=color, linewidth=1.5, alpha=0.7)

    ax_ca.set_xlim(ppm_range[1], ppm_range[0])
    ax_ca.set_xlabel('13C Chemical Shift (ppm)', fontsize=11)
    ax_ca.set_ylabel('Cα positions', fontsize=9)
    ax_ca.set_yticks([])
    ax_ca.grid(True, alpha=0.15)

    # Mini legend for CA panel
    for ss, color in [('helix', '#E74C3C'), ('strand', '#3498DB'), ('coil', '#95A5A6')]:
        ax_ca.axvline(0, color=color, linewidth=2, label=ss, alpha=0)  # invisible, for legend
    ca_legend_patches = [mpatches.Patch(color=c, label=s) for s, c in [('helix', '#E74C3C'), ('strand', '#3498DB'), ('coil', '#95A5A6')]]
    ax_ca.legend(handles=ca_legend_patches, fontsize=8, loc='upper left', ncol=3)

    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"[VIZ Phase 4] Saved 1D spectrum → {output_path}")

    return fig


# ---------------------------------------------------------------------------
# 6. Experimental vs predicted overlay
# ---------------------------------------------------------------------------

def plot_experimental_vs_predicted(
    experimental_shifts: List[float],
    predictions_df: pd.DataFrame,
    atom: str = 'CA',
    ppm_range: Tuple[float, float] = (40, 75),
    sigma: float = 0.2,
    gamma: float = 0.15,
    output_path: Optional[str] = None,
    title: str = None,
) -> plt.Figure:
    """
    Overlay experimental peak positions on the predicted spectrum.
    Shows how well the predicted shifts match experimental data.

    experimental_shifts: list of ppm values from the actual NMR spectrum
    predictions_df: output of ShiftPredictor.predict()
    """
    # Build predicted spectrum for this atom
    atom_preds = predictions_df[predictions_df['atom'] == atom]
    peaks = []
    for _, row in atom_preds.iterrows():
        shift = row.get('predicted_shift')
        if shift is not None and not np.isnan(float(shift)):
            peaks.append({
                'shift': float(shift),
                'residue': row['residue'],
                'atom': atom,
                'ss_class': row.get('ss_class', 'coil'),
            })

    x, pred_spectrum = simulate_1d_spectrum(
        peaks, ppm_range=ppm_range, sigma=sigma, gamma=gamma, lineshape='voigt'
    )

    # Build experimental spectrum
    exp_peaks_list = [{'shift': s, 'atom': atom} for s in experimental_shifts
                      if ppm_range[0] <= s <= ppm_range[1]]
    _, exp_spectrum = simulate_1d_spectrum(
        exp_peaks_list, ppm_range=ppm_range, sigma=sigma * 1.5, gamma=gamma, lineshape='voigt'
    )

    fig, axes = plt.subplots(3, 1, figsize=(14, 9), sharex=True,
                             gridspec_kw={'height_ratios': [2, 2, 1]})

    # Predicted
    axes[0].plot(x, pred_spectrum, color='#3498DB', linewidth=1.2)
    axes[0].fill_between(x, pred_spectrum, alpha=0.2, color='#3498DB')
    axes[0].set_ylabel('Predicted', fontsize=10)
    axes[0].set_title(title or f'{atom} Shift: Predicted vs Experimental', fontsize=12, fontweight='bold')

    # Experimental
    axes[1].plot(x, exp_spectrum, color='#E74C3C', linewidth=1.2)
    axes[1].fill_between(x, exp_spectrum, alpha=0.2, color='#E74C3C')
    axes[1].set_ylabel('Experimental', fontsize=10)

    # Difference
    diff = pred_spectrum - exp_spectrum
    axes[2].plot(x, diff, color='#2C3E50', linewidth=0.8)
    axes[2].axhline(0, color='black', linewidth=0.5, linestyle='--')
    axes[2].fill_between(x, diff, 0, where=(diff > 0), alpha=0.3, color='#3498DB', label='Predicted > Exp')
    axes[2].fill_between(x, diff, 0, where=(diff < 0), alpha=0.3, color='#E74C3C', label='Exp > Predicted')
    axes[2].set_ylabel('Difference', fontsize=10)
    axes[2].set_xlabel('13C Chemical Shift (ppm)', fontsize=11)
    axes[2].legend(fontsize=8)

    for ax in axes:
        ax.set_xlim(ppm_range[1], ppm_range[0])  # NMR convention
        ax.grid(True, alpha=0.15)
        ax.set_yticks([])

    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"[VIZ Phase 4] Saved overlay plot → {output_path}")
    return fig


# ---------------------------------------------------------------------------
# 7. NMRPipe-compatible export
# ---------------------------------------------------------------------------

def export_nmrpipe_tab(
    predictions_df: pd.DataFrame,
    output_path: str,
    atoms: List[str] = None,
    sequence: str = None,
) -> pd.DataFrame:
    """
    Export predicted shifts to NMRPipe .tab format for use with nmrDraw / Sparky.

    Format:
        VARS   INDEX RESID RESNAME ATOMNAME SHIFT
        FORMAT %5d   %4d  %4s     %6s       %8.3f

    Returns the exported DataFrame for inspection.
    """
    if atoms is None:
        atoms = ['CA', 'CB', 'C', 'N']

    df = predictions_df[predictions_df['atom'].isin(atoms)].copy()
    df = df.dropna(subset=['predicted_shift'])
    df = df.sort_values(['seq_id', 'atom']).reset_index(drop=True)

    lines = [
        "REMARK  Predicted chemical shifts — Paravastu Pipeline (Phase 4)",
        "REMARK  Generated by viz_module_phase4.py",
        "REMARK",
        "VARS   INDEX RESID RESNAME ATOMNAME SHIFT",
        "FORMAT %5d   %4d  %4s     %6s       %8.3f",
        "",
    ]

    for i, row in df.iterrows():
        lines.append(
            f"{i+1:5d}  {int(row['seq_id']):4d}  {row['residue']:4s}  {row['atom']:6s}  {row['predicted_shift']:8.3f}"
        )

    Path(output_path).write_text('\n'.join(lines))
    print(f"[VIZ Phase 4] Exported NMRPipe .tab → {output_path}")
    return df


# ---------------------------------------------------------------------------
# 8. Phase 4 main runner — generates all new plots
# ---------------------------------------------------------------------------

def generate_phase4_plots(results: dict, output_dir: Path, bmrb_id: str = ''):
    """
    Generate all Phase 4 outputs from a pipeline results dict.

    Adds to results/:
      - spectrum_1d_voigt_{bmrb_id}.png     (enhanced 1D with Voigt + regions)
      - spectrum_2d_cacb_{bmrb_id}.png      (2D CA-CB correlation)
      - predictions_{bmrb_id}.tab           (NMRPipe export)
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    predictions = results.get('predictions')
    sequence    = results.get('sequence')
    label       = bmrb_id or results.get('bmrb_id', 'unknown')

    if predictions is None or predictions.empty:
        print("[VIZ Phase 4] No predictions available — run pipeline first.")
        return

    print(f"\n[VIZ Phase 4] Generating Phase 4 plots → {output_dir}")

    # 1D Voigt spectrum
    try:
        plot_1d_spectrum_phase4(
            predictions,
            sequence=sequence,
            ppm_range=(0, 200),
            sigma=0.15,
            gamma=0.15,
            include_sidechains=(sequence is not None),
            experiment='direct',
            output_path=str(output_dir / f"spectrum_1d_voigt_{label}.png"),
            title=f"BMRB {label} — Simulated 13C Solid-State NMR Spectrum",
        )
    except Exception as e:
        print(f"  [1D] Failed: {e}")

    # 2D CA-CB correlation
    try:
        simulate_2d_cacb(
            predictions,
            output_path=str(output_dir / f"spectrum_2d_cacb_{label}.png"),
            title=f"BMRB {label} — Simulated CA-CB 2D Correlation",
        )
    except Exception as e:
        print(f"  [2D] Failed: {e}")

    # NMRPipe export
    try:
        export_nmrpipe_tab(
            predictions,
            output_path=str(output_dir / f"predictions_{label}.tab"),
            sequence=sequence,
        )
    except Exception as e:
        print(f"  [TAB] Failed: {e}")

    n_png = len(list(output_dir.glob("spectrum_*phase4*.png"))) + \
            len(list(output_dir.glob("spectrum_1d_voigt*.png"))) + \
            len(list(output_dir.glob("spectrum_2d_cacb*.png")))
    print(f"[VIZ Phase 4] Done. New plots: 1D Voigt, 2D CA-CB; NMRPipe .tab exported.")
