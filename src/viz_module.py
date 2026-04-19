"""
viz_module.py — Visualization tools for NMR chemical shift analysis.

Generates:
- Residue-level shift plots with SS coloring
- Histograms of shift distributions by SS class
- Simulated 1D NMR spectra
- Random coil deviation plots (secondary structure propensity)
- Coverage heatmaps
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # non-interactive backend (safe for all environments)
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
from typing import Optional, Dict, List

# Color scheme for SS classes
SS_COLORS = {
    'helix':   '#E74C3C',   # red
    'strand':  '#3498DB',   # blue
    'coil':    '#95A5A6',   # grey
    'unknown': '#BDC3C7',   # light grey
    'H':       '#E74C3C',
    'E':       '#3498DB',
    'C':       '#95A5A6',
}


# ---------------------------------------------------------------------------
# 1. Per-residue shift plot with SS coloring
# ---------------------------------------------------------------------------

def plot_shifts_by_residue(
    merged_df: pd.DataFrame,
    atom: str = 'CA',
    sequence: str = None,
    output_path: Optional[str] = None,
    title: str = None,
) -> plt.Figure:
    """
    Plot chemical shift vs residue number, colored by secondary structure.
    """
    df = merged_df[merged_df['atom'] == atom].copy()
    if df.empty:
        print(f"No data for atom {atom}")
        return None

    fig, ax = plt.subplots(figsize=(14, 5))

    for ss_class in df['ss_class'].unique():
        sub = df[df['ss_class'] == ss_class]
        color = SS_COLORS.get(ss_class, 'black')
        ax.scatter(sub['seq_id'], sub['shift'], c=color, label=ss_class,
                   s=40, alpha=0.8, zorder=3)

    ax.set_xlabel('Residue Number', fontsize=12)
    ax.set_ylabel(f'{atom} Chemical Shift (ppm)', fontsize=12)
    ax.set_title(title or f'{atom} Shifts by Residue', fontsize=13, fontweight='bold')
    ax.legend(title='SS Class', loc='upper right')
    ax.grid(True, alpha=0.3)
    ax.invert_yaxis()  # NMR convention: downfield on left

    # Annotate sequence if provided
    if sequence:
        ax.set_xticks(range(1, len(sequence)+1, 5))
        ax.set_xticklabels(
            [f"{i}\n{sequence[i-1]}" for i in range(1, len(sequence)+1, 5)],
            fontsize=7
        )

    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Saved: {output_path}")
    return fig


# ---------------------------------------------------------------------------
# 2. Shift distribution histograms by SS class
# ---------------------------------------------------------------------------

def plot_shift_distributions(
    merged_df: pd.DataFrame,
    residue: str,
    atom: str,
    output_path: Optional[str] = None,
    bins: int = 30,
) -> plt.Figure:
    """
    Histogram of shift values for a specific (residue, atom) split by SS class.
    """
    df = merged_df[
        (merged_df['residue'] == residue) &
        (merged_df['atom'] == atom)
    ]
    if df.empty:
        print(f"No data for {residue} {atom}")
        return None

    ss_classes = [s for s in ['helix', 'strand', 'coil', 'unknown'] if s in df['ss_class'].values]

    fig, ax = plt.subplots(figsize=(9, 5))

    for ss_class in ss_classes:
        sub = df[df['ss_class'] == ss_class]['shift']
        if sub.empty:
            continue
        color = SS_COLORS.get(ss_class, 'black')
        ax.hist(sub, bins=bins, alpha=0.6, color=color, label=f"{ss_class} (n={len(sub)})", density=True)

    ax.set_xlabel(f'{atom} Chemical Shift (ppm)', fontsize=12)
    ax.set_ylabel('Density', fontsize=12)
    ax.set_title(f'{residue} — {atom} shift distribution by SS class', fontsize=13, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Saved: {output_path}")
    return fig


# ---------------------------------------------------------------------------
# 3. Simulated 1D NMR spectrum
# ---------------------------------------------------------------------------

def simulate_spectrum(
    shifts: List[float],
    linewidth: float = 0.3,
    ppm_range: tuple = (10, 200),
    n_points: int = 8000,
    output_path: Optional[str] = None,
    title: str = "Simulated 13C Spectrum",
    labels: List[str] = None,
    label_threshold: float = None,
) -> plt.Figure:
    """
    Simulate a 1D NMR spectrum as a sum of Lorentzian peaks.

    shifts: list of chemical shift values (ppm)
    linewidth: peak FWHM in ppm
    """
    x = np.linspace(ppm_range[0], ppm_range[1], n_points)
    spectrum = np.zeros_like(x)

    for shift in shifts:
        if shift is None:
            continue
        # Lorentzian lineshape
        spectrum += (linewidth / 2) ** 2 / ((x - shift) ** 2 + (linewidth / 2) ** 2)

    fig, ax = plt.subplots(figsize=(14, 4))
    ax.plot(x, spectrum, color='#2C3E50', linewidth=1.0)
    ax.fill_between(x, spectrum, alpha=0.15, color='#3498DB')

    # Mark peak positions
    for i, shift in enumerate(shifts):
        if shift is None:
            continue
        ax.axvline(shift, color='#E74C3C', alpha=0.3, linewidth=0.5)
        if labels and i < len(labels) and label_threshold is not None:
            idx = np.argmin(np.abs(x - shift))
            if spectrum[idx] > label_threshold:
                ax.text(shift, spectrum[idx] + 0.01, labels[i], fontsize=6, ha='center')

    ax.set_xlabel('13C Chemical Shift (ppm)', fontsize=12)
    ax.set_ylabel('Intensity (a.u.)', fontsize=12)
    ax.set_title(title, fontsize=13, fontweight='bold')
    ax.invert_xaxis()  # NMR convention
    ax.set_xlim(ppm_range[1], ppm_range[0])
    ax.grid(True, alpha=0.2)
    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Saved: {output_path}")
    return fig


# ---------------------------------------------------------------------------
# 4. Random coil deviation (secondary structure propensity) plot
# ---------------------------------------------------------------------------

def plot_rc_deviation(
    predictions_df: pd.DataFrame,
    atom: str = 'CA',
    sequence: str = None,
    output_path: Optional[str] = None,
) -> plt.Figure:
    """
    Plot CA (or CB) random coil deviation per residue.

    CA deviation > 0 → helix propensity
    CA deviation < 0 → beta-sheet propensity
    Near 0           → coil / random
    """
    df = predictions_df[predictions_df['atom'] == atom].copy()
    df = df.dropna(subset=['rc_deviation'])
    if df.empty:
        print(f"No RC deviation data for {atom}")
        return None

    fig, ax = plt.subplots(figsize=(14, 4))

    colors = ['#E74C3C' if v > 0 else '#3498DB' for v in df['rc_deviation']]
    ax.bar(df['seq_id'], df['rc_deviation'], color=colors, alpha=0.8, width=0.8)
    ax.axhline(0, color='black', linewidth=0.8, linestyle='--')

    # Reference lines
    ax.axhline(2.0, color='#E74C3C', linewidth=0.5, linestyle=':', alpha=0.6)
    ax.axhline(-2.0, color='#3498DB', linewidth=0.5, linestyle=':', alpha=0.6)

    ax.set_xlabel('Residue Number', fontsize=12)
    ax.set_ylabel(f'Δδ {atom} from Random Coil (ppm)', fontsize=12)
    ax.set_title(f'{atom} Secondary Chemical Shift (Δδ > 0 = helix, Δδ < 0 = sheet)',
                 fontsize=12, fontweight='bold')

    # Legend
    helix_patch  = mpatches.Patch(color='#E74C3C', label='Helix propensity (Δδ > 0)')
    sheet_patch  = mpatches.Patch(color='#3498DB', label='Sheet propensity (Δδ < 0)')
    ax.legend(handles=[helix_patch, sheet_patch])
    ax.grid(True, alpha=0.2)
    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Saved: {output_path}")
    return fig


# ---------------------------------------------------------------------------
# 5. Coverage heatmap
# ---------------------------------------------------------------------------

def plot_coverage_heatmap(
    coverage_df: pd.DataFrame,
    sequence: str = None,
    output_path: Optional[str] = None,
) -> plt.Figure:
    """
    Heatmap showing which atoms are observed per residue.
    Rows = atom types, Columns = residue number.
    """
    # Build matrix
    atom_types = sorted({
        atom
        for atoms in coverage_df['atoms_observed']
        for atom in atoms
    })
    max_res = coverage_df['seq_id'].max()

    matrix = np.zeros((len(atom_types), max_res))
    atom_idx = {a: i for i, a in enumerate(atom_types)}

    for _, row in coverage_df.iterrows():
        resnum = int(row['seq_id']) - 1
        if resnum < 0 or resnum >= max_res:
            continue
        for atom in row['atoms_observed']:
            if atom in atom_idx:
                matrix[atom_idx[atom], resnum] = 1

    fig, ax = plt.subplots(figsize=(max(10, max_res // 4), max(3, len(atom_types) // 2)))
    im = ax.imshow(matrix, aspect='auto', cmap='Blues', interpolation='nearest')

    ax.set_yticks(range(len(atom_types)))
    ax.set_yticklabels(atom_types, fontsize=9)
    ax.set_xlabel('Residue Number', fontsize=11)
    ax.set_title('Chemical Shift Coverage (blue = observed)', fontsize=12, fontweight='bold')

    if sequence and len(sequence) == max_res:
        tick_positions = list(range(0, max_res, 5))
        ax.set_xticks(tick_positions)
        ax.set_xticklabels([f"{i+1}\n{sequence[i]}" for i in tick_positions], fontsize=7)

    plt.colorbar(im, ax=ax, fraction=0.02, label='Observed')
    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Saved: {output_path}")
    return fig


# ---------------------------------------------------------------------------
# 6. Convenience: generate all plots for a pipeline result
# ---------------------------------------------------------------------------

def generate_all_plots(results: dict, output_dir: Path):
    """Generate standard set of plots from pipeline results."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    merged_df  = results.get('merged_df')
    stats_df   = results.get('stats_df')
    predictions = results.get('predictions')
    coverage_df = results.get('coverage_df')
    sequence    = results.get('sequence')
    bmrb_id     = results.get('bmrb_id', 'unknown')

    print(f"\n[VIZ] Generating plots → {output_dir}")

    if merged_df is not None:
        for atom in ['CA', 'CB', 'N']:
            if atom in merged_df['atom'].values:
                plot_shifts_by_residue(
                    merged_df, atom=atom, sequence=sequence,
                    output_path=str(output_dir / f"shifts_{atom}_{bmrb_id}.png"),
                    title=f"BMRB {bmrb_id} — {atom} shifts by residue"
                )

    if merged_df is not None:
        # Simulated spectrum from CA + CB shifts
        ca_shifts = merged_df[merged_df['atom'].isin(['CA', 'CB'])]['shift'].dropna().tolist()
        if ca_shifts:
            simulate_spectrum(
                ca_shifts,
                ppm_range=(0, 80),
                output_path=str(output_dir / f"spectrum_13C_{bmrb_id}.png"),
                title=f"BMRB {bmrb_id} — Simulated 13C Spectrum (Cα + Cβ)"
            )

    if predictions is not None and 'rc_deviation' in predictions.columns:
        plot_rc_deviation(
            predictions, atom='CA',
            output_path=str(output_dir / f"rc_deviation_CA_{bmrb_id}.png")
        )

    if coverage_df is not None:
        plot_coverage_heatmap(
            coverage_df, sequence=sequence,
            output_path=str(output_dir / f"coverage_heatmap_{bmrb_id}.png")
        )

    print(f"[VIZ] Done. {len(list(output_dir.glob('*.png')))} plots saved.")
