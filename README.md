# Paravastu Lab — NMR Structural Annotation Pipeline

A computational pipeline that links experimental NMR chemical shift data with 3D protein structural information, enabling shift prediction, structure validation, and automated spectral assignment. Built for the Paravastu Lab to accelerate protein characterization from solid-state NMR data.

---

## Table of Contents

1. [What This Does](#what-this-does)
2. [Why It Matters](#why-it-matters)
3. [How It Works — The Science](#how-it-works--the-science)
4. [Project Structure](#project-structure)
5. [Setup](#setup)
6. [Running the Pipeline](#running-the-pipeline)
7. [Understanding the Outputs](#understanding-the-outputs)
8. [Using the Python API](#using-the-python-api)
9. [Module Reference](#module-reference)
10. [Troubleshooting](#troubleshooting)
11. [Current Status and Roadmap](#current-status-and-roadmap)

---

## What This Does

This pipeline takes two inputs:

1. **NMR chemical shift data** from the [BMRB](https://bmrb.io) (Biological Magnetic Resonance Bank)
2. **3D protein structure** from the [PDB](https://rcsb.org) (Protein Data Bank)

And produces:

- A **reference database** of chemical shift statistics split by secondary structure (helix / strand / coil)
- **Shift predictions**: given a PDB structure, predict what its NMR spectrum should look like
- **Peak assignment candidates**: given an experimental peak position, rank which residues it likely belongs to
- **ML models** that predict secondary structure class from chemical shift patterns
- **Publication-quality plots**: shift scatter plots, simulated spectra, secondary shift charts, coverage heatmaps

### Current dataset

- BMRB entry: **17561** (85-residue cysteine-rich peptide, neurophysin family)
- PDB structure: **2LBH** (NMR structure, 2 chains x 92 residues)
- Shift records parsed: **948** (CA, CB, N, H, HA, HB atoms)
- Stats rows generated: **270** (split by residue x atom x SS class)

---

## Why It Matters

Proteins control nearly everything in biology — brain signaling, hormones, immune function, cell repair. When proteins misfold, the result is disease: Alzheimer's, Parkinson's, Huntington's, prion disease, amyloidosis.

NMR chemical shifts are one of the fastest ways to observe protein structure in solution, but a raw chemical shift number (e.g., CA = 55.2 ppm) means little in isolation.

**Linked to structure, it becomes informative:**

| Observation | Meaning |
|---|---|
| CA shift above random coil | Helical tendency |
| CA shift below random coil | Beta-sheet tendency |
| Large CB change vs random coil | Structural transition |
| Amide (N, H) shifts | Hydrogen bonding environment |

This pipeline builds the systematic link between shift numbers and structural meaning — enabling faster protein characterization, structure validation, and drug discovery support at lower cost than solving full structures every time.

---

## How It Works — The Science

### Step 1: Parse chemical shifts from BMRB

The BMRB stores NMR-STAR formatted files containing all measured chemical shifts for a protein. We parse these into a table: `seq_id | residue | atom | shift`.

### Step 2: Assign secondary structure with DSSP

DSSP (Define Secondary Structure of Proteins) reads a 3D PDB structure and assigns each residue one of:

| DSSP Code | Meaning | Coarse class |
|---|---|---|
| H | alpha-helix | helix |
| G | 3-10 helix | helix |
| I | pi-helix | helix |
| E | beta-strand | strand |
| B | beta-bridge | strand |
| T | Turn | coil |
| S | Bend | coil |
| C | Random coil | coil |

These become the **ground-truth structural labels** for every residue. If mkdssp cannot run, the pipeline automatically reads HELIX and SHEET records directly from the CIF file — these deposited annotations provide equivalent helix/strand/coil assignments.

### Step 3: Align BMRB sequence to PDB sequence

The BMRB sample and the deposited PDB structure may differ slightly — signal peptides may be cleaved, constructs may be longer. We use pairwise global sequence alignment to map BMRB residue numbers to PDB residue numbers so shifts land on the correct structural labels.

For 2LBH: BMRB has 85 residues, PDB has 92 residues. The alignment correctly maps BMRB residues 1-85 to their corresponding PDB positions, identifying that the PDB has a 7-residue N-terminal extension.

### Step 4: Build a reference database

After merging, we compute statistics for every (residue type, atom, SS class) combination:

```
residue | atom | ss_class | count | mean | median | std | min | max | rc_deviation
    V   |  CA  |   helix  |   2   | 62.96 | 62.96 | 0.5 | ...
    V   |  CA  |  strand  |   1   | 60.22 | 60.22 | 0.0 | ...
    V   |  CA  |   coil   |   1   | 62.05 | 62.05 | 0.0 | ...
```

The `rc_deviation` column is the **secondary chemical shift** — the deviation from the theoretical random coil value. This is the core structural signal:

- **rc_deviation > 0 for CA** means helical tendency
- **rc_deviation < 0 for CA** means beta-sheet tendency
- **Near 0** means coil / disordered

### Step 5: Predict shifts from structure

Given a PDB model and its DSSP labels, predict the expected chemical shift at every position by looking up the mean for that (residue, atom, SS class) combination. This lets you ask: "Does this PDB model agree with my experimental spectrum?"

### Step 6: Assign peaks to residues

Given an experimental peak at, say, 55.3 ppm for CA, score all possible residue assignments using a Gaussian probability model. The residue with the highest log-probability is the most likely candidate.

### Step 7: Train ML models

Build feature matrices from the shift data (raw shifts, RC deviations, context window of neighboring shifts) and train RandomForest + XGBoost classifiers to predict helix / strand / coil from shifts alone. The most informative features discovered from real data:

1. `dev_CA` — Ca deviation from random coil (strongest signal, importance 0.16)
2. `dev_N` — 15N deviation from random coil (importance 0.10)
3. `dev_CB` — Cb deviation from random coil (importance 0.08)
4. `shift_N` — raw 15N shift value (importance 0.08)
5. Context window shifts (neighboring residues Ca and Cb, importance 0.07 each)

---

## Project Structure

```
paravastu/
├── src/
│   ├── pipeline.py       <- Main orchestrator — run this
│   ├── bmrb_module.py    <- BMRB fetching and NMR-STAR parsing
│   ├── dssp_module.py    <- PDB download, DSSP extraction, sequence alignment
│   ├── stats_module.py   <- Statistics, ShiftPredictor, peak assignment
│   ├── viz_module.py     <- All plots and spectrum simulation
│   ├── ml_module.py      <- Feature engineering, RandomForest, XGBoost
│   └── models.py         <- Dataclasses
├── data/
│   ├── bmr17561_3.str              <- BMRB NMR-STAR file (auto-downloaded)
│   ├── 2LBH.cif                   <- PDB structure in mmCIF format (auto-downloaded)
│   └── bmrb_shift_stats_by_residue_atom_ss.csv  <- Pre-existing stats CSV
├── results/                        <- All outputs (auto-generated each run)
│   ├── merged_shifts_17561.csv     <- Shifts + SS labels merged
│   ├── stats_17561.csv             <- Reference database (SS-split statistics)
│   ├── predictions_17561.csv       <- Predicted shifts + RC deviations
│   ├── coverage_17561.csv          <- Which atoms observed per residue
│   ├── dssp_2LBH.csv               <- DSSP assignments
│   ├── model_rf.joblib             <- Trained RandomForest model
│   ├── model_xgb.joblib            <- Trained XGBoost model
│   ├── shifts_CA_17561.png         <- Ca shift plot colored by SS
│   ├── shifts_CB_17561.png         <- Cb shift plot
│   ├── shifts_N_17561.png          <- 15N shift plot
│   ├── spectrum_13C_17561.png      <- Simulated 13C NMR spectrum
│   ├── rc_deviation_CA_17561.png   <- Secondary chemical shift chart
│   └── coverage_heatmap_17561.png  <- Atom coverage per residue
├── notebooks/
│   └── analysis.ipynb              <- Interactive Jupyter notebook
├── requirements.txt
└── README.md
```

---

## Setup

### Prerequisites

- Windows, macOS, or Linux
- [Anaconda](https://www.anaconda.com/download) or Miniconda

### 1. Create a dedicated conda environment

Open Anaconda Prompt and run:

```bash
conda create -n nmr python=3.10
conda activate nmr
```

### 2. Install DSSP

```bash
conda install -c salilab dssp
```

Verify it installed:

```bash
mkdssp --version
```

Expected output: `mkdssp version 4.x.x`

If mkdssp fails with `mmcif_pdbx dictionary not found` when you run the pipeline, this does not block the pipeline — it automatically falls back to reading HELIX/SHEET records from the CIF file, which gives equivalent helix/strand/coil labels. To fix mkdssp fully, download the Windows binary from https://github.com/PDB-REDO/dssp/releases (it bundles all required libraries).

### 3. Install Python packages

```bash
pip install biopython pynmrstar pandas numpy scikit-learn xgboost matplotlib joblib requests
```

### 4. Set up the project folder

Your folder should look like this before running:

```
paravastu/
├── src/          <- all .py files go here
├── data/         <- create this folder; pipeline downloads files into it
└── results/      <- create this folder; all outputs go here
```

Create the folders:

```bash
mkdir data
mkdir results
```

---

## Running the Pipeline

**Always run from the `paravastu/` root directory** (not from inside `src/`).

```bash
conda activate nmr
cd C:\Users\YOUR_NAME\paravastu      # adjust path to your machine
```

### Option A — Full pipeline from BMRB (recommended)

Downloads raw shifts from BMRB, downloads PDB structure, runs DSSP, computes SS-aware statistics, generates all 6 plots, trains ML models. Everything automated.

```bash
python src/pipeline.py --bmrb 17561 --pdb 2LBH
```

This runs all 9 steps and saves everything to `results/`.

### Option B — From an existing shifts CSV

If you already have a CSV file of shift data:

```bash
python src/pipeline.py --csv data/bmrb_shift_stats_by_residue_atom_ss.csv --pdb 2LBH --seq VLDLDVRTCPCGPGGKGRCFGICCGDELGCFVGTAEALRCQEENYLPCQSGQKPCGSGGRCAAAGICCPDGCHEDPACDPEAAFS
```

The pipeline auto-detects whether your CSV is a raw shifts CSV (one row per observation, has `seq_id` and `shift` columns) or a stats CSV (pre-aggregated, has `mean` and `std` columns). It handles both correctly.

### Option C — Skip DSSP (shifts and stats only)

```bash
python src/pipeline.py --bmrb 17561 --no-dssp
```

All statistics will have `ss_class=unknown`. Useful for a quick data quality check or when no PDB structure is available.

### Option D — Batch mode (multiple BMRB entries)

To build a larger reference database, run on multiple entries at once:

```python
# In Python or in a script
import sys
sys.path.insert(0, 'src')
from pipeline import run_batch

combined_stats = run_batch(
    bmrb_ids=[17561, 17562, 17563, 17564],
    pdb_map={17561: '2LBH'},  # optional: override PDB for specific entries
)
# Saves combined_stats.csv to results/
```

### Command-line options

| Flag | Description | Default |
|---|---|---|
| `--bmrb N` | BMRB entry ID to fetch and parse | 17561 |
| `--pdb XXXX` | PDB ID to download and run DSSP on | 2LBH |
| `--csv path` | Use a local CSV instead of fetching BMRB | None |
| `--seq STRING` | Protein sequence, one-letter codes (used with --csv) | None |
| `--no-dssp` | Skip DSSP extraction entirely | False |
| `--no-save` | Run pipeline but do not write any output files | False |

---

## Understanding the Outputs

After a successful run, your `results/` folder contains the following.

### CSV files

**`merged_shifts_17561.csv`** — The core dataset. One row per shift observation, with its SS class label attached. This is what you train ML models on and what you analyze.

Columns: `seq_id, residue_name, residue, atom, shift, ss_class`

```
seq_id  residue  atom  shift    ss_class
37      A        CA    52.3     helix
37      A        N     124.2    helix
21      C        CA    55.7     strand
```

**`stats_17561.csv`** — The reference database. For each (residue type, atom, SS class) combination, the mean, std, count, and random coil deviation. This is what the ShiftPredictor uses for predictions and what you use for peak assignment.

Columns: `residue, atom, ss_class, count, mean, median, std, min, max, rc_deviation`

**`predictions_17561.csv`** — Predicted chemical shifts for every residue in the sequence, based on its SS class from the PDB structure. Includes `rc_deviation` (secondary chemical shift) for each prediction.

Columns: `seq_id, residue, atom, ss_class, predicted_shift, rc_ref, rc_deviation`

**`coverage_17561.csv`** — Shows how many and which atom types are measured per residue. Useful for identifying gaps in the data.

Columns: `seq_id, residue, atoms_observed, n_atoms`

**`dssp_2LBH.csv`** — One row per residue with its SS assignment. If full DSSP ran: also includes phi, psi angles and solvent accessibility. If fallback was used: SS class only.

### Plots

**`rc_deviation_CA_17561.png`** — The most scientifically important plot. Shows the secondary chemical shift (observed Cα minus random coil reference) for every residue. Bars pointing up (red) indicate helical regions; bars pointing down (blue) indicate beta-sheet regions. This is how you identify structural regions from shifts alone.

**`shifts_CA_17561.png`** — Cα chemical shift value for every residue, colored by secondary structure class. Red dots = helix, blue = strand, grey = coil. You can see that helical residues cluster at higher ppm and strand residues at lower ppm.

**`shifts_CB_17561.png`** — Same as above for Cβ. Cβ shifts are particularly diagnostic for amino acid type and secondary structure.

**`shifts_N_17561.png`** — Same for 15N backbone amide shifts. These reflect hydrogen bonding environments.

**`spectrum_13C_17561.png`** — A simulated 1D 13C NMR spectrum generated by placing Lorentzian peaks at all observed Cα and Cβ shift positions. This lets you visually compare expected vs measured spectrum.

**`coverage_heatmap_17561.png`** — A grid showing which atoms are observed (blue) or missing (white) for each residue. Shows where your data is sparse.

### Trained models

**`model_rf.joblib`** and **`model_xgb.joblib`** — Trained scikit-learn and XGBoost models. Load with `joblib.load('results/model_rf.joblib')`. The dict contains the model object, CV scores, and feature importances. Accuracy improves substantially with more training data from additional BMRB entries.

---

## Using the Python API

### Run the full pipeline from Python

```python
import sys
sys.path.insert(0, 'src')
from pipeline import run_pipeline

results = run_pipeline(bmrb_id=17561, pdb_id='2LBH')

# Access results
merged_df   = results['merged_df']    # shifts + SS labels
stats_df    = results['stats_df']     # reference database
predictions = results['predictions']  # predicted shifts
dssp_df     = results['dssp_df']      # DSSP per residue
ml          = results['ml']           # trained models
```

### Predict shifts for a new protein

```python
import sys
sys.path.insert(0, 'src')
import pandas as pd
from stats_module import ShiftPredictor, add_random_coil_deviation

# Load the reference database built from 17561 + any other entries you've run
stats_df = pd.read_csv('results/stats_17561.csv')

# Define your protein sequence and known/predicted secondary structure
sequence = 'VLDLDVRTCPCGPGGKGRCFGICCGDELGCFVGTAEALRCQEENYLPCQSGQKPCGSGGRCAAAGICCPDGCHEDPACDPEAAFS'

# ss_map: {residue_number (1-based) -> 'helix', 'strand', or 'coil'}
ss_map = {i+1: 'coil' for i in range(len(sequence))}
# Mark known helix region (residues 37-46 in BMRB numbering)
for i in range(37, 47):
    ss_map[i] = 'helix'

# Predict
predictor = ShiftPredictor(stats_df)
predictions = predictor.predict(sequence, ss_map, atoms=['CA', 'CB', 'N'])
print(predictions.head(20))
```

### Assign an experimental peak to a residue

```python
from stats_module import rank_assignments

# You observe a peak at 55.3 ppm for CA — which residue is most likely?
candidates = rank_assignments(
    peak_shift=55.3,
    atom='CA',
    sequence=sequence,
    ss_map=ss_map,
    stats_df=stats_df,
    top_n=5,
)
print(candidates)
# Output: seq_id | residue | ss_class | log_prob | rank
# Higher log_prob = more likely assignment
```

### Generate plots manually

```python
from viz_module import plot_shifts_by_residue, simulate_spectrum, plot_rc_deviation
import pandas as pd

merged = pd.read_csv('results/merged_shifts_17561.csv')

# CA shift colored by SS class
plot_shifts_by_residue(
    merged, atom='CA',
    sequence='VLDLDVRTCPCGPGGKGRCFGICCGDELGCFVGTAEALRCQEENYLPCQSGQKPCGSGGRCAAAGICCPDGCHEDPACDPEAAFS',
    output_path='results/my_ca_plot.png'
)

# Simulated 13C spectrum
shifts = merged[merged['atom'].isin(['CA','CB'])]['shift'].tolist()
simulate_spectrum(shifts, ppm_range=(0, 80), output_path='results/my_spectrum.png')
```

### Load and use a trained model

```python
import joblib
import pandas as pd
from ml_module import build_feature_matrix

# Load model
rf_data = joblib.load('results/model_rf.joblib')
model = rf_data['model']

# Build features from new merged data
merged = pd.read_csv('results/merged_shifts_17561.csv')
X, y = build_feature_matrix(merged, window_size=1)
X_filled = X.fillna(X.mean())

# Predict secondary structure from shifts
predictions = model.predict(X_filled)
print(predictions)  # array of 'helix', 'strand', 'coil'
```

---

## Module Reference

### `pipeline.py` — Main orchestrator

| Function | Description |
|---|---|
| `run_pipeline(bmrb_id, pdb_id, run_dssp, save_results)` | Full 9-step pipeline from BMRB ID. Returns results dict. |
| `run_pipeline_from_csv(csv_path, pdb_id, sequence, ...)` | Same pipeline starting from a local CSV. |
| `run_batch(bmrb_ids, pdb_map, run_dssp)` | Run on multiple entries, combine stats. |
| `merge_shifts_with_dssp(shifts_df, dssp_map, seq_mapping)` | Join shift data with DSSP structural labels. |

### `bmrb_module.py` — BMRB access

| Function | Description |
|---|---|
| `fetch_bmrb_text(bmr_id, cache_dir)` | Download NMR-STAR file, cache locally for reuse |
| `load_bmrb_from_file(path)` | Load a locally saved .str file |
| `parse_shifts(nmrstar_text)` | Extract all chemical shifts into a DataFrame |
| `extract_sequence(nmrstar_text)` | Infer amino acid sequence from shift assignments |
| `extract_pdb_id(nmrstar_text)` | Find associated PDB ID from BMRB metadata |
| `compute_coverage(df)` | Count unique atoms observed per residue |
| `compute_shift_stats(df)` | Compute mean/std/count grouped by (residue, atom, ss_class) |

### `dssp_module.py` — Structure and alignment

| Function | Description |
|---|---|
| `download_pdb(pdb_id, data_dir)` | Download .cif (preferred) or .pdb from RCSB |
| `extract_dssp_full(pdb_path)` | Run mkdssp, return DataFrame with phi/psi/accessibility |
| `extract_ss_from_pdb_records(pdb_path)` | Fallback: parse HELIX/SHEET from CIF or PDB file |
| `build_ss_segments(ss_map)` | Convert per-residue SS map into contiguous segments |
| `align_and_map(bmrb_seq, pdb_seq)` | Global alignment, returns {bmrb_id: pdb_position} mapping |
| `build_prediction_ss_map(dssp_df, seq_mapping, length)` | Build BMRB-indexed SS map for use in predictions |
| `get_pdb_sequence_from_file(pdb_path)` | Extract amino acid sequence from PDB or CIF file |
| `find_mkdssp()` | Locate mkdssp executable (searches PATH and conda envs) |

### `stats_module.py` — Analysis and prediction

| Function / Class | Description |
|---|---|
| `compute_shift_stats(df)` | Aggregate statistics by (residue, atom, ss_class) |
| `add_random_coil_deviation(df)` | Add rc_deviation column (observed minus random coil reference) |
| `RANDOM_COIL` | Dict of random coil reference values (Wishart 1995 / BMRB averages) |
| `ShiftPredictor(stats_df)` | Class: predicts shifts for a full protein sequence |
| `.predict(seq, ss_map, atoms, method)` | Predict shifts for entire sequence given SS assignments |
| `.predict_uncertainty(seq, ss_map, n_samples)` | Monte Carlo uncertainty estimates |
| `rank_assignments(peak, atom, seq, ss_map, stats, top_n)` | Rank all residues by assignment likelihood for a peak |
| `score_peak_assignment(peak, residue, atom, ss, stats)` | Gaussian log-probability score for one assignment |

### `viz_module.py` — Visualization

| Function | Output |
|---|---|
| `plot_shifts_by_residue(df, atom, sequence, output_path)` | Scatter plot of shifts vs position, colored by SS class |
| `plot_shift_distributions(df, residue, atom)` | Histogram of shift values split by SS class |
| `simulate_spectrum(shifts, linewidth, ppm_range)` | Simulated 1D NMR spectrum as Lorentzian peaks |
| `plot_rc_deviation(df, atom, sequence)` | Bar chart of secondary chemical shifts |
| `plot_coverage_heatmap(df, sequence)` | Heatmap of atom observation per residue |
| `generate_all_plots(results, output_dir)` | Generates all 5 standard plots in one call |

### `ml_module.py` — Machine learning

| Function | Description |
|---|---|
| `build_feature_matrix(df, window_size, include_rc_deviation)` | Build X (features) and y (labels) from merged shifts |
| `train_random_forest(X, y, n_estimators, cv_folds)` | Train RF with stratified cross-validation |
| `train_xgboost(X, y, cv_folds)` | Train XGBoost with stratified cross-validation |
| `evaluate_model(model, X, y)` | Held-out test set evaluation with confusion matrix |
| `save_model(model_dict, path)` | Save model to .joblib file |
| `load_model(path)` | Load model from .joblib file |
| `run_ml_pipeline(merged_df, results_dir, window_size)` | End-to-end: features, train, evaluate, save both models |

---

## Troubleshooting

### BMRB download fails or times out

```
ERROR: Failed to fetch BMRB 17561
```

Download the file manually from https://bmrb.io, search for entry 17561, click Download, choose NMR-STAR v3 format. Save the file as `data/bmr17561_3.str`. The pipeline detects and uses it automatically on the next run.

### mkdssp fails with dictionary error

```
Error while loading dictionary mmcif_pdbx
```

This does not stop the pipeline. It prints this message and automatically falls back to reading HELIX/SHEET records directly from the CIF file. You still get correct helix/strand/coil labels for all residues.

To fix mkdssp permanently: download the pre-built Windows binary from https://github.com/PDB-REDO/dssp/releases — the release binary bundles all required libraries. Replace the `mkdssp.exe` in your conda env's `bin/` folder with this one.

### ModuleNotFoundError

```
ModuleNotFoundError: No module named 'pynmrstar'
```

You are either not in the `nmr` conda environment or packages were not installed. Run:

```bash
conda activate nmr
pip install pynmrstar biopython pandas numpy scikit-learn xgboost matplotlib joblib requests
```

### KeyError: 'seq_id'

You passed a stats CSV (pre-aggregated, no seq_id column) to a function that expects raw shift rows. The pipeline's `run_pipeline_from_csv` handles this automatically, but if calling functions directly, use the correct input format.

Raw shifts CSV columns: `seq_id, residue, atom, shift, ss_class`
Stats CSV columns: `residue, atom, ss_class, count, mean, median, std, min, max`

### FileNotFoundError when running pipeline.py

```
FileNotFoundError: No CSV file found in data/ directory
```

You are running from the wrong directory. Always run from `paravastu/`:

```bash
cd C:\Users\YOUR_NAME\paravastu    # correct - paravastu is the working dir
python src/pipeline.py ...         # src/ is a subdirectory
```

Not from inside src/:

```bash
cd src
python pipeline.py ...             # WRONG - relative paths break
```

### Predictions all show ss_class = coil or unknown

This happens if the saved predictions CSV is from a previous run before the SS alignment fix. Delete `results/predictions_17561.csv` and re-run. The pipeline always overwrites results.

---

## Current Status and Roadmap

### What is complete and working

- BMRB fetch and NMR-STAR parsing (handles solid-state NMR format quirks)
- Sequence extraction inferred from shift assignments (works even without Entity_poly_seq in BMRB)
- Residue-level atom coverage analysis
- PDB/mmCIF download (automatically prefers CIF for mkdssp 4.x compatibility)
- DSSP extraction with full automatic fallback to PDB HELIX/SHEET records (CIF and PDB formats both supported)
- BMRB to PDB sequence alignment using global pairwise alignment, handles N-terminal extensions and length mismatches
- Chemical shift statistics split by secondary structure class (helix/strand/coil)
- Random coil deviation (secondary chemical shift) analysis with Wishart 1995 reference values
- ShiftPredictor: given a PDB structure and DSSP, predict expected NMR spectrum
- Peak assignment ranking using Gaussian probability model
- RandomForest and XGBoost secondary structure classifiers with cross-validation
- Full visualization suite (6 plot types), automatically generated on every pipeline run
- ML training automatically triggered on every pipeline run when enough labeled data exists
- Batch processing mode for building larger reference databases from multiple BMRB entries
- Interactive Jupyter notebook for exploration

### Planned next steps

**Expand the reference database** — the single most impactful improvement. Currently trained on one 85-residue protein (8 helix, 11 strand, 66 coil examples). Run `run_batch()` on 20-50 more solid-state NMR entries from BMRB to dramatically improve ML accuracy and shift prediction quality.

**2D correlation assignment** — use CA-CB pairs jointly for much more discriminating peak assignment. CA-CB chemical shift pairs are highly specific to amino acid type and SS context.

**ESMFold integration** — predict structure from sequence alone for BMRB entries that have no deposited PDB structure. Run DSSP on the ESMFold model to get approximate SS labels.

**Bayesian assignment engine** — full probabilistic framework that uses sequence composition constraints (if there are 6 serines, assign at most 6 serine peaks), peak intensities, and 2D correlations jointly.

**Disorder detection** — flag residues with anomalous shift distributions (large std, outlier values) as potentially disordered or dynamically flexible.

**Disease protein analysis** — apply pipeline to amyloid-forming and aggregation-prone sequences relevant to Alzheimer's, Parkinson's, and related diseases.

---

## Data Sources

| Resource | URL | What we use |
|---|---|---|
| BMRB | https://bmrb.io | NMR chemical shifts in NMR-STAR format |
| RCSB PDB | https://rcsb.org | 3D protein structures in mmCIF format |
| DSSP | https://github.com/PDB-REDO/dssp | Secondary structure assignment algorithm |

---

## Dependencies

| Package | Version | Purpose |
|---|---|---|
| pynmrstar | >= 3.3 | NMR-STAR file parsing |
| biopython | >= 1.81 | PDB/CIF parsing, sequence alignment |
| pandas | >= 2.0 | Data processing and statistics |
| numpy | >= 1.24 | Numerical computation |
| scikit-learn | >= 1.3 | RandomForest, cross-validation, preprocessing |
| xgboost | >= 2.0 | Gradient boosting classifier |
| matplotlib | >= 3.7 | All plots and spectrum simulation |
| requests | >= 2.28 | HTTP downloads from BMRB and RCSB |
| joblib | >= 1.3 | Model serialization and loading |
