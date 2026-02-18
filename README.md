# Paravastu Structural Annotation Pipeline

## Objective

This project builds a structural annotation pipeline for BMRB entries.
The goal is to:

1. Extract protein sequences from BMRB records
2. Analyze residue-level atom coverage
3. Build structural segment annotations (e.g., DSSP-like segments)
4. Create reproducible infrastructure for large-scale analysis

---

## Current State

- Successfully extracted sequence for BMRB 17561
- Computed residue-level atom coverage
- Calculated average atoms per residue
- Set up Git repository with clean structure

Next step: build DSSP segment classification pipeline.

---

## Dependencies

- Python 3.x
- pandas
- numpy
- requests (if used for BMRB fetching)

See requirements.txt for exact versions.

---

## Planned Roadmap

1. Modularize pipeline into src/
2. Build DSSP segment extraction
3. Add visualization layer
4. Generalize to batch processing of multiple BMRB entries
5. Add testing framework

Long-term goal: scalable structural annotation engine.
