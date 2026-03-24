from dataclasses import dataclass
from typing import List

dataclass class Molecule:
    name: str
    formula: str
    atoms: List[str]

dataclass class ShiftPrediction:
    molecule: Molecule
    shift_values: List[float]
    conditions: str  # e.g., temperature, solvent

dataclass class ExperimentalData:
    molecule: Molecule
    measured_shifts: List[float]
    experimental_conditions: str
