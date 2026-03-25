from __future__ import annotations
from dataclasses import dataclass
from enum import Enum
from typing import Dict, List, Optional

class SecondaryStructure(Enum):
    """DSSP-based secondary structure classification."""
    HELIX = "helix"
    STRAND = "strand"
    COIL = "coil"
    UNKNOWN = "unknown"

@dataclass
class ChemicalShift:
    """Represents a single chemical shift measurement."""
    bmrb_id: Optional[int]
    residue_number: int
    residue_name: str
    atom_type: str
    shift_value: float
    secondary_structure: SecondaryStructure = SecondaryStructure.UNKNOWN

@dataclass
class ResidueAssignment:
    """Residue-level mapping between BMRB/PDB and DSSP assignment."""
    residue_number: int
    residue_name: str
    pdb_chain: Optional[str]
    dssp_class: Optional[str]
    secondary_structure: SecondaryStructure = SecondaryStructure.UNKNOWN

@dataclass
class ProteinStructure:
    """Lightweight container for structure metadata and per-residue SS assignments."""
    bmrb_id: Optional[int]
    pdb_id: Optional[str]
    uniprot_id: Optional[str]
    sequence: Optional[str]
    chain_id: Optional[str]
    ss_assignments: Dict[int, ResidueAssignment]

@dataclass
class ShiftStatistics:
    """Aggregated statistics for a residue-atom-secondary-structure combination."""
    residue_one_letter: str
    atom_type: str
    secondary_structure: SecondaryStructure
    count: int
    mean: float
    std: float
    min: float
    max: float
