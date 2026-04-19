"""
models.py — Core dataclasses for the Paravastu pipeline.
"""

from dataclasses import dataclass, field
from typing import List, Optional, Dict


@dataclass
class ShiftRecord:
    """A single chemical shift observation."""
    seq_id: int
    residue: str          # one-letter code
    residue_name: str     # three-letter code
    atom: str             # e.g. CA, CB, N, H
    shift: float          # ppm value
    ss_class: str = "unknown"  # H, E, C, or unknown


@dataclass
class ResidueStats:
    """Aggregated shift statistics for a (residue, atom, ss_class) triple."""
    residue: str
    atom: str
    ss_class: str
    count: int
    mean: float
    median: float
    std: float
    min: float
    max: float


@dataclass
class DSSPRecord:
    """Secondary structure assignment for one residue."""
    residue_number: int
    chain_id: str
    ss_raw: str       # raw DSSP code: H G I E B T S C
    ss_class: str     # coarse: helix / strand / coil


@dataclass
class StructureSegment:
    """A contiguous run of the same secondary structure."""
    start: int
    end: int
    ss_class: str

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclass
class PipelineResult:
    """Full result of running the pipeline on one BMRB entry."""
    bmrb_id: int
    pdb_id: Optional[str]
    sequence: Optional[str]
    shifts_df: object        # pandas DataFrame
    stats_df: object         # pandas DataFrame
    dssp_map: Dict[int, str] = field(default_factory=dict)
    segments: List[StructureSegment] = field(default_factory=list)
    merged_df: object = None  # shifts + ss_class joined
