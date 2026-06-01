"""Common type definitions for the uONT package."""

from typing import NewType, Optional
from dataclasses import dataclass

FullPath = NewType("FullPath", str)
FullFuturePath = NewType("FullFuturePath", str)


class ReferenceSequence:
    def __init__(self, name: str, sequence: str):
        self.name = name
        self.sequence = sequence


@dataclass
class QCMetrics:
    reads_total_number: Optional[int] = None
    reads_n50: Optional[int] = None
    reads_total_bases: Optional[int] = None
    contigs_total_length: Optional[int] = None
    contigs_total_number: Optional[int] = None
    contigs_n50: Optional[int] = None
    contigs_l90: Optional[int] = None
    contigs_Ns_per_kb: Optional[float] = None
    contigs_gc_content: Optional[float] = None
    genome_depth_estimate: Optional[float] = None
    def to_dict(self):
        return {
            "reads_total_number": self.reads_total_number,
            "reads_n50": self.reads_n50,
            "reads_total_bases": self.reads_total_bases,
            "contigs_total_length": self.contigs_total_length,
            "contigs_total_number": self.contigs_total_number,
            "contigs_n50": self.contigs_n50,
            "contigs_l90": self.contigs_l90,
            "contigs_Ns_per_kb": self.contigs_Ns_per_kb,
            "contigs_gc_content": self.contigs_gc_content,
            "genome_depth_estimate": self.genome_depth_estimate
            
        }