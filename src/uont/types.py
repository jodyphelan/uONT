"""Common type definitions for the uONT package."""

from typing import NewType

FullPath = NewType("FullPath", str)
FullFuturePath = NewType("FullFuturePath", str)


class ReferenceSequence:
    def __init__(self, name: str, sequence: str):
        self.name = name
        self.sequence = sequence
