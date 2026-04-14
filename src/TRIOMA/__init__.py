"""
TRIOMA: TRItium Object-oriented and Modular Analysis

A Python-based tool for designing outer fuel cycles in fusion reactors.
Provides object-oriented models for tritium transport analysis including
extractors, heat exchangers, breeding blankets, and circuit modeling.
"""

__version__ = "0.1.8"
__author__ = "Gabriele Ferrero, Samuele Meschini"
__email__ = "gabriele.ferrero@polito.it"

# Import the tools module which handles all submodule imports
from TRIOMA import tools

# Re-export for convenience
Circuit = tools.Circuit
BreedingBlanket = tools.BreedingBlanket
TriomaClass = tools.TriomaClass

__all__ = [
    "Circuit",
    "BreedingBlanket",
    "TriomaClass",
]
