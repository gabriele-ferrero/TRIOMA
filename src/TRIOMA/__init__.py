"""
TRIOMA: TRItium Object-oriented and Modular Analysis

A Python-based tool for designing outer fuel cycles in fusion reactors.
Provides object-oriented models for tritium transport analysis including
extractors, heat exchangers, breeding blankets, and circuit modeling.
"""

__version__ = "0.1.8"
__author__ = "Gabriele Ferrero, Samuele Meschini"
__email__ = "gabriele.ferrero@polito.it"

# Import main classes for convenience
from TRIOMA.tools.Circuit import Circuit
from TRIOMA.tools.BreedingBlanket import BreedingBlanket

__all__ = [
    "Circuit",
    "BreedingBlanket",
]
