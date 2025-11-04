"""
PyCBS - Python Complex Band Structure Calculator

A Python package for calculating Complex Band Structure (CBS) from Quantum ESPRESSO
simulation data, replacing the PWCOND Fortran code for CBS-only calculations.
"""

__version__ = "0.1.0"
__author__ = "QE Developers"
__license__ = "GPL-2.0-or-later"

from .calculator import CBSCalculator, read_kpoints_file
from .reader import QEDataReader
from .compbs import ComplexBandStructure
from .wfc_reader import WavefunctionReader, read_wavefunction_metadata, GVectorGrid
from .hamiltonian import HamiltonianBuilder, build_test_hamiltonian
from .potential import (
    LocalPotentialReader,
    PseudopotentialManager,
    integrate_local_potential
)

__all__ = [
    "CBSCalculator",
    "QEDataReader", 
    "ComplexBandStructure",
    "read_kpoints_file",
    "WavefunctionReader",
    "read_wavefunction_metadata",
    "GVectorGrid",
    "HamiltonianBuilder",
    "build_test_hamiltonian",
    "LocalPotentialReader",
    "PseudopotentialManager",
    "integrate_local_potential",
]
