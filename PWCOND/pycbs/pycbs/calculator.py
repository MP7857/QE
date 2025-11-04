"""
Main calculator module for CBS calculations.

This module provides the high-level CBSCalculator class that orchestrates
the entire Complex Band Structure calculation workflow.
"""

import numpy as np
from pathlib import Path
from typing import Optional, Dict, Any, Tuple
import sys

from .reader import QEDataReader
from .compbs import ComplexBandStructure
from .writer import CBSWriter
from .kgrid import KPointGrid


class CBSCalculator:
    """
    Main calculator for Complex Band Structure (CBS) calculations.
    
    This class reads Quantum ESPRESSO output files and calculates the
    complex band structure for a given system.
    
    Parameters
    ----------
    outdir : str or Path
        Temporary directory containing QE output files
    prefix : str
        Prefix for the QE calculation files
    band_file : str, optional
        Output file name for band structure data (without extension)
    ecut2d : float, optional
        2D energy cutoff in Ry (default: 0.0)
    ewind : float, optional
        Energy window parameter in Ry (default: 1.0)
    epsproj : float, optional
        Projection threshold (default: 1e-6)
    """
    
    def __init__(
        self,
        outdir: str,
        prefix: str,
        band_file: Optional[str] = None,
        ecut2d: float = 0.0,
        ewind: float = 1.0,
        epsproj: float = 1e-6,
    ):
        self.outdir = Path(outdir)
        self.prefix = prefix
        self.band_file = band_file
        self.ecut2d = ecut2d
        self.ewind = ewind
        self.epsproj = epsproj
        
        # Data structures
        self.reader = None
        self.compbs = None
        self.writer = None
        self.kgrid = None
        
        # Results storage
        self.results = {}
        
        # Energy and k-point settings
        self.energy0 = None
        self.denergy = None
        self.nenergy = None
        self.earr = None
        
        self.nk1 = 1
        self.nk2 = 1
        self.k1 = 0.0
        self.k2 = 0.0
        
    def set_energy_range(
        self,
        energy0: float,
        denergy: float,
        nenergy: int
    ):
        """
        Set the energy range for CBS calculation.
        
        Parameters
        ----------
        energy0 : float
            Starting energy relative to Fermi level in eV
        denergy : float
            Energy step in eV
        nenergy : int
            Number of energy points
        """
        self.energy0 = energy0
        self.denergy = denergy
        self.nenergy = nenergy
        
        # Create energy array
        self.earr = np.array([
            energy0 + i * denergy for i in range(nenergy)
        ])
        
    def set_kpoints_grid(
        self,
        nk1: int = 1,
        nk2: int = 1,
        k1: float = 0.0,
        k2: float = 0.0
    ):
        """
        Set the 2D k-point grid for CBS calculation.
        
        Parameters
        ----------
        nk1 : int
            Number of k-points along first direction
        nk2 : int
            Number of k-points along second direction
        k1 : float
            Shift along first direction
        k2 : float
            Shift along second direction
        """
        self.nk1 = nk1
        self.nk2 = nk2
        self.k1 = k1
        self.k2 = k2
        
    def read_data(self):
        """Read Quantum ESPRESSO data files."""
        print(f"Reading QE data from {self.outdir}/{self.prefix}")
        self.reader = QEDataReader(self.outdir, self.prefix)
        self.reader.read()
        print(f"  Number of atoms: {self.reader.nat}")
        print(f"  Number of k-points: {self.reader.nk}")
        print(f"  Number of bands: {self.reader.nbnd}")
        
    def setup_kpoints(self):
        """Setup the 2D k-point grid."""
        print("Setting up k-point grid...")
        self.kgrid = KPointGrid(
            self.nk1, self.nk2,
            self.k1, self.k2,
            self.reader.b1, self.reader.b2
        )
        print(f"  Total k-points: {self.kgrid.nkpts}")
        
    def run(self):
        """
        Run the complete CBS calculation workflow.
        
        Returns
        -------
        dict
            Dictionary containing calculation results
        """
        # Validate input
        if self.energy0 is None:
            raise ValueError("Energy range not set. Call set_energy_range() first.")
            
        # Read data
        self.read_data()
        
        # Setup k-points
        self.setup_kpoints()
        
        # Initialize CBS calculator
        self.compbs = ComplexBandStructure(
            self.reader,
            ecut2d=self.ecut2d,
            ewind=self.ewind,
            epsproj=self.epsproj
        )
        
        # Initialize writer if output requested
        if self.band_file:
            self.writer = CBSWriter(self.band_file)
            self.writer.open()
        
        # Main calculation loop
        print("\nStarting CBS calculation...")
        for ik in range(self.kgrid.nkpts):
            kpoint = self.kgrid.get_kpoint(ik)
            print(f"\nK-point {ik+1}/{self.kgrid.nkpts}: ({kpoint[0]:.6f}, {kpoint[1]:.6f})")
            
            for ien in range(self.nenergy):
                energy = self.earr[ien]
                print(f"  Energy {ien+1}/{self.nenergy}: E-Ef = {energy:.4f} eV")
                
                # Calculate CBS at this (k, E) point
                result = self.compbs.calculate(kpoint, energy)
                
                # Store results
                key = (ik, ien)
                self.results[key] = result
                
                # Write output if requested
                if self.writer:
                    self.writer.write_point(
                        ik, ien,
                        result['kvals'],
                        result['nchan'],
                        energy,
                        self.kgrid.nkpts,
                        self.nenergy
                    )
                
                # Print summary
                print(f"    Channels: {result['nchan']}")
                
        # Close output files
        if self.writer:
            self.writer.close()
            
        print("\nCBS calculation completed successfully!")
        return self.results
        
    def get_results(self) -> Dict[Tuple[int, int], Dict[str, Any]]:
        """
        Get the calculation results.
        
        Returns
        -------
        dict
            Dictionary with (ik, ien) tuples as keys and result dicts as values
        """
        return self.results
