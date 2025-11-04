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


def read_kpoints_file(filename: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    Read k-points from a file in PWCOND format.
    
    File format:
    ```
    nkpts
    kx(1)  ky(1)  weight(1)
    kx(2)  ky(2)  weight(2)
    ...
    kx(nkpts)  ky(nkpts)  weight(nkpts)
    ```
    
    Parameters
    ----------
    filename : str
        Path to k-points file
        
    Returns
    -------
    kpoints : np.ndarray
        Array of k-point coordinates with shape (nkpts, 2)
    weights : np.ndarray
        Array of k-point weights with shape (nkpts,)
        
    Examples
    --------
    >>> kpoints, weights = read_kpoints_file('kpoints.dat')
    >>> calc.set_custom_kpoints(kpoints, weights)
    """
    with open(filename, 'r') as f:
        # Read number of k-points
        nkpts = int(f.readline().strip())
        
        kpoints = []
        weights = []
        
        for _ in range(nkpts):
            line = f.readline().strip()
            if not line:
                break
            parts = line.split()
            if len(parts) >= 3:
                kx, ky, w = float(parts[0]), float(parts[1]), float(parts[2])
                kpoints.append([kx, ky])
                weights.append(w)
            elif len(parts) == 2:
                # No weight specified, will use equal weights
                kx, ky = float(parts[0]), float(parts[1])
                kpoints.append([kx, ky])
                weights.append(1.0)
    
    kpoints = np.array(kpoints)
    weights = np.array(weights)
    
    # Normalize weights if they don't sum to 1
    if len(weights) > 0 and not np.isclose(weights.sum(), 1.0):
        weights = weights / weights.sum()
    
    return kpoints, weights


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
        
        # Manual k-point specification
        self.custom_kpoints = None
        self.custom_weights = None
        
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
        
    def set_custom_kpoints(
        self,
        kpoints: np.ndarray,
        weights: Optional[np.ndarray] = None
    ):
        """
        Set custom k-points for CBS calculation.
        
        This method allows you to manually specify individual k-points
        instead of using an automatic grid, similar to PWCOND input format.
        
        Parameters
        ----------
        kpoints : np.ndarray
            Array of k-point coordinates with shape (nkpts, 2) or (nkpts, 3).
            For 2D arrays, each row is [kx, ky].
            For 3D arrays, each row is [kx, ky, kz] where kz is ignored.
            Coordinates should be in crystal units (fractional).
        weights : np.ndarray, optional
            Array of k-point weights with shape (nkpts,).
            If not provided, equal weights summing to 1.0 are assigned.
            
        Examples
        --------
        Define 4 specific k-points:
        
        >>> kpts = np.array([
        ...     [0.0, 0.0],
        ...     [0.5, 0.0],
        ...     [0.0, 0.5],
        ...     [0.5, 0.5]
        ... ])
        >>> calc.set_custom_kpoints(kpts)
        
        Or with custom weights:
        
        >>> weights = np.array([0.25, 0.25, 0.25, 0.25])
        >>> calc.set_custom_kpoints(kpts, weights)
        
        Read from file format (like PWCOND):
        
        >>> # File contains:
        >>> # 4
        >>> # 0.0000  0.0000  0.25
        >>> # 0.5000  0.0000  0.25
        >>> # 0.0000  0.5000  0.25
        >>> # 0.5000  0.5000  0.25
        >>> data = np.loadtxt('kpoints.dat', skiprows=1)
        >>> calc.set_custom_kpoints(data[:, :2], data[:, 2])
        """
        kpoints = np.asarray(kpoints)
        
        # Handle 3D k-points by taking only first 2 columns
        if kpoints.ndim == 2 and kpoints.shape[1] >= 2:
            self.custom_kpoints = kpoints[:, :2].copy()
        elif kpoints.ndim == 2 and kpoints.shape[1] == 2:
            self.custom_kpoints = kpoints.copy()
        else:
            raise ValueError(
                f"kpoints must have shape (nkpts, 2) or (nkpts, 3), got {kpoints.shape}"
            )
        
        nkpts = self.custom_kpoints.shape[0]
        
        # Set weights
        if weights is None:
            # Equal weights summing to 1.0
            self.custom_weights = np.ones(nkpts) / nkpts
        else:
            weights = np.asarray(weights)
            if weights.shape != (nkpts,):
                raise ValueError(
                    f"weights must have shape ({nkpts},), got {weights.shape}"
                )
            self.custom_weights = weights.copy()
        
        print(f"Set {nkpts} custom k-points")
        
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
        
        if self.custom_kpoints is not None:
            # Use custom k-points
            self.kgrid = KPointGrid.from_custom_kpoints(
                self.custom_kpoints,
                self.custom_weights,
                self.reader.b1,
                self.reader.b2
            )
            print(f"  Using {self.kgrid.nkpts} custom k-points")
        else:
            # Use automatic grid
            self.kgrid = KPointGrid(
                self.nk1, self.nk2,
                self.k1, self.k2,
                self.reader.b1, self.reader.b2
            )
            print(f"  Generated {self.kgrid.nkpts} k-points on {self.nk1}x{self.nk2} grid")
        
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
        
        # Display important warning about current implementation
        print("\n" + "="*70)
        print("WARNING: SIMPLIFIED MODEL IN USE")
        print("="*70)
        print("The current implementation uses a toy tight-binding model for")
        print("demonstration purposes. It does NOT use the actual Hamiltonian")
        print("and potential from your QE calculation.")
        print()
        print("For accurate results, the following needs to be implemented:")
        print("  - Read QE wavefunction data from binary files")
        print("  - Construct proper Hamiltonian and overlap matrices")
        print("  - Calculate actual 2D problem dimensions from system")
        print()
        print("Current results are qualitative only and should NOT be used")
        print("for publication or production work.")
        print("="*70 + "\n")
        
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
