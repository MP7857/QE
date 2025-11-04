"""
Core Complex Band Structure calculation module.

This module implements the main CBS algorithm, following the PWCOND approach
from compbs_2.f90 and related routines.
"""

import numpy as np
from scipy import linalg
from typing import Dict, Any, Tuple
import warnings


class ComplexBandStructure:
    """
    Complex Band Structure calculator.
    
    This class implements the core CBS calculation algorithm, solving the
    generalized eigenvalue problem to find complex k-vectors for given
    energies and 2D k-points.
    
    Parameters
    ----------
    reader : QEDataReader
        QE data reader with loaded system information
    ecut2d : float
        2D energy cutoff in Ry
    ewind : float
        Energy window parameter in Ry
    epsproj : float
        Projection threshold
    gvec_grid : GVectorGrid, optional
        Pre-constructed G-vector grid (Phase 1 integration)
    """
    
    def __init__(
        self,
        reader,
        ecut2d=0.0,
        ewind=1.0,
        epsproj=1e-6,
        gvec_grid=None
    ):
        self.reader = reader
        self.ecut2d = ecut2d
        self.ewind = ewind
        self.epsproj = epsproj
        self.gvec_grid = gvec_grid  # NEW: G-vector grid from Phase 1
        
        # Constants
        self.RYTOEV = 13.6056980659  # Ry to eV conversion
        self.TPIBA = 2.0 * np.pi / reader.alat  # 2π/alat
        
        # Dimensions
        self.n2d = 0  # 2D dimension
        self.norb = 0  # total orbitals
        self.nocros = 0  # crossing orbitals
        self.noins = 0  # inside orbitals
        self.ntot = 0  # total dimension for eigenvalue problem
        
    def calculate(self, kpoint: np.ndarray, energy: float) -> Dict[str, Any]:
        """
        Calculate CBS for a given 2D k-point and energy.
        
        Parameters
        ----------
        kpoint : np.ndarray
            2D k-point in crystal coordinates [k1, k2]
        energy : float
            Energy relative to Fermi level in eV
            
        Returns
        -------
        dict
            Dictionary containing:
            - 'kvals': complex k-values
            - 'kvecs': corresponding eigenvectors
            - 'nchan': number of propagating channels
        """
        # Convert energy from eV to Ry and add Fermi energy
        eryd = energy / self.RYTOEV + self.reader.ef / self.RYTOEV
        
        # Setup the problem dimensions
        # For a simplified implementation, we'll use a model system
        # In a full implementation, this would involve reading/building
        # the Hamiltonian and overlap matrices from QE data
        
        self._setup_problem_dimensions(kpoint)
        
        # Build matrices for the generalized eigenvalue problem
        amat, bmat = self._build_matrices(kpoint, eryd)
        
        # Solve generalized eigenvalue problem: A v = k B v
        kvals, kvecs = self._solve_gep(amat, bmat)
        
        # Classify solutions into propagating and evanescent states
        nchan = self._count_channels(kvals)
        
        return {
            'kvals': kvals,
            'kvecs': kvecs,
            'nchan': nchan,
            'energy': energy,
            'kpoint': kpoint,
        }
        
    def _setup_problem_dimensions(self, kpoint: np.ndarray):
        """
        Setup problem dimensions based on system and k-point.
        
        NOTE: Now using G-vector grid for proper dimensions.
        
        For production use, additional components need to be added:
        - nocros depends on the atomic orbitals crossing the interface
        - noins depends on the interior orbitals
        
        These need to be extracted from:
        - Pseudopotential information
        - Atomic structure and orbital projections
        """
        # Try to use G-vector grid if available
        if hasattr(self, 'gvec_grid') and self.gvec_grid is not None:
            self.n2d = self.gvec_grid.ngper  # Number of 2D plane waves from G-vectors
            print(f"Using G-vector grid: n2d = {self.n2d}")
        else:
            # Fallback to simplified model
            self.n2d = 10
            print("Warning: Using simplified model for n2d (no G-vector grid)")
        
        # For demonstration, no orbitals yet
        self.nocros = 0  # crossing orbitals
        self.noins = 0  # inside orbitals
        self.norb = 2 * self.nocros + self.noins
        self.ntot = 2 * (self.n2d + self.nocros)
        
        print(f"Problem dimensions: n2d={self.n2d}, nocros={self.nocros}, ntot={self.ntot}")
        
    def _build_matrices(
        self,
        kpoint: np.ndarray,
        eryd: float
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Build A and B matrices for the generalized eigenvalue problem.
        
        This is a simplified model implementation. A full implementation
        would construct these from the Hamiltonian, overlap matrix, and
        potential from QE output.
        
        Parameters
        ----------
        kpoint : np.ndarray
            2D k-point
        eryd : float
            Energy in Ry
            
        Returns
        -------
        amat, bmat : np.ndarray
            Matrices for GEP A v = k B v
        """
        n = self.ntot
        
        # Initialize matrices
        amat = np.zeros((n, n), dtype=complex)
        bmat = np.zeros((n, n), dtype=complex)
        
        # Build simplified model matrices
        # In reality, these come from the tight-binding representation
        # of the Hamiltonian in the basis of Bloch functions
        
        # Model: simple tight-binding chain
        # H = -t (|n><n+1| + |n+1><n|) + epsilon_0 |n><n|
        
        t = 1.0  # hopping parameter
        eps0 = eryd  # on-site energy (epsilon_0)
        
        # Build block structure similar to PWCOND
        # The matrices have a structure: [[A11, A12], [A21, A22]]
        # where the blocks relate to left/right moving states
        
        for i in range(n//2):
            # Diagonal blocks
            amat[i, i] = eps0
            amat[i + n//2, i + n//2] = eps0
            
            # Off-diagonal blocks (coupling)
            if i < n//2 - 1:
                amat[i, i+1] = -t
                amat[i+1, i] = -t
                amat[i + n//2, i + n//2 + 1] = -t
                amat[i + n//2 + 1, i + n//2] = -t
        
        # Overlap matrix (identity for orthogonal basis)
        bmat = np.eye(n, dtype=complex)
        
        # Add k-point dependence
        kx, ky = kpoint[0], kpoint[1]
        phase = np.exp(1j * 2.0 * np.pi * kx)
        
        # Boundary conditions with k-point phase
        for i in range(n//2):
            if i == 0:
                amat[i, n//2-1] = -t * phase
                amat[n//2-1, i] = -t * np.conj(phase)
        
        return amat, bmat
        
    def _solve_gep(
        self,
        amat: np.ndarray,
        bmat: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Solve the generalized eigenvalue problem.
        
        Solves: A v = λ B v
        
        Parameters
        ----------
        amat : np.ndarray
            Matrix A
        bmat : np.ndarray
            Matrix B
            
        Returns
        -------
        kvals : np.ndarray
            Complex k-values (eigenvalues)
        kvecs : np.ndarray
            Eigenvectors
        """
        try:
            # Solve generalized eigenvalue problem
            eigenvalues, eigenvectors = linalg.eig(amat, bmat)
            
            # Extract k-values from eigenvalues
            # The relationship depends on the specific formulation
            # For the transfer matrix method: k is directly the eigenvalue
            kvals = eigenvalues
            
            # Sort by real part, then imaginary part
            idx = np.lexsort((np.imag(kvals), np.real(kvals)))
            kvals = kvals[idx]
            kvecs = eigenvectors[:, idx]
            
            return kvals, kvecs
            
        except linalg.LinAlgError as e:
            warnings.warn(f"Failed to solve eigenvalue problem: {e}")
            # Return empty results
            n = amat.shape[0]
            return np.zeros(n, dtype=complex), np.zeros((n, n), dtype=complex)
            
    def _count_channels(self, kvals: np.ndarray) -> int:
        """
        Count the number of propagating (transmission) channels.
        
        Propagating states have small imaginary parts (real k).
        
        Parameters
        ----------
        kvals : np.ndarray
            Complex k-values
            
        Returns
        -------
        int
            Number of propagating channels
        """
        # A state is propagating if Im(k) is small
        threshold = 1e-3
        
        nchan = 0
        for kval in kvals:
            if abs(np.imag(kval)) < threshold:
                nchan += 1
                
        # Count only right-moving states (positive group velocity)
        # In the simplest case, this is half of the propagating states
        nchan = nchan // 2
        
        return max(nchan, 1)  # At least one channel
