"""
Hamiltonian construction module for CBS calculations.

This module implements the Hamiltonian and overlap matrix construction
following PWCOND's approach, integrating with the G-vector grid.

This is Phase 1.5 / Phase 2 of the implementation roadmap.
"""

import numpy as np
from typing import Tuple, Optional, Dict, Any
from .wfc_reader import GVectorGrid


class HamiltonianBuilder:
    """
    Constructs Hamiltonian and overlap matrices for CBS calculations.
    
    Following PWCOND's structure, this builds the matrices needed for
    the generalized eigenvalue problem A v = k B v.
    
    Parameters
    ----------
    gvec_grid : GVectorGrid
        2D G-vector grid for the plane wave basis
    energy : float
        Energy in Rydberg
    """
    
    def __init__(self, gvec_grid: GVectorGrid, energy: float):
        self.gvec_grid = gvec_grid
        self.energy = energy
        self.n2d = gvec_grid.ngper  # Number of 2D plane waves
        
    def build_kinetic_matrix(self, kz: complex) -> np.ndarray:
        """
        Build kinetic energy matrix in the 2D plane wave basis.
        
        T_ij = δ_ij * [(G_i + k_perp)² + k_z²] / 2
        
        where k_perp is the 2D k-point and k_z is the propagation direction.
        
        Parameters
        ----------
        kz : complex
            Complex wave vector in z-direction
            
        Returns
        -------
        np.ndarray
            Kinetic energy matrix (n2d, n2d)
        """
        n = self.n2d
        T = np.zeros((n, n), dtype=complex)
        
        # Get reciprocal lattice and k-point
        bg = self.gvec_grid.bg
        kpoint = self.gvec_grid.kpoint
        
        # Construct kinetic energy: T = (k + G)² / 2
        for i in range(n):
            # Get G-vector
            gx, gy = self.gvec_grid.gper[0, i], self.gvec_grid.gper[1, i]
            
            # k_parallel = k_perp + G_perp (in 2D)
            kx = kpoint[0] + gx
            ky = kpoint[1] + gy
            
            # Convert to Cartesian coordinates using reciprocal lattice
            kx_cart = kx * bg[0, 0] + ky * bg[0, 1]
            ky_cart = kx * bg[1, 0] + ky * bg[1, 1]
            
            # Total kinetic energy: |k_parallel|² + |k_z|²
            k_perp_sq = kx_cart**2 + ky_cart**2
            k_z_sq = kz**2
            
            # T = (|k|²) / 2 in Hartree units, convert to Ry
            T[i, i] = (k_perp_sq + k_z_sq)  # Already in Ry units
            
        return T
    
    def build_simple_hamiltonian(
        self,
        kz: complex,
        potential: Optional[np.ndarray] = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Build simplified Hamiltonian and overlap matrices.
        
        H = T + V
        S = I (orthogonal basis for now)
        
        This is a simplified version. Full implementation would include:
        - Local pseudopotential from QE charge density
        - Non-local pseudopotential projections
        - Proper overlap matrix for ultrasoft pseudopotentials
        
        Parameters
        ----------
        kz : complex
            Complex wave vector in z-direction
        potential : np.ndarray, optional
            Local potential matrix (if None, V = 0)
            
        Returns
        -------
        H : np.ndarray
            Hamiltonian matrix (n2d, n2d)
        S : np.ndarray
            Overlap matrix (n2d, n2d)
        """
        # Build kinetic energy
        T = self.build_kinetic_matrix(kz)
        
        # Add potential (if provided)
        if potential is not None:
            H = T + potential
        else:
            # Free-electron Hamiltonian for demonstration
            H = T
        
        # Overlap matrix (identity for orthogonal basis)
        S = np.eye(self.n2d, dtype=complex)
        
        return H, S
    
    def build_cbs_matrices(
        self,
        kz: complex,
        nocros: int = 0,
        noins: int = 0
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Build A and B matrices for CBS generalized eigenvalue problem.
        
        Following PWCOND's compbs_2.f90, this constructs the matrices:
        A v = k B v
        
        The structure is:
        - Upper left: 2D Hamiltonian blocks
        - Orbital blocks for crossing/interior orbitals
        - Matching conditions
        
        Parameters
        ----------
        kz : complex
            Complex wave vector in z-direction
        nocros : int
            Number of crossing orbitals
        noins : int
            Number of interior orbitals
            
        Returns
        -------
        amat : np.ndarray
            A matrix for GEP
        bmat : np.ndarray
            B matrix for GEP
        """
        # For now, simplified version without orbitals
        # Just use plane waves
        H, S = self.build_simple_hamiltonian(kz)
        
        # Construct A and B matrices following PWCOND structure
        # For simplified case: A = H - E*S, B = S
        amat = H - self.energy * S
        bmat = S
        
        return amat, bmat


def build_test_hamiltonian(
    ngper: int,
    energy: float,
    kpoint: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Build a test Hamiltonian for demonstration.
    
    This creates a simple 1D tight-binding model for testing
    the CBS algorithm before full integration.
    
    Parameters
    ----------
    ngper : int
        Number of basis functions
    energy : float
        Energy in Ry
    kpoint : np.ndarray
        2D k-point
        
    Returns
    -------
    amat, bmat : np.ndarray
        Test matrices for GEP
    """
    n = 2 * ngper  # ntot = 2 * n2d for no orbitals
    
    # Simple tight-binding Hamiltonian: H = -t Σ(|n⟩⟨n+1| + h.c.) + ε₀ I
    t = 1.0  # Hopping parameter
    eps0 = 0.0  # On-site energy
    
    H = np.diag([eps0] * n).astype(complex)
    for i in range(n - 1):
        H[i, i+1] = -t
        H[i+1, i] = -t
    
    # Overlap matrix (identity for orthogonal basis)
    S = np.eye(n, dtype=complex)
    
    # A = H - E*S, B = S
    amat = H - energy * S
    bmat = S
    
    return amat, bmat
