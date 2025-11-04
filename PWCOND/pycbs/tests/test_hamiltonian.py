"""
Tests for Hamiltonian construction module (Phase 1.5/2).
"""

import pytest
import numpy as np

from pycbs.wfc_reader import GVectorGrid
from pycbs.hamiltonian import HamiltonianBuilder, build_test_hamiltonian


class TestHamiltonianBuilder:
    """Tests for Hamiltonian builder."""
    
    def test_initialization(self):
        """Test Hamiltonian builder initialization."""
        # Create a simple G-vector grid
        bg = np.eye(3) * 2 * np.pi
        kpoint = np.array([0.0, 0.0])
        grid = GVectorGrid(bg, kpoint, ecut2d=50.0, nrx=16, nry=16)
        
        builder = HamiltonianBuilder(grid, energy=10.0)
        
        assert builder.n2d == grid.ngper
        assert builder.energy == 10.0
        
    def test_kinetic_matrix(self):
        """Test kinetic energy matrix construction."""
        bg = np.eye(3) * 2 * np.pi
        kpoint = np.array([0.0, 0.0])
        grid = GVectorGrid(bg, kpoint, ecut2d=50.0, nrx=16, nry=16)
        
        builder = HamiltonianBuilder(grid, energy=10.0)
        
        # Build kinetic matrix for kz = 0
        T = builder.build_kinetic_matrix(kz=0.0)
        
        # Check dimensions
        assert T.shape == (grid.ngper, grid.ngper)
        
        # Check it's diagonal (kinetic energy is diagonal in plane wave basis)
        off_diag = T - np.diag(np.diag(T))
        assert np.allclose(off_diag, 0.0)
        
        # Check all eigenvalues are real and positive
        eigs = np.linalg.eigvalsh(T)
        assert np.all(eigs >= 0)
        
    def test_kinetic_matrix_complex_kz(self):
        """Test kinetic energy matrix with complex kz."""
        bg = np.eye(3) * 2 * np.pi
        kpoint = np.array([0.0, 0.0])
        grid = GVectorGrid(bg, kpoint, ecut2d=50.0, nrx=16, nry=16)
        
        builder = HamiltonianBuilder(grid, energy=10.0)
        
        # Build kinetic matrix for complex kz
        T = builder.build_kinetic_matrix(kz=0.5 + 0.1j)
        
        # Check dimensions
        assert T.shape == (grid.ngper, grid.ngper)
        
        # Still diagonal
        off_diag = T - np.diag(np.diag(T))
        assert np.allclose(off_diag, 0.0)
        
    def test_simple_hamiltonian(self):
        """Test simple Hamiltonian construction."""
        bg = np.eye(3) * 2 * np.pi
        kpoint = np.array([0.0, 0.0])
        grid = GVectorGrid(bg, kpoint, ecut2d=50.0, nrx=16, nry=16)
        
        builder = HamiltonianBuilder(grid, energy=10.0)
        
        H, S = builder.build_simple_hamiltonian(kz=0.0)
        
        # Check dimensions
        assert H.shape == (grid.ngper, grid.ngper)
        assert S.shape == (grid.ngper, grid.ngper)
        
        # Overlap should be identity for orthogonal basis
        assert np.allclose(S, np.eye(grid.ngper))
        
        # H should be Hermitian
        assert np.allclose(H, H.conj().T)
        
    def test_cbs_matrices(self):
        """Test CBS matrix construction."""
        bg = np.eye(3) * 2 * np.pi
        kpoint = np.array([0.0, 0.0])
        grid = GVectorGrid(bg, kpoint, ecut2d=50.0, nrx=16, nry=16)
        
        builder = HamiltonianBuilder(grid, energy=10.0)
        
        amat, bmat = builder.build_cbs_matrices(kz=0.0)
        
        # Check dimensions
        assert amat.shape == (grid.ngper, grid.ngper)
        assert bmat.shape == (grid.ngper, grid.ngper)
        
        # B should be identity for orthogonal basis
        assert np.allclose(bmat, np.eye(grid.ngper))
        
    def test_nonzero_kpoint(self):
        """Test with non-zero k-point."""
        bg = np.eye(3) * 2 * np.pi
        kpoint = np.array([0.25, 0.25])
        grid = GVectorGrid(bg, kpoint, ecut2d=50.0, nrx=16, nry=16)
        
        builder = HamiltonianBuilder(grid, energy=10.0)
        
        T = builder.build_kinetic_matrix(kz=0.0)
        
        # Check dimensions
        assert T.shape == (grid.ngper, grid.ngper)
        
        # Should still be diagonal
        off_diag = T - np.diag(np.diag(T))
        assert np.allclose(off_diag, 0.0)


class TestBuildTestHamiltonian:
    """Tests for test Hamiltonian builder."""
    
    def test_basic_construction(self):
        """Test basic test Hamiltonian construction."""
        ngper = 5
        energy = 10.0
        kpoint = np.array([0.0, 0.0])
        
        amat, bmat = build_test_hamiltonian(ngper, energy, kpoint)
        
        # Check dimensions
        n = 2 * ngper
        assert amat.shape == (n, n)
        assert bmat.shape == (n, n)
        
        # B should be identity
        assert np.allclose(bmat, np.eye(n))
        
    def test_hermitian(self):
        """Test that test Hamiltonian is Hermitian."""
        ngper = 5
        energy = 10.0
        kpoint = np.array([0.0, 0.0])
        
        amat, bmat = build_test_hamiltonian(ngper, energy, kpoint)
        
        # A + E*I should be Hermitian (since A = H - E*I)
        H_reconstructed = amat + energy * bmat
        assert np.allclose(H_reconstructed, H_reconstructed.conj().T)
