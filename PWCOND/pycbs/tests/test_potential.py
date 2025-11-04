"""
Tests for local potential module.
"""

import pytest
import numpy as np
import os
import tempfile
from pycbs.potential import (
    LocalPotentialReader,
    PseudopotentialManager,
    integrate_local_potential
)
from pycbs.wfc_reader import GVectorGrid


class TestLocalPotentialReader:
    """Test LocalPotentialReader class."""
    
    def test_init(self):
        """Test initialization."""
        reader = LocalPotentialReader('/tmp', 'test')
        assert reader.outdir == '/tmp'
        assert reader.prefix == 'test'
        assert reader.save_dir == '/tmp/test.save'
    
    def test_read_charge_density_no_file(self):
        """Test reading charge density when file doesn't exist."""
        reader = LocalPotentialReader('/tmp/nonexistent', 'test')
        rho = reader.read_charge_density()
        assert rho is None
    
    def test_construct_local_potential_2d_placeholder(self):
        """Test local potential construction with placeholder."""
        # Create G-vector grid
        bg = np.eye(3) * 2 * np.pi
        kpoint = np.array([0.0, 0.0])
        grid = GVectorGrid(bg, kpoint, ecut2d=10.0, nrx=8, nry=8)
        
        reader = LocalPotentialReader('/tmp', 'test')
        v_loc = reader.construct_local_potential_2d(None, grid)
        
        # Should return zero matrix for now
        assert v_loc.shape == (grid.ngper, grid.ngper)
        assert np.allclose(v_loc, 0.0)
    
    def test_construct_local_potential_2d_with_density(self):
        """Test local potential construction with dummy density."""
        # Create G-vector grid
        bg = np.eye(3) * 2 * np.pi
        kpoint = np.array([0.0, 0.0])
        grid = GVectorGrid(bg, kpoint, ecut2d=10.0, nrx=8, nry=8)
        
        # Dummy charge density
        rho = np.random.rand(8, 8, 8)
        
        reader = LocalPotentialReader('/tmp', 'test')
        v_loc = reader.construct_local_potential_2d(rho, grid)
        
        # Should return zero matrix for now (not implemented)
        assert v_loc.shape == (grid.ngper, grid.ngper)
        assert np.allclose(v_loc, 0.0)


class TestPseudopotentialManager:
    """Test PseudopotentialManager class."""
    
    def test_init_empty(self):
        """Test initialization without PP files."""
        manager = PseudopotentialManager()
        assert manager.pp_files == []
        assert manager.pp_data == {}
    
    def test_init_with_files(self):
        """Test initialization with PP files."""
        files = ['Al.UPF', 'Si.UPF']
        manager = PseudopotentialManager(files)
        assert manager.pp_files == files
    
    def test_get_local_potential_not_loaded(self):
        """Test getting local potential when not loaded."""
        manager = PseudopotentialManager()
        vloc = manager.get_local_potential('Al')
        assert vloc is None
    
    def test_get_nonlocal_projectors_not_loaded(self):
        """Test getting non-local projectors when not loaded."""
        manager = PseudopotentialManager()
        projectors = manager.get_nonlocal_projectors('Al')
        assert projectors is None


class TestIntegration:
    """Test integration functions."""
    
    def test_integrate_local_potential(self):
        """Test integrate_local_potential function."""
        from pycbs.hamiltonian import HamiltonianBuilder
        
        # Create G-vector grid
        bg = np.eye(3) * 2 * np.pi
        kpoint = np.array([0.0, 0.0])
        grid = GVectorGrid(bg, kpoint, ecut2d=10.0, nrx=8, nry=8)
        
        # Create HamiltonianBuilder
        builder = HamiltonianBuilder(grid, energy=10.0)
        
        # Create potential reader
        reader = LocalPotentialReader('/tmp', 'test')
        
        # Integrate potential
        v_loc = integrate_local_potential(builder, reader, grid)
        
        # Should return zero matrix for now
        assert v_loc.shape == (grid.ngper, grid.ngper)
        assert np.allclose(v_loc, 0.0)
