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
    
    def test_read_charge_density_binary(self):
        """Test reading charge density from binary file."""
        # Create a temporary binary file with mock charge density
        with tempfile.TemporaryDirectory() as tmpdir:
            save_dir = os.path.join(tmpdir, 'test.save')
            os.makedirs(save_dir)
            
            # Create mock charge density file in QE binary format
            charge_file = os.path.join(save_dir, 'charge-density.dat')
            
            nr1, nr2, nr3, nspin = 10, 10, 10, 1
            rho_data = np.random.rand(nr1, nr2, nr3)
            
            with open(charge_file, 'wb') as f:
                # Write header
                header = np.array([nr1, nr2, nr3, nspin], dtype=np.int32)
                rec_len = header.nbytes
                np.array([rec_len], dtype=np.int32).tofile(f)
                header.tofile(f)
                np.array([rec_len], dtype=np.int32).tofile(f)
                
                # Write data
                rho_flat = rho_data.flatten(order='F')
                rec_len = rho_flat.nbytes
                np.array([rec_len], dtype=np.int32).tofile(f)
                rho_flat.tofile(f)
                np.array([rec_len], dtype=np.int32).tofile(f)
            
            # Read it back
            reader = LocalPotentialReader(tmpdir, 'test')
            result = reader.read_charge_density()
            
            assert result is not None
            rho, metadata = result
            assert metadata['nr1'] == nr1
            assert metadata['nr2'] == nr2
            assert metadata['nr3'] == nr3
            assert metadata['nspin'] == nspin
            assert rho.shape == (nr1, nr2, nr3)
            assert np.allclose(rho, rho_data)
    
    def test_read_charge_density_binary_spin_polarized(self):
        """Test reading spin-polarized charge density."""
        with tempfile.TemporaryDirectory() as tmpdir:
            save_dir = os.path.join(tmpdir, 'test.save')
            os.makedirs(save_dir)
            
            charge_file = os.path.join(save_dir, 'charge-density.dat')
            
            nr1, nr2, nr3, nspin = 8, 8, 8, 2
            rho_data = np.random.rand(nr1, nr2, nr3, nspin)
            
            with open(charge_file, 'wb') as f:
                # Write header
                header = np.array([nr1, nr2, nr3, nspin], dtype=np.int32)
                rec_len = header.nbytes
                np.array([rec_len], dtype=np.int32).tofile(f)
                header.tofile(f)
                np.array([rec_len], dtype=np.int32).tofile(f)
                
                # Write data
                rho_flat = rho_data.flatten(order='F')
                rec_len = rho_flat.nbytes
                np.array([rec_len], dtype=np.int32).tofile(f)
                rho_flat.tofile(f)
                np.array([rec_len], dtype=np.int32).tofile(f)
            
            # Read it back
            reader = LocalPotentialReader(tmpdir, 'test')
            result = reader.read_charge_density()
            
            assert result is not None
            rho, metadata = result
            assert metadata['nspin'] == nspin
            assert rho.shape == (nr1, nr2, nr3, nspin)
            assert np.allclose(rho, rho_data)
    
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
