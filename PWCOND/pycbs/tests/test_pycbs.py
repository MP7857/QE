"""
Tests for PyCBS package.

These tests verify the basic functionality of the CBS calculator.
"""

import pytest
import numpy as np
from pathlib import Path
import tempfile
import shutil

from pycbs import CBSCalculator
from pycbs.reader import QEDataReader
from pycbs.kgrid import KPointGrid
from pycbs.compbs import ComplexBandStructure


class TestKPointGrid:
    """Tests for k-point grid generation."""
    
    def test_grid_generation(self):
        """Test k-point grid is generated correctly."""
        # Simple 2x2 grid
        b1 = np.array([1.0, 0.0, 0.0])
        b2 = np.array([0.0, 1.0, 0.0])
        
        grid = KPointGrid(nk1=2, nk2=2, k1=0.0, k2=0.0, b1=b1, b2=b2)
        
        assert grid.nkpts == 4
        assert len(grid.xyk) == 4
        assert len(grid.wkpt) == 4
        
        # Check weights sum to 1
        assert np.isclose(np.sum(grid.wkpt), 1.0)
        
    def test_grid_with_shift(self):
        """Test k-point grid with shift."""
        b1 = np.array([1.0, 0.0, 0.0])
        b2 = np.array([0.0, 1.0, 0.0])
        
        grid = KPointGrid(nk1=2, nk2=2, k1=0.5, k2=0.5, b1=b1, b2=b2)
        
        assert grid.nkpts == 4
        # First point should be shifted
        kpt = grid.get_kpoint(0)
        assert np.isclose(kpt[0], 0.25)  # (0 + 0.5) / 2
        assert np.isclose(kpt[1], 0.25)
        

class TestComplexBandStructure:
    """Tests for CBS calculation."""
    
    def test_initialization(self):
        """Test CBS calculator initialization."""
        # Create a mock reader
        class MockReader:
            alat = 10.0
            ef = 5.0
            
        reader = MockReader()
        cbs = ComplexBandStructure(reader, ecut2d=0.0, ewind=1.0, epsproj=1e-6)
        
        assert cbs.ecut2d == 0.0
        assert cbs.ewind == 1.0
        assert cbs.epsproj == 1e-6
        
    def test_calculation_structure(self):
        """Test that calculation returns expected structure."""
        # Create a mock reader
        class MockReader:
            alat = 10.0
            ef = 5.0
            b1 = np.array([1.0, 0.0, 0.0])
            b2 = np.array([0.0, 1.0, 0.0])
            
        reader = MockReader()
        cbs = ComplexBandStructure(reader)
        
        kpoint = np.array([0.0, 0.0])
        energy = 1.0
        
        result = cbs.calculate(kpoint, energy)
        
        # Check result structure
        assert 'kvals' in result
        assert 'kvecs' in result
        assert 'nchan' in result
        assert 'energy' in result
        assert 'kpoint' in result
        
        # Check types
        assert isinstance(result['kvals'], np.ndarray)
        assert isinstance(result['nchan'], int)
        assert result['nchan'] >= 0


class TestCBSCalculator:
    """Tests for the main CBS calculator."""
    
    def test_initialization(self):
        """Test calculator initialization."""
        calc = CBSCalculator(
            outdir='/tmp/test',
            prefix='test',
            band_file='bands.test'
        )
        
        assert calc.outdir == Path('/tmp/test')
        assert calc.prefix == 'test'
        assert calc.band_file == 'bands.test'
        
    def test_energy_range_setting(self):
        """Test setting energy range."""
        calc = CBSCalculator(outdir='/tmp', prefix='test')
        
        calc.set_energy_range(energy0=10.0, denergy=-0.4, nenergy=5)
        
        assert calc.energy0 == 10.0
        assert calc.denergy == -0.4
        assert calc.nenergy == 5
        assert len(calc.earr) == 5
        
        # Check energy array values
        expected = np.array([10.0, 9.6, 9.2, 8.8, 8.4])
        np.testing.assert_allclose(calc.earr, expected)
        
    def test_kpoint_grid_setting(self):
        """Test setting k-point grid."""
        calc = CBSCalculator(outdir='/tmp', prefix='test')
        
        calc.set_kpoints_grid(nk1=3, nk2=3, k1=0.5, k2=0.5)
        
        assert calc.nk1 == 3
        assert calc.nk2 == 3
        assert calc.k1 == 0.5
        assert calc.k2 == 0.5
    
    def test_custom_kpoints_setting(self):
        """Test setting custom k-points."""
        calc = CBSCalculator(outdir='/tmp', prefix='test')
        
        # Define custom k-points
        kpoints = np.array([
            [0.0, 0.0],
            [0.5, 0.0],
            [0.0, 0.5],
            [0.5, 0.5]
        ])
        weights = np.array([0.25, 0.25, 0.25, 0.25])
        
        calc.set_custom_kpoints(kpoints, weights)
        
        assert calc.custom_kpoints is not None
        assert calc.custom_weights is not None
        assert calc.custom_kpoints.shape == (4, 2)
        assert calc.custom_weights.shape == (4,)
        np.testing.assert_allclose(calc.custom_weights, weights)
    
    def test_custom_kpoints_with_3d_input(self):
        """Test that 3D k-points are handled correctly."""
        calc = CBSCalculator(outdir='/tmp', prefix='test')
        
        # Input with 3 columns (kx, ky, kz) - kz should be ignored
        kpoints_3d = np.array([
            [0.0, 0.0, 0.25],
            [0.5, 0.0, 0.25],
            [0.0, 0.5, 0.25],
            [0.5, 0.5, 0.25]
        ])
        
        calc.set_custom_kpoints(kpoints_3d)
        
        # Should extract only first 2 columns
        assert calc.custom_kpoints.shape == (4, 2)
        expected = kpoints_3d[:, :2]
        np.testing.assert_allclose(calc.custom_kpoints, expected)
        
        # Weights should be auto-generated
        assert np.isclose(calc.custom_weights.sum(), 1.0)
        assert len(calc.custom_weights) == 4


class TestKPointGridCustom:
    """Tests for custom k-point functionality in KPointGrid."""
    
    def test_from_custom_kpoints(self):
        """Test creating grid from custom k-points."""
        kpoints = np.array([
            [0.0, 0.0],
            [0.5, 0.0],
            [0.0, 0.5]
        ])
        weights = np.array([0.33, 0.33, 0.34])
        
        b1 = np.array([1.0, 0.0, 0.0])
        b2 = np.array([0.0, 1.0, 0.0])
        
        grid = KPointGrid.from_custom_kpoints(kpoints, weights, b1, b2)
        
        assert grid.nkpts == 3
        assert grid.nk1 == 0  # Indicates custom k-points
        np.testing.assert_allclose(grid.xyk, kpoints)
        np.testing.assert_allclose(grid.wkpt, weights)


class TestReadKpointsFile:
    """Tests for reading k-points from file."""
    
    def test_read_kpoints_file(self):
        """Test reading k-points from PWCOND format file."""
        import tempfile
        from pycbs import read_kpoints_file
        
        # Create temporary file
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.dat') as f:
            f.write("4\n")
            f.write("0.0000  0.0000  0.25\n")
            f.write("0.5000  0.0000  0.25\n")
            f.write("0.0000  0.5000  0.25\n")
            f.write("0.5000  0.5000  0.25\n")
            fname = f.name
        
        try:
            kpoints, weights = read_kpoints_file(fname)
            
            # Check dimensions
            assert kpoints.shape == (4, 2)
            assert weights.shape == (4,)
            
            # Check values
            expected_kpts = np.array([
                [0.0, 0.0],
                [0.5, 0.0],
                [0.0, 0.5],
                [0.5, 0.5]
            ])
            np.testing.assert_allclose(kpoints, expected_kpts)
            
            # Check weights sum to 1
            assert np.isclose(weights.sum(), 1.0)
            
        finally:
            import os
            os.remove(fname)


def test_package_version():
    """Test that package has version info."""
    import pycbs
    assert hasattr(pycbs, '__version__')
    assert isinstance(pycbs.__version__, str)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
