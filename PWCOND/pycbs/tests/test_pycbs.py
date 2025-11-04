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


def test_package_version():
    """Test that package has version info."""
    import pycbs
    assert hasattr(pycbs, '__version__')
    assert isinstance(pycbs.__version__, str)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
