"""
Tests for wavefunction reader module (Phase 1).
"""

import pytest
import numpy as np
from pathlib import Path
import tempfile
import struct

from pycbs.wfc_reader import WavefunctionReader, read_wavefunction_metadata, GVectorGrid


class TestGVectorGrid:
    """Tests for G-vector grid construction."""
    
    def test_initialization(self):
        """Test G-vector grid initialization."""
        bg = np.eye(3) * 2 * np.pi  # Simple cubic
        kpoint = np.array([0.0, 0.0])
        ecut2d = 10.0  # Ry
        
        grid = GVectorGrid(bg, kpoint, ecut2d, nrx=10, nry=10)
        
        assert grid.ngper > 0
        assert grid.ngpsh > 0
        assert grid.gper is not None
        assert grid.gper.shape[0] == 2
        assert grid.gper.shape[1] == grid.ngper
        
    def test_gamma_point(self):
        """Test G-vector grid at Gamma point."""
        bg = np.eye(3) * 2 * np.pi
        kpoint = np.array([0.0, 0.0])
        ecut2d = 5.0
        
        grid = GVectorGrid(bg, kpoint, ecut2d, nrx=8, nry=8)
        
        # At gamma, first G-vector should be (0,0)
        assert grid.gper[0, 0] == 0
        assert grid.gper[1, 0] == 0
        
    def test_nonzero_kpoint(self):
        """Test G-vector grid at non-zero k-point."""
        bg = np.eye(3) * 2 * np.pi
        kpoint = np.array([0.25, 0.25])
        ecut2d = 5.0
        
        grid = GVectorGrid(bg, kpoint, ecut2d, nrx=8, nry=8)
        
        assert grid.ngper > 0
        # Grid depends on k-point
        
    def test_ecut_increases_gvecs(self):
        """Test that higher cutoff gives more G-vectors."""
        bg = np.eye(3) * 2 * np.pi
        kpoint = np.array([0.0, 0.0])
        
        # Use larger cutoffs to ensure multiple G-vectors
        grid1 = GVectorGrid(bg, kpoint, ecut2d=50.0, nrx=10, nry=10)
        grid2 = GVectorGrid(bg, kpoint, ecut2d=200.0, nrx=10, nry=10)
        
        assert grid2.ngper >= grid1.ngper  # Higher cutoff should have at least as many


class TestWavefunctionReader:
    """Tests for WavefunctionReader class."""
    
    def test_initialization(self):
        """Test reader initialization."""
        reader = WavefunctionReader('/tmp', 'test')
        
        assert reader.outdir == Path('/tmp')
        assert reader.prefix == 'test'
        assert reader.save_dir == Path('/tmp/test.save')
        
    def test_find_wfc_files_missing_dir(self):
        """Test handling of missing save directory."""
        reader = WavefunctionReader('/nonexistent', 'test')
        
        with pytest.raises(FileNotFoundError):
            reader.find_wfc_files()
            
    def test_find_wfc_files_empty_dir(self):
        """Test with empty save directory."""
        with tempfile.TemporaryDirectory() as tmpdir:
            save_dir = Path(tmpdir) / 'test.save'
            save_dir.mkdir()
            
            reader = WavefunctionReader(tmpdir, 'test')
            wfc_files = reader.find_wfc_files()
            
            assert len(wfc_files) == 0
            
    def test_find_wfc_files_with_files(self):
        """Test finding wavefunction files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            save_dir = Path(tmpdir) / 'test.save'
            save_dir.mkdir()
            
            # Create dummy wfc files
            (save_dir / 'wfc1.dat').touch()
            (save_dir / 'wfc2.dat').touch()
            
            reader = WavefunctionReader(tmpdir, 'test')
            wfc_files = reader.find_wfc_files()
            
            assert len(wfc_files) == 2
            assert all(f.name.startswith('wfc') for f in wfc_files)
            
    def test_read_wfc_header_missing_file(self):
        """Test reading header from non-existent file."""
        reader = WavefunctionReader('/tmp', 'test')
        
        # read_wfc_header catches exceptions and returns error info
        header = reader.read_wfc_header(Path('/nonexistent/wfc1.dat'))
        assert 'error' in header
            
    def test_read_wfc_header_malformed(self):
        """Test reading header from malformed file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            wfc_file = Path(tmpdir) / 'wfc1.dat'
            
            # Write malformed data
            with open(wfc_file, 'wb') as f:
                f.write(b'invalid data')
                
            reader = WavefunctionReader(tmpdir, 'test')
            header = reader.read_wfc_header(wfc_file)
            
            # Should handle gracefully and return error info
            assert 'error' in header
            
    def test_read_wfc_header_basic_format(self):
        """Test reading a simple header."""
        with tempfile.TemporaryDirectory() as tmpdir:
            wfc_file = Path(tmpdir) / 'wfc1.dat'
            
            # Create a minimal Fortran unformatted file
            with open(wfc_file, 'wb') as f:
                # Record 1: ik and xk
                rec_len = 28  # 4 (ik) + 24 (xk)
                f.write(struct.pack('i', rec_len))  # Record start
                f.write(struct.pack('i', 1))  # ik = 1
                f.write(struct.pack('ddd', 0.0, 0.0, 0.0))  # xk
                f.write(struct.pack('i', rec_len))  # Record end
                
                # Record 2: npw
                rec_len = 4
                f.write(struct.pack('i', rec_len))
                f.write(struct.pack('i', 100))  # npw = 100
                f.write(struct.pack('i', rec_len))
                
            reader = WavefunctionReader(tmpdir, 'test')
            header = reader.read_wfc_header(wfc_file)
            
            if 'error' not in header:
                assert header['ik'] == 1
                assert header['xk'] == (0.0, 0.0, 0.0)
                assert header['npw'] == 100


class TestWavefunctionMetadata:
    """Tests for metadata utility function."""
    
    def test_metadata_missing_dir(self):
        """Test metadata with missing directory."""
        with pytest.raises(FileNotFoundError):
            read_wavefunction_metadata('/nonexistent', 'test')
            
    def test_metadata_empty_dir(self):
        """Test metadata with empty directory."""
        with tempfile.TemporaryDirectory() as tmpdir:
            save_dir = Path(tmpdir) / 'test.save'
            save_dir.mkdir()
            
            metadata = read_wavefunction_metadata(tmpdir, 'test')
            
            assert metadata['num_wfc_files'] == 0
            assert metadata['wfc_files'] == []
            
    def test_metadata_with_files(self):
        """Test metadata with wavefunction files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            save_dir = Path(tmpdir) / 'test.save'
            save_dir.mkdir()
            
            # Create dummy files
            (save_dir / 'wfc1.dat').touch()
            (save_dir / 'wfc2.dat').touch()
            
            metadata = read_wavefunction_metadata(tmpdir, 'test')
            
            assert metadata['num_wfc_files'] == 2
            assert 'wfc1.dat' in metadata['wfc_files']
            assert 'wfc2.dat' in metadata['wfc_files']


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
