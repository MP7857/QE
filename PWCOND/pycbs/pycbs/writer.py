"""
Output writer for CBS calculation results.

This module handles writing CBS results to files in formats compatible
with PWCOND output.
"""

import numpy as np
from pathlib import Path
from typing import Optional


class CBSWriter:
    """
    Writer for CBS calculation results.
    
    Writes output in the same format as PWCOND for compatibility.
    
    Parameters
    ----------
    band_file : str
        Base name for output files (without extension)
    cutplot : float, optional
        Cutoff for plotting evanescent states (default: 2.0)
    """
    
    def __init__(self, band_file: str, cutplot: float = 2.0):
        self.band_file = band_file
        self.cutplot = cutplot
        
        # File handles
        self.f_re = None
        self.f_im = None
        self.f_co_re = None
        self.f_co_im = None
        self.f_3d = None
        
        self.is_open = False
        
    def open(self):
        """Open all output files and write headers."""
        try:
            self.f_re = open(f"{self.band_file}.re", 'w')
            self.f_im = open(f"{self.band_file}.im", 'w')
            self.f_co_re = open(f"{self.band_file}.co_re", 'w')
            self.f_co_im = open(f"{self.band_file}.co_im", 'w')
            self.f_3d = open(f"{self.band_file}.3d", 'w')
            
            # Write headers
            self.f_re.write("# Re(k), E-Ef\n")
            self.f_im.write("# Im(k), E-Ef\n")
            self.f_co_re.write("# Re(k), E-Ef\n")
            self.f_co_im.write("# Im(k), E-Ef\n")
            self.f_3d.write("# Re(k), Im(k), E-Ef\n")
            
            self.is_open = True
        except IOError as e:
            # Clean up any opened files on error
            self.close()
            raise IOError(f"Failed to open output files: {e}")
        
    def write_point(
        self,
        ik: int,
        ien: int,
        kvals: np.ndarray,
        nchan: int,
        energy: float,
        nkpts: int,
        nenergy: int
    ):
        """
        Write CBS data for a single (k-point, energy) combination.
        
        Parameters
        ----------
        ik : int
            K-point index
        ien : int
            Energy index
        kvals : np.ndarray
            Complex k-values
        nchan : int
            Number of propagating channels
        energy : float
            Energy value in eV
        nkpts : int
            Total number of k-points
        nenergy : int
            Total number of energy points
        """
        if not self.is_open:
            raise RuntimeError("Files not open. Call open() first.")
            
        # Write k-point marker at start of each k-point
        if ien == 0:
            for f in [self.f_re, self.f_im, self.f_co_re, self.f_co_im, self.f_3d]:
                f.write(f"# k-point {ik+1:5d}\n")
        
        # Sort k-values for output
        n = len(kvals)
        # Number of states in the left-moving block
        # Assumes symmetric structure: half states move right, half move left
        n_states_left = n // 2
        
        # Write propagating states (real k)
        for i in range(nchan):
            if i < len(kvals):
                # Right-moving state
                kval = kvals[i]
                self._write_propagating_state(kval, energy)
                
                # Left-moving state (if exists)
                if n_states_left + i < len(kvals):
                    kval = kvals[n_states_left + i]
                    self._write_propagating_state(kval, energy)
        
        # Write evanescent states (complex k)
        for k in range(nchan, n_states_left):
            if k < len(kvals):
                # Forward evanescent
                kval = kvals[k]
                self._write_evanescent_state(kval, energy)
                
                # Backward evanescent
                if n_states_left + k < len(kvals):
                    kval = kvals[n_states_left + k]
                    self._write_evanescent_state(kval, energy)
        
        # Flush buffers periodically
        if (ik * nenergy + ien) % 10 == 0:
            self._flush()
            
    def _write_propagating_state(self, kval: complex, energy: float):
        """Write a propagating state to output files."""
        re_k = np.real(kval)
        im_k = np.imag(kval)
        
        # Real part file
        self.f_re.write(f"{re_k:10.4f}{energy:10.4f}\n")
        
        # 3D file
        self.f_3d.write(f"{re_k:10.4f}{im_k:10.4f}{energy:10.4f}\n")
        
    def _write_evanescent_state(self, kval: complex, energy: float):
        """Write an evanescent state to output files."""
        re_k = np.real(kval)
        im_k = np.imag(kval)
        dim = abs(im_k)
        
        # Only write if within cutoff
        if dim <= self.cutplot:
            # Classification based on real part
            if abs(re_k) <= 1e-3:
                # Pure imaginary k
                self.f_im.write(f"{-dim:10.4f}{energy:10.4f}\n")
            elif abs(re_k - 0.5) <= 1e-3 or abs(re_k + 0.5) <= 1e-3:
                # Near zone boundary
                self.f_im.write(f"{0.5 + dim:10.4f}{energy:10.4f}\n")
            else:
                # Complex k with significant real part
                self.f_co_re.write(f"{re_k:10.4f}{energy:10.4f}\n")
                self.f_co_im.write(f"{-0.5 - dim:10.4f}{energy:10.4f}\n")
            
            # Always write to 3D file
            self.f_3d.write(f"{re_k:10.4f}{im_k:10.4f}{energy:10.4f}\n")
    
    def _flush(self):
        """Flush all output buffers."""
        if self.is_open:
            for f in [self.f_re, self.f_im, self.f_co_re, self.f_co_im, self.f_3d]:
                if f:
                    f.flush()
        
    def close(self):
        """Close all output files."""
        if self.is_open:
            for f in [self.f_re, self.f_im, self.f_co_re, self.f_co_im, self.f_3d]:
                if f:
                    f.close()
            self.is_open = False
            
    def __enter__(self):
        """Context manager entry."""
        self.open()
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.close()
