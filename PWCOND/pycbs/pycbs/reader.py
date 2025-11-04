"""
Module for reading Quantum ESPRESSO data files.

This module provides functionality to read QE XML output files and extract
the necessary data for CBS calculations.
"""

import numpy as np
from pathlib import Path
import xml.etree.ElementTree as ET
from typing import Optional, Tuple
import struct


class QEDataReader:
    """
    Reader for Quantum ESPRESSO data files.
    
    This class reads QE output in XML format and extracts structural,
    electronic, and potential information needed for CBS calculations.
    
    Parameters
    ----------
    outdir : Path or str
        Directory containing QE output files
    prefix : str
        Prefix of the QE calculation
    """
    
    def __init__(self, outdir, prefix):
        self.outdir = Path(outdir)
        self.prefix = prefix
        self.data_file_schema = self.outdir / f"{prefix}.save" / "data-file-schema.xml"
        
        # System properties
        self.nat = 0  # number of atoms
        self.ntyp = 0  # number of atom types
        self.nbnd = 0  # number of bands
        self.nk = 0  # number of k-points
        self.nspin = 1  # number of spin components
        self.npol = 1  # number of spin polarizations
        
        # Lattice vectors and reciprocal lattice
        self.alat = 0.0  # lattice parameter in Bohr
        self.a1 = np.zeros(3)
        self.a2 = np.zeros(3)
        self.a3 = np.zeros(3)
        self.b1 = np.zeros(3)
        self.b2 = np.zeros(3)
        self.b3 = np.zeros(3)
        self.omega = 0.0  # cell volume
        
        # Atomic positions and types
        self.tau = []  # atomic positions
        self.ityp = []  # atomic types
        self.atm = []  # atom names
        
        # Electronic structure
        self.ef = 0.0  # Fermi energy in eV
        self.nelec = 0.0  # number of electrons
        self.ecutwfc = 0.0  # wavefunction cutoff in Ry
        self.ecutrho = 0.0  # charge density cutoff in Ry
        
        # Wavefunction grid
        self.nr1 = 0
        self.nr2 = 0
        self.nr3 = 0
        
        # K-points
        self.xk = []  # k-point coordinates
        self.wk = []  # k-point weights
        
        # Additional flags
        self.noncolin = False
        self.lspinorb = False
        
    def read(self):
        """Read all necessary data from QE output files."""
        if not self.data_file_schema.exists():
            # Try old XML format
            old_xml = self.outdir / f"{self.prefix}.xml"
            if old_xml.exists():
                self._read_old_xml(old_xml)
            else:
                raise FileNotFoundError(
                    f"Could not find QE data file: {self.data_file_schema} or {old_xml}"
                )
        else:
            self._read_xml_schema()
            
    def _read_xml_schema(self):
        """Read data from new XML schema format."""
        tree = ET.parse(self.data_file_schema)
        root = tree.getroot()
        
        # Read lattice information
        cell = root.find(".//atomic_structure")
        if cell is not None:
            self.nat = int(cell.get('nat', 0))
            alat_attr = cell.get('alat')
            if alat_attr:
                self.alat = float(alat_attr)
            
            # Read lattice vectors
            cell_vecs = cell.find('cell')
            if cell_vecs is not None:
                a1_elem = cell_vecs.find('a1')
                a2_elem = cell_vecs.find('a2')
                a3_elem = cell_vecs.find('a3')
                if a1_elem is not None and a1_elem.text:
                    self.a1 = np.array([float(x) for x in a1_elem.text.split()])
                if a2_elem is not None and a2_elem.text:
                    self.a2 = np.array([float(x) for x in a2_elem.text.split()])
                if a3_elem is not None and a3_elem.text:
                    self.a3 = np.array([float(x) for x in a3_elem.text.split()])
            
            # Read atomic positions
            atoms = cell.findall('.//atom')
            self.tau = []
            self.atm = []
            for atom in atoms:
                name = atom.get('name', '')
                pos_text = atom.text
                if pos_text:
                    pos = np.array([float(x) for x in pos_text.split()])
                    self.tau.append(pos)
                    self.atm.append(name)
        
        # Read basis information
        basis = root.find(".//basis")
        if basis is not None:
            # FFT grid
            fft_grid = basis.find(".//fft_grid")
            if fft_grid is not None:
                self.nr1 = int(fft_grid.get('nr1', 0))
                self.nr2 = int(fft_grid.get('nr2', 0))
                self.nr3 = int(fft_grid.get('nr3', 0))
        
        # Read band structure information
        band_structure = root.find(".//band_structure")
        if band_structure is not None:
            # Number of spin components
            nspin_elem = band_structure.find('noncolin')
            if nspin_elem is not None and nspin_elem.text:
                self.noncolin = nspin_elem.text.strip().lower() == 'true'
                if self.noncolin:
                    self.nspin = 4
                    self.npol = 2
            
            # Spin-orbit
            spinorb_elem = band_structure.find('spinorbit')
            if spinorb_elem is not None and spinorb_elem.text:
                self.lspinorb = spinorb_elem.text.strip().lower() == 'true'
            
            # Number of bands
            nbnd_elem = band_structure.find('nbnd')
            if nbnd_elem is not None and nbnd_elem.text:
                self.nbnd = int(nbnd_elem.text)
            
            # Number of electrons
            nelec_elem = band_structure.find('nelec')
            if nelec_elem is not None and nelec_elem.text:
                self.nelec = float(nelec_elem.text)
            
            # Fermi energy
            ef_elem = band_structure.find('fermi_energy')
            if ef_elem is not None and ef_elem.text:
                self.ef = float(ef_elem.text)
            
            # K-points
            nks_elem = band_structure.find('nks')
            if nks_elem is not None and nks_elem.text:
                self.nk = int(nks_elem.text)
            
            # Read individual k-points
            ks_energies = band_structure.findall('.//ks_energies')
            self.xk = []
            self.wk = []
            for ks in ks_energies:
                k_point = ks.find('k_point')
                if k_point is not None and k_point.text:
                    kpt = np.array([float(x) for x in k_point.text.split()])
                    self.xk.append(kpt)
                    weight = float(k_point.get('weight', 1.0))
                    self.wk.append(weight)
        
        # Calculate reciprocal lattice vectors
        self._calculate_reciprocal_lattice()
        
    def _read_old_xml(self, xml_file):
        """Read data from old XML format (QE < 6.5)."""
        # Simplified reader for old format
        # This is a placeholder - full implementation would parse the old XML
        raise NotImplementedError(
            "Old XML format reading not yet implemented. "
            "Please use QE version 6.5 or later."
        )
        
    def _calculate_reciprocal_lattice(self):
        """Calculate reciprocal lattice vectors from direct lattice."""
        # Volume of unit cell
        self.omega = np.dot(self.a1, np.cross(self.a2, self.a3))
        
        if abs(self.omega) < 1e-10:
            raise ValueError("Cell volume is too small or lattice vectors are coplanar")
        
        # Reciprocal lattice vectors (in units of 2π/alat)
        tpi = 2.0 * np.pi
        self.b1 = tpi * np.cross(self.a2, self.a3) / self.omega
        self.b2 = tpi * np.cross(self.a3, self.a1) / self.omega
        self.b3 = tpi * np.cross(self.a1, self.a2) / self.omega
        
    def get_lattice_parameter(self) -> float:
        """Get the lattice parameter in Bohr."""
        return self.alat
        
    def get_reciprocal_vectors(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Get reciprocal lattice vectors."""
        return self.b1, self.b2, self.b3
        
    def get_fermi_energy(self) -> float:
        """Get Fermi energy in eV."""
        return self.ef
