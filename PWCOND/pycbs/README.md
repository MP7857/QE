# PyCBS - Python Complex Band Structure Calculator

PyCBS is a Python package that provides Complex Band Structure (CBS) calculation functionality, serving as a Python-based replacement for the PWCOND Fortran code in Quantum ESPRESSO. It reads simulation data files from Quantum ESPRESSO and calculates the complex band structure without performing transmission calculations.

## ⚠️ Current Status - Phase 2 CONTINUING

**IMPORTANT**: Significant progress on physics-based implementation! Local potential module framework added.

**Phase 1 & 2 Progress - Wavefunction Reader, Hamiltonian & Potential:**
- ✅ **Completed**: Basic wavefunction file reader (`wfc_reader.py`)
- ✅ File discovery and metadata extraction
- ✅ Basic binary header parsing
- ✅ 2D G-vector grid construction (`GVectorGrid` class - following init_gper.f90)
- ✅ Kinetic energy matrix construction (`HamiltonianBuilder` class)
- ✅ Free-electron Hamiltonian in 2D plane wave basis
- ✅ CBS matrix formulation for generalized eigenvalue problem
- ✅ Integration with CBS calculator (uses HamiltonianBuilder when G-vector grid available)
- ✅ **NEW**: Local potential module framework (`potential.py`)
- ✅ **NEW**: `LocalPotentialReader` and `PseudopotentialManager` classes
- ⏳ **In Progress**: Binary charge density file reader
- ⏳ **In Progress**: 3D → 2D potential projection algorithm
- ⏳ **In Progress**: UPF pseudopotential file parser
- ⏳ **Next**: Non-local pseudopotential projections
- ⏳ **Next**: Full coefficient reading from binary files

**What works NOW:**
- Free-electron complex band structure calculations with proper plane wave basis
- Physics-based kinetic energy matrices T = (G+k)²/(2m)
- Proper generalized eigenvalue problem formulation
- Automatic fallback to toy model when G-vector grid not available
- Local potential integration framework (placeholder implementations)

**For accurate results with pseudopotentials**, the following still needs implementation:
- Binary charge density file reading and potential construction
- UPF pseudopotential file parsing
- Non-local pseudopotential projections  
- Complete wavefunction coefficient reading
- Validation against PWCOND results

**Current use cases:**
- ✅ Learning the CBS calculation workflow
- ✅ Testing the package structure and API
- ✅ Understanding output formats
- ✅ Exploring QE wavefunction file structure
- ✅ Constructing 2D G-vector grids (following PWCOND init_gper.f90)
- ✅ Building kinetic energy and Hamiltonian matrices
- ✅ Understanding CBS generalized eigenvalue problem
- ✅ Free-electron CBS calculations
- ✅ **NEW**: Understanding local potential integration framework
- ❌ Production calculations with pseudopotentials (Phase 2 continuation needed)
- ❌ Publication-quality data

For production use with accurate physics and pseudopotentials, please use the original PWCOND Fortran code.

See `ROADMAP.md` for detailed implementation status, `examples/example_wfc_reader.py` for G-vector functionality, `examples/example_hamiltonian.py` for Hamiltonian construction, `examples/example_integrated_cbs.py` for integrated physics-based CBS calculation, and **NEW** `examples/example_potential.py` for local potential framework.

## Features

- Read Quantum ESPRESSO output files (XML and binary formats)
- Calculate Complex Band Structure (CBS) for materials
- Support for both collinear and non-collinear calculations
- Output in multiple formats (.re, .im, .3d files)
- Pure Python implementation with NumPy/SciPy for numerical operations

## Installation

Install from the source directory:

```bash
cd PWCOND/pycbs
pip install -e .
```

## Usage

### As a Python library

```python
from pycbs import CBSCalculator

# Initialize calculator with QE data directory and prefix
calc = CBSCalculator(
    outdir='./tmp',
    prefix='al',
    band_file='bands.al'
)

# Set energy range and k-points
calc.set_energy_range(energy0=10.0, denergy=-0.4, nenergy=26)
calc.set_kpoints_grid(nk1=1, nk2=1)

# Run CBS calculation
calc.run()

# Access results
results = calc.get_results()
print(f"Number of channels: {results['nchan']}")
print(f"K-values: {results['kvals']}")
```

### Using custom k-points (like PWCOND)

You can specify custom transverse k-points instead of using an automatic grid:

```python
from pycbs import CBSCalculator, read_kpoints_file
import numpy as np

# Method 1: Define k-points directly
calc = CBSCalculator(outdir='./tmp', prefix='al', band_file='bands.al')

kpoints = np.array([
    [0.0000, 0.0000],  # Gamma
    [0.5000, 0.0000],  # X
    [0.0000, 0.5000],  # Y
    [0.5000, 0.5000],  # M
])
weights = np.array([0.25, 0.25, 0.25, 0.25])

calc.set_custom_kpoints(kpoints, weights)
calc.set_energy_range(energy0=10.0, denergy=-0.4, nenergy=26)
calc.run()

# Method 2: Read from file (PWCOND format)
# File format:
# 4
# 0.0000  0.0000  0.25
# 0.5000  0.0000  0.25
# 0.0000  0.5000  0.25
# 0.5000  0.5000  0.25

kpoints, weights = read_kpoints_file('kpoints.dat')
calc.set_custom_kpoints(kpoints, weights)
```

### Command-line interface

```bash
pycbs --outdir ./tmp --prefix al --band-file bands.al \
      --energy0 10.0 --denergy -0.4 --nenergy 26
```

## Input Files

PyCBS requires the following input files generated by Quantum ESPRESSO's pw.x:

1. XML data file: `{prefix}.xml` in the `{outdir}` directory
2. Binary data files in `{outdir}/{prefix}.save/`

## Output Files

When `band_file` is specified, PyCBS generates:

- `{band_file}.re` - Real part of k-vectors vs energy
- `{band_file}.im` - Imaginary part of k-vectors vs energy  
- `{band_file}.3d` - Combined real and imaginary k-vectors
- `{band_file}.co_re` - Real part (complementary representation)
- `{band_file}.co_im` - Imaginary part (complementary representation)

## Comparison with PWCOND

PyCBS implements the Complex Band Structure calculation (ikind=0) from PWCOND:

| Feature | PWCOND | PyCBS |
|---------|--------|-------|
| CBS Calculation | ✓ | ✓ |
| Transmission Calculation | ✓ | ✗ |
| Language | Fortran | Python |
| Dependencies | QE libraries | NumPy, SciPy |

## Theory

Complex band structure describes the behavior of electronic states in materials, including both propagating (real k) and evanescent (imaginary k) states. This is essential for understanding:

- Quantum transport in nanostructures
- Tunneling phenomena
- Contact resistance in devices
- Electron injection at interfaces

## Contributing

Contributions are welcome! Please see the main Quantum ESPRESSO repository for contribution guidelines.

## License

This package is distributed under the GNU General Public License v2 or later, consistent with Quantum ESPRESSO.

## Citation

If you use PyCBS in your research, please cite:
- Quantum ESPRESSO: P. Giannozzi et al., J. Phys.: Condens. Matter 21, 395502 (2009)
- PWCOND method: A. Smogunov et al., Phys. Rev. B 70, 045417 (2004)
