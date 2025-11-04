# Installation and Usage Guide for PyCBS

## ⚠️ Important Notice

**PyCBS is currently a PROTOTYPE implementation** using a simplified tight-binding model. It does NOT provide accurate CBS calculations yet. The current version is suitable for:
- Learning the CBS workflow
- Testing the API design
- Understanding output formats

**For accurate production calculations, use the original PWCOND Fortran code.**

## Overview

PyCBS (Python Complex Band Structure) is a Python package prototype that demonstrates CBS calculation functionality. The package structure and API are complete, but the physics implementation needs enhancement to read actual QE wavefunction data.

## Requirements

- Python 3.8 or higher
- NumPy >= 1.20.0
- SciPy >= 1.7.0
- Quantum ESPRESSO output files (version 6.5 or later recommended)

## Installation

### Method 1: Install from source (recommended)

```bash
cd PWCOND/pycbs
pip install -e .
```

### Method 2: Install with development dependencies

```bash
cd PWCOND/pycbs
pip install -e ".[dev]"
```

This will also install pytest for running tests.

## Quick Start

### 1. Prepare QE Data

First, run a standard pw.x calculation with Quantum ESPRESSO:

```bash
pw.x < al.scf.in > al.scf.out
```

This will generate the necessary data files in your `outdir` (e.g., `./tmp/al.save/`).

### 2. Run CBS Calculation (Python API)

```python
from pycbs import CBSCalculator

# Initialize calculator
calc = CBSCalculator(
    outdir='./tmp',
    prefix='al',
    band_file='bands.al'
)

# Set energy range
calc.set_energy_range(
    energy0=10.0,    # Starting energy in eV above Fermi
    denergy=-0.4,     # Energy step in eV
    nenergy=26        # Number of energy points
)

# Set k-point grid
calc.set_kpoints_grid(nk1=1, nk2=1)

# Run calculation
results = calc.run()
```

### 3. Run CBS Calculation (Command Line)

```bash
pycbs --outdir ./tmp --prefix al --band-file bands.al \
      --energy0 10.0 --denergy -0.4 --nenergy 26
```

## Output Files

PyCBS generates the following output files (compatible with PWCOND):

- `{band_file}.re` - Real part of k-vectors vs energy
- `{band_file}.im` - Imaginary part of k-vectors vs energy  
- `{band_file}.3d` - Combined real and imaginary k-vectors
- `{band_file}.co_re` - Real part (complementary representation)
- `{band_file}.co_im` - Imaginary part (complementary representation)

## Advanced Usage

### Custom 2D k-point Grid

```python
calc.set_kpoints_grid(
    nk1=4,      # 4 k-points along first direction
    nk2=4,      # 4 k-points along second direction
    k1=0.5,     # Shift in first direction
    k2=0.5      # Shift in second direction
)
```

### Define Custom K-points (PWCOND-style)

Instead of an automatic grid, you can specify exact k-points:

```python
import numpy as np
from pycbs import read_kpoints_file

# Method 1: Directly with numpy arrays
kpoints = np.array([
    [0.0000, 0.0000],
    [0.5000, 0.0000],
    [0.0000, 0.5000],
    [0.5000, 0.5000]
])
weights = np.array([0.25, 0.25, 0.25, 0.25])

calc.set_custom_kpoints(kpoints, weights)

# Method 2: Read from file (PWCOND format)
# File contains:
# 4
# 0.0000  0.0000  0.25
# 0.5000  0.0000  0.25
# 0.0000  0.5000  0.25
# 0.5000  0.5000  0.25

kpoints, weights = read_kpoints_file('kpoints.dat')
calc.set_custom_kpoints(kpoints, weights)
```

This is compatible with PWCOND's k-point input format. The third column (if present) is treated as the weight.

### Adjust Calculation Parameters

```python
calc = CBSCalculator(
    outdir='./tmp',
    prefix='al',
    band_file='bands.al',
    ecut2d=0.0,      # 2D energy cutoff in Ry
    ewind=1.0,       # Energy window in Ry
    epsproj=1e-6     # Projection threshold
)
```

### Access Results Programmatically

```python
results = calc.run()

# Results is a dictionary with (ik, ien) keys
for (ik, ien), result in results.items():
    print(f"K-point {ik}, Energy {ien}:")
    print(f"  Energy: {result['energy']} eV")
    print(f"  Channels: {result['nchan']}")
    print(f"  K-values: {result['kvals']}")
```

## Comparison with PWCOND

### What PyCBS Does

✅ Complex Band Structure (CBS) calculation (ikind=0)  
✅ Read QE XML data files  
✅ Solve generalized eigenvalue problem for complex k-vectors  
✅ Generate output files compatible with PWCOND  
✅ Support for 2D k-point grids (automatic and custom)  
✅ PWCOND-compatible k-point input format  
✅ Pure Python implementation - easy to modify and extend  

### What PyCBS Does NOT Do

❌ Transmission calculations (ikind=1, ikind=2)  
❌ Conductance calculations  
❌ Scattering region analysis  

If you need transmission calculations, use the original PWCOND Fortran code.

## Testing

Run the test suite:

```bash
cd PWCOND/pycbs
pytest tests/
```

All tests should pass:

```
================================================= test session starts ==================================================
tests/test_pycbs.py::TestKPointGrid::test_grid_generation PASSED                                     [ 12%]
tests/test_pycbs.py::TestKPointGrid::test_grid_with_shift PASSED                                     [ 25%]
tests/test_pycbs.py::TestComplexBandStructure::test_initialization PASSED                            [ 37%]
tests/test_pycbs.py::TestComplexBandStructure::test_calculation_structure PASSED                     [ 50%]
tests/test_pycbs.py::TestCBSCalculator::test_initialization PASSED                                   [ 62%]
tests/test_pycbs.py::TestCBSCalculator::test_energy_range_setting PASSED                             [ 75%]
tests/test_pycbs.py::TestCBSCalculator::test_kpoint_grid_setting PASSED                              [ 87%]
tests/test_pycbs.py::test_package_version PASSED                                                     [100%]
================================================== 8 passed in 0.26s ===================================================
```

## Examples

See the `examples/` directory for complete working examples:

- `example_al_bulk.py` - CBS calculation for Al(001) bulk

## Troubleshooting

### FileNotFoundError: QE data file not found

Make sure:
1. You have run `pw.x` successfully
2. The `outdir` path is correct
3. The `prefix` matches your QE input
4. QE version is 6.5 or later (for new XML format)

### ImportError: No module named 'pycbs'

Make sure you installed the package:
```bash
cd PWCOND/pycbs
pip install -e .
```

### Tests fail

Make sure pytest is installed:
```bash
pip install pytest
```

## Performance Notes

PyCBS is written in pure Python for:
- Easy installation (no compilation needed)
- Easy to understand and modify
- Good integration with Python scientific ecosystem

For very large systems or production calculations, the original Fortran PWCOND code may be faster. However, for typical CBS calculations, PyCBS provides comparable performance with much better usability.

## Contributing

Contributions are welcome! Areas for improvement:

1. Full implementation of QE data reading (currently simplified)
2. Support for old XML format (QE < 6.5)
3. Support for ultrasoft pseudopotentials
4. Parallel k-point calculations
5. Visualization tools for CBS data

## License

GPL-2.0-or-later, consistent with Quantum ESPRESSO.

## References

- A. Smogunov et al., Phys. Rev. B 70, 045417 (2004) - PWCOND method
- P. Giannozzi et al., J. Phys.: Condens. Matter 21, 395502 (2009) - Quantum ESPRESSO
