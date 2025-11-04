# Quick Guide: Using Custom K-points in PyCBS

This guide shows you how to define custom transverse k-points for CBS calculations, similar to PWCOND.

## Method 1: Define k-points directly in Python

```python
import numpy as np
from pycbs import CBSCalculator

# Initialize calculator
calc = CBSCalculator(
    outdir='./tmp',
    prefix='al',
    band_file='bands.al'
)

# Define your custom k-points (transverse k-points)
kpoints = np.array([
    [0.0000, 0.0000],  # Gamma point
    [0.5000, 0.0000],  # X point
    [0.0000, 0.5000],  # Y point
    [0.5000, 0.5000],  # M point
])

# Optional: define weights (if not provided, equal weights are used)
weights = np.array([0.25, 0.25, 0.25, 0.25])

# Set the custom k-points
calc.set_custom_kpoints(kpoints, weights)

# Set energy range and run
calc.set_energy_range(energy0=10.0, denergy=-0.4, nenergy=26)
calc.run()
```

## Method 2: Read from file (PWCOND format)

### Step 1: Create a k-points file

Create a file named `kpoints.dat` with the following format:

```
4
0.0000  0.0000  0.25
0.5000  0.0000  0.25
0.0000  0.5000  0.25
0.5000  0.5000  0.25
```

**Format explanation:**
- First line: number of k-points
- Following lines: `kx  ky  weight` for each k-point
- `kx, ky` are the transverse k-point coordinates (in crystal units)
- `weight` is the k-point weight

### Step 2: Use in Python

```python
from pycbs import CBSCalculator, read_kpoints_file

# Read k-points from file
kpoints, weights = read_kpoints_file('kpoints.dat')

# Initialize calculator and set k-points
calc = CBSCalculator(outdir='./tmp', prefix='al', band_file='bands.al')
calc.set_custom_kpoints(kpoints, weights)

# Set energy range and run
calc.set_energy_range(energy0=10.0, denergy=-0.4, nenergy=26)
calc.run()
```

## Method 3: Converting from 3D format

If your k-points include the third component (kz), you can still use them:

```python
import numpy as np

# K-points with 3 components (kx, ky, kz)
kpoints_3d = np.array([
    [0.0000, 0.0000, 0.25],
    [0.5000, 0.0000, 0.25],
    [0.0000, 0.5000, 0.25],
    [0.5000, 0.5000, 0.25]
])

# PyCBS will automatically use only the first two columns
calc.set_custom_kpoints(kpoints_3d)
# The third column is interpreted as weights
```

## Advantages of Custom K-points

1. **Exact control**: Specify exactly which k-points to calculate
2. **Non-uniform grids**: Use different spacing in different regions
3. **K-paths**: Calculate along specific high-symmetry paths
4. **PWCOND compatible**: Use the same input format as PWCOND
5. **Flexible**: Mix different k-point densities as needed

## Comparison: Grid vs Custom

### Automatic Grid
```python
# Generates regular Monkhorst-Pack grid
calc.set_kpoints_grid(nk1=4, nk2=4)
# Results in 16 evenly-spaced k-points
```

### Custom K-points
```python
# Specify exactly the k-points you want
kpoints = np.array([
    [0.00, 0.00],
    [0.25, 0.00],
    [0.50, 0.00],
    [0.50, 0.50],
])
calc.set_custom_kpoints(kpoints)
# Results in 4 k-points at your chosen positions
```

## Complete Example

```python
#!/usr/bin/env python3
import numpy as np
from pycbs import CBSCalculator

# Initialize
calc = CBSCalculator(
    outdir='./tmp',
    prefix='al',
    band_file='bands.al_custom'
)

# Define 4 custom k-points
kpoints = np.array([
    [0.0, 0.0],  # Gamma
    [0.5, 0.0],  # X
    [0.0, 0.5],  # Y
    [0.5, 0.5],  # M
])

calc.set_custom_kpoints(kpoints)

# Energy range: 10 eV to 0 eV in steps of -0.4 eV
calc.set_energy_range(
    energy0=10.0,
    denergy=-0.4,
    nenergy=26
)

# Run CBS calculation
print("Running CBS calculation with custom k-points...")
results = calc.run()

print(f"\nCalculation complete!")
print(f"Calculated {len(results)} (k, E) combinations")
print(f"Output files: bands.al_custom.re, .im, .3d, etc.")
```

## See Also

- Full example: `examples/example_custom_kpoints.py`
- API documentation: `README.md`
- Installation guide: `INSTALL.md`
