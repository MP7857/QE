# PWCOND Tools

This directory contains utility scripts for working with PWCOND calculations.

## check_pwcond_phases.py

A Python script for checking the phase structure of PWCOND's `four.f90` projectors.

### Prerequisites

- Python 3.x
- NumPy library (`pip install numpy`)

### Usage

**Note:** As of this PR, debug output is built into `four.f90` and can be enabled without code modifications.

1. **Enable debug output** by setting an environment variable before running PWCOND:

   ```bash
   export PWCOND_DEBUG_PHASES=1
   ```

2. **Run PWCOND** - it will automatically generate `w0_debug.dat` in the working directory

3. **Run the phase checker**:

   ```bash
   python3 tools/check_pwcond_phases.py w0_debug.dat --gauge canonical
   ```

   or

   ```bash
   python3 tools/check_pwcond_phases.py w0_debug.dat --gauge pwcond
   ```

### Debug File Format

The `w0_debug.dat` file contains:
- Global header (written once when file is created) with format information
- For each call to `four()`: 
  - Block separator and metadata comments (lb, nz1, ngper, nm)
  - Data lines: `kz ig m Re(w0) Im(w0)` (zero-based indices for Python)
  - A subset of the w0 array (first 3 kz values, first 5 ig values, all m channels)
- **Data is appended** across multiple iterations, allowing analysis of multiple orbitals

To disable debug output, simply unset the environment variable:
```bash
unset PWCOND_DEBUG_PHASES
```

**Note on file size:** For long-running calculations with many orbitals, the debug file can grow large. To manage this:
- Delete or rename `w0_debug.dat` before starting a new calculation series
- Use the subset output (only first 3 kz, first 5 ig values per orbital)
- Disable debug output when not actively analyzing phase patterns

### Output

The script will print a table showing:
- Estimated phase from numerical data
- Expected phase from theory
- Phase difference in degrees

Large phase differences (e.g., ~180°) indicate sign/gauge issues, while small differences indicate numerical noise only.

### Example Output

```
Phase analysis per w0(:,:,m):
 m  |  phi_est (num)      |  phi_ref (theory)   |  Δangle [deg]
----+----------------------+---------------------+-------------
  0 | +1.000+0.000i | +1.000+0.000i |       0.000
  1 | +0.000-1.000i | -0.000-1.000i |       0.000
  2 | +0.000-1.000i | -0.000-1.000i |       0.000
```

### Extending the Script

You can extend this script to:
- Parse QE binary formats
- Compute reference phases directly from UPF files
- Add automated tests with tolerance checks
- Support more complex m-channel mappings for d and f orbitals
