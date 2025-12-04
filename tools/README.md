# PWCOND Tools

This directory contains utility scripts for working with PWCOND calculations.

## check_pwcond_phases.py

A Python script for checking the phase structure of PWCOND's `four.f90` projectors.

### Prerequisites

- Python 3.x
- NumPy library (`pip install numpy`)

### Usage

1. **Instrument `four.f90`** to write a debug file. Add the following code after `w0` is fully assembled:

   ```fortran
   ! Debug output for phase checking
   integer :: debug_unit
   if (ig == 1 .and. kz == 1) then
      call find_free_unit(debug_unit)
      open(unit=debug_unit, file='w0_debug.dat', status='unknown')
      do m = 1, 7
         write(debug_unit, '(3I4,2ES20.10)') kz-1, ig-1, m-1, &
              real(w0(kz,ig,m)), aimag(w0(kz,ig,m))
      enddo
      close(debug_unit)
   endif
   ```

2. **Run PWCOND** to generate the debug file (e.g., `w0_debug.dat`)

3. **Run the phase checker**:

   ```bash
   python3 tools/check_pwcond_phases.py w0_debug.dat --gauge canonical
   ```

   or

   ```bash
   python3 tools/check_pwcond_phases.py w0_debug.dat --gauge pwcond
   ```

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
