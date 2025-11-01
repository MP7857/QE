# F-Orbital Implementation in PWCOND - Summary

## Changes Made

This implementation adds support for f-orbitals (l=3) to the PWCOND module of Quantum ESPRESSO for complex band structure calculations in ballistic transport.

### Modified Files

1. **PWCOND/src/four.f90** - Core implementation
   - Extended array dimensions from 5 to 7 for all 7 f-orbital components (m=-3 to +3)
   - Added Bessel functions J₀, J₁, J₂, J₃ computations
   - Implemented boundary condition handling for l=3
   - Added normalization factors s1, s2, s3 for f-orbitals
   - Included z-dependent polynomial factors (z³, z², z)
   - Added auxiliary arrays x5, x6, fx5, fx6, wadd2

2. **PWCOND/src/scatter_forw.f90** - Calling routine
   - Updated w0 array allocation from dimension 5 to 7

### Documentation Created

3. **PWCOND/Doc/f_orbital_theory.md** - Comprehensive theoretical framework
   - Mathematical background and formulation
   - Spherical harmonics conventions and ordering
   - Bessel function integration theory
   - Normalization derivations
   - Complete code mapping and explanation

4. **PWCOND/Doc/theory_code_comparison.md** - Detailed verification
   - Side-by-side comparison of theory and implementation
   - Tables showing all orbital types (s, p, d, f)
   - Sign convention verification
   - Normalization factor validation

## Implementation Details

### Orbital Ordering
Following QE convention, the 7 f-orbital components are ordered as:
1. f_{z³} (m=0)
2. f_{xz²} (m=+1, real)
3. f_{yz²} (m=+1, imaginary)
4. f_{z(x²-y²)} (m=+2, real)
5. f_{xyz} (m=+2, imaginary)
6. f_{x(x²-3y²)} (m=+3, real)
7. f_{y(3x²-y²)} (m=+3, imaginary)

### Key Features
- **Bessel Functions**: Uses J₀ through J₃ with appropriate radial weight factors
- **Z-Dependence**: Implements up to cubic terms (5z³ - 3z for m=0)
- **Normalization**: Three separate factors (s1, s2, s3) for different m-value groups
- **Boundary Treatment**: Linear interpolation at slice boundaries with proper 1/|z|³ scaling
- **Memory Management**: All arrays properly allocated and deallocated

### Normalization Factors (l=3)
```fortran
s1 = (2π/4S) * sqrt(105/(4π))     ! For m=±3 and in-plane m=±2 components
s2 = (2π/2S) * sqrt(21/(2π))      ! For z-dependent m=±2,±1 components  
s3 = (2π/S) * sqrt(7/(32π))       ! For m=0 component with z³
```

## Testing and Validation

### Completed
- ✅ Syntax validation (no compilation errors in isolation)
- ✅ Array allocation balance verified
- ✅ Conditional logic completeness checked
- ✅ Pattern consistency with s, p, d implementations confirmed
- ✅ Normalization factors derived and validated
- ✅ Sign conventions verified

### Recommended Testing
- [ ] Full compilation with QE build system
- [ ] Test with f-electron pseudopotentials (rare earth elements, actinides)
- [ ] Verify transmission coefficients are physical (0 ≤ T ≤ 1)
- [ ] Check complex band structure consistency
- [ ] Performance benchmarking

## Usage

To use f-orbitals in PWCOND calculations:

1. Ensure pseudopotentials include f-orbitals (e.g., for lanthanides or actinides)
2. Run standard PWCOND workflow with systems containing f-electrons
3. The code will automatically detect l=3 and use the new implementation
4. No changes to input files required

## Mathematical Background

The implementation computes 2D Fourier transforms:
```
w₀(z,g,m) = (1/S) ∫ β(r) Y_{lm}(θ,φ) exp(-ig·r_⊥) dr_⊥
```

Using cylindrical coordinates and Bessel function expansions:
```
exp(-ig·r_⊥) = Σₙ (-i)ⁿ Jₙ(g·ρ) exp(in(φ-φ_g))
```

For l=3, this requires integrals of:
- β(r) × J₃(g·ρ) × (ρ³/r³) for m=±3 components
- β(r) × J₂(g·ρ) × (ρ²/r³) for m=±2 components  
- β(r) × J₁(g·ρ) × (ρ/r³) for m=±1 components
- β(r) × J₀(g·ρ) × (1/r³) for m=0 component

Combined with angular factors cos(mφ), sin(mφ) and z-polynomial factors.

## Consistency with QE Conventions

The implementation:
- Uses real spherical harmonics as defined in upflib/ylmr2.f90
- Follows the same orbital ordering convention as other QE modules
- Maintains sign conventions from existing s, p, d implementations
- Uses the same Bessel function routines (bessj.f90)
- Applies consistent normalization based on spherical harmonic theory

## References

See the detailed documentation files for:
- Complete mathematical derivations
- Verification of existing implementations
- Detailed code-to-theory mapping
- Comprehensive testing recommendations

## Authors

Implementation and documentation: GitHub Copilot
Based on original PWCOND code by A. Smogunov (2003)

## License

GNU General Public License v2 (same as Quantum ESPRESSO)
