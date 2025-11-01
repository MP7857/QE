# Quick Reference: f-Orbital Implementation

## At a Glance

**Purpose**: Add f-orbital (l=3) support to PWCOND for complex band structure calculations

**Files Modified**: 
- `PWCOND/src/four.f90` (core implementation)
- `PWCOND/src/scatter_forw.f90` (array dimensions)

**Documentation**: 3 comprehensive markdown files in `PWCOND/Doc/`

## Key Implementation Points

### Array Dimension Change
```fortran
! Before:
complex(DP) :: w0(nz1, ngper, 5)

! After:
complex(DP) :: w0(nz1, ngper, 7)
```

### Bessel Functions Used for l=3
- J₀(g·ρ) - zeroth order
- J₁(g·ρ) - first order  
- J₂(g·ρ) - second order
- J₃(g·ρ) - third order (NEW for f-orbitals)

### The 7 f-Orbital Components
1. f_{z³} - along transport axis
2. f_{xz²} - mixed, x-component
3. f_{yz²} - mixed, y-component
4. f_{z(x²-y²)} - mixed, symmetric
5. f_{xyz} - mixed, antisymmetric
6. f_{x(x²-3y²)} - in-plane, x-dominant
7. f_{y(3x²-y²)} - in-plane, y-dominant

### Normalization Factors
```fortran
s1 = 2π/(4S) × √(105/4π)    ! For high-m components
s2 = 2π/(2S) × √(21/2π)     ! For mid-m components
s3 = 2π/S × √(7/32π)        ! For m=0 component
```

### Z-Dependence
- **m=0**: 5z³ - 3z (cubic)
- **m=±1**: z² (quadratic)
- **m=±2**: z (linear)
- **m=±3**: constant (in-plane)

## Code Pattern

Each angular momentum l follows this pattern:

| l | Orbitals | Components | Bessel Orders | Z-power | Arrays |
|---|----------|------------|---------------|---------|--------|
| 0 | s        | 1          | J₀            | z⁰      | 1      |
| 1 | p        | 3          | J₀,J₁         | z¹      | 2      |
| 2 | d        | 5          | J₀,J₁,J₂      | z²      | 4      |
| 3 | f        | 7          | J₀,J₁,J₂,J₃   | z³      | 6      |

## Usage

No input file changes needed! Just:
1. Use f-electron pseudopotentials (lmax=3)
2. Run pwcond.x normally
3. Code auto-detects l=3 and uses new implementation

## Testing Checklist

- [ ] Compile with full QE build system
- [ ] Test with Ce, Gd, or Yb pseudopotentials
- [ ] Verify 0 ≤ T ≤ 1 for transmission
- [ ] Check complex band structure
- [ ] Verify current conservation

## For More Details

See comprehensive documentation in `PWCOND/Doc/`:
- `f_orbital_theory.md` - Full mathematical theory
- `theory_code_comparison.md` - Verification tables
- `F_ORBITAL_README.md` - Implementation summary
