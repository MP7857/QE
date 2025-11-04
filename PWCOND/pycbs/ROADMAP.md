# PyCBS Implementation Roadmap

## Current Status

PyCBS is currently a **prototype** that demonstrates the CBS calculation workflow using a simplified tight-binding model. It is **NOT suitable for production use** or accurate physics calculations.

## What's Implemented

✅ **Package Infrastructure**
- Clean Python package structure with setuptools
- Command-line interface (CLI)
- Comprehensive test suite (pytest)
- Documentation (README, INSTALL, guides)

✅ **API Design**
- `CBSCalculator` - High-level interface
- `QEDataReader` - XML parsing for system metadata
- `KPointGrid` - K-point management (automatic and custom)
- `CBSWriter` - PWCOND-compatible output format
- `read_kpoints_file()` - PWCOND format compatibility

✅ **Basic Functionality**
- Reading QE XML schema (lattice, atoms, reciprocal vectors)
- K-point grid generation
- Custom k-point specification
- Energy range configuration
- Output file generation (.re, .im, .3d files)

## What's NOT Implemented (Critical for Accuracy)

The following components use placeholder/simplified implementations and need to be replaced with real physics:

### 1. Wavefunction Data Reading

**Current:** Not implemented  
**Needed:** Read QE binary wavefunction files

**Files to read:**
- `${outdir}/${prefix}.save/wfc*.dat` - Wavefunction coefficients
- `${outdir}/${prefix}.save/charge-density.dat` - Charge density
- `${outdir}/${prefix}.save/data-file-schema.xml` - Full system description

**Implementation tasks:**
- Parse QE binary format (Fortran unformatted files)
- Extract plane wave coefficients
- Handle both norm-conserving and ultrasoft pseudopotentials
- Support spin-polarized and non-collinear cases

### 2. Hamiltonian Construction

**Current:** Simple tight-binding model with fixed parameters  
**Needed:** Construct actual Hamiltonian from QE data

**In `compbs.py`, function `_build_matrices()`:**

Current (lines 123-190):
```python
# Simple tight-binding: H = -t(|n><n+1| + h.c.) + ε₀|n><n|
t = 1.0  # Fixed hopping
eps0 = eryd  # Energy-dependent diagonal
```

Should be replaced with:
```python
# Build from QE potential and kinetic energy
# H = T + V_local + V_nl
# where T is kinetic energy operator in plane wave basis
# V_local is local pseudopotential
# V_nl is non-local pseudopotential
```

**Key tasks:**
- Compute kinetic energy matrix in 2D plane wave basis
- Add local potential (from charge density + pseudopotentials)
- Add non-local pseudopotential contributions
- Handle k-point dependency correctly
- Include spin-orbit coupling if present

### 3. Problem Dimensioning

**Current:** Fixed dimensions (n2d=10, nocros=2, noins=0)  
**Needed:** Calculate from actual system properties

**In `compbs.py`, function `_setup_problem_dimensions()`:**

Current (lines 97-121):
```python
self.n2d = 10      # Fixed
self.nocros = 2    # Fixed
self.noins = 0     # Fixed
```

Should be:
```python
# Calculate n2d from 2D cutoff
self.n2d = count_2d_plane_waves(ecut2d, b1, b2, kpoint)

# Calculate orbital counts from atomic structure
self.nocros, self.noins = calculate_orbital_counts(
    atoms, positions, orbitals, z_layers
)
```

**Key tasks:**
- Determine 2D plane wave count from cutoff energy
- Identify crossing orbitals based on atomic positions
- Classify interior vs crossing orbitals
- Handle different pseudopotential types

### 4. Overlap Matrix

**Current:** Identity matrix (orthogonal basis assumed)  
**Needed:** Actual overlap for ultrasoft pseudopotentials

For ultrasoft pseudopotentials, the overlap matrix is non-trivial:
```python
# S = I + Σ Q_ij |β_i><β_j|
# where Q_ij are augmentation charges
# and β_i are projector functions
```

### 5. Boundary Conditions

**Current:** Simple periodic boundary with fixed phase  
**Needed:** Proper scattering boundary conditions

The CBS calculation requires:
- Matching conditions at interfaces
- Evanescent wave decay conditions
- Proper normalization for bound/resonant states

## Implementation Priority

To make PyCBS production-ready, implement in this order:

### Phase 1: Core Physics (Essential)
1. **Wavefunction reader** - Parse QE binary files
2. **Kinetic energy matrix** - T = (k + G)²/(2m) in plane wave basis
3. **Local potential** - Read and interpolate from charge density
4. **Problem dimensions** - Calculate from system geometry

### Phase 2: Pseudopotentials (Important)
5. **Norm-conserving PP** - Simpler case, good for testing
6. **Non-local potential** - V_nl projections
7. **Ultrasoft PP** - Augmentation charges and overlap matrix
8. **Spin-orbit coupling** - For systems with heavy elements

### Phase 3: Advanced Features (Nice to have)
9. **Optimization** - Use sparse matrices, parallel k-points
10. **Validation** - Compare with PWCOND results systematically
11. **Additional output** - Decay lengths, channel decomposition
12. **GUI/visualization** - Plot complex bands interactively

## Validation Strategy

Once implemented, validate against PWCOND:

1. **Test case 1:** Al bulk (simple metal, norm-conserving PP)
2. **Test case 2:** Si bulk (semiconductor, norm-conserving PP)
3. **Test case 3:** Au bulk (metal with spin-orbit)
4. **Test case 4:** Al wire (1D system, matching PWCOND example01)

For each case:
- Run both PWCOND and PyCBS with identical inputs
- Compare complex k-values at multiple energies
- Check number of propagating channels
- Verify decay rates for evanescent states

## Resources Needed

### Technical Skills
- Fortran binary file I/O in Python
- Plane wave DFT theory
- Pseudopotential formalism
- Linear algebra (sparse matrices, eigensolvers)

### Reference Materials
- PWCOND source code (`PWCOND/src/compbs_2.f90`, `form_zk.f90`)
- QE documentation on file formats
- Smogunov et al., PRB 70, 045417 (2004) - Theory paper
- Quantum ESPRESSO user guide

### Estimated Effort
- Phase 1: 2-3 weeks (core physics)
- Phase 2: 2-3 weeks (pseudopotentials)
- Phase 3: 1-2 weeks (polish)
- **Total: ~6-8 weeks** for full implementation

## Getting Started with Implementation

If you want to contribute to making PyCBS production-ready:

1. **Start with Phase 1, item 1**: Implement the wavefunction reader
   - Look at `iotk` routines in QE for binary format
   - Use `numpy.fromfile()` or `struct` module
   - Test by reading a simple Al calculation

2. **Reference PWCOND code**: Study these key files:
   - `local_set.f90` - Sets up local potential
   - `init_gper.f90` - Initializes 2D problem
   - `form_zk.f90` - Constructs k-vectors
   - `compbs_2.f90` - Solves eigenvalue problem

3. **Incremental validation**: After each component:
   - Add unit tests comparing to known values
   - Test with PWCOND reference data
   - Document assumptions and limitations

## Conclusion

PyCBS provides an excellent foundation with clean architecture, comprehensive tests, and good documentation. The main gap is the physics implementation in `compbs.py`. 

With the roadmap above, PyCBS can become a fully functional CBS calculator. Until then, **use PWCOND for any calculations requiring accurate results**.

For questions or contributions, see the main repository.
