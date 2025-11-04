#!/usr/bin/env python3
"""
Example: Using the Hamiltonian Builder (Phase 1.5/2)

This demonstrates the Hamiltonian construction capability,
which bridges the G-vector grid to actual CBS calculations.
"""

from pycbs import GVectorGrid, HamiltonianBuilder
import numpy as np


def example_kinetic_energy():
    """Example: Building kinetic energy matrix."""
    print("=" * 70)
    print("Example 1: Kinetic Energy Matrix Construction")
    print("=" * 70)
    
    # Define reciprocal lattice vectors (simple cubic)
    bg = np.eye(3) * 2 * np.pi
    
    # Define 2D k-point
    kpoint = np.array([0.0, 0.0])  # Gamma point
    
    # Construct G-vector grid
    print(f"\nConstructing G-vector grid at k = ({kpoint[0]}, {kpoint[1]})")
    grid = GVectorGrid(bg, kpoint, ecut2d=50.0, nrx=16, nry=16)
    
    # Build Hamiltonian
    energy = 10.0  # Ry
    builder = HamiltonianBuilder(grid, energy)
    
    print(f"\nBuilding kinetic energy matrix...")
    print(f"  Basis size: {builder.n2d} plane waves")
    
    # Build T matrix for different kz values
    for kz in [0.0, 0.5, 1.0, 0.5+0.1j]:
        T = builder.build_kinetic_matrix(kz)
        
        # Get eigenvalues
        eigs = np.linalg.eigvalsh(T.real)  # T is diagonal, so eigenvalues are diagonal elements
        
        print(f"\n  kz = {kz}:")
        print(f"    Matrix shape: {T.shape}")
        print(f"    Min kinetic energy: {np.min(eigs):.3f} Ry")
        print(f"    Max kinetic energy: {np.max(eigs):.3f} Ry")


def example_full_hamiltonian():
    """Example: Building full Hamiltonian matrix."""
    print("\n" + "=" * 70)
    print("Example 2: Full Hamiltonian Construction")
    print("=" * 70)
    
    # Setup
    bg = np.eye(3) * 2 * np.pi
    kpoint = np.array([0.25, 0.25])  # Non-zero k-point
    grid = GVectorGrid(bg, kpoint, ecut2d=50.0, nrx=16, nry=16)
    
    energy = 10.0  # Ry
    builder = HamiltonianBuilder(grid, energy)
    
    print(f"\nBuilding Hamiltonian at k = ({kpoint[0]}, {kpoint[1]})")
    print(f"  Energy: {energy} Ry")
    
    # Build simple Hamiltonian (free electron for now)
    H, S = builder.build_simple_hamiltonian(kz=0.5)
    
    print(f"\n  Hamiltonian H:")
    print(f"    Shape: {H.shape}")
    print(f"    Is Hermitian: {np.allclose(H, H.conj().T)}")
    
    print(f"\n  Overlap S:")
    print(f"    Shape: {S.shape}")
    print(f"    Is identity: {np.allclose(S, np.eye(S.shape[0]))}")
    
    # Solve eigenvalue problem H ψ = E S ψ
    from scipy.linalg import eigh
    eigenvalues, eigenvectors = eigh(H, S)
    
    print(f"\n  Eigenvalue spectrum:")
    print(f"    Min: {np.min(eigenvalues.real):.3f} Ry")
    print(f"    Max: {np.max(eigenvalues.real):.3f} Ry")
    print(f"    Number of states: {len(eigenvalues)}")


def example_cbs_matrices():
    """Example: Building CBS matrices for complex k-vector problem."""
    print("\n" + "=" * 70)
    print("Example 3: CBS Matrix Construction")
    print("=" * 70)
    
    # Setup
    bg = np.eye(3) * 2 * np.pi
    kpoint = np.array([0.0, 0.0])
    grid = GVectorGrid(bg, kpoint, ecut2d=50.0, nrx=16, nry=16)
    
    energy = 10.0  # Ry
    builder = HamiltonianBuilder(grid, energy)
    
    print(f"\nConstructing CBS matrices for generalized eigenvalue problem")
    print(f"  A v = k B v")
    print(f"  where k is the complex propagation constant")
    
    # Build CBS matrices
    amat, bmat = builder.build_cbs_matrices(kz=0.5+0.1j)
    
    print(f"\n  Matrix A:")
    print(f"    Shape: {amat.shape}")
    print(f"    Type: complex" if np.iscomplexobj(amat) else "    Type: real")
    
    print(f"\n  Matrix B:")
    print(f"    Shape: {bmat.shape}")
    print(f"    Is positive definite: {np.all(np.linalg.eigvalsh(bmat.real) > 0)}")
    
    # Solve generalized eigenvalue problem
    from scipy.linalg import eig
    eigenvalues, eigenvectors = eig(amat, bmat)
    
    # Complex k-values
    kvals = eigenvalues
    print(f"\n  Complex k-vectors found: {len(kvals)}")
    print(f"    Real parts: {np.min(kvals.real):.3f} to {np.max(kvals.real):.3f}")
    print(f"    Imag parts: {np.min(kvals.imag):.3f} to {np.max(kvals.imag):.3f}")
    
    # Count propagating vs evanescent modes
    propagating = np.sum(np.abs(kvals.imag) < 1e-6)
    evanescent = len(kvals) - propagating
    print(f"\n  Mode analysis:")
    print(f"    Propagating modes: {propagating}")
    print(f"    Evanescent modes: {evanescent}")


def example_integration_status():
    """Show integration status."""
    print("\n" + "=" * 70)
    print("Example 4: Integration Status & Next Steps")
    print("=" * 70)
    
    print("\nPhase 1.5/2 Progress: Hamiltonian Construction")
    print("  ✅ Kinetic energy matrix T = (G+k)²/(2m)")
    print("  ✅ Simple Hamiltonian construction")
    print("  ✅ CBS matrix formulation (A, B for GEP)")
    print("  ✅ Integration with G-vector grid")
    print("  ⏳ Local pseudopotential (next)")
    print("  ⏳ Non-local pseudopotential (next)")
    
    print("\nCurrent capabilities:")
    print("  - Free-electron Hamiltonian in 2D plane wave basis")
    print("  - Proper kinetic energy with complex k-vectors")
    print("  - CBS generalized eigenvalue problem structure")
    
    print("\nNext integration steps:")
    print("  1. Read local potential from QE charge density")
    print("  2. Add non-local pseudopotential projections")
    print("  3. Integrate with compbs.py for full CBS calculation")
    print("  4. Replace toy tight-binding model with real physics")
    
    print("\nThe Hamiltonian builder provides:")
    print("  - Proper basis (G-vector grid)")
    print("  - Kinetic energy operator")
    print("  - Framework for adding potentials")
    print("  - CBS matrix construction")


if __name__ == '__main__':
    print("PyCBS Hamiltonian Builder Examples")
    print("Phase 1.5/2: Hamiltonian Construction")
    print()
    
    example_kinetic_energy()
    example_full_hamiltonian()
    example_cbs_matrices()
    example_integration_status()
    
    print("\n" + "=" * 70)
    print("For more information, see:")
    print("  - ROADMAP.md - Full implementation plan")
    print("  - pycbs/hamiltonian.py - Source code")
    print("  - tests/test_hamiltonian.py - Test suite")
    print("=" * 70)
