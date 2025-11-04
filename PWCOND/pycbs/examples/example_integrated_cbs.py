#!/usr/bin/env python3
"""
Example: Integrated CBS calculation with physics-based Hamiltonian.

This example demonstrates how to use the integrated CBS calculator with
the HamiltonianBuilder for physics-based calculations.

This shows Phase 1/2 integration where:
- G-vector grid provides proper 2D plane wave basis
- HamiltonianBuilder constructs kinetic energy matrices
- CBS solver finds complex k-vectors using real physics

Note: This still uses free-electron Hamiltonian (no pseudopotential yet).
"""

import numpy as np
from pycbs import GVectorGrid, HamiltonianBuilder, ComplexBandStructure
from pycbs.reader import QEDataReader

def main():
    print("=" * 70)
    print("Integrated CBS Calculation with Physics-Based Hamiltonian")
    print("=" * 70)
    print()
    
    # ========================================================================
    # Step 1: Define system parameters (simulating QE output)
    # ========================================================================
    print("Step 1: Setting up system parameters")
    print("-" * 70)
    
    # Reciprocal lattice vectors (example: cubic cell with a=5 Bohr)
    alat = 5.0  # Bohr
    a = 2.0 * np.pi / alat
    bg = np.array([
        [a, 0, 0],
        [0, a, 0],
        [0, 0, a]
    ])
    
    print(f"Reciprocal lattice vectors (2π/a):")
    print(f"  bg[0] = [{bg[0,0]:.4f}, {bg[0,1]:.4f}, {bg[0,2]:.4f}]")
    print(f"  bg[1] = [{bg[1,0]:.4f}, {bg[1,1]:.4f}, {bg[1,2]:.4f}]")
    print(f"  bg[2] = [{bg[2,0]:.4f}, {bg[2,1]:.4f}, {bg[2,2]:.4f}]")
    print()
    
    # ========================================================================
    # Step 2: Construct 2D G-vector grid
    # ========================================================================
    print("Step 2: Constructing 2D G-vector grid")
    print("-" * 70)
    
    kpoint = np.array([0.0, 0.0])  # Gamma point in 2D
    ecut2d = 20.0  # Rydberg
    nrx, nry = 12, 12  # FFT grid
    
    grid = GVectorGrid(bg, kpoint, ecut2d, nrx, nry)
    
    print(f"2D k-point: ({kpoint[0]:.4f}, {kpoint[1]:.4f})")
    print(f"Energy cutoff: {ecut2d:.2f} Ry")
    print(f"FFT grid: {nrx} x {nry}")
    print(f"Number of G-vectors: {grid.ngper}")
    print(f"Number of shells: {grid.ngpsh}")
    print()
    
    # ========================================================================
    # Step 3: Build Hamiltonian matrices
    # ========================================================================
    print("Step 3: Building Hamiltonian matrices")
    print("-" * 70)
    
    energy = 10.0  # Rydberg
    builder = HamiltonianBuilder(grid, energy)
    
    # Test with different k_z values
    kz_values = [0.0, 0.5, 1.0, 0.5 + 0.1j]
    
    print(f"Energy: {energy:.2f} Ry")
    print(f"Testing kinetic energy matrices at different k_z:")
    print()
    
    for kz in kz_values:
        T = builder.build_kinetic_matrix(kz)
        
        # Check if Hermitian (for real k_z)
        is_hermitian = np.allclose(T, T.conj().T)
        
        # Compute eigenvalues to see energy spectrum
        eigvals = np.linalg.eigvalsh(T.real) if is_hermitian else np.linalg.eigvals(T)
        
        print(f"  k_z = {kz}:")
        print(f"    Matrix shape: {T.shape}")
        print(f"    Hermitian: {is_hermitian}")
        print(f"    Energy range: [{eigvals.min():.3f}, {eigvals.max():.3f}] Ry")
        print()
    
    # ========================================================================
    # Step 4: Build CBS matrices
    # ========================================================================
    print("Step 4: Building CBS matrices for generalized eigenvalue problem")
    print("-" * 70)
    
    kz_guess = 0.5
    amat, bmat = builder.build_cbs_matrices(kz_guess, nocros=0, noins=0)
    
    print(f"k_z guess: {kz_guess}")
    print(f"A matrix shape: {amat.shape}")
    print(f"B matrix shape: {bmat.shape}")
    print()
    
    # ========================================================================
    # Step 5: Solve generalized eigenvalue problem
    # ========================================================================
    print("Step 5: Solving generalized eigenvalue problem: A v = k B v")
    print("-" * 70)
    
    try:
        from scipy.linalg import eig
        
        # Solve GEP
        eigvals, eigvecs = eig(amat, bmat)
        
        # Sort by real part of eigenvalue
        idx = np.argsort(eigvals.real)
        eigvals = eigvals[idx]
        eigvecs = eigvecs[:, idx]
        
        print(f"Found {len(eigvals)} eigenvalues (k_z values)")
        print()
        
        # Analyze propagating vs evanescent modes
        tol = 1e-6
        propagating = np.abs(eigvals.imag) < tol
        evanescent = ~propagating
        
        n_prop = np.sum(propagating)
        n_evan = np.sum(evanescent)
        
        print(f"Propagating states (Im(k) ≈ 0): {n_prop}")
        print(f"Evanescent states (Im(k) ≠ 0): {n_evan}")
        print()
        
        if n_prop > 0:
            print("First 5 propagating states (k_z):")
            for i, kval in enumerate(eigvals[propagating][:5]):
                print(f"  {i+1}. k_z = {kval.real:.6f} + {kval.imag:.6e}j")
        
        print()
        
        if n_evan > 0:
            print("First 5 evanescent states (k_z):")
            for i, kval in enumerate(eigvals[evanescent][:5]):
                print(f"  {i+1}. k_z = {kval.real:.6f} + {kval.imag:.6f}j")
        
    except Exception as e:
        print(f"Error solving GEP: {e}")
    
    print()
    print("=" * 70)
    print("Integration Status:")
    print("  ✓ G-vector grid construction (Phase 1)")
    print("  ✓ Kinetic energy matrix (Phase 2)")
    print("  ✓ Free-electron Hamiltonian (Phase 2)")
    print("  ✓ CBS matrix formulation (Phase 2)")
    print("  ⏳ Pseudopotential integration (next)")
    print("  ⏳ Full CBS calculator integration (next)")
    print("=" * 70)


if __name__ == '__main__':
    main()
