#!/usr/bin/env python3
"""
Example: Local potential integration (Phase 2)

This example demonstrates the local pseudopotential module structure.
Currently returns placeholder values as full implementation requires
reading QE binary charge density files.

This is Phase 2 of the PWCOND implementation roadmap.
"""

import numpy as np
from pycbs import (
    GVectorGrid,
    HamiltonianBuilder,
    LocalPotentialReader,
    PseudopotentialManager,
    integrate_local_potential
)

print("=" * 70)
print("Local Potential Integration Example (Phase 2)")
print("=" * 70)
print()

# Setup: Create G-vector grid for 2D plane wave basis
print("1. Setting up 2D G-vector grid...")
bg = np.eye(3) * 2 * np.pi  # Simple cubic reciprocal lattice
kpoint = np.array([0.0, 0.0])  # Gamma point in 2D
ecut2d = 50.0  # 2D energy cutoff in Rydberg

grid = GVectorGrid(bg, kpoint, ecut2d=ecut2d, nrx=16, nry=16)
print(f"   Found {grid.ngper} G-vectors in {grid.ngpsh} shells")
print()

# Create Hamiltonian builder
print("2. Creating Hamiltonian builder...")
energy = 10.0  # Energy in Rydberg
builder = HamiltonianBuilder(grid, energy=energy)
print(f"   Energy: {energy} Ry")
print(f"   Basis size: {builder.n2d}")
print()

# Create local potential reader
print("3. Creating local potential reader...")
outdir = './tmp'
prefix = 'al'
pot_reader = LocalPotentialReader(outdir, prefix)
print(f"   Output directory: {pot_reader.outdir}")
print(f"   Save directory: {pot_reader.save_dir}")
print()

# Try to read charge density
print("4. Attempting to read charge density...")
rho = pot_reader.read_charge_density()
if rho is None:
    print("   ⚠️  Charge density file not found (expected)")
    print("   📝 NOTE: This requires QE calculation output")
else:
    print(f"   ✅ Charge density shape: {rho.shape}")
print()

# Construct local potential matrix
print("5. Constructing local potential matrix...")
v_loc = pot_reader.construct_local_potential_2d(rho, grid)
print(f"   V_loc shape: {v_loc.shape}")
print(f"   V_loc dtype: {v_loc.dtype}")
print(f"   V_loc norm: {np.linalg.norm(v_loc):.6e}")
print()

# Integrate with Hamiltonian
print("6. Integrating local potential into Hamiltonian...")
v_integrated = integrate_local_potential(builder, pot_reader, grid)
print(f"   Integrated potential shape: {v_integrated.shape}")
print(f"   Integrated potential norm: {np.linalg.norm(v_integrated):.6e}")
print()

# Build full Hamiltonian (kinetic + potential)
print("7. Building full Hamiltonian...")
kz = 0.5 + 0.0j  # Real k_z for propagating state
T = builder.build_kinetic_matrix(kz)
H_full = T + v_integrated

print(f"   Kinetic energy matrix shape: {T.shape}")
print(f"   Full Hamiltonian shape: {H_full.shape}")
print(f"   Kinetic energy norm: {np.linalg.norm(T):.6f}")
print(f"   Potential energy norm: {np.linalg.norm(v_integrated):.6e}")
print(f"   Total Hamiltonian norm: {np.linalg.norm(H_full):.6f}")
print()

# Pseudopotential manager
print("8. Pseudopotential manager (framework)...")
pp_manager = PseudopotentialManager(['Al.UPF'])
print(f"   PP files: {pp_manager.pp_files}")
print(f"   📝 NOTE: UPF file parsing not yet implemented")
print()

# Summary
print("=" * 70)
print("Summary: Local Potential Module Structure")
print("=" * 70)
print()
print("✅ Completed:")
print("   - LocalPotentialReader class structure")
print("   - PseudopotentialManager class structure")
print("   - Integration framework with HamiltonianBuilder")
print("   - Module exports and tests")
print()
print("⏳ Remaining for Phase 2:")
print("   - Binary charge density file reader")
print("   - 3D → 2D potential projection algorithm")
print("   - UPF pseudopotential file parser")
print("   - Non-local pseudopotential projections")
print("   - Ultrasoft PP augmentation charges")
print()
print("📖 See ROADMAP.md for complete implementation plan")
print()
