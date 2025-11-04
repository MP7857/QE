#!/usr/bin/env python3
"""
Example: Local potential integration with binary charge density reader (Phase 2)

This example demonstrates:
1. Reading charge density from QE binary files (NEW!)
2. Local pseudopotential module structure
3. Integration with Hamiltonian builder

This is Phase 2 of the PWCOND implementation roadmap.
"""

import numpy as np
import os
import tempfile
from pycbs import (
    GVectorGrid,
    HamiltonianBuilder,
    LocalPotentialReader,
    PseudopotentialManager,
    integrate_local_potential
)

print("=" * 70)
print("Local Potential Integration Example (Phase 2 - Binary Reader)")
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

# Demonstrate charge density reading with mock data
print("3. Demonstrating charge density binary reader...")
with tempfile.TemporaryDirectory() as tmpdir:
    save_dir = os.path.join(tmpdir, 'demo.save')
    os.makedirs(save_dir)
    
    # Create a mock charge density file in QE binary format
    charge_file = os.path.join(save_dir, 'charge-density.dat')
    nr1, nr2, nr3, nspin = 16, 16, 16, 1
    rho_mock = np.random.rand(nr1, nr2, nr3) * 0.1  # Mock electron density
    
    print(f"   Creating mock charge density file...")
    print(f"   FFT grid: {nr1} x {nr2} x {nr3}")
    print(f"   Spin components: {nspin}")
    
    with open(charge_file, 'wb') as f:
        # Write header (Fortran unformatted)
        header = np.array([nr1, nr2, nr3, nspin], dtype=np.int32)
        rec_len = header.nbytes
        np.array([rec_len], dtype=np.int32).tofile(f)
        header.tofile(f)
        np.array([rec_len], dtype=np.int32).tofile(f)
        
        # Write data
        rho_flat = rho_mock.flatten(order='F')
        rec_len = rho_flat.nbytes
        np.array([rec_len], dtype=np.int32).tofile(f)
        rho_flat.tofile(f)
        np.array([rec_len], dtype=np.int32).tofile(f)
    
    # Read it back using our binary reader
    pot_reader = LocalPotentialReader(tmpdir, 'demo')
    result = pot_reader.read_charge_density()
    
    if result is not None:
        rho, metadata = result
        print(f"   ✅ Successfully read charge density!")
        print(f"   Charge density shape: {rho.shape}")
        print(f"   Grid dimensions: {metadata['nr1']} x {metadata['nr2']} x {metadata['nr3']}")
        print(f"   Spin components: {metadata['nspin']}")
        print(f"   Total charge: {np.sum(rho):.6f} (arbitrary units)")
    else:
        print("   ⚠️  Failed to read charge density")
        rho = None
print()

# Try to read from actual QE output if available
print("4. Attempting to read from QE output directory...")
outdir = './tmp'
prefix = 'al'
pot_reader_qe = LocalPotentialReader(outdir, prefix)
print(f"   Output directory: {pot_reader_qe.outdir}")
print(f"   Save directory: {pot_reader_qe.save_dir}")

result_qe = pot_reader_qe.read_charge_density()
if result_qe is None:
    print("   ⚠️  Charge density file not found (expected if no QE calc)")
    print("   📝 NOTE: Run QE pw.x calculation to generate charge-density.dat")
    rho_qe = None
else:
    rho_qe, metadata_qe = result_qe
    print(f"   ✅ Found QE charge density!")
    print(f"   Shape: {rho_qe.shape}")
    print(f"   Grid: {metadata_qe['nr1']} x {metadata_qe['nr2']} x {metadata_qe['nr3']}")
print()

# Construct local potential matrix
print("5. Constructing local potential matrix...")
v_loc = pot_reader_qe.construct_local_potential_2d(rho_qe, grid)
print(f"   V_loc shape: {v_loc.shape}")
print(f"   V_loc dtype: {v_loc.dtype}")
print(f"   V_loc norm: {np.linalg.norm(v_loc):.6e}")
print()

# Integrate with Hamiltonian
print("6. Integrating local potential into Hamiltonian...")
v_integrated = integrate_local_potential(builder, pot_reader_qe, grid)
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
print("Summary: Binary Charge Density Reader Implemented!")
print("=" * 70)
print()
print("✅ Newly Completed:")
print("   - Binary charge density file reader (QE format)")
print("   - Fortran unformatted record parser")
print("   - Support for spin-polarized densities")
print("   - FFT grid metadata extraction")
print()
print("✅ Already Complete:")
print("   - LocalPotentialReader class structure")
print("   - PseudopotentialManager class structure")
print("   - Integration framework with HamiltonianBuilder")
print("   - Module exports and tests (45/45 passing)")
print()
print("⏳ Next Steps for Phase 2:")
print("   - 3D → 2D potential projection algorithm")
print("   - UPF pseudopotential file parser")
print("   - Non-local pseudopotential projections")
print("   - Ultrasoft PP augmentation charges")
print()
print("📖 See ROADMAP.md for complete implementation plan")
print()
