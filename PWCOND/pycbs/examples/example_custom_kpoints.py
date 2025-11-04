#!/usr/bin/env python3
"""
Example showing how to use custom k-points in PyCBS.

This demonstrates the two ways to specify custom k-points:
1. Using numpy arrays directly
2. Reading from a PWCOND-format file
"""

import numpy as np
from pycbs import CBSCalculator, read_kpoints_file


def example_with_arrays():
    """Example using numpy arrays to define custom k-points."""
    print("=" * 70)
    print("Example 1: Using numpy arrays for custom k-points")
    print("=" * 70)
    
    # Initialize calculator
    calc = CBSCalculator(
        outdir='./tmp',
        prefix='al',
        band_file='bands.al_custom'
    )
    
    # Define custom k-points (similar to PWCOND format)
    # These are the transverse k-points
    kpoints = np.array([
        [0.0000, 0.0000],  # Gamma point
        [0.5000, 0.0000],  # X point
        [0.0000, 0.5000],  # Y point
        [0.5000, 0.5000],  # M point
    ])
    
    # Optional: define custom weights
    # If not provided, equal weights will be used
    weights = np.array([0.25, 0.25, 0.25, 0.25])
    
    # Set the custom k-points
    calc.set_custom_kpoints(kpoints, weights)
    
    # Set energy range
    calc.set_energy_range(
        energy0=10.0,
        denergy=-0.4,
        nenergy=26
    )
    
    # Run calculation
    print("\nRunning CBS calculation with custom k-points...")
    # results = calc.run()
    
    print("\nCustom k-points set successfully!")
    print(f"Number of k-points: {len(kpoints)}")
    for i, (kpt, w) in enumerate(zip(kpoints, weights)):
        print(f"  k[{i}] = ({kpt[0]:.4f}, {kpt[1]:.4f}), weight = {w:.4f}")


def example_with_file():
    """Example reading k-points from a file."""
    print("\n" + "=" * 70)
    print("Example 2: Reading k-points from file")
    print("=" * 70)
    
    # Create a sample k-points file (PWCOND format)
    kpoints_file = '/tmp/kpoints.dat'
    with open(kpoints_file, 'w') as f:
        f.write("4\n")
        f.write("0.0000  0.0000  0.25\n")
        f.write("0.5000  0.0000  0.25\n")
        f.write("0.0000  0.5000  0.25\n")
        f.write("0.5000  0.5000  0.25\n")
    
    print(f"\nCreated k-points file: {kpoints_file}")
    print("File contents:")
    with open(kpoints_file, 'r') as f:
        print(f.read())
    
    # Read k-points from file
    kpoints, weights = read_kpoints_file(kpoints_file)
    
    print(f"Read {len(kpoints)} k-points from file:")
    for i, (kpt, w) in enumerate(zip(kpoints, weights)):
        print(f"  k[{i}] = ({kpt[0]:.4f}, {kpt[1]:.4f}), weight = {w:.4f}")
    
    # Initialize calculator and set k-points
    calc = CBSCalculator(
        outdir='./tmp',
        prefix='al',
        band_file='bands.al_file'
    )
    
    calc.set_custom_kpoints(kpoints, weights)
    
    # Set energy range
    calc.set_energy_range(
        energy0=10.0,
        denergy=-0.4,
        nenergy=26
    )
    
    print("\nK-points loaded and ready for CBS calculation!")


def example_comparison():
    """Show the difference between grid and custom k-points."""
    print("\n" + "=" * 70)
    print("Example 3: Comparison of grid vs custom k-points")
    print("=" * 70)
    
    # Method 1: Automatic grid
    print("\nMethod 1: Automatic 2x2 grid")
    calc1 = CBSCalculator(outdir='./tmp', prefix='al')
    calc1.set_kpoints_grid(nk1=2, nk2=2)
    # This generates: (0,0), (0,0.5), (0.5,0), (0.5,0.5)
    print("  Generates regular 2x2 Monkhorst-Pack grid")
    
    # Method 2: Custom k-points
    print("\nMethod 2: Custom k-points")
    calc2 = CBSCalculator(outdir='./tmp', prefix='al')
    custom_kpts = np.array([
        [0.0, 0.0],
        [0.5, 0.0],
        [0.0, 0.5],
        [0.5, 0.5]
    ])
    calc2.set_custom_kpoints(custom_kpts)
    print("  Uses exactly specified k-points")
    print("  K-points:")
    for i, kpt in enumerate(custom_kpts):
        print(f"    k[{i}] = ({kpt[0]}, {kpt[1]})")
    
    print("\nAdvantages of custom k-points:")
    print("  - Specify exact k-points you want")
    print("  - Non-uniform grids possible")
    print("  - Can follow a specific k-path")
    print("  - Compatible with PWCOND input format")


if __name__ == '__main__':
    example_with_arrays()
    example_with_file()
    example_comparison()
    
    print("\n" + "=" * 70)
    print("All examples completed successfully!")
    print("=" * 70)
