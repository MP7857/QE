#!/usr/bin/env python3
"""
Example script for CBS calculation using PyCBS.

This example shows how to use PyCBS to calculate the complex band structure
of aluminum bulk along the (001) direction.
"""

from pycbs import CBSCalculator


def main():
    """Run CBS calculation example."""
    
    # Initialize calculator
    # Assumes you have run pw.x with prefix='al' in outdir='./tmp'
    calc = CBSCalculator(
        outdir='./tmp',
        prefix='al',
        band_file='bands.al',
        ecut2d=0.0,
        ewind=1.0,
        epsproj=1e-6
    )
    
    # Set energy range
    # Start at 10 eV above Fermi level, step down by 0.4 eV, 26 points
    calc.set_energy_range(
        energy0=10.0,
        denergy=-0.4,
        nenergy=26
    )
    
    # Set k-point grid (single k-point at Gamma for this example)
    calc.set_kpoints_grid(
        nk1=1,
        nk2=1,
        k1=0.0,
        k2=0.0
    )
    
    # Run calculation
    print("Starting CBS calculation for Al(001)...")
    print("=" * 70)
    
    results = calc.run()
    
    print("\n" + "=" * 70)
    print("Calculation completed!")
    print(f"Total k-points calculated: {len(results)}")
    
    # Print some results
    print("\nSample results:")
    for (ik, ien), result in list(results.items())[:3]:
        print(f"\nK-point {ik}, Energy point {ien}:")
        print(f"  Energy: {result['energy']:.4f} eV")
        print(f"  Number of channels: {result['nchan']}")
        print(f"  First few k-values:")
        for i, kval in enumerate(result['kvals'][:5]):
            print(f"    k[{i}] = {kval.real:.6f} + {kval.imag:.6f}i")
    
    print("\nOutput files generated:")
    print("  bands.al.re - Real part of k-vectors")
    print("  bands.al.im - Imaginary part of k-vectors")
    print("  bands.al.3d - Combined 3D data")
    print("  bands.al.co_re - Complex real part")
    print("  bands.al.co_im - Complex imaginary part")
    

if __name__ == '__main__':
    main()
