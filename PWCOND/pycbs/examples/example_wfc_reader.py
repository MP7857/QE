#!/usr/bin/env python3
"""
Example: Using the Wavefunction Reader (Phase 1)

This demonstrates the new wavefunction reading capability,
which is the first step toward full PWCOND implementation.
"""

from pycbs import WavefunctionReader, read_wavefunction_metadata


def example_metadata():
    """Example: Get metadata about wavefunction files."""
    print("=" * 70)
    print("Example 1: Reading Wavefunction Metadata")
    print("=" * 70)
    
    # This assumes you have run pw.x and have a prefix.save directory
    outdir = './tmp'
    prefix = 'al'
    
    print(f"\nLooking for wavefunction data in {outdir}/{prefix}.save/")
    
    try:
        metadata = read_wavefunction_metadata(outdir, prefix)
        
        print(f"\nFound {metadata['num_wfc_files']} wavefunction files:")
        for wfc_file in metadata['wfc_files']:
            print(f"  - {wfc_file}")
            
        if metadata['num_wfc_files'] > 0:
            print(f"\nSample header from first file:")
            for key, value in metadata.get('sample_header', {}).items():
                print(f"  {key}: {value}")
        else:
            print("\nNo wavefunction files found.")
            print("Make sure you have run pw.x with this prefix.")
            
    except FileNotFoundError as e:
        print(f"\nError: {e}")
        print("\nTo use this example:")
        print("1. Run pw.x calculation: pw.x < input.in > output.out")
        print("2. Make sure outdir and prefix match your calculation")


def example_detailed_reader():
    """Example: Using WavefunctionReader directly."""
    print("\n" + "=" * 70)
    print("Example 2: Using WavefunctionReader Class")
    print("=" * 70)
    
    outdir = './tmp'
    prefix = 'al'
    
    # Initialize reader
    reader = WavefunctionReader(outdir, prefix)
    print(f"\nInitialized reader for {prefix}")
    print(f"Save directory: {reader.save_dir}")
    
    try:
        # Find all wavefunction files
        wfc_files = reader.find_wfc_files()
        print(f"\nFound {len(wfc_files)} wavefunction files")
        
        # Read header from first file
        if wfc_files:
            print(f"\nReading header from {wfc_files[0].name}...")
            header = reader.read_wfc_header(wfc_files[0])
            
            print("Header contents:")
            for key, value in header.items():
                if key == 'error':
                    print(f"  Note: {value}")
                else:
                    print(f"  {key}: {value}")
                    
            # Note about full implementation
            print("\n" + "-" * 70)
            print("NOTE: Full wavefunction coefficient reading is not yet implemented.")
            print("This requires:")
            print("  - Proper Fortran unformatted file handling")
            print("  - QE version-specific format support")
            print("  - Handling of different data types (gamma_only, spin, etc.)")
            print("\nSee ROADMAP.md for full implementation plan.")
            print("-" * 70)
            
    except FileNotFoundError as e:
        print(f"\nError: {e}")


def example_integration_status():
    """Show integration status with main CBS calculator."""
    print("\n" + "=" * 70)
    print("Example 3: Integration Status")
    print("=" * 70)
    
    print("\nPhase 1 Progress: Wavefunction Reader")
    print("  ✅ File discovery (finding wfc*.dat files)")
    print("  ✅ Basic header reading (k-point index, npw)")
    print("  ✅ Metadata extraction")
    print("  ⏳ Full coefficient reading (in progress)")
    print("  ⏳ G-vector reading (in progress)")
    
    print("\nNext Steps (Phase 1 remaining):")
    print("  1. Complete binary format parser for all QE versions")
    print("  2. Read G-vector indices from gvectors.dat")
    print("  3. Read plane wave coefficients for all bands")
    print("  4. Handle gamma_only, spin-polarized, non-collinear cases")
    
    print("\nPhase 2: Hamiltonian Construction (not started)")
    print("  ⏳ Build kinetic energy matrix")
    print("  ⏳ Add local pseudopotential")
    print("  ⏳ Add non-local pseudopotential")
    print("  ⏳ Construct overlap matrix")
    
    print("\nTo use the wavefunction reader in CBS calculations:")
    print("  - It will be integrated into compbs.py")
    print("  - Will replace the toy tight-binding model")
    print("  - Provides actual QE data to CBS algorithm")


if __name__ == '__main__':
    print("PyCBS Wavefunction Reader Examples")
    print("Phase 1: Reading QE Binary Data")
    print()
    
    example_metadata()
    example_detailed_reader()
    example_integration_status()
    
    print("\n" + "=" * 70)
    print("For more information, see:")
    print("  - ROADMAP.md - Full implementation plan")
    print("  - pycbs/wfc_reader.py - Source code")
    print("  - tests/test_wfc_reader.py - Test suite")
    print("=" * 70)
