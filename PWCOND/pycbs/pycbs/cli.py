"""
Command-line interface for PyCBS.

This module provides a command-line interface for running CBS calculations.
"""

import argparse
import sys
from pathlib import Path

from . import CBSCalculator


def main():
    """Main entry point for command-line interface."""
    parser = argparse.ArgumentParser(
        description='PyCBS - Python Complex Band Structure Calculator',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic CBS calculation
  pycbs --outdir ./tmp --prefix al --band-file bands.al \\
        --energy0 10.0 --denergy -0.4 --nenergy 26

  # With custom k-point grid
  pycbs --outdir ./tmp --prefix al --band-file bands.al \\
        --energy0 10.0 --denergy -0.4 --nenergy 26 \\
        --nk1 4 --nk2 4
        """
    )
    
    # Required arguments
    parser.add_argument(
        '--outdir',
        type=str,
        required=True,
        help='Temporary directory containing QE output files'
    )
    parser.add_argument(
        '--prefix',
        type=str,
        required=True,
        help='Prefix for QE calculation files'
    )
    
    # Energy range arguments
    parser.add_argument(
        '--energy0',
        type=float,
        required=True,
        help='Starting energy relative to Fermi level (eV)'
    )
    parser.add_argument(
        '--denergy',
        type=float,
        required=True,
        help='Energy step (eV)'
    )
    parser.add_argument(
        '--nenergy',
        type=int,
        required=True,
        help='Number of energy points'
    )
    
    # Optional arguments
    parser.add_argument(
        '--band-file',
        type=str,
        default=None,
        help='Output file name for band structure (without extension)'
    )
    parser.add_argument(
        '--nk1',
        type=int,
        default=1,
        help='Number of k-points along first direction (default: 1)'
    )
    parser.add_argument(
        '--nk2',
        type=int,
        default=1,
        help='Number of k-points along second direction (default: 1)'
    )
    parser.add_argument(
        '--k1',
        type=float,
        default=0.0,
        help='Shift along first direction (default: 0.0)'
    )
    parser.add_argument(
        '--k2',
        type=float,
        default=0.0,
        help='Shift along second direction (default: 0.0)'
    )
    parser.add_argument(
        '--ecut2d',
        type=float,
        default=0.0,
        help='2D energy cutoff in Ry (default: 0.0)'
    )
    parser.add_argument(
        '--ewind',
        type=float,
        default=1.0,
        help='Energy window parameter in Ry (default: 1.0)'
    )
    parser.add_argument(
        '--epsproj',
        type=float,
        default=1e-6,
        help='Projection threshold (default: 1e-6)'
    )
    
    args = parser.parse_args()
    
    # Validate inputs
    if not Path(args.outdir).exists():
        print(f"Error: Directory {args.outdir} does not exist", file=sys.stderr)
        return 2  # Exit code 2 for missing files/directories
    
    try:
        # Initialize calculator
        calc = CBSCalculator(
            outdir=args.outdir,
            prefix=args.prefix,
            band_file=args.band_file,
            ecut2d=args.ecut2d,
            ewind=args.ewind,
            epsproj=args.epsproj
        )
        
        # Set energy range
        calc.set_energy_range(
            energy0=args.energy0,
            denergy=args.denergy,
            nenergy=args.nenergy
        )
        
        # Set k-point grid
        calc.set_kpoints_grid(
            nk1=args.nk1,
            nk2=args.nk2,
            k1=args.k1,
            k2=args.k2
        )
        
        # Run calculation
        print("=" * 70)
        print("PyCBS - Python Complex Band Structure Calculator")
        print("=" * 70)
        calc.run()
        
        print("\n" + "=" * 70)
        print("Calculation completed successfully!")
        print("=" * 70)
        
        if args.band_file:
            print(f"\nOutput files:")
            for ext in ['.re', '.im', '.3d', '.co_re', '.co_im']:
                filename = f"{args.band_file}{ext}"
                if Path(filename).exists():
                    print(f"  {filename}")
        
        return 0
        
    except FileNotFoundError as e:
        print(f"Error: Required file not found - {e}", file=sys.stderr)
        return 2  # Exit code 2 for missing files
    except ValueError as e:
        print(f"Error: Invalid input - {e}", file=sys.stderr)
        return 3  # Exit code 3 for invalid input
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1  # Exit code 1 for general errors


if __name__ == '__main__':
    sys.exit(main())
