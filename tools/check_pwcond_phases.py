#!/usr/bin/env python3
"""
check_pwcond_phases.py

Scaffold for checking the phase structure of PWCOND's four.f90 projectors.

Intended usage:
1. Modify four.f90 to optionally write a debug file (e.g. 'w0_debug.dat')
   containing w0(kz,ig,m) for one or a few (kz,ig).
2. Run this script to:
   - load the debug data,
   - group channels by (l, m),
   - estimate the complex phase of each channel,
   - compare against the expected phase pattern.

You (or GitHub Copilot) can extend this file by:
- Parsing real QE binary formats,
- Computing reference "canonical" phases directly from the UPF,
- Adding numerical tolerances and automated tests.
"""

import numpy as np
from pathlib import Path
from typing import Dict, Tuple

# ----------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------

# Expected canonical phase pattern for real harmonics per |m|
# (up to an overall ± sign per orbital)
EXPECTED_PHASES = {
    0: 1.0 + 0.0j,   # m = 0 → +1 (real)
    1: -1.0j,        # m = 1 → -i
    2: -1.0 + 0.0j,  # m = 2 → -1
    3: 1.0j,         # m = 3 → +i
}

# If you want to check the current PWCOND gauge (+1,+i,-1,-i), set this instead:
PWCOND_GAUGE_PHASES = {
    0: 1.0 + 0.0j,
    1: 1.0j,
    2: -1.0 + 0.0j,
    3: -1.0j,
}

# Choose which pattern you want to check against:
REFERENCE_PHASES = EXPECTED_PHASES  # or PWCOND_GAUGE_PHASES


# ----------------------------------------------------------------------
# Debug file parsing
# ----------------------------------------------------------------------

def load_w0_debug_ascii(path: Path) -> np.ndarray:
    """
    Example loader for a simple ASCII debug file written by four.f90.

    Expected format (you can adapt four.f90 accordingly):
        kz ig m  Re(w0)  Im(w0)

    Returns:
        w0[kz, ig, m] as a complex numpy array,
        with zero-based indices for kz, ig, m.
    """
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data[None, :]

    kz = data[:, 0].astype(int)
    ig = data[:, 1].astype(int)
    m  = data[:, 2].astype(int)
    re = data[:, 3]
    im = data[:, 4]

    nkz = kz.max() + 1
    nig = ig.max() + 1
    nm  = m.max() + 1

    w0 = np.zeros((nkz, nig, nm), dtype=np.complex128)
    for kk, gg, mm, rr, ii in zip(kz, ig, m, re, im):
        w0[kk, gg, mm] = rr + 1j * ii

    return w0


# ----------------------------------------------------------------------
# Phase analysis
# ----------------------------------------------------------------------

def estimate_phase(z: np.ndarray) -> complex:
    """
    Estimate the dominant complex phase of a vector of complex numbers
    (e.g. all w0(kz,ig,m) for various kz,ig) by least-squares.

    Returns a complex phase of unit modulus (approximately).
    """
    # Avoid zeros
    mask = np.abs(z) > 1e-12
    if not np.any(mask):
        return 0.0 + 0.0j

    z_nonzero = z[mask]
    # Normalize amplitudes and average
    z_norm = z_nonzero / np.abs(z_nonzero)
    phase = z_norm.mean()
    # Normalize to unit modulus
    if np.abs(phase) > 0:
        phase /= np.abs(phase)
    return phase


def phase_difference(phi_num: complex, phi_ref: complex) -> float:
    """
    Compute the angle (in radians) between two complex phases.
    """
    if phi_num == 0 or phi_ref == 0:
        return np.nan
    ratio = phi_num / phi_ref
    return np.angle(ratio)


def analyze_w0_phases(
    w0: np.ndarray,
    m_map: Dict[int, int],
    reference_phases: Dict[int, complex] = REFERENCE_PHASES,
) -> Dict[int, Tuple[complex, complex, float]]:
    """
    Analyze the phase of each m-block.

    Args:
        w0: complex array w0[kz,ig,m]
        m_map: mapping from code m-index → |m| (e.g. {0:0, 1:1, 2:1, 3:2, ...})
        reference_phases: expected phase per |m| (canonical or PWCOND gauge)

    Returns:
        result[m_index] = (phi_est, phi_ref, delta_angle)
    """
    nkz, nig, nm = w0.shape
    result = {}

    for m_idx in range(nm):
        if m_idx not in m_map:
            continue
        abs_m = m_map[m_idx]
        if abs_m not in reference_phases:
            continue

        phi_ref = reference_phases[abs_m]
        # Flatten over kz,ig for this m
        z = w0[:, :, m_idx].ravel()
        phi_est = estimate_phase(z)
        delta = phase_difference(phi_est, phi_ref)
        result[m_idx] = (phi_est, phi_ref, delta)

    return result


def print_phase_report(result: Dict[int, Tuple[complex, complex, float]]) -> None:
    """
    Print a human-readable report of phase differences per m-index.
    """
    print("Phase analysis per w0(:,:,m):")
    print(" m  |  phi_est (num)      |  phi_ref (theory)   |  Δangle [deg]")
    print("----+----------------------+---------------------+-------------")

    for m_idx in sorted(result.keys()):
        phi_est, phi_ref, delta = result[m_idx]
        if np.isnan(delta):
            deg = "NaN"
        else:
            deg = f"{np.degrees(delta):7.3f}"

        print(
            f" {m_idx:2d} | "
            f"{phi_est.real:+6.3f}{phi_est.imag:+6.3f}i | "
            f"{phi_ref.real:+6.3f}{phi_ref.imag:+6.3f}i | "
            f"{deg:>11}"
        )


# ----------------------------------------------------------------------
# Main CLI
# ----------------------------------------------------------------------

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Check PWCOND four.f90 projector phases against theory."
    )
    parser.add_argument(
        "debug_file",
        type=str,
        help="Path to ASCII debug file written by four.f90 (kz ig m Re Im).",
    )
    parser.add_argument(
        "--gauge",
        choices=["canonical", "pwcond"],
        default="canonical",
        help="Which reference phase pattern to use.",
    )
    args = parser.parse_args()

    if args.gauge == "canonical":
        ref = EXPECTED_PHASES
    else:
        ref = PWCOND_GAUGE_PHASES

    path = Path(args.debug_file)
    w0 = load_w0_debug_ascii(path)

    # Example: mapping code m-index → |m|
    # For s,p,d,f with PWCOND-like ordering:
    #   m=0: s or z-like,
    #   m=1,2: p_-x, p_-y (|m|=1),
    #   m=3..5: d-states,
    #   etc.
    #
    # You MUST adapt this mapping to your actual w0 layout.
    m_map_example = {
        0: 0,  # s or z
        1: 1,  # cos φ  (p_-x)
        2: 1,  # sin φ  (p_-y)
        # Add entries for d and f as needed:
        # 3: 0, # d_z^2-1
        # 4: 1, # d_-xz
        # 5: 1, # d_-yz
        # 6: 2, # d_{x^2-y^2}
        # 7: 2, # d_{xy}
        # ...
    }

    result = analyze_w0_phases(w0, m_map_example, reference_phases=ref)
    print_phase_report(result)


if __name__ == "__main__":
    main()
