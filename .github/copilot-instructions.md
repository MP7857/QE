# Copilot Instructions for PWCOND Ylm / four.f90

This repository contains a modified `four.f90` (in Quantum ESPRESSO's PWCOND)
with extended support for **f-orbitals (l = 3)** and a carefully defined
phase convention for the 2D Fourier transform of nonlocal projectors.

The key physics invariants that **MUST NOT** be violated by any automatic edits
(Copilot, refactoring, etc.) are:

1. **UPF Ylm definition (baseline)**  
   - UPF pseudopotentials define the angular part of each projector as the
     standard complex spherical harmonics \( Y_l^m(\theta,\phi) \) in the
     Condon–Shortley convention.
   - The radial functions stored in UPF are \( \chi_{li}(r) \), and the
     full projector is
     \[
     \beta_{lmi}(\mathbf{r}) = \frac{\chi_{li}(r)}{r}\,Y_l^m(\hat{\mathbf{r}}).
     \]
   - **Do not** introduce ad-hoc phase factors or change this assumption.

2. **Real spherical harmonics**  
   - Internally, QE/PWCOND uses a **real** spherical-harmonics basis,
     constructed from complex \( Y_l^m \) as cos/sin combinations.
   - The canonical real-harmonic phase pattern for the 2D Fourier transform
     is (for each m-block, ignoring overall orbital sign):
     \[
     m=0:\ +1,\quad
     m=1:\ -i,\quad
     m=2:\ -1,\quad
     m=3:\ +i.
     \]
   - This comes from the exact plane-wave expansion
     \[
     e^{-i \mathbf{g}\cdot\mathbf{r}_\perp}
       = \sum_m (-i)^m J_m(g r_\perp)\,e^{i m(\phi-\phi_g)}.
     \]

3. **Separation of concerns: radial vs. angular**  
   - The radial part is handled by Simpson integration over `r` via arrays
     `x1..x6` → `fx1..fx6`.
   - **Do not mix physics and numerics**: phase corrections and Ylm-related
     sign changes are only allowed in the **final angular assembly**, not in
     the radial integrals themselves.

4. **PWCOND gauge vs. canonical theory**  
   - The canonical real-harmonic pattern above may be implemented directly
     or with a systematic "gauge" (e.g. global `-1` factor for all odd-m
     orbitals).
   - Any change to `four.f90` MUST preserve a clearly documented mapping
     between:
       - UPF complex \( Y_l^m \),
       - real spherical harmonics (cos/sin),
       - cubic orbitals used in PWCOND (e.g., `p_z`, `p_-x`, `d_-xz`,
         `f_{z(x^2-y^2)}`).
   - Copilot should **not** introduce new minus signs or swap orbitals in
     the `w0(kz,ig,m)` ordering unless explicitly requested and documented.

5. **m-block consistency**  
   - All orbitals sharing the same \(|m|\) must share the same Bessel order
     \(J_m\), the same phase factor \((-i)^m\), and differ only by cos/sin
     or polynomial prefactors in \(z\) and \(r_\perp\).
   - For example, for a given \(l\):
     - `w0(:,:,2)` and `w0(:,:,3)` (cos/sin) in the m=1 block must have the
       **same complex phase** up to a sign.
     - Analogous for m=2 and m=3.

6. **Testing / diagnostics**  
   - Use the script `tools/check_pwcond_phases.py` as a basis to check that
     the numerical output from `four.f90` is consistent with the theoretical
     phase pattern and with UPF radial functions.
   - If Copilot adds tests or diagnostics, they must:
     - Distinguish clearly between **phase** errors (Ylm, angular gauge)
       and **quadrature** errors (integration/endpoints).
     - Never silently "fix" phase by adjusting the radial integrals.

## What Copilot should do

- When editing `four.f90`, preserve:
  - The mapping between `lb` (l) and the number/order of m-channels.
  - The interpretation of `x1..x6`, `fx1..fx6` as distinct radial integrals.
  - The separation of:
    1. radial integration over `r`,
    2. angular assembly via `cs`, `sn`, `cs2`, `sn2`, `cs3`, `sn3`,
    3. final normalization via `s1..s4`.
- When asked for "improvements", Copilot should suggest:
  - Better diagnostics (printing out phases, writing debug files).
  - Additional Python checks around `w0(kz,ig,m)`.
  - Cleaner code comments documenting the Ylm and phase conventions.
- Do **not**:
  - Randomly flip signs,
  - Change the order of `w0(:,:,m)` channels,
  - Merge different m-channels into a single array.

## Reference goals

The main goals for this repository are:

1. Maintain a **rigorous, documented mapping** between UPF \(Y_l^m\) and PWCOND
   real projectors, including f-orbitals.
2. Provide automated checks (via `tools/check_pwcond_phases.py` and/or CI)
   for:
   - Phase correctness per m-block.
   - Stability of the computed projectors vs. `nz1`, `ewind`, and other
     numerical parameters.
3. Enable safe iteration with Copilot without breaking the physics.
