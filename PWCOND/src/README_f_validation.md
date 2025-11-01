# PWCOND f-orbital (l=3) validation plan

Goal: verify that adding l=3 (f) in PWCOND yields complex band structures (CBS) consistent with pw.x band structures when using the same pseudopotentials and k-paths.

## 1) Quick checks
- Rebuild QE with PWCOND.
- Run an existing s/p/d PWCOND example; ensure identical results to reference.
- Run a case with f-electron pseudopotentials and check PWCOND completes without runtime errors.

## 2) Mapping sanity (optional)
- If your UPF generator or pw.x uses a different f-order than the analytical order used in four.f90, set:
  `map_f_in(k_in) = k_target` in `PWCOND/src/four.f90` (search for `map_f_in`). Default is identity.

Analytical internal order used by four.f90:
1. m=0 → z(5z^2-3r^2)
2. cos φ → x(5z^2-r^2)
3. sin φ → y(5z^2-r^2)
4. cos 2φ → z(x^2-y^2)
5. sin 2φ → xyz
6. cos 3φ → x(x^2-3y^2)
7. sin 3φ → y(3x^2-y^2)

## 3) Numerical regression
- Run PWCOND CBS for the lead with k∥ on a symmetry line (e.g., Γ–M).
- Compare:
  - s/p/d CBS: unchanged vs. baseline.
  - f CBS: continuous dispersion and symmetry-consistent degeneracies.

## 4) pw.x vs. pwcond alignment
- Do an NSCF band calculation in pw.x along the same in-plane k-path and compare at the same energies:
  - Propagating branches: real kz from pw.x should agree with real branches of PWCOND CBS.
  - Evanescent branches: check onset thresholds in gaps.
- If a discrepancy is isolated to a single |m| family, adjust `map_f_in` and re-test.

## 5) Diagnostics (optional)
- Temporarily add prints in the f-block of `four.f90` to log (l,m), selected slot, and normalization for a few G-shells.
- Sweep `ecut2d` and z-slab resolution slightly to verify numerical stability.

## Troubleshooting checklist
- Confirm f pseudopotential channels and (l,m) indexing with `ld1.x` output (or UPF details).
- Ensure lead boundaries (bd1/bd2) do not truncate within beta radii.