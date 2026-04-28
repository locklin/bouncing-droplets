# Bouncing Droplet Born Rule & Khrennikov Correction Notes

## Date: 2026-04-21

## Context

We have a working Level 0 bouncing droplet simulator in
`claude-general/bouncing-droplets/` (C + raylib). It shows a droplet
bouncing inside a corral (circle, Bunimovich stadium, or rectangle),
with a wave field display and position histogram. The question: what
is the analog of the Born rule in this system, and can we detect
higher-order corrections predicted by Khrennikov's prequantum theory?

---

## The Born Rule Analog (Brady & Anderson, Section 5.7)

The chain of reasoning from the paper:

1. The wave field near a slow droplet: h = R cos(theta - omega_0 t),
   where R = amplitude envelope, theta = phase.

2. Droplet velocity: v = (c^2 / omega_0) * grad(theta). The droplet
   follows the phase gradient — identical to Bohmian mechanics where
   v = grad(S)/m and psi = R * exp(iS/hbar).

3. The imaginary part of the Schrodinger-like equation gives:
   dR^2/dt + div(R^2 * v) = 0
   This is a continuity equation: R^2 is a conserved density flowing
   with the droplet.

4. Therefore if initial position is distributed as R^2, it stays
   distributed as R^2. This is the Born rule: P(x) = |psi|^2 = R^2.

## What the histogram shows

The histogram counts where the droplet *has been* — a time-averaged
occupation density. The Born rule prediction is R^2(x,y) at any
instant. These are the same if the system is ergodic (time average =
ensemble average).

In chaotic regimes (high memory, especially in corrals), experiments
show the histogram does converge to something matching quantum
eigenmode |psi|^2. This is the famous Couder/Harris result.

## What we could compute

Since our wave field (the stroboscopic Bessel sum) is effectively the
amplitude envelope R, wave^2 at each grid point is the instantaneous
Born rule prediction. Time-averaging wave^2 gives the Born rule
prediction for the histogram.

---

## The Khrennikov Correction

### Background

Khrennikov's PCSFT (prequantum classical statistical field theory)
says quantum mechanics emerges from a classical random field, and the
Born rule is only a 2nd-order approximation:

    P(x) = |psi(x)|^2 + epsilon * |psi(x)|^4 + ...

Key predictions:
- epsilon is small
- Only even-order corrections (no h^3, h^5, etc.)
- The correction is detectable only in asymmetric potentials

### Why this system is a natural testbed

The wave field h(x,y,t) is deterministic (superposition of Bessel
functions from known bounce positions). But after many bounces in a
chaotic regime, it *looks* like a random field because memory of
distant past bounces creates an effectively stochastic background.
This is deterministic chaos masquerading as randomness — exactly the
type of "sub-quantum" random field Khrennikov envisions.

### Why only even orders

The wave equation is linear and J_0 is even. The field h has no
preferred sign. Over long times in the chaotic regime, the
distribution of h at any point is symmetric around zero:

- <h>   = 0  (zero mean)
- <h^2> != 0 (Born rule — the "energy density")
- <h^3> = 0  (odd moment of symmetric distribution)
- <h^4> != 0 (Khrennikov correction)
- <h^5> = 0  (odd again)

The even-only structure falls out of the wave equation symmetry
without needing to impose it.

### The Gaussianity subtlety

For a Gaussian random field, <h^4> = 3 * <h^2>^2 (kurtosis = 3).
So a truly Gaussian background makes h^4 degenerate with h^2 squared.

The interesting signal is departures from Gaussianity, which happen
in confined geometries where eigenmodes create structure. The excess
kurtosis at each point measures this:

    kappa(x) = <h^4>(x) - 3 * <h^2>(x)^2

So the proper test is:

    P(x) ~ a * <h^2>(x) + b * (<h^4>(x) - 3 * <h^2>(x)^2)

where:
- a * <h^2> is the Born rule
- b * (excess kurtosis) is the Khrennikov correction
- b should be ~0 for symmetric geometries (circle)
- b should be nonzero for asymmetric geometries (stadium, asymmetric corrals)

---

## Implementation Plan

### What to accumulate (per grid point, per bounce)

```c
h2_sum[iy][ix] += wave[iy][ix] * wave[iy][ix];
h4_sum[iy][ix] += wave[iy][ix] * wave[iy][ix] * wave[iy][ix] * wave[iy][ix];
n_samples++;     // same for all grid points (= total bounces)
```

Then:
- <h^2>(x) = h2_sum / n_samples
- <h^4>(x) = h4_sum / n_samples
- excess_kurtosis(x) = <h^4> - 3 * <h^2>^2

### What to display

Option A: Third panel showing <h^2> or the excess kurtosis
Option B: Overlay comparing histogram vs Born rule prediction
Option C: Residual plot: histogram - a * <h^2>, see if it correlates
with excess kurtosis

### Statistical test

After N bounces:
1. Normalize histogram and <h^2> to sum to 1 over interior points
2. Fit: histogram ~ a * <h^2> (linear regression, single coefficient)
3. Compute residual: r = histogram - a * <h^2>
4. Fit: r ~ b * excess_kurtosis
5. Is b significantly nonzero? Is R^2 improved?
6. Compare b across geometries (circle vs stadium vs asymmetric)

### Geometries to test

- Circle (symmetric) — expect b ≈ 0
- Bunimovich stadium (symmetric but chaotic) — might see b != 0
- Asymmetric corral (e.g., quarter-circle, D-shape) — best chance for b != 0
- Rectangle (integrable, not chaotic) — control case

Consider adding a D-shaped (half-circle) or asymmetric stadium geometry
to maximize the chance of seeing the correction.

---

## Questions to think about

1. How many bounces are needed for the statistics to converge? The
   histogram needs O(GRID^2) ~ 40,000 visits for decent coverage.
   At 3 bounces/frame * 30 fps, that's ~20 minutes of wall time, or
   ~7 minutes at speed=10. Probably want speed=50 and let it run.

2. Is the wave field the right object, or should we use the wave
   amplitude at the droplet position specifically? The field at the
   droplet position determines the kick, so it's more physically
   relevant than the field far from the droplet.

3. The memory parameter mu controls how "random-looking" the field
   becomes. Higher mu = more memory = more wave structure = more
   non-Gaussian. There might be an optimal mu for seeing the
   correction.

4. The Bessel-kernel model (Level 0) doesn't enforce boundary
   conditions on the wave field. The eigenmode model (Gilet) does.
   The correction might be more visible with proper eigenmodes since
   the boundary conditions create the asymmetry. Worth implementing
   the eigenmode version for a controlled test.

---

## References

- Brady & Anderson (2014), arXiv:1401.4356, Sections 5.6-5.7
- Khrennikov, "Beyond Quantum" (2014), Cambridge University Press
- Khrennikov, "Prequantum classical statistical field theory" series
- Gilet (2014), Phys. Rev. E 90, 052917 (eigenmode corral model)
- Harris, Bush et al. (2013) — corral histogram matching eigenmodes

## Code location

- Simulator: claude-general/bouncing-droplets/main.c
- Feasibility report: claude-general/docs/bouncing-droplet-simulation.md
- These notes: claude-general/docs/bouncing-droplet-born-rule-notes.md
