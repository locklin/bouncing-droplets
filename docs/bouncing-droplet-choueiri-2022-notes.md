# Choueiri et al. (2022) — Crises and Chaotic Scattering in Pilot-Wave Experiments

## Citation
Choueiri, G., Suri, B., Merrin, J., Serbyn, M., Hof, B. & Budanur, N.B.
(2022). "Crises and chaotic scattering in hydrodynamic pilot-wave
experiments." Chaos 32, 093138. arXiv:2206.01531.

Affiliations: IST Austria (1,3,4), U. Toledo (2), Indian Inst. Sci. (3),
Max Planck Institute for Physics of Complex Systems, Dresden (4).
Budanur is corresponding author — same as in the 2019 theory paper.

---

## What This Paper Adds Beyond Budanur & Fleury (2019)

The 2019 paper was **purely numerical** (simulations of the Oza model
with harmonic potential). This 2022 paper is **primarily experimental**
with the same symmetry-reduction analysis applied to real droplet data.
The key advance: they observe **crisis bifurcations** and **chaotic
scattering** in the experiment, phenomena that were hinted at but not
fully characterized in the 2019 simulations.

| Aspect | Budanur & Fleury 2019 | Choueiri et al. 2022 |
|--------|----------------------|---------------------|
| Data | Simulation (Oza model) | Experiment (real droplets) |
| Geometry | Harmonic potential | Circular corral (concentric cylinders) |
| Me range | 10-19 | 20-34 |
| Symmetry reduction | Yes (first Fourier mode slice) | Yes (same method, applied to measured trajectories) |
| Main finding | Bifurcation sequence to chaos | Crisis bifurcations + chaotic scattering |
| Periodic orbits | Found by Newton's method | Inferred from trajectory classification |
| Wave field | Fully resolved (Fourier-Bessel) | Not measured (only droplet position) |
| State space dim. | 55D (complete) | 4D projection (x, y, vx, vy) |

## Experimental Setup

- Circular corral: diameter 19.95 mm, depth 6 mm, in silicone oil
- Concentric overflow region creates effective radial confining force
  (not exactly harmonic, but qualitatively similar)
- Vibration frequency f_0 = 75 Hz
- Droplet diameter D = 0.85 mm
- Faraday wavelength lambda_F = 5.3 mm
- Camera records at ~57 ms intervals, cubic spline interpolation for velocity
- Memory Me = (1 - gamma/gamma_F)^{-1}, varied from 20 to 34

## Symmetry Reduction

Same method as Budanur & Fleury (2019): rotate the coordinate system so
the velocity vector always points in the +x direction:

    x_tilde = R(-theta_hat) * x,  v_tilde = R(-theta_hat) * v
    where theta_hat = arg(vx + i*vy)

This maps all rotation-equivalent trajectories to the same curve.
After reduction, the circle becomes a point, the lemniscate becomes a
figure-8 in the reduced plane, etc.

Poincare section: intersections with the half-plane x_tilde * [1,0] = 0,
v_tilde * [1,0] > 0 (velocity in +x direction when crossing the y-axis).

## Key Experimental Results

### Orbit Types (Fig. 2)
Same zoology as 2019 theory but at different Me values:
- Me = 20: Circle (stable) AND Lemniscate (stable, coexisting!)
- Me = 26.48: Chaotic lemniscate
- Me = 31.48: Chaotic oval
- Me = 31.86: Full chaos (intermittent lemniscate + oval)

### Bifurcation Sequence (Fig. 3)
Orbit diagram created by sweeping Me up and down in steps of Delta_Me ≈ 0.36:

1. **Me ∈ [20, 22.5]**: Stable lemniscate (two fixed points on Poincare section)
2. **Me ≈ 22.5**: Lemniscate destabilizes → chaotic lemniscate attractor
3. **Me ≈ 26.6**: **Boundary crisis** — chaotic lemniscate loses stability,
   system jumps to circular orbit
4. **Me ∈ [26.6, 29.5]**: Stable circular orbit (returns to regular dynamics!)
5. **Me ≈ 29.5**: Circular orbit destabilizes → chaotic oval attractor
6. **Me ≈ 31.65**: **Symmetry-increasing crisis** — oval chaos merges with
   lemniscate repeller → full chaotic attractor with angular momentum
   sign changes

**Hysteresis**: Sweeping Me down from 33 to 20, the system stays on
circular orbits through the entire range [20, 26.6] where lemniscates
coexist — classic multistability.

### Crisis Bifurcations (Fig. 4)
Two distinct crises identified:

1. **Boundary crisis (Me ≈ 26.6)**: Chaotic lemniscate attractor intersects
   its basin boundary (the stable manifold of the circular orbit). The
   lemniscate chaos disappears and the system falls onto the circle.

2. **Symmetry-increasing crisis (Me ≈ 31.65)**: The chaotic oval attractor
   (confined to L > 0, angular momentum positive) suddenly expands to
   include L < 0 (negative angular momentum). The lemniscate repeller,
   which separates the two oval basins, becomes part of the attractor.
   This is a MERGER of two reflection-symmetric chaotic sets.

   The authors note this is DIFFERENT from standard crisis bifurcations
   (where a chaotic attractor collides with an unstable periodic orbit).
   Here, it collides with a chaotic repeller (the lemniscate set).
   This opens new theoretical questions for high-dimensional chaos.

### Chaotic Scattering (Fig. 5-6)
After the symmetry-increasing crisis (Me > 31.65):

- The droplet alternates between oval-shaped and lemniscate-shaped
  trajectory segments
- Angular momentum L switches sign during lemniscate episodes
- **Lifetime distributions are exponential** — the hallmark of chaotic
  scattering from nonattracting chaotic sets (repellers)
- Lemniscate lifetimes: T_w ≈ 105 T_F, roughly constant across Me values
- Oval lifetimes: T_w decreases from 105 to 38 T_F as Me increases from
  32.1 to 33.6 (ovals become more transient farther from crisis)

The dynamics can be understood as **consecutive scatterings between
two chaotic repellers** (oval and lemniscate), connected by transition
trajectories that pass through the lemniscate set.

## Connection to Turbulence

The authors explicitly draw the parallel to transitional turbulence:
- In pipe flow, turbulence appears via formation of chaotic repellers
  (transient turbulent puffs) that eventually merge into sustained
  turbulence via crisis-like bifurcations
- The pilot-wave system provides a low-dimensional experimental analog
  where the same phenomena can be studied systematically
- References to Budanur's turbulence work (refs 22, 23, 27, 29, 55-61)

## What's Missing / Open Questions

1. No wave field data — they only track the droplet position. The full
   infinite-dimensional state (including the wave field) is not measured.
   Their 4D Poincare section is a projection of the full state space.

2. No periodic orbit computation from experiment. The orbits are
   classified by trajectory shape (circle, oval, lemniscate) but not
   found precisely by Newton's method as in the 2019 theory paper.

3. The confining potential is not exactly harmonic — the concentric
   cylinder geometry creates a more complex effective potential.

4. No Lyapunov exponent measurements reported.

5. No code or data publicly available.

---

## Implications for Our Simulator

### What we can already reproduce
Our `harmonic` binary implements the same Oza model as the 2019 paper
with the Fourier-Bessel decomposition. By varying Me from 10 to 19, we
should see the same bifurcation sequence (circles → ovals → lemniscates
→ chaos). The 2022 paper's experimental Me range (20-34) is higher,
partly because the corral geometry is different from a harmonic potential.

### What we should add

1. **Orbit diagram**: Sweep Me automatically, record Poincare section
   intersections, plot y_tilde_P vs Me. This directly reproduces Fig. 3.
   Could be a new mode in `harmonic` or a separate diagnostic tool.

2. **Symmetry reduction**: Implement the Budanur-Fleury velocity-angle
   slice in `harmonic`. Rotate (x, y, vx, vy) so vy = 0 at each display
   step. This maps circles to points and reveals the reduced dynamics.

3. **Angular momentum tracking**: Compute L = m(x*vy - y*vx) at each
   step, plot time series. This is the key observable for identifying
   lemniscate vs oval episodes and crisis transitions.

4. **Lifetime distributions**: In the chaotic regime, classify trajectory
   segments as "oval" (|L| > threshold) or "lemniscate" (|L| < threshold),
   measure episode durations, plot survival functions. Exponential tails
   confirm chaotic scattering.

5. **Corral confinement model**: To match the 2022 experiment more
   closely, replace the harmonic potential with the concentric-cylinder
   effective potential. This might shift the bifurcation values of Me
   to match the experimental range better.

### For the Khrennikov test
The crisis bifurcation at Me ≈ 31.65 creates a symmetry-breaking →
symmetry-restoring transition. Before the crisis, the attractor has
broken reflection symmetry (L > 0 only). After, it's symmetric.
This is exactly the type of symmetry change that the Khrennikov
correction depends on. Testing whether the excess kurtosis changes
character at the crisis point would be a clean experiment.

---

## Key References from the Paper

- Ref 18: Budanur & Fleury (2019) — the theory companion paper
- Ref 14: Bush & Oza (2020) — comprehensive review
- Ref 9: Harris et al. (2013) — corral wavelike statistics
- Ref 17: Perrard & Labousse (2018) — chaos in harmonic well
- Ref 42: Grebogi, Ott & Yorke (1982) — crisis bifurcations theory
- Ref 47: Chossat & Golubitsky (1988) — symmetry-increasing bifurcation
- Refs 22-29: Budanur's turbulence papers (pipe flow, Kolmogorov flow)

---

## File Locations

- This file: `docs/bouncing-droplet-choueiri-2022-notes.md`
- 2019 theory notes: `docs/bouncing-droplet-budanur-fleury-notes.md`
- Born rule notes: `docs/bouncing-droplet-born-rule-notes.md`
- Feasibility report: `docs/bouncing-droplet-simulation.md`
- Simulator code: `bouncing-droplets/` (walker, oza, harmonic, billiard, schrodinger)
