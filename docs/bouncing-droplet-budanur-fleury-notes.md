# Budanur & Fleury (2019) — State Space Geometry of Chaotic Pilot-Wave

## Citation
Budanur, N.B. & Fleury, M. (2019). "State space geometry of the chaotic
pilot-wave hydrodynamics." Chaos 29, 013122. arXiv:1812.09011.

Budanur is from Cvitanovic's group (IST Austria / Georgia Tech lineage).

---

## The Model (Oza + Harmonic Force)

The trajectory equation (eq 27):

    m x_ddot + D x_dot = -m*g * grad h(r, tau) - k*r

where k*r is the harmonic confining force. The wave field (eq 28):

    h(r, tau) = A * sum_{n} J_0(kF |r - r(n*TF)|) * exp(-(tau - n*TF) / (TF*Me))

Key: they decompose h into a Fourier-Bessel series (eq 30):

    h(r, tau) = A * sum_{n=0}^{inf} (2 - delta_{n,0}) * J_n(kF*r)
                * [C_n(tau) cos(n*theta) + S_n(tau) sin(n*theta)]

where the mode amplitudes evolve as ODEs (eqs 33-34):

    C_n_dot = -C_n / (TF*Me) + J_n(kF*r) * cos(n*theta) / TF
    S_n_dot = -S_n / (TF*Me) + J_n(kF*r) * sin(n*theta) / TF

**This is the key insight we were missing.** Instead of eigenmodes of the
cavity, the wave field is decomposed into Bessel-Fourier modes in polar
coordinates centered at the origin. No boundary needed — the harmonic
potential confines the droplet. Each mode evolves as a simple linear ODE
with exponential decay + source from the droplet position.

## Dimensionless Form

After scaling lengths by kF and time by TF (eqs 44-47):

    x_ddot = -eta * x_dot - chi * x
             - mu * sum (2-delta) [J_{n-1}(r) K_n cos(theta) - n*J_n(r)/r * L_n]

    y_ddot = -eta * y_dot - chi * y
             - mu * sum (2-delta) [J_{n-1}(r) K_n sin(theta) + n*J_n(r)/r * M_n]

    C_n_dot = -Me^{-1} * C_n + J_n(r) cos(n*theta)
    S_n_dot = -Me^{-1} * S_n + J_n(r) sin(n*theta)

where:
    K_n = C_n cos(n*theta) + S_n sin(n*theta)
    L_n = C_n cos((n+1)*theta) + S_n sin((n+1)*theta)
    M_n = C_n cos((n+1)*theta) - S_n sin((n+1)*theta)

Parameters from experiment (eq 48, 55):
    eta = D*TF/m = 0.2        (drag)
    chi = k*TF^2/m = 0.008    (spring constant)
    mu = g*A*kF^2*TF^2 = 0.0375...  (wave-droplet coupling)
    Me in [10, 19]             (memory parameter)

## State Space

The full state vector (eq 49):

    a = (x, y, x_dot, y_dot, C_0, C_1, S_1, C_2, S_2, ...)

Truncating at N=25 Fourier modes gives 2N+5 = 55 dimensions.
They use scipy odeint with adaptive stepping, delta_tau = 0.01.

## Symmetry Reduction (The Cvitanovic Slice Method)

The system has continuous rotation symmetry (SO(2)). They use the
"first Fourier mode slice" method (Budanur et al. 2015):

1. Define a "slice template" a' = (x=0, y=0, vx=1, vy=0, h=0)
2. Define rotation angle: phi_hat = Arg(vx + i*vy) (velocity direction)
3. The symmetry-reduced state: a_hat = g(-phi_hat) * a
   (rotate so velocity always points in +x direction)

In reduced coordinates, v_y = 0 always, so the system drops from
55D to 54D (or effectively 53D since v_y is trivially zero).

**Relative equilibria** (REQ): droplets on circular orbits become fixed
points in the reduced space.

**Relative periodic orbits** (RPO): ovals, lemniscates, trefoils become
periodic orbits in the reduced space.

## Invariant Solutions Found (Table I)

| Name    | Type | Shape       | Me  | Notes |
|---------|------|-------------|-----|-------|
| REQ0c   | REQ  | Circle      | 10  | Stable for Me < 13.8 |
| RPO1o   | RPO  | Oval        | 15  | Period ~40, stable until Me~16.9 |
| RPO4os  | RPO  | Oval (p4)   | 18  | Period-4, stable |
| RPO4ou  | RPO  | Oval (p4)   | 18  | Period-4, unstable |
| PO2l    | PO   | Lemniscate  | 14  | Figure-8, angular momentum reversal |
| RPO2l   | RPO  | Lemniscate  | 16  | |
| RPO1t   | RPO  | Trefoil     | 14  | |

## Bifurcation Sequence

1. Me ~ 10: Stable circular orbit (REQ0c)
2. Me ~ 13.8: Hopf bifurcation → circular orbit destabilizes
3. Me ~ 14-15: Oval orbits (RPO1o) appear via supercritical Hopf
4. Me ~ 16.9: Neimark-Sacker bifurcation → invariant 2-torus
5. Me ~ 17.1-17.2: Saddle-node bifurcation on torus → period-4 orbits
6. Me ~ 17.8-17.9: Period doubling cascade → chaos
7. Me ~ 19: Global attractor, angular momentum sign changes

## Numerical Methods

- Integration: scipy.integrate.odeint (adaptive), error < 10^{-8}
- Periodic orbit finding: Newton's method on the Poincare section
- Continuation: pseudo-arclength continuation for tracking bifurcations
- Floquet analysis: eigenvalues of reduced Jacobian for stability
- All in Python

## Connection to Periodic Orbit Theory / Gutzwiller

The paper doesn't explicitly use the Gutzwiller trace formula, but the
entire framework is set up for it:

1. The symmetry reduction gives a proper quotient space
2. Periodic orbits in the reduced space are cleanly identified
3. Floquet multipliers (stability eigenvalues) are computed
4. The chaotic attractor is described via its embedded periodic orbits

The next step (which they don't take but Cvitanovic would) is to use
the periodic orbit expansion / cycle expansion to compute:
- Lyapunov exponents as weighted averages over periodic orbits
- Escape rates
- Long-time statistical averages (the "dynamical zeta function")

This is exactly where the Gutzwiller trace formula enters for the
quantum/semiclassical analog.

---

## Implications for Our Simulator

### For oza.c (Level 1):
The Budanur-Fleury model (eqs 44-47) is the RIGHT implementation of
the Oza model. Key differences from what we have:

1. **Fourier-Bessel mode decomposition** instead of summing over past
   bounce positions. This is equivalent but much more efficient:
   - Our approach: O(N_memory) Bessel evaluations per time step
   - Their approach: O(N_modes) simple ODE updates per time step
   - With N_modes = 25 and N_memory = 512, theirs is ~20x faster

2. **Harmonic confining potential** instead of a hard corral. This gives
   rotation symmetry and enables the full Cvitanovic toolbox.

3. **Continuous-time ODE** instead of discrete bounces. More physical
   and enables adaptive time stepping.

### New binary ideas:

1. **`harmonic.c`** — Oza model with harmonic confinement, Fourier-Bessel
   decomposition. This is the "proper" Level 1 that can reproduce all
   of Budanur & Fleury's results. The confining potential replaces the
   corral walls — no boundary issues.

2. **`billiard.c`** — Classical billiard (hard ball bouncing in corral,
   no waves). Compare trajectory statistics with the pilot-wave case.
   Use DynamicalBilliards.jl approach or write from scratch in C.

3. **`schrodinger.c`** — Solve Schrodinger equation in the same geometry.
   Compare quantum eigenstate |psi|^2 with the pilot-wave histogram.

4. **`gutzwiller.c`** — Semiclassical quantization via the trace formula.
   Use periodic orbits from the billiard to approximate quantum spectrum.

### For the Khrennikov test:
The harmonic potential is radially symmetric, so the Khrennikov
correction should vanish by symmetry. But adding an asymmetric
perturbation (e.g., chi_x != chi_y, making an elliptical potential)
breaks the symmetry while keeping the mode decomposition tractable.

---

## DynamicalBilliards.jl

- Julia package by George Datseris (JuliaDynamics)
- v4.1.0, last updated June 2024, 993 commits
- Supports arbitrary 2D billiard geometries
- Computes Lyapunov exponents, boundary maps
- No quantum/semiclassical features
- Could be used for the classical billiard comparison

GitHub: https://github.com/JuliaDynamics/DynamicalBilliards.jl
JOSS paper: https://joss.theoj.org/papers/10.21105/joss.00458

## ChaosBook Resources

- Chapter 9: Billiards (classical dynamics, stability)
- Chapter 37-39: WKB, semiclassical evolution, Gutzwiller trace formula
- Chapter 42: Helium atom via periodic orbit theory
- Computational tools in the "extras" archive

URL: https://chaosbook.org/paper.shtml#billiards

---

## Code / File Locations

- Level 0 simulator: claude-general/bouncing-droplets/main.c (walker)
- Level 1 Oza model: claude-general/bouncing-droplets/oza.c
- Feasibility report: claude-general/docs/bouncing-droplet-simulation.md
- Born rule / Khrennikov notes: claude-general/docs/bouncing-droplet-born-rule-notes.md
- This file: claude-general/docs/bouncing-droplet-budanur-fleury-notes.md
