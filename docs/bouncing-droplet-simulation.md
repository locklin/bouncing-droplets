# Simulating Bouncing Droplet Pilot-Wave Hydrodynamics

## Feasibility Report and Implementation Plan

### 1. What Is the System?

A millimetric silicone oil droplet bounces on the surface of a vibrating bath of
the same liquid. A thin air film prevents coalescence. When the bath vibration
amplitude is tuned just below the Faraday instability threshold, each bounce
creates standing waves that persist between impacts — the "memory" regime.
Past a critical driving acceleration the droplet begins to self-propel laterally:
it becomes a "walker," surfing its own wave field.

The Brady & Anderson (2014) paper in this repository shows that this system
is, to a remarkable approximation, *Lorentz covariant* with respect to the
surface wave speed *c*, and that many quantum-mechanical phenomena emerge:
single- and double-slit diffraction, tunnelling, quantised orbits, spin-half
behaviour, Pauli exclusion, and an analogue of the Schrodinger equation where
Planck's constant is replaced by b = mc^2/omega.

The question: **can we simulate this from first principles, faithfully enough to
reproduce the quantum-like phenomenology, and at what computational cost?**

---

### 2. Hierarchy of Models (Cheapest to Most Expensive)

The literature has converged on a clear hierarchy. Each level trades physical
fidelity for computational tractability.

#### Level 0: Discrete Dynamical Maps

**What it is.** An iterated map x_{n+1} = f(x_n, wave_field_n), updated once per
bounce. The wave field is a superposition of eigenmodes of the container
geometry, excited at each impact point and decaying between bounces.

**Key references.**
- Gilet (2014, 2016, Phys. Rev. E) — circular corral
- Rahman & Blackmore (2020) — Neimark-Sacker bifurcations
- Rahman, Ghanem & Blackmore (2023) — elliptical corral, stochastic extensions

**Cost.** Seconds to minutes on a laptop per 10^6 bounces.

**Fidelity.** Captures long-time statistics, quantised orbits, chaotic dynamics in
corrals. Does not resolve the actual wave field shape or bouncing dynamics.
Container geometry must be chosen to have known eigenmodes.

**Code availability.** No public repositories found. Equations are simple enough
to implement from the papers (< 200 lines of Python/Julia).

#### Level 1: Stroboscopic Integro-Differential Trajectory Equation

**What it is.** The Oza-Rosales-Bush (2013) model. The vertical bouncing
dynamics are averaged out; the droplet moves continuously with an
integro-differential equation of motion:

  m x_ddot + D x_dot = -F grad h(x, t)

where h(x,t) is a superposition of Bessel-function wave kernels deposited at
every past position, exponentially decaying with a memory timescale T_M:

  h(x, t) = integral_0^t A J_0(k_F |x(t) - x(s)|) exp(-(t-s)/T_M) ds

**Key references.**
- Oza, Rosales & Bush, JFM 737 (2013) — foundational paper
- Turton, Couchman & Bush, Chaos 28 (2018) — comprehensive review
- Durey & Milewski, JFM 891 (2020) — spectral methods for corrals

**Cost.** Real-time to 10x real-time on a single core. A million bounces in
minutes. The integral must be truncated at ~50-100 memory times.

**Fidelity.** Captures walking, orbital quantisation, chaotic statistics in corrals,
tunnelling probabilities, single-slit diffraction. Does NOT capture:
- The actual bounce dynamics (vertical motion)
- Droplet deformation
- Air-layer dynamics
- Wave transients during impact
- Doppler effect details
- Multi-droplet interactions (requires extensions)

This is the "workhorse" model used for most quantum-analogy investigations.

**Code availability.** No public repo found, but the equations are well-specified.
Implementation is ~500 lines. NJIT (Oza's group) and MIT (Bush's group) have
internal MATLAB codes.

#### Level 2: Quasi-Potential Coupled Drop-Wave Model

**What it is.** The Milewski-Galeano-Rios-Nachbin-Bush (2015) model. Solves the
weakly viscous quasi-potential flow equations for the bath surface, coupled to a
vertical bouncing model (logarithmic spring or kinematic match) for the
droplet. The droplet is treated as a rigid sphere; the air layer is parameterised.

**Key references.**
- Milewski, Galeano-Rios, Nachbin & Bush, JFM 778 (2015) — the foundational paper
- Galeano-Rios, Milewski & Vanden-Broeck, JFM 826 (2017) — kinematic match
- Couchman, Turton & Bush, JFM 871 (2019) — improved wave damping

**Cost.** 10^2 to 10^4 times real-time per bounce. A typical simulation of ~100
bounces takes minutes to hours depending on resolution.

**Fidelity.** First model to capture:
- Rapid transient waves during impact
- Doppler effect in the wave field
- Walker-walker interactions
- Correct transition thresholds (bouncing -> walking -> chaos)

Does NOT capture droplet deformation or the air-layer dynamics.

**Code availability.** No public repo. Internal MATLAB/Fortran codes at Bath
(Milewski) and IMPA (Nachbin/Galeano-Rios).

#### Level 3: Lubrication-Mediated Reduced Model (2D)

**What it is.** The Phillips-Cimpeanu-Milewski (2025) model. Decomposes the
problem into three coupled regions — drop, air film, bath — each solved with
appropriate reduced equations. The air layer is resolved via thin-film
lubrication theory. Currently implemented in 2D (cylindrical droplet geometry).

**Key references.**
- Phillips, Cimpeanu & Milewski, Proc. R. Soc. A (2025) — 2D rebound model
- Galeano-Rios, Cimpeanu et al., JFM 958 (2023) — kinematic match DNS comparison

**Cost.** Several hours on a standard desktop (MATLAB) for multi-bounce
simulations with Faraday forcing. This is described as "prohibitively expensive
for DNS" in the equivalent geometry.

**Fidelity.** Captures:
- Air layer evolution (the physics that prevents coalescence)
- Correct rebound coefficient of restitution
- Multiple rebounds on a vibrating bath
- Low-Weber-number capillary dynamics

Does NOT capture 3D geometry, lateral self-propulsion (walking), or far-field
wave memory effects. Currently only models repeated vertical bouncing, not
walking.

**Code availability.** The BouncingDropletsLiquid2D repository on GitHub
(rcsc-group) contains both the MATLAB reduced model and Basilisk DNS code
for the 2D problem, including a Faraday forcing variant.

#### Level 4: Full DNS (Navier-Stokes with VOF)

**What it is.** Direct numerical simulation of the full two-phase (or three-phase:
drop + air + bath) Navier-Stokes equations using Volume-of-Fluid interface
tracking on adaptive meshes. The state of the art uses Basilisk (successor to
Gerris), developed by Stephane Popinet.

**Key references.**
- Alventosa, Cimpeanu & Harris, JFM 958 (2023) — inertio-capillary rebound
- rcsc-group/BouncingDroplets (GitHub) — axisymmetric DNS for drop-bath impact
- rcsc-group/BouncingDropletsMovingPool3D (GitHub) — 3D extension (2026)

**Cost.** Each individual bounce takes approximately **24 CPU-hours on 4-8
cores** for a related axisymmetric problem. A 3D simulation would be orders of
magnitude more expensive. Simulating 100 bounces of a walker in 3D with
full wave memory would require on the order of **10^4 to 10^5 CPU-hours** —
weeks on a cluster.

**Fidelity.** In principle, captures everything:
- Full Navier-Stokes flow in both fluids
- Interface deformation of both drop and bath
- Air layer dynamics
- Faraday wave excitation and memory
- Viscous dissipation in boundary layers

In practice, nobody has yet performed a full DNS of a *walking* droplet —
only individual bounce events. The multi-scale challenge (micron-scale air
layer vs centimetre-scale wave field vs seconds-long memory) makes this
extremely difficult.

**Code availability.**
- [rcsc-group/BouncingDroplets](https://github.com/rcsc-group/BouncingDroplets) — Basilisk C, axisymmetric, single bounce
- [rcsc-group/BouncingDropletsLiquid2D](https://github.com/rcsc-group/BouncingDropletsLiquid2D) — Basilisk C + MATLAB, 2D multi-bounce with Faraday forcing
- [rcsc-group/BouncingDropletsMovingPool3D](https://github.com/rcsc-group/BouncingDropletsMovingPool3D) — Basilisk C, 3D moving pool

#### Level 5: Full DNS with Resolved Air Layer (SPH)

**What it is.** Smoothed Particle Hydrodynamics applied to the full drop-bath
system including the air layer and walking dynamics.

**Key references.**
- Molteni & Lupo (2017, Computers & Fluids) — SPH simulation of walking droplet

**Cost.** Not reported in the paper. SPH is typically slower than grid-based VOF
for equivalent resolution but handles free surfaces more naturally.

**Fidelity.** Molteni & Lupo reported successful wave-droplet coupling and
constant-velocity walking, with "relevant simplifying assumptions." The
method is promising but has not been reproduced or extended.

**Code availability.** No public code.

---

### 3. Existing Open-Source Code

| Repository | Method | What it does | Language |
|---|---|---|---|
| [rcsc-group/BouncingDroplets](https://github.com/rcsc-group/BouncingDroplets) | Basilisk DNS (VOF) | Single axisymmetric bounce, non-coalescing | C |
| [rcsc-group/BouncingDropletsLiquid2D](https://github.com/rcsc-group/BouncingDropletsLiquid2D) | Basilisk DNS + MATLAB reduced model | 2D multi-bounce with Faraday forcing | C, MATLAB |
| [rcsc-group/BouncingDropletsMovingPool3D](https://github.com/rcsc-group/BouncingDropletsMovingPool3D) | Basilisk DNS (VOF) | 3D drop impact on moving pool | C |
| [simrandahia/Quantum-Walking-Droplet-Analogy](https://github.com/simrandahia/Quantum-Walking-Droplet-Analogy) | LBM (D3Q19) + MD | Student project (IPT2021), walking droplet analogy | C++ |
| [erkara/TrackingWalkers-YOLOv8](https://github.com/erkara/TrackingWalkers-YOLOv8) | Deep learning | Tracking walkers in experimental videos | Python |

**Notable absence:** There is no public implementation of the Oza-Bush
stroboscopic model or the Milewski quasi-potential model, despite these being
the most widely used in the literature. This represents a significant gap and
an opportunity.

---

### 4. Can We Simplify the Physics?

**Can the droplet be treated as rigid?** Yes — this is standard practice in
Levels 1-3. Molacek & Bush (2013) showed that the droplet deformation
contributes only a logarithmic correction to the impact force. The
"kinematic match" approach (Galeano-Rios et al. 2017) treats the droplet as
a rigid sphere that kinematically matches the bath surface during contact.
This works well for We < 1 (Weber number), which covers the walking regime.

**Can the liquid model be simplified?** Yes, substantially:
- The bath flow is well-approximated by quasi-potential flow (Level 2),
  saving the cost of resolving viscous boundary layers everywhere.
- The Faraday waves are nearly monochromatic (wavelength ~5 mm), so
  spectral methods work well.
- Viscous damping enters only through a surface damping term, not through
  resolving the full boundary layer.
- The air layer can be parameterised rather than resolved (Levels 1-2).

**What cannot be simplified?** The wave *memory* — the integral over all past
bouncing positions — is the essential physics that produces quantum-like
behaviour. Any model must track this. In the stroboscopic model, this is an
integral over O(50-100) Faraday periods. In DNS, this means simulating for
much longer than a single bounce.

---

### 5. Feasibility Assessment

#### Can we reproduce the quantum-like phenomena numerically?

| Phenomenon | Minimum model level | Status | Difficulty |
|---|---|---|---|
| Walking (self-propulsion) | Level 1 | Done (Oza 2013) | Easy |
| Orbital quantisation | Level 1 | Done (Oza 2014) | Easy |
| Single-slit diffraction | Level 1 | Done (Pucci 2018) | Medium |
| Double-slit "interference" | Level 1 | Disputed — see below | Hard |
| Tunnelling | Level 1 | Done (Nachbin 2017) | Medium |
| Corral statistics | Level 0-1 | Done (Gilet 2016, Durey 2020) | Medium |
| Spin-half behaviour | Level 2+ | Not yet simulated | Hard |
| Pauli exclusion analogue | Level 2+ | Not yet simulated | Very hard |
| Anderson localisation | Level 1 | Partially done | Hard |

#### The Double-Slit Controversy

The original Couder & Fort (2006) double-slit result — showing
interference-like statistics — has **not been reproduced**. In 2015, three
independent teams (Bohr/Andersen in Denmark, Bush at MIT, Batelaan at
Nebraska) attempted to replicate it and failed: droplets went through slits
in nearly straight lines with no interference pattern. Ellegaard & Levinsen
(Phys. Rev. E, 2024) performed a comprehensive phase-space study and
concluded the original result was an overinterpretation of data.

However, the Pucci et al. (2018, JFM) numerical simulations using the
stroboscopic model *do* show slit-dependent deflection statistics that
qualitatively match single-slit diffraction. The situation for double-slit
remains unclear and controversial.

A simulation of the double-slit experiment would therefore be scientifically
valuable, but should be understood as testing what the *model* predicts,
not necessarily what the physical system does.

---

### 6. Proposed Implementation Plan

#### Stage 0: Feasibility Demonstration (2-3 weeks, 1 person)

**Goal.** Implement the Oza-Bush stroboscopic model (Level 1) and reproduce
published results for walking, orbital quantisation, and slit diffraction.

**Deliverables.**
1. Python/Julia implementation of the integro-differential trajectory equation
2. Reproduce walking speed vs. memory parameter (Oza 2013, Fig. 3)
3. Reproduce quantised orbits in a harmonic potential (Oza 2014)
4. Single walker through a single slit with deflection histogram

**Cost.** One developer, 2-3 weeks. Zero hardware cost (runs on a laptop).

**Why this first.** The stroboscopic model is the entry point to the field. If we
cannot reproduce published results at this level, higher-fidelity simulations
are pointless. This also gives us a working simulator for designing experiments.

#### Stage 1: Double-Slit and Corral Experiments (2-4 weeks)

**Goal.** Use the Level 1 simulator to investigate the double-slit controversy
and corral dynamics systematically.

**Deliverables.**
1. Double-slit deflection statistics over thousands of droplet trajectories
2. Corral statistics (circular and elliptical geometries)
3. Tunnelling probability vs. barrier width
4. Comparison with theoretical predictions from Brady & Anderson

**Cost.** Negligible hardware. Main cost is developer time for setting up boundary
conditions (slits, corrals, barriers) in the wave field calculation.

**Technical challenge.** The wave field near boundaries is handled differently in
different papers. Slits require careful treatment of wave reflection/diffraction
off submerged barriers. The standard approach uses image sources or
Dirichlet-to-Neumann maps for the specific geometry.

#### Stage 2: Coupled Drop-Wave Model (1-2 months)

**Goal.** Implement the Milewski et al. (2015) quasi-potential model to capture
the full wave field and bouncing dynamics.

**Deliverables.**
1. 2D wave field solver (spectral or pseudospectral) coupled to bouncing model
2. Reproduce the Doppler effect in the wave field
3. Two-walker interactions (attraction, repulsion, orbiting)
4. Spin-half behaviour investigation (Brady & Anderson Section 6)

**Cost.** ~500 CPU-hours for systematic parameter sweeps. Requires careful
numerical implementation (stiff ODEs, spectral wave solver).

#### Stage 3: High-Fidelity DNS Benchmarking (2-3 months)

**Goal.** Use the Basilisk-based rcsc-group code to perform DNS of individual
bounces and validate reduced models.

**Deliverables.**
1. DNS of single bounce with Faraday forcing, compare to Level 2/3 models
2. DNS of 2-3 consecutive bounces, assess wave memory accumulation
3. Quantify the error introduced by each simplification level

**Cost.** ~1000-5000 CPU-hours on a cluster. Each bounce is ~24 CPU-hours
in axisymmetric geometry, much more in 3D.

---

### 7. Cost-Benefit Analysis

| Stage | Calendar time | CPU cost | Developer effort | Scientific value |
|---|---|---|---|---|
| 0: Stroboscopic model | 2-3 weeks | ~0 (laptop) | Medium | High — reproduces most quantum analogies |
| 1: Slit/corral experiments | 2-4 weeks | ~0 (laptop) | Medium | Very high — addresses open controversy |
| 2: Quasi-potential model | 1-2 months | ~500 CPU-hr | High | High — captures wave field and multi-droplet |
| 3: DNS benchmarking | 2-3 months | ~5000 CPU-hr | Very high | Medium — validates but unlikely to reveal new physics |

**The overwhelming conclusion is that Stage 0-1 (stroboscopic model +
experiments) provides 80% of the scientific value at 1% of the cost.** The
DNS (Stage 3) is primarily useful for validating the reduced models, not for
discovering new phenomenology.

---

### 8. Is Full DNS Feasible for a Double-Slit Experiment?

**Short answer: No, not with current resources.**

A double-slit experiment requires:
- Walking the droplet across ~50 wavelengths (total domain ~25 cm)
- Maintaining wave memory for ~50-100 Faraday periods (~1-2 seconds)
- Resolving the air layer (~1 micron) and the wave field (~5 mm wavelength)
- Running ~1000 trajectories for statistics

This gives a scale separation of ~10^5 (1 micron to 25 cm), requiring
extremely fine adaptive meshing. Each trajectory would need ~100 bounces.
At ~24 CPU-hours per bounce (axisymmetric, single bounce, no memory), the
cost would be:

  100 bounces x 1000 trajectories x ~100 CPU-hours/bounce (3D with memory)
  = 10^7 CPU-hours
  ~ 1000 GPU-years or $2-5 million in cloud compute

This is clearly infeasible. Even a single walking trajectory in full 3D DNS
with wave memory has never been done.

**Are we better off with tubs of vibrating oil?** For the double-slit question
specifically, **yes, absolutely.** A physical experiment costs perhaps $5000 in
equipment and produces thousands of trajectories in a few days. The Ellegaard
& Levinsen (2024) and Pucci et al. (2018) experimental campaigns demonstrate
this is practical.

However, **the stroboscopic model (Level 1) is competitive with experiment**
for many questions and has the advantage of perfect repeatability, easy
parameter variation, and the ability to test geometries that are hard to build
physically (e.g., exotic corral shapes, non-physical memory parameters, etc.).

---

### 9. Recommendation

**Start with the stroboscopic model (Level 1).** This is where the field lives.
It is the model used in the majority of published walking-droplet papers.
Implementation is straightforward, cost is negligible, and it provides access
to all the key quantum-analogy phenomenology.

The specific investigations from Brady & Anderson that would be most
interesting to simulate numerically:

1. **Diffraction patterns** (Section 5) — does the model reproduce the correct
   wavelength and angular distribution?
2. **Lorentz covariance** (Section 3) — can we verify the gamma^2(v^2/c^2 + 1/2)
   relationship computationally?
3. **Orbiting pairs and spin** (Section 6) — do orbiting droplets in the model
   show the predicted spin-half behaviour and Pauli exclusion?
4. **Inverse-square force** (Section 4) — can we measure the force law between
   walkers and boundaries numerically?

For items 1-2, the stroboscopic model suffices. For items 3-4, the
quasi-potential model (Level 2) is needed to capture the detailed wave field
that mediates the interactions.

---

### 10. Key References

**Foundational experimental papers:**
- Couder, Protiere, Fort & Boudaoud, Nature 437 (2005)
- Couder & Fort, Phys. Rev. Lett. 97 (2006)
- Eddi, Sultan, Moukhtar, Fort, Rossi & Couder, JFM 674 (2011)

**The paper in this repository:**
- Brady & Anderson, arXiv:1401.4356 (2014) — "Why bouncing droplets are a pretty good model of quantum mechanics"

**Stroboscopic model (Level 1):**
- Oza, Rosales & Bush, JFM 737 (2013)
- Molacek & Bush, JFM 727 (2013)
- Turton, Couchman & Bush, Chaos 28 (2018) — review

**Quasi-potential model (Level 2):**
- Milewski, Galeano-Rios, Nachbin & Bush, JFM 778 (2015)
- Couchman, Turton & Bush, JFM 871 (2019)

**DNS and reduced models (Levels 3-4):**
- Alventosa, Cimpeanu & Harris, JFM 958 (2023)
- Phillips, Cimpeanu & Milewski, Proc. R. Soc. A (2025)

**Discrete dynamical models (Level 0):**
- Gilet, Phys. Rev. E 90 (2014) and 93 (2016)
- Rahman & Blackmore, Chaos Solitons Fractals (2020)

**Double-slit controversy:**
- Pucci, Harris, Faria & Bush, JFM 835 (2018) — simulations + experiments
- Ellegaard & Levinsen, Phys. Rev. E 109 (2024) — comprehensive negative result
- Batelaan group (Nebraska, 2015) — failed replication

**Recent energetics and theory:**
- Durey & Bush, JFM (2025) — energetics of pilot-wave hydrodynamics
- Sankaran (MIT, 2024) — generalised tunnelling framework
- Darrow & Bush (2024) — Lagrangian pilot-wave framework

**Open-source code:**
- [rcsc-group/BouncingDroplets](https://github.com/rcsc-group/BouncingDroplets)
- [rcsc-group/BouncingDropletsLiquid2D](https://github.com/rcsc-group/BouncingDropletsLiquid2D)
- [rcsc-group/BouncingDropletsMovingPool3D](https://github.com/rcsc-group/BouncingDropletsMovingPool3D)
- [simrandahia/Quantum-Walking-Droplet-Analogy](https://github.com/simrandahia/Quantum-Walking-Droplet-Analogy)

**Resource pages:**
- [MIT Bush group](https://thales.mit.edu/bush/)
- [DotWave.org](http://dotwave.org/) — comprehensive walking droplet resource
- [UNC Physical Mathematics Lab](https://www.pml.unc.edu/hqas)
- [NJIT Walking Droplets](https://cfsm.njit.edu/walking_droplets.php)
- [Aminur Rahman (pilotwaves)](https://pilotwaves.github.io/)
