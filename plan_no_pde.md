# Undeniably-Superfluid Orbit Engine

**Compute orbits using only superfluid hydrodynamics; compare to GR-1PN in a separate module. No (G), no Newtonian ansatz.**

## 0) Scope & success criteria

* Integrate (N) bodies using **forces from superfluid surface integrals** (momentum & pressure flux) around each body.
* Optionally evolve mass via **intake flux** across the same surface.
* Provide an **incompressible baseline** (fast, robust) and a **small-Mach compressible extension** with proper **near-field renormalization**.
* Side-by-side **GR-1PN comparator** exists in a separate module for plotting/validation, never used to compute slab forces.
* Pass the acceptance tests in §10 (emergent (1/r^2); coefficient match at <0.5%; stable orbits; compressible precession finite & (c_s^{-2})-scaled).

---

## 1) Model & parameters

**Medium**

* (\rho_0) (background density)
* (c_s) (slab sound speed) — controls small-Mach corrections
* (\beta_0) → (Q/M = 1/\beta_0) (Solar-System: treat as constant)
* (Optional) (\gamma_\beta) for rarification: (\beta(\rho)=\beta_0(\rho/\rho_0)^{-\gamma_\beta}) (default (0))

**Bodies (for each (a))**

* (M_a), (\mathbf x_a), (\mathbf v_a)
* (Q_a = M_a/\beta_0) (negative for sinks)
* Control-surface radius (R_a) (choose small vs separations; fixed in time)

**Never import (G) or kg**. If you need “gravitational parameters”, compute/fit (\mu_a \equiv K M_a) with (K=\rho_0/(4\pi\beta_0^2)) only for reporting; do not use (\mu) to compute forces.

---

## 2) Fluid field (no grid)

Potential flow of point sinks:
[
\phi(\mathbf x)= -\sum_{b}\frac{s_b}{4\pi,|\mathbf x-\mathbf x_b|},\quad
\mathbf v(\mathbf x)=\nabla\phi=\sum_b \frac{s_b}{4\pi}\frac{\mathbf r_b}{r_b^3},\quad s_b\equiv Q_b/\rho_0.
]
This is the **only** field sampled.

---

## 3) Control-surface mechanics (the force you integrate)

### 3.1 Master surface integral

For a small sphere (\partial B_a) of radius (R_a) about body (a):
[
\boxed{\ \mathbf F_a
=\oint_{\partial B_a}!\Big[\rho,\mathbf v,(\mathbf v!\cdot!\hat{\mathbf n})
-\big(P+\tfrac12\rho v^2\big)\hat{\mathbf n}\Big],dA\ } \tag{1}
]

**Incompressible baseline** (fast & stable): set (\rho=\rho_0), drop the bracket ( -\big(P+\tfrac12\rho v^2\big)\hat{\mathbf n}) (it cancels for steady potential flow), and compute
[
\mathbf F_a^{(\text{inc})}=\rho_0\oint_{\partial B_a}\mathbf v(\mathbf v!\cdot!\hat{\mathbf n}),dA. \tag{2}
]

### 3.2 Decomposition and near-field identities

On (\partial B_a), split (\mathbf v=\mathbf v_{\rm self}+\mathbf v_{\rm ext}) where:

* (\mathbf v_{\rm self} = \dfrac{s_a}{4\pi R_a^2},\hat{\mathbf n}) (radial monopole),
* (\mathbf v_{\rm ext}) is the sum of all (b\neq a) contributions, smooth on the small sphere.

Useful integrals on a sphere (with area (4\pi R_a^2)):
[
\oint \hat{\mathbf n},dA=0,\qquad
\oint (\hat{\mathbf n}!\cdot!\mathbf A),\hat{\mathbf n},dA=\frac{4\pi R_a^2}{3},\mathbf A\ \ (\mathbf A\ \text{constant over the sphere}). \tag{3}
]

**Resulting emergent inverse-square (incompressible):**
[
\boxed{\ \mathbf F_a^{(\text{inc})}=\frac{4}{3},\frac{Q_a}{4\pi},\mathbf v_{\rm ext}(\mathbf x_a)\ } \tag{4}
]
and since (\mathbf v_{\rm ext}(\mathbf x_a)=\sum_{b\neq a}\dfrac{Q_b}{4\pi\rho_0}\dfrac{\mathbf r_{ab}}{r_{ab}^3}),
[
\boxed{\ \mathbf F_a^{(\text{inc})}
=\sum_{b\neq a}\frac{4}{3},\frac{Q_aQ_b}{4\pi\rho_0}\frac{\mathbf r_{ab}}{r_{ab}^3}\ } \tag{5}
]
which is exactly the coefficient we measured numerically (to (4\times10^{-4})).

**Implementation note:** you can evaluate (2) by quadrature (works), or **use (4) directly** (same physics, ~100× faster). Keep both paths behind a switch for auditability.

---

## 4) Small-Mach compressible extension (stable, renormalized)

Naïvely inserting (\rho=\rho_0\big(1-\tfrac{v^2}{2c_s^2}\big)), (P=c_s^2(\rho-\rho_0)) into (1) blows up because (\mathbf v_{\rm self}\propto R_a^{-2}). The fix is **near-field renormalization**:

* **Rule A (momentum term):** keep the **full** (\mathbf v) in (\rho,\mathbf v(\mathbf v!\cdot!\hat{\mathbf n})) — this contains the physical (v_{\rm self}!\times!v_{\rm ext}) cross-term that generates the force.
* **Rule B (thermodynamic bracket):** evaluate (\rho(v), P(v)) using **external velocity only**, (\mathbf v_{\rm ext}), i.e.
  [
  \rho^\star=\rho_0\Big(1-\frac{v_{\rm ext}^2}{2c_s^2}\Big),\quad
  P^\star=c_s^2(\rho^\star-\rho_0)=-\tfrac12\rho_0 v_{\rm ext}^2.\tag{6}
  ]
  This **subtracts the singular self-contribution** and implements the standard throat counterterm.

Then compute:
[
\boxed{\ \mathbf F_a^{(\text{comp})}
=\oint_{\partial B_a}!\Big[\rho^\star,\mathbf v,(\mathbf v!\cdot!\hat{\mathbf n})
-\big(P^\star+\tfrac12\rho^\star v_{\rm ext}^2\big)\hat{\mathbf n}\Big],dA\ .} \tag{7}
]

* In the (c_s!\to!\infty) limit, (\rho^\star!\to!\rho_0) and the bracket cancels, recovering (2)–(5).
* For finite (c_s), you get a small, finite correction (\propto v_{\rm ext}^2/c_s^2), hence perihelion precession (\propto c_s^{-2}).

**Fast path:** expand (7) to leading order in (\mathrm{Ma}^2) and evaluate analytically with identities (3); code both the **quadrature** and the **closed-form** correction; they must agree within noise.

---

## 5) Mass intake (optional)

Two interchangeable modes (both referenced to the same surface):

* **Parametric:** (\dot M_a^{\rm param}=M_a/\beta_0) (Solar-System: negligible).
* **Flux:** (\dot M_a^{\rm flux}=-\oint \rho^\star,\mathbf v!\cdot d\mathbf A) with (\rho^\star) from (6) (subtracts self divergence). Keep (M_a(t)) and (Q_a(t)=M_a/\beta_0) in sync.

---

## 6) Time integration

Use **velocity-Verlet** (leapfrog). For each step:

1. (If enabled) update (Q_a(t)) from (M_a(t)).
2. Compute (\mathbf v_{\rm ext}(\mathbf x_a)) from other bodies; (optional) compute full quadrature for audit snapshots only.
3. Get (\mathbf F_a) from (5) for incompressible, plus analytic compressible correction or evaluate (7) by quadrature (skip-steps recommended).
4. (\mathbf a_a=\mathbf F_a/M_a); leapfrog update (\mathbf v_a,\mathbf x_a).
5. (If enabled) update mass by intake flux every (N_{\text{intake}}) steps (tiny effect).
6. Diagnostics/logging.

**Performance knobs**

* Precompute Fibonacci normals for each (R_a).
* Use **analytic (5)** for 99% of steps; run a **quadrature audit** every (N_{\text{audit}}) steps (e.g., 100–1000) to confirm.
* Vectorize pairwise sums; cache (\mathbf v_{\rm ext}).

---

## 7) GR-1PN comparator (separate module)

* Implement EIH 1PN accelerations (harmonic gauge) with physical (c).
* Same initial conditions; integrate with the same leapfrog (or RK4).
* Post-process **perihelion precession** and **elements**.
* Never feed GR internals back into slab; used only for plots/validation.

---

## 8) File/Module layout

```
slab/
  geometry.py          # fibonacci_sphere(), sphere identities
  medium.py            # Medium(rho0, cs, beta0, gamma_beta)
  bodies.py            # Body dataclass; Q=M/beta; sync M↔Q
  field.py             # v_total(), v_ext(), v_self()
  surface.py           # force_incompressible(), force_compressible_renorm(), intake_flux()
  dynamics.py          # step_verlet(), assemble_forces()
  diagnostics.py       # osculating_elements(), peri_finder(), energy-like checks
  gr1pn.py             # eih_1pn_accel(); separate
  io_cfg.py            # load/validate config; unit checks
  run.py               # CLI entry
```

**Key APIs**

```python
def v_ext_at(xa, bodies, a_idx, rho0): ...
def force_incompressible(a_idx, bodies, medium): ...         # uses eq. (5)
def force_incompressible_quad(a_idx, bodies, medium, npts): ...  # audit via eq. (2)
def force_compressible_renorm(a_idx, bodies, medium): ...     # eq. (7), analytic O(Ma^2)
def force_compressible_quad(a_idx, bodies, medium, npts): ... # eq. (7), quadrature
def intake_flux(a_idx, bodies, medium, npts): ...             # with ρ* from (6)
def step_verlet(bodies, medium, dt, opts): ...
def eih_1pn_accel(bodies, c_light): ...
```

---

## 9) Configuration (YAML)

```yaml
medium:
  rho0: 1.0
  cs: 1.0e4
  beta0: 1.0e10
  gamma_beta: 0.0

bodies:
  - name: Sun
    M: 1.0
    x: [0,0,0]
    v: [0,0,0]
    R: 1.0e-3
  - name: Mercury
    M: 3.3e-7
    x: [a0, 0, 0]
    v: [0, v0, 0]
    R: 5.0e-4

numerics:
  dt: 2.0e-3
  steps: 200000
  audit_every: 500
  npts_audit: 512
  use_compressible: true
  use_flux_mass: false
  intake_every: 2000

compare_gr_1pn:
  enable: true
  c_light: 63239.7263    # AU/yr
  measure_peri: true

outputs:
  save_every: 1000
  write_csv: true
  plots: [orbit, precession_vs_time, force_decomp]
```

---

## 10) Tests & acceptance criteria

### 10.1 Emergent inverse-square (slab, incompressible)

* Place two equal sinks at separation (r), measure (|\mathbf F|) on one via **quadrature** and fit (C) in (F \approx C/r^2).
* **Pass if** (\big|C/C_{\rm theory}-1\big|<5\times10^{-3}), with
  [
  C_{\rm theory} = \frac{4}{3}\frac{|Q_1Q_2|}{4\pi\rho_0}.
  ]

### 10.2 Orbit sanity

* Two-body, low-(e): using **analytic (5)** forces, keep (a,e) constant over (>50) orbits (energy-like drift small).
* **Pass if** (|\Delta a|/a<10^{-5}) over 50 orbits at default (dt).

### 10.3 Compressible correction: finite & (c_s^{-2})

* Enable compressible renormalized force (7); run slightly eccentric orbit & extract perihelion precession per orbit.
* **Pass if** (a) precession is finite and (b) scales (\propto c_s^{-2}) when you vary (c_s) by a factor of 2 (slope within 10%).

### 10.4 Quadrature audit

* Every `audit_every` steps, compare analytic force to quadrature force (incompressible and compressible modes separately).
* **Pass if** relative difference (<10^{-3}) (incompressible) and (<5\times10^{-3}) (compressible).

### 10.5 GR-1PN comparator

* Independent module reproduces Mercury’s (\sim 43''/)century within a few %.
* Reports slab vs GR precession on the same initial orbit without sharing parameters.

---

## 11) Numerical & stability notes

* **Quadrature:** Fibonacci sphere (N\in[256,1024]) is smooth; cache normals per body.
* **Skip strategy:** do analytic forces every step; do quadrature audits every few hundred steps.
* **Self-field:** never evaluate (\rho,P) with the self velocity; use (\mathbf v_{\rm ext}) for those.
* **Units:** choose convenient code units (e.g., AU, yr, code mass); keep (K=\rho_0/(4\pi\beta_0^2)) as a reported diagnostic only.

---

## 12) What to log (for “it’s fluid, not GR”)

* For each body & audit step:

  * (\mathbf v_{\rm ext}(\mathbf x_a)), (\mathbf F^{\text{inc}}_a) (analytic & quad), (\mathbf F^{\text{comp}}_a) (analytic & quad).
  * Decompose the momentum integrand into **self–ext cross-term** and **ext–ext** piece.
  * Intake flux (\dot M_a) (both parametric and surface estimate).
* Plots: (F) vs (r) (log–log slope (-2)), coefficient (C), precession vs (1/c_s^2).

---

## 13) Pseudocode (core loop)

```python
for step in range(steps):
    # 1) Q from M
    for a in bodies: a.Q = a.M / medium.beta0

    # 2) External velocities at each body (vectorized pairwise sum)
    for a in bodies:
        v_ext[a] = sum_{b!=a} (a.Q, b.Q) -> -(Q_b/4π) * r_ab / |r_ab|^3

    # 3) Forces
    for a in bodies:
        F_inc = medium.rho0 * a.Q * v_ext[a]                  # control-surface lemma
        if use_compressible:
            F_comp = compressible_correction_analytic(a, v_ext[a], bodies, medium)  # O(Ma^2), from eq. (7)
        else:
            F_comp = 0
        Fa = F_inc + F_comp
        if audit_step:
            Fa_quad_inc  = force_incompressible_quad(a, bodies, medium, npts)
            Fa_quad_comp = force_compressible_quad(a, bodies, medium, npts)
            log_audit(a, Fa, Fa_quad_inc, Fa_quad_comp)

    # 4) Leapfrog update
    v += 0.5 * (F/M) * dt
    x += v * dt
    recompute F at new x
    v += 0.5 * (F/M) * dt

    # 5) Optional mass intake every intake_every steps
    if step % intake_every == 0 and use_flux_mass:
        for a in bodies:
            Mdot = intake_flux(a, bodies, medium, npts_intake)   # uses ρ* from v_ext
            a.M += Mdot * (intake_every * dt)

    # 6) Diagnostics & output
```

---

## 14) Deliverables

* CLI app (`run.py`) that reads YAML, runs slab simulation, optional GR comparator, writes:

  * `state.csv` (time, x,v, M, Q), `forces.csv` (analytic & quad), `elements.csv`, `diagnostics.json`
  * quick plots (orbit tracks, (F)–(r), precession vs (1/c_s^2))
* Unit tests for the identities (3), coefficient (C), and compressible scaling.

---

This is everything the AI needs to implement the full simulation with the right numerics and the “no-way-it’s-GR” evidentiary trail.

