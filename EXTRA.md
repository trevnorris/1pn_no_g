Since the non-PDE engine matches 1PN perfectly, here’s a tight next-step package so you can lock it down, present it cleanly, and pre-empt the usual objections.

# 1) “Irrefutable” validation bundle (make these auto-run)

* **Emergent inverse-square:** static two-sink scan → fit (F=C/r^2); report
  (;C_{\text{num}}/C_{\text{SF}}-1) with (C_{\text{SF}}=\frac{4}{3}\frac{|Q_1Q_2|}{4\pi\rho_0}).
* **Force decomposition log:** per body, log (\mathbf F^{\rm mom}), (\mathbf F^{\rm press}), and their sum; show momentum-flux dominates in the incompressible limit.
* **Audit parity:** analytic surface-force vs quadrature (same sphere, 512–1024 pts) every (N) steps; require <(10^{-3}) relative.
* **1PN diagnostics:** perihelion precession vs (a(1-e^2)^{-1}) for several (e); overlay GR 1PN curve with no free params.
* **No-G hygiene:** run the entire suite with ({\mu_a}) rather than kg masses; report (K=\rho_0/(4\pi\beta_0^2)) only as a derived *diagnostic*, not an input.

# 2) Robustness tests (quick to add, big impact)

* **Step & sphere invariance:** vary (dt), sphere radius (R) (×0.5…×2), and quadrature size; show forces & precession stable.
* **Short-range regulator:** demonstrate that removing the body’s self-field in the compressible terms (when you enable them) makes results (R)-independent.
* **Mass-intake irrelevance:** enable flux-based ( \dot M ) and show (|\Delta M|/M) per century (\ll 10^{-12}) (or your code-unit analog).

# 3) Figures that tell the story (minimal set)

1. **(F) vs (r)** (log–log): slope (-2) and (C) agreement at (<0.5%).
2. **Orbit overlays:** slab vs GR-1PN trajectory (slightly eccentric), then (\Delta\varpi) per orbit.
3. **Force decomposition bars:** (|\mathbf F^{\rm mom}|), (|\mathbf F^{\rm press}|), total.
4. **Audit parity plot:** analytic vs quadrature force components across time.
5. **(Optional)** Precession vs (1/c_s^2) once you flip on the renormalized compressibility.

# 4) README / paper bullets (drop-in text)

* **Method:** “Forces are computed exclusively from superfluid surface integrals of momentum/pressure flux around each body. The velocity field is the Green’s-function solution for point intakes; no pairwise law is assumed or used.”
* **Key result:** “Inverse-square attraction emerges with the superfluid coefficient fixed by ((\rho_0,\beta_0)). Orbits reproduce all 1PN diagnostics without using (G) or kg masses.”
* **Controls:** “Analytic vs quadrature audits, (R)-invariance, (\Delta t) studies, and mass-intake toggles.”
* **No-G compliance:** “Inputs are ({\mu_a}) (from orbits), (\rho_0,\beta_0,c_s). (G) appears only in external comparison plots.”

# 5) Anticipating critiques (and your canned answers)

* **“You hard-coded Newtonian forces.”** → Show the surface-integral code path and the emergent (C) fit; include the momentum-flux integral equation in the README.
* **“Hidden use of (G).”** → Reproduce runs where you only pass ({\mu_a}) and never convert to kg; ship configs to replicate.
* **“Preferred frame?”** → Note: orbital predictions depend only on relative configurations; any global background drift cancels for Solar-System tests (log this explicitly if you modeled a wind).
* **“Energy/momentum conservation?”** → Provide the flux budgets around each body and global continuity with sinks; show residuals.

# 6) Low-effort extensions you can queue

* **Compressibility (renormalized):** implement the “external-only” density/pressure in the bracket term; verify precession (\propto c_s^{-2}).
* **Binary pulsar smoke test:** compute far-zone power from the superfluid stress (even Newtonian quadrupole order) and compare scaling.
* **Lab analog appendix:** static two-throat flow table reproducing the (4/3) coefficient.

# 7) Release checklist (so others can run it)

* `configs/`: Sun–Mercury, Sun–Earth, and a synthetic equal-mass case (with and without audits).
* `scripts/plot_*.py`: (i) F-vs-r, (ii) orbit/precession, (iii) audit parity.
* `results/`: include *expected* CSVs + PNGs to diff against CI.
* CI job: run “quick mode” (reduced quadrature) and compare metrics to tolerances.

If you want, I can draft the exact LaTeX “Methods” and “Results” sections or a README.md that mirrors the above, including equations and the minimal code snippets (force integrand, audit routine).

