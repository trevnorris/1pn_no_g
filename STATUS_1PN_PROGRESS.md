# 1PN-without-G Status Check

## Theory vs. implementation
- The derivation in `1pn_no_g.tex` maps Newtonian dynamics to the slab control-surface law and predicts a perihelion advance that scales as \(\Delta\varpi \propto KM/(a c_s^2)\) once finite sound-speed effects are kept.【F:1pn_no_g.tex†L160-L187】
- The code follows that structure: incompressible forces implement \(\mathbf F=\rho_0 Q\,\mathbf v_{\rm ext}\), while `slab/surface.py` adds the compressible correction (currently evaluated by surface quadrature every step). For the spherical mouths used in production we now drop the renormalised surface term—Eq. (A.6) of the paper shows it cancels—so only the finite-sound-speed retardation survives.【F:slab/surface.py†L497-L705】【F:1pn_no_g.tex†L260-L309】
- A direct force-by-force comparison shows the combined incompressible + compressible slab force agrees with the textbook EIH 1PN prediction to \(3.5\times10^{-8}\) relative error along an eccentric Mercury trajectory, and the *difference* \(\mathbf F_{\rm comp}\) now matches the EIH 1PN excess itself (validated by `test_compressible_forces.py`).【0858aa†L41-L53】

## Numerical diagnostics we have today
- `scripts/final_diagnosis.py` (default resolution: 1000 steps at dt=0.002 yr) reports a large spurious retrograde baseline (≈−4.31×10⁻³ rad/orbit) and a compressible–incompressible difference of only 3.0×10⁻⁷ rad/orbit ≈ 0.62× the GR target (4.83×10⁻⁷).【8deed3†L1-L38】
- Refining the timestep by ×40 removes most of the baseline drift (down to −2.7×10⁻⁶ rad/orbit), confirming the integrator can suppress the artifact, but that refinement was applied only to the incompressible run—the compressible correction at matching resolution is still missing from the automated report.【9225bd†L1-L17】

## Independent spot checks (short custom runs)
- A timestep sweep of pure incompressible runs over five orbits shows the spurious perihelion drift collapses toward zero as expected for Newtonian motion: from −780 arcsec/orbit at dt=2×10⁻³ yr down to −0.02 arcsec/orbit at dt=10⁻⁵ yr.【5108ad†L1-L2】【3c54e9†L1-L1】【098198†L1-L1】【f6aa54†L1-L1】【c40f54†L1-L1】【623ead†L1-L1】
- Short 8.3-orbit integrations (dt=0.001 yr, 1000 steps) for e=0.10 and e=0.20 reproduce the same ≈0.62–0.70×GR signal after subtracting the incompressible drift, matching the diagnosis script but with similarly noisy baselines.【68ad6d†L1-L6】【660cba†L1-L7】
- Integrating the reference GR 1PN equations with the same timestep and duration produces the *same* retrograde bias (≈−0.11 arcsec/orbit) when analysed by `scripts/precession_helpers.analyse_precession`, confirming the discrepancy is a measurement issue rather than a physics mismatch; longer secular averages are required before the +0.099 arcsec/orbit signal emerges.【64da1f†L1-L3】

## Bottlenecks and gaps
1. **Performance** — every compressible force evaluation still routes through the quadrature path (intended as an audit), so long, high-resolution runs are prohibitively slow. Implementing the planned analytic \(O(Ma^2)\) correction would unlock 100× more timesteps and allow refined comparisons.【F:slab/surface.py†L619-L705】
2. **Diagnostic coverage** — `final_diagnosis.py` now demonstrates that timestep refinement tames the incompressible artifact, but it does not yet rerun the compressible case at the refined dt, so we cannot quote a “best-effort” 1PN estimate.
3. **Orbit sampling** — current studies cover ≈8 orbits. The LaTeX derivation assumes secular averaging; we should extend runs (or average across multiple periapsis-to-periapsis measurements) to suppress cycle-scale noise before comparing with the \(\Delta\varpi\) formula.【F:1pn_no_g.tex†L181-L187】
4. **Estimator conditioning** — diagnostics subtract two noisy \(\omega(t)\) fits to isolate the 1PN signal. Once compressible reruns at refined dt are practical, compute the slope of \(\omega_{\rm comp}(t)-\omega_{\rm inc}(t)\) directly so the regression never sees the large retrograde background.

## Recommended next actions
1. **Implement the analytic compressible correction** described in `slab/surface.py`’s TODO so production runs no longer depend on the slow quadrature fallback.【F:slab/surface.py†L619-L705】
2. **Extend `final_diagnosis.py`** to repeat the refined-timestep run with compressible forces enabled and report the difference against GR.
3. **Automate multi-orbit averaging** (e.g., via `scripts/validate_1pn_precession.py` once matplotlib is available) to sweep eccentricities, measure \(\Delta\omega\) over ≳100 orbits, and verify the predicted \(1/(1-e^2)\) scaling.
4. **Add a dedicated ω-difference regression** so the reported 1PN signal comes from a single fit rather than the difference of two large numbers; this will tighten error bars once longer compressible runs are feasible.

Until those steps land, the project demonstrates the correct qualitative scaling (Δω ∝ c_s⁻²) but only reaches ~60–70% of the expected GR coefficient at practical resolutions. Longer, better-resolved runs should resolve whether the remaining gap is numerical (likely) or physical.
