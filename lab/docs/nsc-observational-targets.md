# Cosmological observables: targets and required closure

Checked against primary collaboration papers on **2026-09-08**. This page
specifies the path from the proposed nested boundary mechanism to observations.
It preserves the nested-space hypothesis while distinguishing measured targets,
conservation requirements, and predictions still to be derived. The executable
[observable bridge](nsc-observable-bridge.md) and its
[record](../results/nsc-7-observable-bridge.json) implement the associated
controls; they do not constitute a cosmological fit or a resolved tension.

## Normalize the source before predicting expansion

Use `c=1`, proper cosmological time, `N=ln(a)`, and prime `d/dN`. Let
`rho=rho_b+rho_d`, `p=p_b+p_d`, and `f=rho_d/rho`. Here dark denotes the combined
unresolved response, including effects ordinarily assigned to dark matter and
dark energy. Let `Q_b>0` denote internal transfer into the visible component;
`J_b,J_d` denote external supply to the respective components of the room.

\[
\dot\rho_b+3H(\rho_b+p_b)=Q_b+J_b,\qquad
\dot\rho_d+3H(\rho_d+p_d)=-Q_b+J_d,
\]
\[
\boxed{\dot\rho+3H(\rho+p)=J,\qquad J=J_b+J_d.}
\]

Internal exchange cancels from total continuity. With `w_i=p_i/rho_i`,
`epsilon_b=Q_b/(H rho)`, and `j_i=J_i/(H rho)`, exact differentiation gives

\[
\boxed{f'=3f(1-f)(w_b-w_d)-\epsilon_b+(1-f)j_d-fj_b.}
\]

This identity needs no assumed gravitational field equation. The pressure
difference and the external source partition both affect the fraction.
For visible dust and a closed total ledger it becomes
`f'=-3 w_d f(1-f)-epsilon_b`.

An external source for a room subsystem requires its reservoir, boundary stress,
and gravitational matching to be accounted for. The complete stress in ordinary
constant-G Einstein equations is conserved by the Bianchi identity. Adding `J`
to its continuity equation while retaining every original Einstein equation is
inconsistent. If one retains only the flat Friedmann constraint
`H^2=8 pi G rho/3` with sourced continuity, differentiation instead requires

\[
q_{\rm dec}=\frac{1+3p/\rho-j}{2},\qquad j=J/(H\rho).
\]

This corresponds to `p_eff=p-J/(3H)` and needs its own physical closure; it does
not justify using the unsourced Raychaudhuri equation with the original `p`.

For signature `+---`, unit timelike `u` and outward unit spacelike `s`, the
outward energy flux is `F_out=-T^{mu nu}u_mu s_nu` in `c=1` units. A spatial
two-surface integral of this flux is power; its integral over a timelike
worldtube, including proper time, is transferred energy. Only

\[
\boxed{Q_b=P_{\rm deposited}/V_c}
\]

has the required units of energy per proper volume per proper time. In SI the
physical flux has units `W m^-2`; the appropriate factor of `c` must be restored
when contracting an energy-density-normalized stress tensor. For a static
metric with lapse `mathcal N` and `xi=partial_t=mathcal N u`,

\[
E_\xi=\int\mathcal N\rho\,dV,\qquad
\frac{dE_\xi}{dt}=-\oint\mathcal N^2 F_{\rm out}\,dA
                  +\int\mathcal N^2 Q\,dV.
\]

The receiver's proper power is incoming Killing power divided by its uniform
lapse squared. Moving boundaries and evolving geometry also require their work
terms. The [energy-transfer owner](nsc-energy-transfer.md) provides the general
work identity. Equivalent window sources and surface fluxes must not be counted
twice. Incoming Dirac energy may become radiation or heat; its conversion into
fixed-mass baryonic dust remains a physical dependency, including baryon-number
accounting.

## A plateau constrains the mechanism; it does not supply it

For the closed, flat, constant-G Einstein baseline with visible dust and positive
densities, `q_dec=(1+3 w_d f)/2`. A finite constant fraction `0<f_*<1` requires

\[
\boxed{\epsilon_{b*}=(1-f_*)(1-2q_{{\rm dec}*}),\qquad
\frac{Q_{b*}}{H\rho_b}=1-2q_{{\rm dec}*}.}
\]

An accelerating plateau therefore needs `Q_b/(H rho_b)>1` in this convention.
This necessary balance is not a derived transfer law or a proof of attraction.
The full state and geometry perturbations must establish stability.

For `Q_b=0`, a finite positive limiting density ratio with regular asymptotic
evolution implies `w_d -> 0` and `q_dec -> 1/2`. The conclusion also requires
the preceding gravity, dust, and conservation assumptions. Zero transfer alone
does not imply it: the flat positive-Lambda control has `Q_b=0` and
`q_dec -> -1`, while its dark-to-baryon ratio diverges with finite curvature.
The [plateau owner](nsc-plateau-conditions.md) derives these distinctions.

Even a known function `Q_b(t)` does not select the pressure or expansion.
For any regular chosen expanding `a(t)` in the closed dust baseline,

\[
\rho_b(t)=a^{-3}(t)\left[a^3(t_0)\rho_b(t_0)
                  +\int_{t_0}^t a^3(s)Q_b(s)\,ds\right],
\]
\[
\rho_d=\frac{3H^2}{8\pi G}-\rho_b,\qquad
p_d=-\frac{2\dot H+3H^2}{8\pi G}
\]

satisfy the background equations wherever the densities remain positive.
Different admissible histories can share the same transfer. The common-action
stress and its state evolution must select among them. Likewise a retarded
boundary map does not determine an occupation state: the existing static
vacuum and prepared excitation have the same Hamiltonian but different currents.

## From one state and action to measurable predictions

1. Specify the common action, physical state, boundary and initial data, scale,
   and stress renormalization. Solve their coupled variation to obtain
   `T_mu_nu`, transfer, boundary work, and component deposition in one ledger.
2. Solve the metric and conservation equations with that stress to predict
   `H(z)`, component densities, pressures, and fractions. Expansion distances
   follow from integrals of `1/H(z)`. BAO and CMB comparison also require the
   early thermal history, sound horizon, and relevant perturbation transfer.
3. Derive the linear response of that same state and action: `delta Q`, momentum
   transfer, pressure perturbations, anisotropic stress, and metric response.
   Together with initial perturbations these determine matter clustering
   `P_m(k,z)` and the lensing potential `Phi+Psi`. A background scalar `Q` alone
   does not determine them. Different momentum-transfer completions can alter
   growth; see the primary analysis by
   [Clemson et al.](https://arxiv.org/abs/1109.6234).
4. Propagate those quantities into growth, weak lensing, galaxy clustering, and
   CMB observables, including calibration and nuisance uncertainties in the
   likelihood. `S_8=sigma_8 sqrt(Omega_m/0.3)` uses the matter abundance; the
   combined dark fraction defined above is not `Omega_m`.

Lambda-CDM is a predictive benchmark: its parameters and assumptions jointly
determine expansion and perturbation observables. Its phenomenological inputs
leave questions of physical origin open. CPL, `w(a)=w_0+w_a(1-a)`, is a useful
two-parameter background description; it supplies neither a unique microscopic
mechanism nor a complete perturbation law by itself. NSC's constructive target
is to derive the common mechanism and its response, then compare its predictions
with these benchmarks using matched data and uncertainty accounting.

## Current development targets

- **Expansion:** H0DN reports `H_0=73.50 +/- 0.81 km s^-1 Mpc^-1`. Its quoted
  differences are `7.1 sigma` from Planck+SPT+ACT under flat Lambda-CDM and
  `5.0 sigma` from BBN+DESI2 BAO. A tension significance must name its comparison.
  [H0DN Collaboration, A&A 708 A166 (2026)](https://arxiv.org/abs/2510.23823).
- **Clustering and lensing:** DES Y6 finds `S_8=0.789 +/- 0.012` in Lambda-CDM.
  Its difference from the combined primary-CMB datasets is `2.6 sigma` in the
  `S_8` projection and `1.8 sigma` in the full parameter space.
  [DES Collaboration, January 2026](https://arxiv.org/abs/2601.14559).
- **A complementary lensing result:** KiDS-Legacy finds
  `S_8=0.815 (+0.016,-0.021)`, agreeing with Planck at `0.73 sigma`.
  [Wright et al., A&A 703 A158 (2025)](https://arxiv.org/abs/2503.19441).

These are measured comparison targets, not NSC predictions or evidence that an
unspecified transfer resolves a discrepancy. Their differing outcomes require
dataset-specific likelihood tests. The quoted analyses use their stated model
assumptions; testing NSC requires predicting the underlying observables rather
than importing a fitted parameter value as a model-independent measurement.

Because these data have now informed development, they are development targets.
Confirmation requires an explicitly reserved, previously unseen comparison or
future independent data, with the operator, state prescription, parameters, and
prediction procedure frozen before inspecting that result.
