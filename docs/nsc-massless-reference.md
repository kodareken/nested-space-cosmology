# Absolute massless reference and its regulator match

A massless Hadamard reference can be constructed on the actual child tail
through its smooth conformal extension. The standard anomalous conformal
stress law then fixes its leading absolute stress in a declared
renormalization convention. Its leading radial null stress is negative in
that convention. The calculation retains all angular modes in the field
definition; it does not use a lowest-Landau-level vacuum approximation.

Using this renormalized reference in the finite-cutoff theory requires a
specific conversion functional. The converted reference and complement
sum to the original determinant. Their reallocation does not produce a
new vacuum source or close the full self-sourcing equation.

## DHP applicability audit

| Imported formulation | Application here | NSC-specific work still required |
|---|---|---|
| Dappiaggi–Hack–Pinamonti, Hadamard Dirac stress | A conserved renormalized stress with stated local ambiguities; differences use the smooth bispinor remainder | Actual parent state, compact/transmission domain and finite-cutoff matching |
| Hollands, ultrastatic and general globally hyperbolic Dirac states | A Hadamard state propagated through a smooth auxiliary metric | Match its reference covariance to the parent and preserve flux |
| General anomalous conformal-stress law | Transport the massless reference from the smooth conformal geometry; retain nonzero Weyl terms | Evaluate its finite-radius state-dependent part and its regulator conversion |
| Stored sphere spectrum and heat coefficients | Normalize the invariant limit and check light/complement allocation | Full angular stress on the actual time-dependent child |
| Euclidean-to-in-in response methods | A route to causal nonlocal kernels under their state assumptions | Establish the appropriate state kernel; do not assume the asymptotically flat in-vacuum rule automatically applies |

Primary owners are [DHP, Definition 4.1 and Theorem 4.1](https://arxiv.org/html/0904.0612),
[Hollands](https://arxiv.org/abs/gr-qc/9901069), and
[Barvinsky–Wachowski, equations 4.17, 4.32 and 4.34](https://arxiv.org/html/2306.03780v3).
The last gives a conformal-stress relation beyond conformally flat metrics.
These established methods are applied here; their general theorems are
not rederived or claimed as new.

## A definite auxiliary reference specification

Use the already verified child conformal metric, with u=-1/rho:

\[
\bar g=g/r^2=\frac{du^2}{(1+u^2)^2W(u)}-W(u)dz^2-d\Omega_2^2,
\quad W(u)=3(\pi-\arctan u)+\frac{3u-u^2}{1+u^2}.
\]

W(0)=3 pi and W is smooth. Continue this metric to u<0. The implementation
agrees with W down to u=-d/2, smoothly interpolates to W=3 pi on
[-d,-d/2], and is constant below -d, with d=0.1. The interpolation uses
the standard smooth switch exp(-1/t)/[exp(-1/t)+exp(-1/(1-t))]. Its width
is an auxiliary reference choice, not a parameter fitted to stress.
The continued W is positive: on the interpolation interval the analytic
W is at least 3 pi-d^2, and the interpolation is a convex combination.

Convert the lapse to conformal proper time eta. The constant region is
ultrastatic and can be extended freely in eta. Its spatial sections are
the complete unwrapped R_z x S2. The untwisted massless Dirac Hamiltonian
has omega^2=k^2/(3 pi)+kappa^2 with |kappa|>=1, so there is no zero-mode
choice. Let C_out be its negative-energy ground-state projector and set

\[
\bar C_{\rm ref}(\eta)=U_{\bar D}(\eta,\eta_{\rm out})
 C_{\rm out}U_{\bar D}^{\dagger}(\eta,\eta_{\rm out}),
\qquad \psi=r^{-3/2}\bar\psi.
\]

Hadamard propagation and conformal covariance supply a reference on the
tail. A projector imposed at an arbitrary instant of the actual changing
metric would not establish this property. The construction is reflection
symmetric in z and carries no net longitudinal reference flux. It is not
yet a parent-selected state. That requires the matched parent covariance,
cross-boundary correlations and flux, not just the retarded boundary map.

Different smooth auxiliary continuations give bounded conformal stress
differences at u=0; the existing r^-4 result controls their physical
difference. Their finite-radius stress can differ. The numerical full
covariance transport and finite-radius angular stress have not been
evaluated in this record. The construction is massless and untwisted;
massive compact modes acquire an unbounded conformal mass m r at u=0
and need their own state construction.

## Absolute reference stress on the actual asymptotic tail

Work in +---, hbar=c=L_throat=1. Assign the renormalized light field the
covariant proper-time subtraction convention inherited from the project:

\[
T^\mu{}_\mu=\frac{1}{16\pi^2}
 [-C^2/20+11E_4/360-\Box R/30].
\]

No independent dimensionful vacuum or Einstein term is assigned to that
renormalized light contribution. Their existing finite-cutoff contribution
is retained by the conversion below. A Hadamard implementation must match
the finite R-squared convention as well; its removable box-R coefficient
cannot be imported with an incompatible sign or convention. No complete
physical NSC coefficient is set to zero by this choice of allocation.

The general conformal law contains r^-4 times the barred stress, a barred
Weyl term with two derivatives of log(r) times the barred Weyl tensor,
and local curvature tensors H1,H3 and Ricci contracted with Weyl. On the
actual geometry, the already established conformal expansion gives

\[
R_L=-36\pi+18\pi u^2+O(u^4),\qquad
R^\mu{}_\nu\longrightarrow-9\pi\delta^\mu{}_\nu,
\qquad C_{\rm orth}=O(u^3).
\]

Thus the barred Weyl tensor is O(u), and its differentiated log term,
after multiplication by r^-4, is O(u^3). Bounded barred stress and local
curvature terms contribute O(u^4). H1 tends to zero; H3 tends to
3 H_child^4 times the identity. The leading reference stress is therefore

\[
\lim\langle T^\mu{}_\nu\rangle_{\rm ref}
=\frac{11H_{\rm child}^4}{960\pi^2}\delta^\mu{}_\nu,
\qquad H_{\rm child}=\sqrt{3\pi}.
\]

This establishes the limit for the specified conformally regular reference
class, rather than assigning exact de Sitter stress at finite radius. In
the declared heat convention the new curvature-tensor application also gives

\[
\rho_{\rm ref}=\frac{33}{320}-\frac{3u^2}{32}+O(u^3),\qquad
p_{\parallel,\rm ref}=p_{\perp,\rm ref}
=-\frac{33}{320}+\frac{u^2}{32}+O(u^3),
\]
\[
\rho_{\rm ref}+p_{\parallel,\rm ref}
=-\frac{u^2}{16}+O(u^3).
\]

The trace matches the inherited heat anomaly and the expansion obeys
the energy/work identity through this order. A finite R-squared change
of convention can change the u^2 coefficient; the complementary conversion
must then change with it. This negative reference contribution alone
does not establish the sign of the complete finite-cutoff source. The
expansion is not valid at the neck, where u is not small.

## Exact conversion to the adopted finite cutoff

The compact matching already defines Gamma5=Gamma_light,nu+Gamma_H,nu.
Introduce the canonical renormalized light functional at fixed normalization
mu, and define

\[
\mathcal C_{\nu,\mu}[g]
=\Gamma_{\rm light,\nu}[g]-\Gamma_{\rm light,ren}[g;\mu].
\]

If K(s)=Tr exp(-s D4^2) and A0,A2,A4 are its integrated heat coefficients,
the same ultraviolet subtraction gives

\[
\mathcal C_{\nu,\mu}=I_{\rm div}(\nu,\mu)
-\frac12\int_0^{\nu^{-2}}\frac{ds}{s}
\left[K(s)-\frac{A_0+A_2s+A_4s^2}{(4\pi s)^2}\right],
\]
\[
I_{\rm div}=\frac{A_0\nu^4/2+A_2\nu^2+A_4\log(\nu^2/\mu^2)}{32\pi^2}.
\]

This expression requires the same boundary and infrared prescription.
Its metric variation retains both the local terms and the finite remainder.
The latter cannot be discarded at the stored child curvature. The correct
allocation is

\[
\boxed{\Gamma_5=\Gamma_{\rm light,ren}
 +\left(\Gamma_{H,\nu}+\mathcal C_{\nu,\mu}\right).}
\]

For the invariant control, let a=1/H_child be the curvature radius; it is
distinct from the increasing areal radius r(T). The existing spectrum gives
partial_log_a Gamma_light,nu=Tr exp[-D4^2/nu^2]. Its power terms are
2[(nu a)^4-(nu a)^2]/3, and its renormalized limit is 11/90, in agreement
with the reference density. The known compact determinant is read from
its record, not recomputed.

At the stored a=0.325735 and nu=1, in the declared units:

| Contribution to the invariant density projection | Value |
|---|---:|
| Renormalized massless reference | 0.103125000 |
| Original finite-cutoff light contribution | 2.8628e-16 |
| Complement after the exact conversion | -0.102079088 |
| Unchanged total Dirac modulus projection | 0.00104591235 |

The pieces depend on normalization; their sum preserves the original
functional. Adding the renormalized reference to the unconverted complement
would change that functional. This conversion is not a mechanism for
adjusting the vacuum energy or creating a stationary solution.

## Full angular stress and source-budget gate

The accepted next calculation still needs the actual state-dependent
barred tensor at finite child radius, the full angular and massive compact
content, and the same-prescription causal conversion. The source to compare
with geometry is the complete matched tensor, including independent lapse,
radial and sphere variations and the boundary flux.

[Barvinsky's Euclidean-to-in-in rule, equation 2.29](https://arxiv.org/html/1408.6112)
is useful prior machinery, but its stated Poincare-invariant in-vacuum at
past asymptotically flat infinity is an assumption to verify. It does not
by itself supply the actual interior state kernel. A spatial Hamiltonian
cutoff cannot replace the spacetime regulator without an equality proof.

The angular source calculation should use the actual Dirac two-point/CTP
kernel, with the conformal reference and parent difference carried together,
then vary the conversion functional with the same state prescription.
The static source and response owners provide reusable variations. A short
curvature expansion is not controlled at the neck. The full finite-radius
angular source and the neck/throat budget have not been completed here.

Evidence: [record](../results/development/massless-reference.json), reproduced
by scripts/check_nsc_massless_reference.py --check. It checks the new
asymptotic tensor/trace/work map, regulator allocation and light-cutoff
controls. Existing compact and endpoint generators are not rerun.
Grok source audits independently assessed Hadamard applicability and the
leading asymptotic limit. A separate attempted review of the u^2 coefficient
timed out without a usable report and is not counted as verification; that
coefficient is checked here against the explicit conformal-stress tensors,
their trace and the energy/work identity.
