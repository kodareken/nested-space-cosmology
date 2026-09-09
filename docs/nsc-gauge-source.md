# A definite interaction correction to the magnetic Dirac source

The retained Dirac and gauge fields supply an explicit interaction in the
projected magnetic sector. Applying known QED2 bosonization gives one charged
collective mode with a radius-dependent screening scale. Its transmission
through the recorded exterior now determines a conditional correction to
the outgoing power and, by current conservation and regularity, the horizon
null source.

This is a contribution from fields already present in the candidate action.
The collective variable is not an inserted fundamental scalar, a resolution
compensator, or a predicted four-dimensional particle mass. The geometry,
full gauge coupling and state remain stipulated inputs to this development
calculation. Full self-sourcing is still open.

## Decision and reused physics

The obstacle was asking the free determinant to supply a different,
unspecified common-functional contribution. Three paths were considered:
reweighting the scale average, quantizing all metric/gauge fluctuations, and
including a definite interaction among the retained light fields. The first
still needs its measure; the second needs the full constrained Hessians and
domain. The projected Dirac–gauge interaction already has an exact known
reduction, so it is the next bounded source calculation.

The distinction in [Andrianov–Kurkov–Lizzi, equations19,25–26 and47–48](https://arxiv.org/html/1106.3263v1)
is retained: restoring Weyl invariance and representing a physical-cutoff
determinant by collective fields are different constructions. Neither
specifies the remaining geometric measure. No protected thesis text or
historical anomaly result is changed here.

Reuse [Maldacena–Milekhin–Popov, Appendix A](https://arxiv.org/html/1807.04726v3)
for the magnetic light sector and Schwinger mechanism. They already derive
one massive gauge/collective mode and the remaining massless flavour sector.
Curved-space bosonization is established, for example in
[Alimohammadi–Mohseni Sadjadi, section2](https://arxiv.org/html/hep-th/0011232v2).
The new work maps that action, its source variations and transmission to
the NSC exterior; it does not rediscover bosonization or Hawking radiation.

## Retain the same gauge and fermion normalization

Use the canonical massless 4D Dirac zero field already supplied by the
paired compact candidate. In a magnetic sector with q=|q_mag|, it supplies
q lowest angular Dirac modes. Retain the spherical two-dimensional gauge
field along time and radius. The four-dimensional Maxwell convention is

\[
S_A=-\frac1{4g_4^2}\int\sqrt{|g_4|}F^2,
\qquad e_2^2(r)=\frac{g_4^2}{4\pi r^2}.
\]

The physical matched g4 is not selected here. In particular, the
[Dirac complement coefficient](nsc-compact-matching.md) alone is not its
complete inverse coupling. That coefficient already includes its 5D Dirac
UV contribution; adding another copy of the same UV Maxwell term would
double-count it. The retained light contribution and other same-action
completion remain necessary.

QED2 bosonization separates the total U(1) current from a neutral sector
of central charge q-1. The canonically normalized charge mode X has

\[
\mu^2(r)=\frac{q e_2^2(r)}{\pi}
        =\frac{q g_4^2}{4\pi^2r^2}.
\]

This is a local collective screening scale. To see its source explicitly,
write the electric part in the zero-electric-flux sector as

\[
\mathcal L_E=\frac{E^2}{2e_2^2}-\sqrt{\frac q\pi}EX.
\]

Eliminating E gives E=e2^2 sqrt(q/pi) X and L_E=-mu^2 X^2/2.
Thus the charged action is the canonical two-dimensional kinetic term
minus this positive potential. Nonzero electric-flux or topological sectors
require their additional data; they are not silently included.

The new determinant contribution replaces the charged part of the free
LLL theory:

\[
\Delta\Gamma_{\rm gauge}
 =\Gamma_{\rm scalar}[\mu(r)]-\Gamma_{\rm scalar}[0]
  +\Gamma_{\rm local/zero\ modes}.
\]

The q-1 neutral CFT units remain. Adding a whole new massive scalar
determinant to the original free LLL action would count the charge mode
twice. Local renormalization and zero-mode terms must be matched to the
same parent action before this becomes an absolute stress tensor.

## Radius variation supplies the Maxwell pressure

The potential is also the eliminated radial Maxwell energy:

\[
\rho_E=\frac{E^2}{2g_4^2}
      =\frac{\mu^2X^2}{8\pi r^2},\qquad
p_{\rm sphere}=\frac{\mu^2X^2}{8\pi r^2}.
\]

The angular pressure follows by varying r in the same reduced action.
It must not be dropped when the mode is represented as a two-dimensional
scalar. On shell,

\[
\nabla_a t^a{}_b=\tfrac12 X^2\partial_b\mu^2.
\]

Lifting t_ab/(4 pi r^2) and including T_theta^theta=-p_sphere cancels that
radius-dependent force in the full spherical Ward identity. These are
operator/action identities. Their quantum composite expectation values
still need one consistent renormalization prescription. They supply the
correct source variation, not an uncomputed absolute vacuum expectation.

## The screening scale becomes a transmission barrier

Reuse the horizon and surface gravity from the
[unwrapped source record](nsc-horizon-source.md). Outside its horizon,

\[
ds_2^2=A(dt^2-dx_*^2),\quad
\frac{dx_*}{d\rho}=\frac1A,\quad r^2=1+\rho^2.
\]

The collective mode obeys

\[
[-\partial_{x_*}^2+V(\rho)]X_\omega=\omega^2X_\omega,
\qquad
V(\rho)=A(\rho)\frac{qg_4^2}{4\pi^2(1+\rho^2)}.
\]

There is no four-dimensional minimally coupled scalar r''/r term: X is
the canonical two-dimensional current/electric-displacement mode. V is
nonnegative outside the horizon and tends to zero at both asymptotic
scattering ends. In particular mu(r) tends to zero in the parent exterior.
A constant-mass Boltzmann suppression is therefore not the correct flux
calculation. This exterior operator has a nonnegative quadratic form; that
does not establish the health of the full coupled metric theory.

With the specified Unruh state and no incoming parent bath, the projected
Killing power is

\[
P_\infty=(q-1)\frac{\kappa_h^2}{48\pi}
 +\frac1{2\pi}\int_0^\infty
       \frac{\omega\,\mathcal T_{\rm charge}(\omega)}
            {e^{\omega/T_H}-1}\,d\omega,
\qquad T_H=\frac{\kappa_h}{2\pi}.
\]

The unit-transmission limit recovers the existing free result q kappa_h^2/(48 pi).
The charged mode uses Bose statistics; the central-charge normalization
agrees with the free complex Dirac channel. The underlying flux normalization
is the established [Iso–Umetsu–Wilczek result](https://arxiv.org/html/hep-th/0602146v2).
The numerical boundary condition is unit ingoing amplitude at the future
horizon. For this real, stationary potential and equal asymptotic wave
speeds, scattering reciprocity gives the same transmission for emission
from the horizon. No electric chemical potential or rotation is included.

At the stored kappa_h=0.238325799634, for q=1 and unselected test couplings:

| g4 input | Power / one free LLL channel | Parent Killing power |
|---:|---:|---:|
| 0.25 | 0.990575 | 0.0003731112 |
| 0.5 | 0.964044 | 0.0003631182 |
| 1 | 0.873478 | 0.0003290055 |

These are synthetic coupling controls, not a selection or fit of the physical
coupling. The ratios describe escaping energy in this projected sector;
they are not cosmological dark fractions.

For a stationary conserved source with future-regular PG components,
J^rho=-T_tau_tau at the horizon. Consequently the same calculated power gives

\[
\langle T(K,K)\rangle_h=-\frac{P_\infty}{4\pi r_h^2}.
\]

This inference uses the current and regularity, not a new full stress-tensor
integration. A nonzero net power still requires geometry evolution or a
compensating flux for a self-consistent solution; it cannot support an
exactly static geometry by itself.

## Numerical and physical scope

The primary solver evolves incoming/outgoing amplitudes. An independent
scalar-field IVP agrees at three frequencies. The exact Wronskian supplies
the stable transmission formula; raw amplitude extraction and current
defects are retained. Relative current defects are below 6e-12 in the three
base coupling calculations. No horizon root or earlier spectrum is rerun.

Use the established [Boonserm–Visser bounds](https://arxiv.org/pdf/0901.0944),
applied separately to the two omitted tails and combined as transfer matrices.
Writing T=sech^2(theta), their added scattering rapidity is bounded by

\[
S_{\rm tail}=\frac{qg_4^2}{8\pi^2\omega}
 [\arctan\rho_L-\arctan\rho_h+\arctan(1/\rho_R)].
\]

Thus T_full lies between sech^2(theta+S_tail) and
sech^2(max(theta-S_tail,0)). Omitted thermal frequencies also have explicit
positive bounds. At g4=0.5, the base charged-power ratio lies between about
0.9640376 and0.9641120 when those tail comparisons are applied. Numerical
ODE and quadrature errors are assessed separately; these are not rigorous
floating-point interval enclosures. Frequency nodes, both spatial ends and
ODE tolerance are changed independently in the record.

This is the specified LLL plus spherical-photon theory. Higher angular
photons, Landau levels, compact excitations, induced interactions and possible
four-dimensional magnetic catalysis are outside the calculation. Angular
gaps shrink toward infinity, so local near-throat scale separation is not
a uniform proof of a full four-dimensional truncation. The q-1 neutral
count is the projected QED2 result, not a theorem about all interacting 4D
channels. The full state, coupling, absolute source and metric backreaction
remain to be matched. No cosmological Q or observational prediction is claimed.

`scripts/check_nsc_gauge_source.py --check` authenticates its dependencies
and compares every field of `results/development/gauge-source.json` under
the existing absolute3e-9/relative3e-8 policy. New records refuse overwrite.
The source, assumptions and error controls are retained together so the
next calculation can use this interaction rather than repeat the free one.

A bounded Grok review checked the projected action and its source references.
An additional pressure/bound review did not return usable evidence and is
not counted as verification. Those identities are explicitly checked in
the runner and against the cited scattering bounds; no prior result was
rerun to compensate for that missing report.
