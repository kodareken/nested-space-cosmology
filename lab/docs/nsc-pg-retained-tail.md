# High-energy tail inputs for the other retained massive groups

Group 13 already has a full-energy common-PG covariance. This owner supplies
only the high-energy tail integrals for the remaining 31 retained massive
groups. Low-energy, subgap and mid-band inputs are assembled separately. Full
retained C1b remains OPEN.

## Method

Reuse the published group-13 endpoint series. For each signed family,

$$
K=\texttt{centered\_integral}(),
\qquad
G=\texttt{integral}((1,1,1)),
$$

from

```text
MassiveTailPackets(order=16, points_per_unit=64, minimum_energy=L)
EndpointTailSeries(lower=L, maximum_power=24, normalization_terms=6)
```

`K` and `G` are positive-energy 8×8 matrices. They are not signed-folded and
are not a covariance. The Fourier/normalization remainder of the series is a
separate representation bound; it is not a global physical mode-error theorem.

Groups 1–32 excluding 13: nonzero angular labels use both signs; angular-zero
group 23 is plus only. That is 61 signed families. Per group the sum of
representation C-bounds must be ≤ 10⁻⁹ and Gram-bounds ≤ 2×10⁻⁹ (one bound
when angular is zero). If `L=160` fails, one numerical split to `L=320` is
used for both angular signs. Failing 160 data are kept as diagnostics. This
`L` is a quadrature split, not a physical refit.

## Commands

```sh
python3 scripts/derive_nsc_pg_retained_tail.py --workers 2
python3 scripts/derive_nsc_pg_retained_tail.py --check
```

The record hashes sources and inputs. It does not change A, q, Ω, ζ or
V_full, substitute an LLL map, compute stress, or take a metric timestep.

The initial worker used a thread fallback after its sandbox denied process
semaphores. Before integration, an independent calculation evaluated every new
L=320 exponential-integral table in separate processes with cloned 50-digit
contexts. All 29 numerical matrices were identical. The
[precision audit](nsc-pg-tail-moment-audit.json) preserves that evidence; the
reproducer now permits separate processes or explicit serial execution only.
No physical tail data were changed to pass the audit.

All 61 signed families satisfy their selected representation-error budgets.
Groups 10, 11, 12, 31 and 32 use the numerical transition at 320; the others
use 160. The largest summed covariance bound is $9.73\times10^{-10}$ and
the largest summed Gram bound is $1.95\times10^{-9}$. The same transition
is used for both angular signs of a group.
