# Prescribed spectral pulse: work and the action-reference gap

The [bounded consumer](../src/recursive_horizons/nsc_discovery_spectral_work.py)
uses the **actual uniform case** of
`nsc-discovery-dynamic-preparation-v2.json/npz`, authenticated against its
producer `67ab570`. It reads the stored six source columns, occupations,
eigenvalues and uniform nodal geometry. It performs no source preparation
or geometry solve. The [CLI](../scripts/derive_nsc_discovery_spectral_work.py)
previews by default; [tests](../tests/test_nsc_discovery_spectral_work.py) use
small Fock/block examples and read-only source forecasts.

## Declared state and operator

For this case $Q_0=0.46394471443118$, $r_0=2.74634407266958$ and
$N_0=r_0Q_0$. On AP momentum $p_n=2\pi(n+1/2)/8$,

$$H_0(p)=p\sigma_2+\kappa Q_0\sigma_1,\qquad G=\partial_QH=\kappa\sigma_1.$$

The sector has $\kappa=1$ and multiplicity $M=4$ **outside** $H$ and the
column evolution law. The stored covariance is the full finite CAR covariance:
six positive states with occupations $(.75,.75,.5,.5,.25,.25)$ and empty
negative states. It is stationary with respect to this frozen field operator,
not a claim that the prepared geometry is stationary.

The consumer derives AP coefficients from the stored columns, checks the
actual eigen-equation and paired populations, and records covariance,
negative-occupation and uniformity roundoff. No dense ambient covariance is
formed. The occupied-channel positive-frequency gaps and residues are:

| Gap | Residue per AP channel, before $M$ | Degeneracy |
|---|---:|---:|
| $1.21565993076$ | $-0.313052244370$ | 2 |
| $2.53231737107$ | $-0.432868524929$ | 2 |
| $4.03512523584$ | $-0.236780396598$ | 2 |

These residues are $(c_--c_+)|G_{-+}|^2$. Only gaps and residues are recorded;
no singular pole is evaluated.

## One imposed pulse and a locked forecast

Prescribe uniform $Q(t)=Q_0+J(t)$, fixed $r=r_0$, $L=Q$ and zero shift:

$$J(t)=A\sin^2(\pi t/T)\cos(\omega t),\quad
 A=0.001,\quad\omega=0.8\Delta_1,\quad T=4\pi/\omega.$$

Here $\omega=0.972527944610$ and $T=12.9213465628$. This pulse has
$J=\dot J=0$ at both endpoints and $\int Jdt=0$. Its actual and reference
proper-clock endpoints therefore agree, though their clocks differ during
the pulse. Positivity follows from $Q\ge Q_0-A>0$.

Before nonlinear evolution, seal the second-order work prediction

$$W^{(2)}=M\sum_{a<b}(c_a-c_b)\Delta\epsilon\,
 |G_{ab}|^2|\widetilde J(\Delta\epsilon)|^2.$$

Stable finite Fourier integrals evaluate $\widetilde J$, including coincident
frequencies, without division by a pole. The analytic preview gives
$W^{(2)}_{\rm full}=-2.26972145\times10^{-5}$ and
$W^{(2)}_{\rm ref}=+7.56923703\times10^{-6}$.
These are forecasts; the nonlinear production pulse is root work.

## Nonlinear measurement and mathematical reference

DOP853 evolves only independent local $2\times2$ AP blocks in interaction
picture. Endpoint transition probabilities give work without subtracting
large initial/final energies. Independently integrated
$M\int\mathrm{Tr}[C(t)G]\dot Jdt$ supplies the energy/work closure indicator.
The initial constant force is analytically removed from the integrand;
its endpoint term is exactly zero. The pulse-area integral is also computed
and checkpointed to check the proper-clock endpoint. Unitarity, CAR,
up/down transition agreement and the positive-Q bound stay explicit.
Only local block endpoints and scalar integrals are retained; no dense global
propagator, covariance or history is stored.

The passive mathematical control is $C_{\rm ref}=P_-+C_+$ on this **same
frozen field pulse**. It is not fed to gravity and is not called a physical
vacuum. A separate $C_{\rm vac}=P_-$ control permits mean subtraction.
Means are linear in $C$, so reference-minus-vacuum work can match the declared
full-source work. Gaussian noise and generating functions are nonlinear in
$C$: for this actual source, the normalized characteristic ratio differs
from the full-source characteristic by about $0.07836$ at angle $0.2$.
Noise and characteristic numbers refer to one representative canonical CAR
sector; $M$ multiplies the mean/action work outside $H$.

Mapping the repository's declared full finite covariance to a physical
filled-sea excitation correction and the same $\Gamma$ remains unresolved.
The present mean-field realization stays valid in its declared domain.
Activity with respect to $M\mathrm{Tr}(CH)$ is not a universal failure of the
model or an antimatter restriction. This prescribed-pulse control introduces
no new force, interaction, autonomous regeneration or cosmological claim.

## Root production and replay

After freezing these producers, create one new directory through these stages:

```sh
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_spectral_work.py
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_spectral_work.py --prepare
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_spectral_work.py --predict
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_spectral_work.py --run
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_spectral_work.py --check
.venv/validation/bin/python scripts/lab.py -m pytest tests/test_nsc_discovery_spectral_work.py -q
```

Default output is `lab/results/development/nsc-discovery-spectral-work-v1`.
`--output` selects a new child prefix there. Prepare binds source/operator;
predict exclusively seals the work forecast; run exclusively creates the
nonlinear endpoint record; check reads and recomputes without evolution.
The aggregate limit is 30 CPU seconds and each JSON/NPZ pair is capped at
64 MiB. Inputs and source hashes are checked, and producing commits are
recorded only when source bytes match their Git blobs. Historical replay
can authenticate old producers through those blobs without healing records.
