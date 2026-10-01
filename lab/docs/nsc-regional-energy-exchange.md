# Regional coordinate ledger and normal-observer balance

Diagnostic only. The global number near zero is the constrained coordinate
Hamiltonian. Renewal, a continuum limit, and a frozen-link embedding are not
claimed. Geometry-derived \(B\) is not used. Six owned tests passed in 0.42 s.

\(F_L\), \(F_Q\), and \(F_\beta\) are nodal partials of \(M\operatorname{Tr}(CH)\),
with \(M=4\kappa\) once. The coordinate generator is the nodal sample
\(e=LF_L+\beta F_\beta\). \(F_Q\dot Q\) is coordinate metric work and includes
\(\dot Q=\partial_x(\beta Q)\). It is not observer pressure work. Proper
densities divide by the node spacing:

\[
\rho=\frac{F_L}{4\pi r^4 Q\,dx},\quad
p_r=-\frac{F_Q}{4\pi r^4 L\,dx},\quad
p_\perp=\frac{L F_L+Q F_Q}{8\pi r^4 LQ\,dx},\quad
j=-\frac{F_\beta}{4\pi r^4 Q^2\,dx}.
\]

The observer balance is proper pressure work plus the spatial lapse gradient
\(-V(j/q)\partial_x N\). The checked residual is that identity on the nodal
shell \(F_L/r=V\rho\,dx\), so the discrete flux and the lapse-gradient sample
both carry \(dx\). Hamiltonian derivatives, including \(F_Q/r\) and the
chain-rule \(\partial E/\partial r\), are not \(p_r\) or \(p_\perp\).

| Path | sha256 |
|---|---|
| `lab/src/recursive_horizons/nsc_regional_energy_exchange.py` | `9339f401ab3db9d0d2f36d2ae43dded15cb9f83d0c9d3463e4c48726eab73ce7` |
| `lab/tests/test_nsc_regional_energy_exchange.py` | `4b688de82022b88c566048299034e833396b5378a849dba57a22d52940fff722` |
| `lab/results/development/nsc-regional-energy-exchange-v1.json` | `19fe2b6194e2be3bb92d59c90e7b0bc067857c860739a813a6d5fda5c3d1da57` |

The record's module hash matches that source. Its `v5_npz_sha256` is
`0c92d3d545afffb47f22c9e8d705dc4d70b2df2ca2fdf5c6abc254f37af35b0d`, the same
final npz digest in the v5 record. Current coupling, Galerkin, feedback-action,
and conformal-source hashes match v5 `hashes_after`. Those files were not edited.

## Two different energies

Nodal \(F_L,F_Q,F_\beta\) are partial derivatives of \(M\operatorname{Tr}(CH)\),
\(M=4\kappa\) once. A density divides by the node spacing \(dx\).
\(N=rL\), \(q=rQ\), and \(V=4\pi q r^2\).

The coordinate generator is the nodal sample. It sums directly to the matter
trace. It is not the proper density \(\rho\):

\[
e=L F_L+\beta F_\beta=M(aK+mS-\beta J),\qquad a=L/Q,\quad m=\kappa L.
\]

\(L F_L=N(F_L/r)\) is that shell energy multiplied by the lapse. The
normal-observer shell sample is \(F_L/r=V\rho\,dx\), and

\[
\rho=\frac{F_L}{4\pi r^4 Q\,dx},\quad
p_r=-\frac{F_Q}{4\pi r^4 L\,dx},\quad
p_\perp=\frac{L F_L+Q F_Q}{8\pi r^4 LQ\,dx},\quad
j=-\frac{F_\beta}{4\pi r^4 Q^2\,dx}.
\]

\(F_Q/r\) and the chain-rule \(\partial E/\partial r\) stay Hamiltonian
derivatives. They are not \(p_r\) and \(p_\perp\).

## Coordinate flux and coordinate work

With \(L\) and \(\beta\) held fixed,

\[
\partial_t e+\partial_x(\Phi_{\mathrm{shift}}+\Phi_{\mathrm{proper}}+\Phi_{\mathrm{cross}})=F_Q\dot Q,
\]

\[
\Phi_{\mathrm{shift}}=-\beta e,\qquad
\Phi_{\mathrm{proper}}=M a^2 J,\qquad
\Phi_{\mathrm{cross}}=-M a\beta K=\beta Q F_Q.
\]

\(\Phi_{\mathrm{cross}}\) is the shift/conformal-pressure cross term. Summed
over columns it is one coordinate flux, not a current between spatial regions.
\(F_Q\dot Q\) is coordinate metric work. It uses the full \(\dot Q\), so
\(\dot Q=\partial_x(\beta Q)\) keeps it nonzero when the proper expansion vanishes.
A comoving window \(\partial_t W=\beta\partial_x W\) cancels only
\(\Phi_{\mathrm{shift}}\).

The geometric SBP density cancels this same \(F_Q\dot Q\) locally. Its canonical
flux is \((\partial\varepsilon/\partial r_x)\dot r+(\partial\varepsilon/\partial\chi_x)\dot\chi+(\partial\varepsilon/\partial p_{Q,x})\dot p_Q\), split into the shift pieces of those velocities and the remainder. Windowed matter plus gravity reproduces the global Hamiltonian. The difference from the windowed constraint smearing is the summation-by-parts term in \(2L(\partial_x F)/Q\).

## Observer balance

\[
K_r=\frac{\dot q-\partial_x(\beta q)}{Nq},\qquad
K_\perp=\frac{\dot r-\beta r_x}{Nr},
\]

\[
\partial_t(V\rho)+\partial_x\big[V(Nj/q-\beta\rho)\big]
=-NV(p_r K_r+2p_\perp K_\perp)-V(j/q)\partial_x N.
\]

The first term on the right is proper pressure work. The second is momentum
work from a spatial lapse gradient. Neither is \(F_Q\dot Q\).

## What was measured

Band-limited \(N=64\) coordinate closure errors are \(9.92\times10^{-10}\)
(matter), \(9.86\times10^{-10}\) (comoving), and \(-3.07\times10^{-8}\) (matter
plus gravity). Those channels were \(-0.1869\), \(1.3559\), \(-0.3245\), and
\(0.1081\). Dropping any one of them leaves an error of that size. The maximum
gap from the previous coordinate draft is \(5.09\times10^{-11}\).

On the same slice with momenta set to zero, \(K_\perp=0\) and \(K_r=3.64\times10^{-7}\),
while coordinate-work \(\ell_1\) is \(3.005\) and proper-pressure-work \(\ell_1\)
is \(9.19\times10^{-7}\).

The proper-energy residual on that slice has maximum \(8.66\times10^{-6}\).
The shell identity \(F_L/r=V\rho\,dx\) holds to \(5.6\times10^{-17}\). The
Hamiltonian chain \(NF_N+qF_q+rF_r\) holds to \(5.6\times10^{-17}\).

Manufactured Galerkin columns, \(n_f=32\), \(n_q=128\): unprojected coordinate
closure \(2.08\times10^{-8}\), lifted closure \(-2.072\), leakage \(-2.072\).
Proper residual maximum \(8.94\times10^{-7}\) unprojected and \(2.472\) on the
lifted rate.

The physical packet's frozen-geometry coordinate flux alias is \(55.39\) at
\(N=32\) and \(8.06\) at \(N=64\). That fall is measured aliasing, not a
continuum proof.

Saved v5 \(n_f=512\), \(n_q=2048\), initial state only, then two steps
\(dt=5\times10^{-4}\) (\(t=0.001\), not \(T=0.005\)). Quadrature radius matches
the npz. Global coordinate energy \(-8.60\times10^{-13}\). Summed matter
\(24.7500000000\), summed gravity the opposite, summed normal-observer shell
energy \(10.0121377379\). Partition gap \(2.2\times10^{-16}\).

| Window | Normal shell \(F_L/r\) | Lapse-weighted \(LF_L\) | Coordinate \(e\) | Coordinate work | Proper pressure work |
|---|---:|---:|---:|---:|---:|
| \([0,2]\) | \(5.9878147633\) | \(11.8361697398\) | \(11.8361697398\) | \(-2.7645533449\) | \(\sim10^{-13}\) |
| \([2,4]\) | \(3.4982391842\) | \(11.7044806298\) | \(11.7044806298\) | \(-4.7446625387\) | \(\sim10^{-13}\) |
| \([4,6]\) | \(0.1415897795\) | \(0.6737247199\) | \(0.6737247199\) | \(-0.3857563048\) | \(\sim10^{-14}\) |
| \([6,8]\) | \(0.3844940108\) | \(0.5356249105\) | \(0.5356249105\) | \(-0.0854555064\) | \(\sim10^{-15}\) |

Initial \(K_r\) and \(K_\perp\) are about \(10^{-10}\), so the large coordinate work
is shift drag of \(Q\). Unprojected proper-window closure stays below
\(2.5\times10^{-9}\); the lifted proper-window gap stays below \(2.2\times10^{-9}\).
Coordinate lifted closure stays below \(2.0\times10^{-8}\). After two steps the
global coordinate energy is \(-7.40\times10^{-10}\) and \(r,Q\) stay positive.
Regional coordinate totals move by \(0.00378\), \(-0.01305\), \(0.00597\), and
\(0.00330\).

## Limits

The packet alias and the manufactured Galerkin lift are explicit extra
exchange. The v5 state is the saved diagnostic initial slice; this note does
not rerun \(T=0.005\) or change its full-\(C\) status. Historical vacuum
branch \(B\) is not a prerequisite. Embedding stays optional. Continuum
renewal is not claimed. The saved \(T=0.05\) episode is the separate
record `results/development/nsc-spherical-feedback-episode-v1.json`. This
note does not rerun it.

```sh
python scripts/lab.py -m pytest tests/test_nsc_regional_energy_exchange.py -q
```
