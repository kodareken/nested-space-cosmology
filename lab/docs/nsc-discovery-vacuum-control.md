# Locked source-free Bertotti–Robinson control

This conditional control uses the existing conformal action and its original
discrete rates. It prepares a homogeneous magnetic-radius state with zero
canonical matter covariance. Its initial geometry differs from the occupied
discovery family, so it does not isolate the effect of removing matter from
that family's saved initial state. The locked induced and magnetic action
terms remain present. “Source-free” here means zero canonical matter source;
it is not an absolute continuum vacuum or a sea subtraction.

Owners:

- [Module](../src/recursive_horizons/nsc_discovery_vacuum_control.py)
- [Creation/check driver](../scripts/derive_nsc_discovery_vacuum_control.py)
- [Focused tests](../tests/test_nsc_discovery_vacuum_control.py)

## Derivation from the locked action

The existing action has
\[
F=-4\pi A r^2+2\alpha\chi,\quad F_\chi=2\alpha,\quad
Z=-24\pi A,\quad
V=8\pi A r^2-\alpha(4\chi+\chi^2)-2\pi C_F f^2,
\qquad \alpha=-4\pi C_W/3.
\]
The audited coefficients select
\[
r_m^2=\frac{C_F f^2}{4A}=1,\qquad Q_0=b_0/a_0\simeq0.2513249333.
\]
In the conformal gauge \(L=Q,\beta=0\), put \(r=r_m,\chi=0,p_Q=0\).
The action's original momenta then give
\[
p_\chi=2F_\chi D=4\alpha D,\qquad
p_r=2F_rD=-16\pi A r_mD,\qquad
\Pi=p_r-(F_r/F_\chi)p_\chi=0,\quad D=\dot Q/Q.
\]
Here \(V=0\). Both pre-gauge constraints and all canonical matter forces
vanish. The original conformal rates give
\[
\dot p_\chi=-4\alpha Q^2,\qquad
\dot p_r=16\pi A r_mQ^2,\qquad
\dot D=-Q^2,\qquad \dot Q=QD.
\]
With all momenta zero at \(t=0\), the solution is
\[
Q=Q_0\operatorname{sech}(Q_0t),\qquad
D=-Q_0\tanh(Q_0t),\qquad D^2+Q^2=Q_0^2.
\]
The total discrete Hamilton energy has the exact target zero. No coefficient,
pressure force, multiplicity, constraint, or radius equation is replaced.
No initial-radius solve or post-step projection is called.

## Metric and source domain

For the actual metric \(g=(rQ)^2(dt^2-dx^2)-r^2d\Omega^2\), signature \(+---\),
the exact control has
\[
R_h=2,\quad R_4=C^2=0,\quad R_{ab}R^{ab}=4/r_m^4,\quad K=8/r_m^4,
\]
and normal covariant tides \(R_{0101}=-1/r_m^2\), \(R_{0202}=R_{0303}=0\).
The normal worldlines are geodesic in this homogeneous chart. Their proper
time from the turnover is \(r_m\arctan[\sinh(Q_0t)]\).
The period and child proper lengths are \(8r_mQ\) and \(2r_mQ\).
They contract although the areal radius and all these four-dimensional
invariants and tides remain constant. Proper-length or \(Q\) contraction
alone therefore does not establish physical instability. The occupied
family's actual radius and tide growth provide additional measurements;
this control does not decide that family's stability.

The original rate API expects six columns. Six zero spinor columns and six
zero occupations give \(C=0\) and the CAR bound \(0\le C\le I\), with every
covariance eigenvalue zero. The occupied rank-six packet's Gram-identity
preparation requirement does not apply to this degenerate source. Packet
localization consumers are not used on zero probability.

Integration uses the native Galerkin CauchyState: coarse nodal momentum
densities \(p_Q,p_r,p_\chi\). The NPZ also labels canonical momenta in the
identity geometry frame, \(\pi_g=\Delta x_g p_g\). These are distinct from
the mechanical combination \(\Pi\) above and are never fed into a nodal
rate without decoding.

## Bounded execution and evidence

At nf=32 and nf=64, direct instantaneous residual checks use the closed
solution at \(T=0,1,3,8\). A separately composed RK4 method calls the full
original discrete rates for 300 steps of \(dt=0.01\) through \(T=3\), retaining
states at 0, 1 and 3. No exact solution is imposed on its stages. Actual
metric jets are reconstructed from each state's projected rates, with
Fourier spatial derivatives and no substitution of \(\chi\) for curvature.
The record compares state and rate errors, constraints, energy, metric
invariants, tides, \(\Pi\), and \(D^2+Q^2\) in their own units.

One numerical thread and a shared thirty CPU-second budget bound the whole
control. The longer \(T=8\) check is instantaneous; it is not an evolved
campaign or an extension of the short control. Source hashes, the unchanged
locked coefficient input, and the NPZ payload hash bind the record.
JSON and NPZ each have a one-MiB size bound.

From the repository root:

~~~sh
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_vacuum_control.py
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_vacuum_control.py --write
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_vacuum_control.py --check
.venv/validation/bin/python scripts/lab.py -m pytest tests/test_nsc_discovery_vacuum_control.py -q
~~~

The default is a preview and writes no files. Explicit creation writes
record.json and payload.npz only in a new directory under
lab/results/development/nsc-discovery-vacuum-control-v1. An optional path
after --write or --check selects a descendant. Existing directories are
never overwritten. Checking reads bindings and arrays without evaluating
rates or evolving a state. Production creation follows the coordinating
agent's source-freezing commit; no production record is created by the tests.
