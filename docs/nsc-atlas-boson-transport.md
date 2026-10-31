# Inherited covariance on the evolving volume--matter parent

The [nonlinear parent](nsc-atlas-parent-gh.md) now supplies a finite evolving
geometry and mean fields. The next connection is its state-independent
quadratic matter propagator, followed by covariance inheritance. The prior
Euclidean preparations, finite energy minimizers and late local limits
already have owners; another vacuum prescription would not supply finite
ancestral data.

## Action and finite question

Keep the canonical pair S=sqrt(84)(s-s_star), chi and the same potential
W=V(s)+mu^2(s)chi^2/2. Its Hessian is

$$
W_{AB}=\begin{pmatrix}
(V''+72\mu^2\chi^2)/84&-12\mu^2\chi/\sqrt{84}\\
-12\mu^2\chi/\sqrt{84}&\mu^2
\end{pmatrix}.
$$

This interaction is inherited from the existing action. Fluctuations are
Gaussian matter probes on its classical metric; metric fluctuations,
quantum backreaction and a continuum Hadamard limit are outside this
finite calculation. Canonical quantization is imported explicitly.

Use the checked epsilon=.25, degree16 parent through its existing proper
endpoint4.364851930252. A passive sampler stores actual DOP853 dense
polynomials. It cannot steer adaptive field steps. Five sparse snapshots
are insufficient for the oscillatory Hessian. One background continuation
with the sampler is the required input; there is no new physical grid
campaign or changed source.

## Canonical projection and composition

For reference-orthonormal real spatial modes, project canonical field and
momentum density. The Hamiltonian is

$$
H_2=\tfrac12q^TKq+p^TBq+\tfrac12p^TAp,\qquad
G=\begin{pmatrix}K&B^T\\B&A\end{pmatrix},\qquad \dot M=JGM.
$$

A integrates N/sqrt(det h), K integrates the spatial-gradient and Hessian
terms with N sqrt(det h), and B integrates beta^i D_i. A is a direct
canonical projection, not the inverse of a truncated Lagrangian mass.
The physical normalization r^(3/2) of harmonics cancels the common r^-3
reference measure. The spatial derivatives retain their factor r.

The field-exchange sector S(++), chi(--) contains the original two source
harmonics. Probe degrees2,4,6 test its finite spectral projection; they
are not complete continuum entropies. Commutator-free fourth-order Magnus
steps preserve the Hamiltonian structure without projecting the answer
back onto the symplectic group.

For an inherited bosonic covariance Gamma_parent,

$$
\Gamma(t)=M(t)\Gamma_{\rm parent}M(t)^T,\qquad M^TJM=J.
$$

The retarded response is state independent. Noise and reduced entropy
depend on the incoming covariance and its cross correlations. The full
retained-sector symplectic spectrum and Gaussian entropy are preserved. An invertible M
cannot map two distinct full incoming covariances to the same finite-time
state; coarse local limits remain a separate question. No instantaneous
vacuum reset or fitted fluctuation width is introduced.

There is an exact onset connection. For an initial chi-momentum kick f,
with both field displacements and S momentum zero, successive canonical
equations give q_S(0)=q_S'(0)=q_S''(0)=0 and

$$
q_S(t)=-\frac{(t-t_0)^3}{6}A_S(t_0)K_{S\chi}(t_0)A_\chi(t_0)f
+O((t-t_0)^4).
$$

The kinetic and shift matrices are field diagonal, so their time
derivatives cannot create an earlier cross-field term. Nonzero inherited
mixing gives a covariance-independent response. The observer's induced
delta log(g_KK) is -7 delta S/sqrt(84), evaluated on the fixed background
geodesic. This excludes metric/observer-displacement perturbations.

The finite harmonic projection is retarded in time; it does not establish
the exact spatial light-cone support of a continuum field. The prescribed
source harmonics have global spatial support.

## Geometry-only local onset

The [independent initial-slice check](../results/summaries/nsc-atlas-boson-onset-v1.json)
uses the field equation directly, without the projected Hamiltonian. For
a unit canonical momentum-density profile f in the physical reference
measure, normal(deltachi)=f/sqrt(det h). The initially normal background
observer therefore measures

$$
\delta S(\tau,X_0)=-\frac{W_{S\chi}(0,X_0)f(X_0)}{6\sqrt{\det h(0,X_0)}}\tau^3
+O(\tau^4).
$$

Scalar displacement and both S initial data vanish. Spatial, expansion
and acceleration terms act on these vanishing data at the first two
orders. Differentiating the potential force supplies the displayed third
derivative. Lapse factors cancel on changing to the initial observer's
proper parameter. The local coupling response is -7 delta S/sqrt(84).

For the existing epsilon=.25, degree16 slice, W_Schi=-10.8875320314 and
sqrt(det h)=.976720138154 at the prescribed observer. The two normalized
physical-reference profiles give local delta-log(g_KK) cubic coefficients
-.0368819999195 and-.0921662836119. These are per unit canonical impulse
in model units, not a selected perturbation amplitude or an SI coupling.
They include background discretization but no finite probe cutoff.

## The same propagator supplies Schur memory

Let x=(q_chi,p_chi) be retained and y=(q_S,p_S) unresolved. Partition
D=JG into these canonical pairs. If U_y solves the homogeneous D_yy
equation, variation of constants gives

$$
\dot x(t)=D_{xx}(t)x(t)+\int_{t_0}^t\mathcal K(t,s)x(s)\,ds+L(t)y_0,
$$

$$
\mathcal K(t,s)=D_{xy}(t)U_y(t,s)D_{yx}(s),\qquad
L(t)=D_{xy}(t)U_y(t,t_0).
$$

This is the [causal reduction programme](schur-memory.md) with vertices
derived from the actual evolving parent. U_y is obtained from the
decoupled eliminated operator, not a subblock of the coupled propagator.
Its force covariance is L(t) Gamma_yy L(s)^T and its correlation with the
initial retained data is L(t) Gamma_yx. The latter cannot be discarded.
Nonstationary memory is not automatically a local vacuum stress or Lambda.

The quadratic-response and covariance machinery is attributed to
[Gaussian canonical-state methods](https://arxiv.org/html/1110.3234v1)
and the [parametric-oscillator environment programme](https://arxiv.org/abs/gr-qc/9312035).
The actual Hamiltonian and vertices here come from this atlas's fields
and metric. Our covariance convention is hbar=1, vacuum nu=1/2; the
review uses hbar=2. Entropy is computed from the symplectic spectrum of
the supplied covariance, not chosen by this map.

## Current finite continuation

The resolved history reproduces the checked degree16 endpoint with phase
gap1.37e-13 radians and density gap6.7e-15 V_star. Seven probe controls at
degrees2,4,6,8,10, a degree6 half-step and a degree8 decoupled control are preserved in the
[response record](../results/summaries/nsc-atlas-boson-response-v1.json).
Time refinement changes the two endpoint response ratios by3.84e-7 and
4.06e-7, much less than spatial refinement. Degree6 leaves22.36%
of the source2 scalar reference-norm power in its last two degrees.
Degree8 reduces that tail to1.578%, while its corresponding volume tail
is3.258%. Its matched decoupled control has zero cross-field response.
The degree8 memory at final/midpoint is nonzero. The final degree10 run
returns reference response ratios .003168467238 and .003023500301,
with degree8--10 spreads .000018088553 and .000065708365. Those changes
contract from .000253824284 and .000235287590 at the preceding comparison.
These are reference norms, not energy fractions. They are empirical
controls; time refinement is at degree6, not degree10. Late pointwise
response remains more sensitive and is not admitted as converged.
The original stronger-amplitude parent and the incoming covariance remain
open. This uses one resolved background, without another physical field grid.

## Implementation and stopping conditions

- `nsc_atlas_gh_history.py`: regenerable dense background, no extrapolation.
- `nsc_atlas_boson_hamiltonian.py`: physical metric, shift and two-field Hessian.
- `nsc_atlas_boson_transport.py`: symplectic propagation and supplied-covariance operations.

The budget is one background history within2000 wall seconds, followed by
finite probe time/degree controls within900 seconds each. The strongest
failure risks are underresolved history, omitted probe harmonics, and an
unacknowledged state choice. Dense-polynomial checks, independent time and
degree comparisons, and separate response/covariance outputs address them.
Failure records its specific cause before any further calculation. A later
independent inherited-state prediction must reproduce the stored retarded
response and its conditional noise without changing the source law.

## Reproduction

The following commands generate one finite response and its matched Schur
control from a clean workspace. Output paths must be fresh. The recorded
history producer is `771c5a78` and the response producer `292eb146`; the compact
records retain exact source hashes and parameters. Repeat the response
command at degrees2,4,6,10 and use `--max-step 0.01` at degree6 to reproduce
the remaining controls.

```powershell
python scripts/run_python.py scripts/evolve_nsc_atlas_parent_gh.py --degree 16 --epsilon 0.25 --threads 8 --budget-seconds 2000 --dense-cache .build/boson-replay-history.npz --output .build/boson-replay-background.json
python scripts/run_python.py scripts/propagate_nsc_atlas_boson_parent.py --history .build/boson-replay-history.npz --background-record .build/boson-replay-background.json --probe-degree 8 --max-step 0.02 --threads 8 --budget-seconds 900 --matrix-cache .build/boson-replay-coupled.npz --output .build/boson-replay-coupled.json
python scripts/run_python.py scripts/propagate_nsc_atlas_boson_parent.py --history .build/boson-replay-history.npz --background-record .build/boson-replay-background.json --probe-degree 8 --max-step 0.02 --threads 8 --budget-seconds 900 --decoupled-control --matrix-cache .build/boson-replay-control.npz --output .build/boson-replay-control.json
python scripts/run_python.py scripts/summarize_nsc_atlas_boson_response.py --input .build/boson-replay-coupled.json --input .build/boson-replay-control.json --history .build/boson-replay-history.npz --coupled-record .build/boson-replay-coupled.json --coupled-matrix .build/boson-replay-coupled.npz --control-record .build/boson-replay-control.json --control-matrix .build/boson-replay-control.npz --output .build/boson-replay-summary.json
python scripts/run_python.py scripts/plot_nsc_atlas_boson_response.py
```

The plot uses the committed small summaries. Complete matrices and dense
polynomials remain local and regenerable. No observation calibrates these
calculations; none of these controls selects the incoming covariance.
