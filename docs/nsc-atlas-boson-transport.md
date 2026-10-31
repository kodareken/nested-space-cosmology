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
symplectic spectrum and Gaussian entropy are preserved. An invertible M
cannot map two distinct full incoming covariances to the same finite-time
state; coarse local limits remain a separate question. No instantaneous
vacuum reset or fitted fluctuation width is introduced.

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
