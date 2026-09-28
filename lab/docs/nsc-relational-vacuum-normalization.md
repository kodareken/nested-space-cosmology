# Recursive zero-tadpole normalization determines the local vacuum term

The missing coefficient can be fixed by a law already native to the framework:
only a difference relative to another physical state has local meaning. Apply
that principle to the metric one-point function of the complete common action.

Let \(\mathcal P_0\) extract the homogeneous, local, zero-derivative metric
tadpole after every fermionic, bosonic, gauge, boundary and measure contribution
has been assembled. Define

\[
\boxed{
\Gamma_{\rm rel}=(1-\mathcal P_0)\Gamma_{\rm one},
\qquad
\left.\frac{\delta\Gamma_{\rm rel}}{\delta g_{\mu\nu}}
\right|_{\rm homogeneous,\,unlinked}=0.}
\]

The same condition is imposed in every room after expressing it in inherited
dimensionless units. This is the recursive normalization law.

## What the projector removes

In the local coefficient basis

\[
\Gamma_{\rm local}=\int\sqrt{-g}\,(A R-CF^2-V+\cdots),
\]

the projector is

\[
\mathcal P_0(V,A,C)=(V,0,0),
\qquad
(1-\mathcal P_0)(V,A,C)=(0,A,C).
\]

It removes the local zero-momentum gravitational tadpole at fixed internal
spectral data. It does not divide the
complete determinant by a disconnected determinant, which would also remove
useful local response. Einstein, gauge, curvature, Casimir, link and nonlocal
terms remain in the action.

When a compact length is dynamical, its topology-dependent finite Casimir
remainder remains a separate function of that length. The projector acts only
on the local \(a_0\) density; it does not erase this finite modulus potential.

The subtraction coefficient is not independently chosen. It is exactly the
homogeneous \(V\) produced by the completed same action:

\[
S_{0,\rm CTP}=-V_{\rm full}
\left(\operatorname{Vol}[g_+]-\operatorname{Vol}[g_-]\right).
\]

Equal histories remain normalized. Direct sums remain additive. When
\(V\mapsto\Omega^4V\) and \(\operatorname{Vol}\mapsto
\Omega^{-4}\operatorname{Vol}\), the subtraction is invariant. The projector
is idempotent and is the unique linear map on \((V,A,C)\) that annihilates the
vacuum basis while preserving the Einstein and gauge bases.

This is a finite gravitational normalization condition allowed by the local
stress-tensor ambiguity. The raw spectral action still has its nonzero
\(a_0\) coefficient; the NSC law specifies which part of that raw action is
the relational gravitating source. The spectral action and its volume term are
standard [Chamseddine–Connes](https://arxiv.org/abs/hep-th/9606001); local
covariance and stress conservation retain finite local gravitational
renormalizations [Hollands–Wald](https://arxiv.org/abs/gr-qc/0404074).

## MMP gate

Applying the law to the retained coefficients gives

\[
V_{\rm full}^{\rm rel}=0,
\qquad
\lambda_4=0,
\qquad
\Xi=\lambda_4 r_e^2=0.
\]

The induced \(G_N\), gauge normalization, charge radius and imported MMP
throat length are unchanged. Thus both the charged seed bound and the
asymptotically flat vacuum condition pass without fitting a cosmological
coefficient.

The physical interpretation is direct: homogeneous unlinked vacuum defines
the zero of relational gravitational energy. Curvature and boundary gradients
remain sources. Consequently the local \(a_0\) term cannot select \(\Omega\),
\(\zeta\) or \(\Phi\); their values must come from the recursive boundary and
gradient response.

## Remaining MMP conditions

This closes the absolute-vacuum coefficient gate. A full physical embedding
still needs the complete induced \(A,C\) normalization, a closed magnetic
return domain with its covariance, and a decision whether the lowest-Landau
sector remains massless after the physical link is included.

## Reproduction

```sh
python3 scripts/derive_nsc_relational_vacuum.py --check
```

The calculation uses the three preceding gate records and a finite coefficient
projector. It does not read the preserved untracked normalization/Weyl files,
rerun MMP gravity, or execute a historical scientific generator.
