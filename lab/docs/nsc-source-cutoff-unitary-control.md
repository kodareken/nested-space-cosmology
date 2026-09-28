# A globally unitary control for the source-cutoff limit

An energy-dependent phase applied separately to each source column does
not by itself define a unitary operator on the spatial Hilbert space.
It therefore cannot prove the absence of a cutoff/coincidence exchange
in the [common-kernel pairing](nsc-common-subtracted-ks-source.md).
The following exact control isolates that issue without modifying NSC.

On the circle use the negative Fourier basis
`e_n(z)=exp(-i*n*z)/sqrt(2*pi)`, `n>=1`. Let the real symmetric Jacobi
operator have entries

\[
B_{n+1,n}=B_{n,n+1}=\frac1{2(n+1/2)},\qquad n\ge1.
\]

Its row sums are bounded, so B is bounded self-adjoint and `U(t)=exp(i*t*B)`
is unitary. Extend U by the identity on the complementary Fourier sector.
It commutes with the full negative-sector projector, so the complete
covariance and its current are unchanged, exactly.

Let `P_N` instead retain only source labels `1,...,N`. At `t=0`,
`d(U P_N U†)/dt=i[B,P_N]` has just two nonzero entries, at `(N+1,N)` and
`(N,N+1)`. Evaluating their separated kernel at `z+eta/2,z-eta/2` gives

\[
\dot C_{N,\eta}(z)=
\frac{\sin z}{2\pi(N+1/2)}e^{-i(N+1/2)\eta},\qquad
-i\partial_\eta\dot C_{N,\eta}(z)=
-\frac{\sin z}{2\pi}e^{-i(N+1/2)\eta}.
\]

The kernel tends uniformly to zero. Its differentiated coincidence value
does not: it is `-sin(z)/(2*pi)` at every N. Thus convergence of the kernel
alone cannot justify differentiating after setting the split to zero.
For completeness, Abel averaging over cutoffs with weights
`(1-r)*r^(N-1)` multiplies the current by

\[
\frac{(1-r)e^{-3i\eta/2}}{1-r e^{-i\eta}}.
\]

At each fixed nonzero split modulo `2*pi`, this tends to zero as `r->1`.
The reverse order gives one. This is an explicit order-of-limits example,
not an assertion that the unaveraged oscillatory sequence converges.

The [executable control](../scripts/check_nsc_source_cutoff_unitary_control.py)
checks the identities exactly and independently differentiates finite
matrix exponentials of the Hermitian generator. No physical source,
geometry, coupling or stress value enters it.

For NSC this proves only that global unitarity does **not** eliminate the
named source-cutoff/coincidence hole on its own. It neither determines nor
licenses any proposed NSC correction. That coefficient and its remainder
must follow from the actual Dirac evolution and the owned subtraction.
No term is added to the incoming constraints or to `Gamma_rest`.

```sh
python scripts/check_nsc_source_cutoff_unitary_control.py --check
```
