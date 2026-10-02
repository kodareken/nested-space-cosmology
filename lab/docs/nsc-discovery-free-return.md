# Saved-field nearly free carrier return

The [consumer](../src/recursive_horizons/nsc_discovery_free_return.py) compares
an exact free matrix transport of the actual saved T3 field with its saved T8
continuation. The [driver](../scripts/derive_nsc_discovery_free_return.py) previews
or creates/checks a new record, and the [tests](../tests/test_nsc_discovery_free_return.py)
check the AP identity, fixed spatial window, saved join, operator, covariance,
and immutable replay. No trajectory, solver, or future field evolution is run.

## Exact operator and saved representation

The same conformal field operator has \(L=Q,\beta=0\), hence
\[
H(t)=\sigma_2P+\kappa\sigma_1 U_f^\dagger\operatorname{diag}(Q(t))U_f,\qquad
P=U_f^\dagger P_{\rm fine}U_f.
\]
The representative generator carries no angular multiplicity factor.
The free comparison is
\[
U_0(\Delta t)=e^{-i\sigma_2P\Delta t}.
\]
The original AP carrier has period \(L=8\) and wavenumbers
\(p_m=2\pi(m+\tfrac12)/L\). Therefore, as an operator on the retained band,
\[
U_0(8)=-I,\qquad C(8)=C(0).
\]
The sign cancels in every quadratic measurement. Every fixed spatial window
probability returns after one period for any source columns and weights;
this identity requires no specially prepared packet or renewal.

The implementation diagonalizes the actual retained Hermitian \(P\), records
its gap from the owned AP spectrum, and applies cosine/sine factors to all six
spinor columns. It retains the complete covariance
\(C=\Phi\,\operatorname{diag}(c)\Phi^\dagger\), including cross terms.
The fixed child interval is \((1,3)\), using the existing observer's trigonometric
window integral of weighted half-density mass. It is not a modal projector.

Saved geometry arrays are coefficients in the frozen W frame. Saved momenta
are canonical pi and additionally require the original spatial-spacing decode.
The consumer calls episode.state_from_arrays, then nested.reconstruct_state
before galerkin.prolong_state. Fermion columns remain in their original AP
band. Prolonging W-encoded Q directly would produce the wrong metric and is
not a valid replay.

## Measured relation and its domain

The primary input is nf256 coupled half-step, episode-v1 ordinal 2 at T3,
identical in every per-array hash to crossing-v1 ordinal 0. The target is
crossing ordinal 1 at T8. The consumer applies \(U_0(5)\) to the actual T3
columns and compares the result with the complete actual T8 columns,
the source covariance, and the fixed child window.

The primary weighted field relative error is approximately
\(2.982\times10^{-6}\), while the fixed-child probability difference is about
\(10^{-12}\). These saved-data comparisons attribute the coupled return to
nearly free carrier transport in its contracting conformal geometry; they do
not establish renewal or a new gravity solution. Optional retained controls
include nf128, both saved timestep caps, and the nf256 frozen-geometry field.
A frozen return alone does not prove the free limit: its complete field is
compared separately with the same exact free operator.

The record checks source Gram matrices and CAR eigenvalues, weighted norms,
the exact period return, a zero-duration identity, and the actual saved
Hamiltonian's free-plus-mass decomposition. An algebraic zero-mass operator
check on the saved positive metric confirms the generator convention; it
does not change any production coupling coefficient or state.

## Conditional Duhamel statement

Along the realized metric, unitary field propagators give
\[
\|\Phi(8)-U_0(5)\Phi(3)\|_c
\le |\kappa|\,\|\Phi(3)\|_c
\int_3^8\|Q(t)\|_\infty\,dt .
\]
No Grönwall factor is needed. The weighted norm retains the six source
occupations. This is a conditional statement about the exact semidiscrete
field equation on that realized metric, while saved RK4 drift is a separate
numerical effect.

The record reports decoded endpoint samples of \(\|Q\|_\infty\). Those samples
are not an upper bound on the integral. The required intervening Q history
and an integral upper bound remain unknown to this consumer, so its certified
integral and certified Duhamel error-bound fields are null.

## Provenance and commands

The crossing observed-run binding supplies its authenticated producing commit
b7f0dab and records that it is a post-run coordinator observation. The original
1c9e770 physics envelope is verified separately. All consumed checkpoint JSON,
NPZ, observation files, the exact AP numerical owners, and the four new consumer
files are bound with repository-relative SHA256 paths. Original records are
preserved. The exact joins also check the frozen W and six weights; no source,
clock, geometry, or covariance reset is introduced.

One numerical thread and a ten CPU-second budget bound the saved-evidence
calculation. Default controls are selected; --no-controls selects the primary
case only. From the repository root, after freezing these sources:

~~~sh
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_free_return.py
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_free_return.py --write
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_free_return.py --check
.venv/validation/bin/python scripts/lab.py -m pytest tests/test_nsc_discovery_free_return.py -q
~~~

Explicit creation writes only a new
lab/results/development/nsc-discovery-free-return-v1.json and refuses an existing
file. Checking authenticates sources and inputs and reapplies the exact free
matrix operator without an ODE, trajectory, or write. The default previews
without creating a production record.
