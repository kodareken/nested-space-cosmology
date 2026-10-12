# Principal coefficients on the compatible incoming surface family

The principal matrix at the unchanged baseline is approximately

$$
\begin{pmatrix}
0.0573992775862764&-0.178743096482006\\
-0.0297051917511208&0.0416979944511213
\end{pmatrix}.
$$

Its directed determinant interval is
`[-0.0029161631968961416,-0.0029161631968961255]`. After eliminating U,
the coefficient of `w'''` lies in
`[-0.05080487628982585,-0.050804876289825454]`. It is nonzero: there is
no principal third-order cancellation in this included constraint system.

Reuse the fixed incoming surface, the established local/reference action,
the closed lapse coefficient c_u, and the full spatial Weyl projector.
The new family is `delta r=T*w(z)+T^3*U(z)/6`, `delta a=0`; at T=0 its
six slots are `w,w_z,w_zz,w_zzz,U,U_z`. No spatial interval, boundary,
length, history or physical initial data is selected.

The missing connection is the principal matrix at the original baseline:

$$
\begin{pmatrix}
\partial_U\mathcal E_N&\partial_{w_{zz}}\mathcal E_N\\
\partial_{U_z}\mathcal E_\beta&\partial_{w_{zzz}}\mathcal E_\beta
\end{pmatrix}.
$$

Its first entry is the already-certified c_u and is reused. The other
entries must retain the spatial Weyl terms and full two-dimensional Euler
derivatives. Stop after an exact/directed principal relation or a precise
mathematical obstruction. Do not rerun c_u/c_v, old reference quadratures
or baseline probe generators. A zero determinant is a structural cancellation,
not a non-existence result; it redirects the surface problem to lower order.

Risks are a mixed-derivative factorial, dropping third spatial Weyl order,
and mistaking a point-germ relation for a compatible surface solution.
The coefficients are action responses, not new couplings or refit parameters.

## New local Euler entries

Write the entries as `[[A,B],[d,e]]`, where A is the saved c_u certificate.
The original two-jet density depends on `N_T,N_z,N_zz` in the lapse
variation and on `beta_Tz` in the highest shift variation. The actual full
Euler operator therefore extracts

$$
B=-L_{N_T,r_{zz}}-L_{N_z,r_{Tz}}+L_{N_{zz},r_T},\quad
d=L_{\beta_{Tz},r_{TT}},\quad e=L_{\beta_{Tz},r_{zz}}.
$$

These are derivatives of the density, before restriction to the background;
the mixed coordinate derivative appears once. A separate total-derivative
jet calculation verifies these formulas. Exact-degree four-corner density
stencils check only these new entries against the unchanged action owners;
they do not rerun the old c_u or c_v responses.

For `Ha=a_T/a`, `Hr=r_T/r`, define

$$K=-\frac{32\pi C_W}{3}+\frac{h_q+\log r}{30\pi}.$$

The summed local entries are

$$
B_{\rm loc}=\frac r a\left[K(2H_a-3H_r)+\frac{2H_a+H_r}{60\pi}\right],
\quad d_{\rm loc}=ar\left[K+\frac1{60\pi}\right],
\quad e_{\rm loc}=-d_{\rm loc}/a^2.
$$

The constant Euler, Einstein, Maxwell, cylinder, LLL geometry, WZ Euler and
WZ gauge channels have zero contribution to these highest entries. Their
lower-order responses are not discarded.

## Full spatial Weyl reference entries

Formal labels `omega` and `q_formal` extract derivatives; they are not chosen
physical frequencies or lengths. For a linear variation proportional to
`exp(-i*omega*T+i*q_formal*z)`, the owned spatial Weyl product shifts the
background Hamiltonian to `H(k+q_formal/2)` on the left and
`H(k-q_formal/2)` on the right. Its exact off-band response satisfies both
linear transport and projector idempotency. Expanding the scalar current
through the `omega^3*q_formal` and `omega*q_formal^3` terms retains the
third spatial Weyl contribution.

With `p=k/a`, `L=ell/r`, `M2=m^2+L^2`, the corresponding point kernels are

$$
d(p)=-\frac{L^2p^2}{16r(p^2+M2)^{7/2}},\qquad
e(p)=\frac{L^2p^2(5M2-2p^2)}{32ra^2(p^2+M2)^{9/2}}.
$$

Full-line integration with the existing `D/(2*pi)`, `dk=a*dp`, gives

$$d_g=-\frac{D a L^2}{120\pi r M2^2},\qquad e_g=-d_g/a^2,$$

where D is the complete copy/degeneracy multiplicity. The independent
[reference lapse helper](../src/recursive_horizons/nsc_incoming_surface_reference_lapse.py)
retains all surviving order-four projector-idempotency traces and gives

$$B_g=\frac{D L^2[(3H_r-2H_a)M2+8H_rm^2]}{120\pi a r M2^3}.$$

Its third-order point term is odd in p and has zero full-line integral.
Both reference calculations use standard finite beta-function moments,
without preparing new momentum quadrature or physical modes.

## Exact relation and directed decision

Combining the new entries with the **reused** closed c_u equations gives
the exact per-channel identities

$$B=-\frac{2A+H_rd}{a^2},\qquad e=-\frac d{a^2},$$

and hence

$$
\det\begin{pmatrix}A&B\\d&e\end{pmatrix}
=\frac{d(A+H_rd)}{a^2},\qquad
e-dB/A=\frac{d(A+H_rd)}{a^2A}.
$$

The first matrix entry is imported from its prior directed certificate.
The new d interval is `[-0.029705191751120803,-0.02970519175112073]`.
Directed evaluation of the new entries and their exact relation gives the
strictly negative determinant and eliminated coefficient quoted above.
No old coefficient proof, reference quadrature or source calculation runs.

This removes the principal obstruction to the local compatible-function
ODE near the baseline. It is not a selected physical root, a global Cauchy
surface, parent preparation, endpoint matching or extended stationarity.
No domain length, boundary condition, metric step or Gamma_rest term is added.

```sh
python3 scripts/derive_nsc_incoming_surface_principal.py --prepare
python3 scripts/derive_nsc_incoming_surface_principal.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_surface_principal.py
```
