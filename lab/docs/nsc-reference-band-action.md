# Reference subtraction from its effective band action

The existing fourth-order projector now supplies an action density, including
its time connection. This completes a formal action-level construction before
finite regulator matching; it does not create a new physical covariance or
select a transmitting geometry.

## One reference band, with its connection

The [spatial symbol](nsc-spatial-reference-symbol.md) already obeys
$P\star P=P$ and $i\varepsilon\partial_TP=[H,P]_\star$ through order four.
Apply the imported effective-band construction of
[Panati, Spohn and Teufel, section 4.4](https://arxiv.org/abs/math-ph/0201055).
For a fixed rank-one reference chart $\Pi$, form

$$
S=\Pi P+(I-\Pi)(I-P),\qquad
U=S\star(S^\dagger\star S)^{-1/2}_{\star}.
$$

The inverse square root is a **star inverse**, evaluated order by order.
An ordinary pointwise inverse would omit spatial derivatives. The code
retains both unitarity identities, the intertwining relation
$U\star P\star U^\dagger=\Pi$, and normalization-commutator residuals.
The two gapped signs of the massless zero-angular channel use their separate
fixed band charts; no zero-overlap denominator is regularized.
Their reference projector is exactly constant because the owned LLL operator
is chiral and diagonal. That identity is used directly to avoid differentiating
roundoff in a generic square-root ratio. Its difference from the stored generic
point calculation is reported; no massive channel or physical state is changed.

Conjugating the already owned operator $H-i\varepsilon\partial_T$ fixes the
time-connection sign:

$$
h=U\star H\star U^\dagger
+i\varepsilon(\partial_TU)\star U^\dagger.
$$

The occupied-band subtraction density is $\mathrm{tr}(\Pi h)$, with positive
sign in $S_{\mathrm{sub}}$. Its zero-order value is the negative-normal-band
energy. The subtraction CTP phase is the branch difference
$S_{\mathrm{sub}}[g_+]-S_{\mathrm{sub}}[g_-]$; the physical Gaussian action
and $C_0$ remain unchanged. Here $U$ is a formal reference-band frame, **not**
$U_{L0}$, $V_c$, or a completed physical `EndpointBranchJets` object.

## Vary the action before discarding boundary terms

For $D=\delta U\star U^\dagger$, the executable variation identity is

$$
\delta h=U\star\delta H\star U^\dagger
+[D,h]_{\star}+i\varepsilon\partial_TD.
$$

The record keeps the action derivative, the reference vertex
$\mathrm{tr}(P\star\delta H)$, the time derivative, and the remaining
star-trace exchange separately. Their sum is checked for all four raw-KS
amplitude directions. The state/time convention is the same one that the
preceding spatial-symbol record matched to the homogeneous reference.

Under a cyclic full phase-space trace and fixed endpoint variations, the
time and star-trace terms integrate to boundary contributions, giving the
required positive reference vertex. They are **not** dropped at a finite
momentum endpoint or under a nonconstant regulator. Compact spacetime support
alone does not make a momentum-weighted star trace cyclic. The reported
local exchange term is not itself an integrated cutoff mismatch or stress.

## Finite matching remains part of the same action

The [vacuum matching owner](nsc-vacuum-matched-ctp.md) specifies the conversion
$B^E=\Gamma^E_{\mathrm{heat}}-\Gamma^E_{\mathrm{canonical,vac}}$.
The locked local coefficients use the compact complement
$H_\mu=h_\Lambda-E_1(y/\mu^2)$, not an arbitrary momentum weight on this band
symbol. Independent lapse, shift, measure and regulator variations must
therefore enter the general conversion before source composition.

For a common channel inventory the target identity is

$$
\left[\Gamma^E_{G,\mathrm{vac}}-
\Gamma^E_{\mathrm{ref,matched}}+
\Gamma^E_{\mathrm{local,locked}}-
\Gamma^E_{\mathrm{heat,full}}\right]_{\le4}=0.
$$

The subtraction phase above is the negative of the canonical reference
contribution; its finite heat-scheme conversion has not been evaluated here.
The one-light-field matching cannot silently be assigned to the positive
compact canonical fields. No coefficient is refitted or presumed zero.

The [record](../results/development/nsc-reference-band-action.json) covers
the inherited local control point and paired control momenta. It keeps the
formal action gate separate from regulated matching, the complete source,
physical endpoint selection and stationarity. No extra
$\Gamma_{\mathrm{rest}}$, physical stress or metric evolution is introduced.

```sh
python3 scripts/derive_nsc_reference_band_action.py --prepare
python3 scripts/derive_nsc_reference_band_action.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_reference_band_action.py
```
