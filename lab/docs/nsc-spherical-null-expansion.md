# Spherical null expansion on the saved \(T=0.05\) chart

Kinematic postprocessing of the stored autonomous episode. The production
step remains `nsc_spherical_galerkin_coupling`. This measurement does not
evolve a state, rebuild the column source, or recompute inheritance. It
does not read \(\chi\). Record
`results/development/nsc-spherical-null-expansion-v1.json`, status
`MEASURED_SAMPLE_TRAPPING_CHARACTER_CHANGES`. Assessment CPU `0.070397` s.
The episode JSON and NPZ bytes are unchanged.

```sh
python scripts/lab.py scripts/derive_nsc_spherical_null_expansion.py
python scripts/lab.py -m pytest tests/test_nsc_spherical_null_expansion.py -q
```

The positive-boost flag compares the expansions returned by the rescaling.
The saved-node check uses the areal-maximum sample and the interpolated
\(r_x\) zeros.

## Equations and sign convention

Signature \(+---\). The saved chart uses the static clock of
`nsc_spherical_coupling` and the normal observer already stored by the
episode. \(L\) is the calibration lapse factor and \(\beta\) is the
calibration shift, both held fixed. The ledger in
`nsc_regional_energy_exchange` uses the same \(N=rL\), \(q=rQ\), and

\[
K_\perp=\frac{\dot r-\beta r_x}{Nr}.
\]

The line element and the null legs are

\[
g=r^2\bigl[L^2\mathrm dt^2-Q^2(\mathrm dx+\beta\mathrm dt)^2-\mathrm d\Omega^2\bigr],
\]
\[
n=\frac{\partial_t-\beta\partial_x}{N},\qquad e_1=\frac{\partial_x}{q},\qquad
k_\pm=n\pm e_1.
\]

Then \(g(n,n)=+1\), \(g(e_1,e_1)=-1\), \(k_\pm\) are null, and
\(g(k_+,k_-)=2\). The legs are not divided by \(\sqrt{2}\). The spherical
expansions are the directional derivatives of the areal radius,

\[
\theta_\pm=\frac{2}{r}\left(\frac{\dot r-\beta r_x}{N}\pm\frac{r_x}{q}\right)
=\frac{2}{r}\,k_\pm(r).
\]

Their product is the inverse-metric square of the areal gradient. Angular
derivatives of \(r\) vanish on each symmetry sphere, so only the \((t,x)\)
block contributes:

\[
\theta_+\theta_-=\frac{4}{r^2}g^{ab}\partial_a r\,\partial_b r.
\]

Under \(+---\), a negative product means a spacelike gradient and opposite
expansion signs. A positive product means a timelike gradient and equal
signs. A zero product is a null gradient. The labels are:

| Sample | Condition |
|---|---|
| Trapped | \(\theta_+<0\) and \(\theta_-<0\) |
| Anti-trapped | \(\theta_+>0\) and \(\theta_->0\) |
| Untrapped | the two signs differ |
| Marginal indicator | an expansion is zero, or it lies inside the declared margin |

A positive reciprocal rescaling \(k_+\to\lambda k_+\),
\(k_-\to k_-/\lambda\) multiplies \(\theta_+\) by \(\lambda\) and \(\theta_-\)
by \(1/\lambda\). The product and each sign stay put. The stored normal
velocity is the episode field `frame_proper`,
\((\dot r_{\mathrm{lifted}}-\beta r_x)/N\), and `frame_K_perp` is that
velocity divided by \(r\). The sum \(\theta_++\theta_-=4K_\perp\) links the
expansions to that observer. The spatial slope \(r_x\) is the periodic
spectral derivative with the Nyquist symbol set to zero, the same symbol as
`periodic_derivative`.

Three manufactured points fix the signs without the saved state.
\(r=2\), \(L=Q=1\), \(\beta=0\):

- \(\dot r=-1/2\), \(r_x=0\) gives \(\theta_\pm=-1/4\), both negative.
- \(\dot r=+1/2\), \(r_x=0\) gives \(\theta_\pm=+1/4\), both positive.
- \(\dot r=0\), \(r_x=1/2\) gives \(\theta_+=+1/4\) and \(\theta_-=-1/4\).

On those points and on eight further positive-chart samples, the
inverse-metric product disagrees with \(\theta_+\theta_-\) by at most
`2.842170943040401e-14`. The null-norm defect is at most
`4.440892098500626e-16`. The boost product defect is at most
`2.842170943040401e-14`. The FFT derivative matches the owner matrix on
64 points, Nyquist included, to `8.362754932988992e-14`.

## Chart and scope

The coordinate \(x\) has period \(8\). Each \(x\) carries a round sphere of
areal radius \(r(t,x)\). The measured object is a scalar pair
\((\theta_+,\theta_-)\) on the quadrature nodes of this finite periodic
\(S^1\times S^2\). A both-negative sample is a local property of the areal
gradient. A global event horizon, an interior, an exterior, black-hole
birth, a child region, and regeneration are outside the measurement.
Eleven frames, spaced by \(0.005\), do not locate a continuum horizon.

The clock on the packet arc is the calibration literal:
\(L(0)=b_0=0.2457850789\), \(\beta(0)=\beta_0=0.325985800803\), and
\(L(2)=b_0\Omega^2\). The bridge on \([4,8]\) is `profile_coordinate`.
\(L\) on the saved nodes runs from `0.19114480326227243` to
`1.599976373372201`, and \(\beta\) from `0.25351616965379487` to
`2.1220555034255035`.

## What the saved samples do

Geometry stays positive on all four runs. On `nf512` at \(dt=0.0005\), \(r\)
starts in \([4.210504072341751,\,4.876819634418227]\). Initial \(Q\) lies in
\([0.25132493326910366,\,0.251324933269177]\), a spread of about \(7\times 10^{-14}\)
about the gauge ratio \(b_0/a_0\). At \(T=0.05\), \(r\) is in
\([4.2123529367075765,\,4.874395995293925]\) and \(Q\) is in
\([0.2195220548716644,\,0.26128339207850754]\). The initial normal velocity
is at most `1.6655044770048675e-10` on `nf512` and
`1.4981042168143584e-08` on `nf256`. The episode metadata matches those
stored endpoint velocities with gap `0`. The historical \(10^{-8}\) flag is
not used as a veto.

At \(T=0\) every quadrature sample is untrapped. The product is negative on
every node; its least negative primary value is `-4.446198031841503e-10`.
The initial areal maximum and minimum are already in the chart. The primary
maximum node has \(\theta_+=3.113501917378642\times 10^{-4}\),
\(\theta_-=-3.113501506382507\times 10^{-4}\), and
\(r_x=9.305250666967035\times 10^{-4}\). Where
\(r_x=0\), both expansions equal \(2v_n/r\). Linear zeros of the stored
slope put those critical values at `8.353813094039343e-12` and
`-1.9636126344543498e-11` on the fine grid, and at about \(10^{-9}\) on the
coarse grid. Both sit inside the margin below, so those interpolated zeros
are marginal within the initial normal-velocity noise. The quadrature node
above is a different sample. The later arcs open on these same extrema.

From the first stored frame \(T=0.005\) through \(T=0.05\), each run has one
both-negative arc containing the areal maximum and one both-positive arc
containing the areal minimum. The complementary nodes stay opposite-sign.
Each later frame has four product crossings, one at each edge of those two
arcs. No stored node has an expansion that is exactly zero.

Primary `nf512`, \(dt=0.0005\), of 2048 nodes:

| \(T\) | Both negative | Both positive | Opposite sign | Both-negative nodes | Both-positive nodes |
|---:|---:|---:|---:|---|---|
| 0 | 0 | 0 | 2048 |  |  |
| 0.005 | 4 | 8 | 2036 | \(x\in[1.5352,1.5469]\) | \(x\in[5.6602,5.6875]\) |
| 0.05 | 45 | 182 | 1821 | \(x\in[1.4296875,1.6015625]\), span `0.171875` | \(x\in[5.09765625,5.8046875]\), span `0.70703125` |

Coarse `nf256` has half the nodes. At \(T=0.05\) it counts 23, 91, and 910.
The both-positive fraction \(91/1024=182/2048\). The both-negative fraction
differs by one coarse-node weight, \(23/1024\) against \(45/2048\). The
coarse both-negative nodes occupy the same endpoint coordinates
\([1.4296875,1.6015625]\). The coarse both-positive arc is
\([5.1015625,5.8046875]\), span `0.703125`, one fine spacing inside the fine
arc. Shared-node signs
disagree at `0` samples on all eleven frames. Halving \(dt\) leaves every
raw count unchanged. The largest \(dt\) movement of \(\theta_+\) is
`2.1634388863711607e-09` on `nf512` and `2.1784517800771397e-09` on
`nf256`. The largest shared-node movement of \(\theta_+\) is
`1.8325977887805045e-07`.

At \(T=0.05\) the interpolated critical values are both negative at the
maximum, `-0.018596390301270596`, and both positive at the minimum,
`0.010981780986768712`. The final primary product runs from
`-0.01586935920781164` to `0.00037204551501067413`.

## Margin

An expansion is strict when its absolute value exceeds ten times the larger
of the measured \(\theta\) comparison gap and the initial floor
\(2|v_n|/r_{\min}\). That margin is `1.8325977887805045e-06`. No sample on
any stored frame falls inside it: strict counts equal raw counts. The
smallest same-sign expansion is `1.6192002416757407e-05`, which is `88.3555`
times the comparison gap. High Fourier modes of stored \(r\), above the odd
geometry band, contribute at most `7.94694272106008e-11` to \(r_x\) on the
fine grid. The stored identity \(K_\perp r=\mathrm{proper}\) misses by
`1.3877787807814457e-17`. Clock null norms on the final fine frame miss
\(g(n,n)=1\) by `6.661338147750939e-16`.

Node counts and four crossings describe the sampled edges. They do not
prove a continuous marginal surface. The sub-grid position of each product
zero is known only to one spacing, `8/2048` or `8/1024`.

## Resolved and open

Resolved on these samples, with \(r>0\) and \(Q>0\): the quadrature
classification changes from entirely opposite-sign at \(T=0\) to one
both-negative neighborhood of the existing areal maximum and one
both-positive neighborhood of the existing areal minimum by \(T=0.005\),
and that split remains through \(T=0.05\). The \(dt\) pairs agree in counts.
The `nf512` and `nf256` shared nodes agree in sign. The product, null norms,
and positive boost identity hold on the manufactured controls. Source hashes
of the coupling, Galerkin owner, regional ledger, feedback action, conformal
source, and episode driver match `hashes_after` in the episode JSON, and
the NPZ hash matches the hash recorded there.

Open: the instant inside \((0,0.005)\) when the product first becomes
positive; the continuum coordinate of each crossing inside one grid spacing;
the sign of \(2v_n/r\) at the initial critical points inside the
normal-velocity noise; any global causal reading; any time after \(T=0.05\);
any identification of these expansions with \(\chi\).

The local consequence for the spherical gradient boundary is the product
itself. \(\theta_+\theta_-=4r^{-2}g^{ab}\partial_a r\,\partial_b r\) is the
causal character of the areal gradient. Its zero is that boundary on this
chart. The saved samples move it from the initial critical points, marginal
within the normal-velocity noise, onto the edges of the two finite arcs.

## Local-boundary successor

`results/development/nsc-local-boundary-review-v1.json` stays immutable.
Later independent checks are not written back into that review. The
successor is `results/development/nsc-local-boundary-review-v2.json`. It
binds the same episode JSON and NPZ, null-expansion JSON, coupled-response
JSON and NPZ, stored observer string, and pre-correction module hash.

That pre-correction module hash is
`1aa47f0b3d680621d23ac9329e0d9e1e38070b7df42fc2ec48848ea367fd82d5`.
Those bytes are the current module with the boost-return edit inverted, so
the flag again compares `factor * theta` and `theta / factor`. The current
module hash is
`5c57b61c30008f994e1e889b93584488a5ea1675580ab6dfa68281bfc833e89a`.
The flag on the current module reads the arrays returned by
`boosted_expansions`. Saved samples and control arrays are unchanged. The
replay is the in-memory measurement. The writer is not used, because it
would replace the v1 JSON.

The independent test named by the v1 review is
`68c2fa1559d705c3aae2d3b64d693af72311802e51cb5914432e1422063b4ea9`.
The current independent test is
`d2d6afb8471481a74b816223d58e33066f112abf1026d15eae74b291222e44fa`.
The later file adds
`test_control_policy_keeps_g_initial_and_proxies_off_the_gate` and
`test_v2_reads_stored_endpoint_without_touching_v1`. Both recovered files
are pinned in `.source-history` as content carriers. The original eighteen
cache objects stay in place. The carrier commit is not an earlier
laboratory commit. The v2 review binds the current bytes of this note, the
null-expansion module, the null-expansion tests, the coupled-response
module, and the coupled-response note.
