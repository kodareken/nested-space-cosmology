# Sharp weighted Fourier sums of the reciprocal-radius perturbation

Root assembles integrated Weyl remainders and global-projector bounds. This
owner supplies the missing **omega-only** weighted sums of the true
reciprocal perturbation on the saved pure-radius control, so that
\(k/\omega\) integrals need not be replaced by a single global high moment
\(A_n(q)\).

It reuses the committed \(p=16\) tables \(W^1,\ldots,W^4\), the high-radius
loader, alias/tail identities and directed denominator checks. It does not
evolve a field or source, fit a profile, or change the action. Canonical
\(K\) is a numerical momentum split, not a source-energy cutoff. The
physical local incoming gate remains OPEN.

## Reciprocal polynomial and its true remainder

On the existing control, \(U\equiv 0\), \(\lvert\chi s\rvert\le\sigma\), and
physical \(W\) already contains \(\alpha=0.001\). That factor is not applied
again. With constant positive \(r_{\mathrm{ref}}\),

\[
q=\frac1{r_g}-\frac1{r_{\mathrm{ref}}},\qquad
q_4=\sum_{p=1}^{4}(-\chi s)^p\frac{W^p}{r_{\mathrm{ref}}^{p+1}}.
\]

The homogeneous reference potential is only the \(\ell=0\) mode and is
absent from \(q_4\). Positive-order transfers skip \(\ell=0\). The mode
envelope of the polynomial is the triangle inequality on the saved balls,

\[
\lvert\widehat{(q_4)}_\ell\rvert
\le\sum_{p=1}^{4}
\frac{\sigma^p}{r_{\mathrm{ref,min}}^{p+1}}
\bigl\lvert\widehat{(W^p)}_\ell\bigr\rvert.
\]

The unknown remainder \(q-q_4\) is **not** a difference of two Wiener
uppers and is not a truncated inverse used as physics. In the Wiener
algebra, \(\eta_j=A_j(\delta r)/r_{\mathrm{ref,min}}\) with
\(A_j(\delta r)=\sigma A_j(W)\) and \(\eta_0<1\). The coefficient-positive
majorant \(f(x)=x^5/(1-x)\) is composed with the exponential generating
function of \(\eta\) by `arb_series` through at least order 4:

\[
A_j(q-q_4)
\le j!\,[t^j]\,\frac{f(\eta(t))}{r_{\mathrm{ref,min}}}.
\]

The denominator \(1-\eta_0\) is a directed lower. Pointwise \(r_{\min}\) is
never used as a Wiener inverse.

Omitted \(W^p\) Fourier bands use their **own** certified \(L^1\) tails
from `alias_and_tail_bounds` at \(p=16\), recomputed from the stored
derivative \(L^1\) dyadics, not JSON binary64 displays.

## Omega-only kernels

Write \(w=\lvert\omega\rvert\). None of the following includes
\(\mathrm{d}k/(2\pi)\), a sum over \(\mathrm{sign}(k)\), or spin/source
multiplicity. Callers apply those factors once. Both \(\pm\omega\) modes
already sit in the stored coefficient inventory; there is no extra factor
two.

### Small transfer, \(\lvert k\rvert\ge\max(K,w)\)

\[
f_{\mathrm{small}}(n,s)
=\frac{(w/2)^n}{(5-s)\max(K,w)^{5-s}}
=\int_{\max(K,w)}^\infty k^{s-6}\,(w/2)^n\,\mathrm{d}k,
\qquad n=1,\ldots,5,\quad s\in\{0,1\}.
\]

The value is zero at \(w=0\). The two branches agree at \(w=K\). Known
coefficients are summed against this kernel; they are not replaced by a
global \(A_n\). Polynomial envelopes:

- \(n\le 5-s\): constant \(K^{n-(5-s)}/(2^n(5-s))\);
- \(n>5-s\) (only \(n=5,s=1\)): \(w/(2^n(5-s))\).

### Large-transfer shifted symbol

\[
f_{\mathrm{shift}}(s)
=1_{w>K}\frac{w^{s+1}-K^{s+1}}{s+1},\qquad s\in\{0,1\}.
\]

Envelope \(w^{s+1}/(s+1)\). The \(n=0,p=0\) term for \(P_\infty\) / global
\(P_0\) uses this kernel separately.

### Large-transfer Taylor terms at unshifted \(k\)

\[
f_{\mathrm{taylor}}(n,p,s)
=1_{w>K}\,(w/2)^n\int_K^w k^{s-p}\,\mathrm{d}k,
\qquad n=0,\ldots,4,\quad p\ge 1,\quad s\in\{0,1\}.
\]

Owned projector calls use \(p=j+1+n\) for \(j=0,\ldots,4\). Their
polynomial envelopes have degree at most 4:

- \(p>s+1\): \(w^n K^{s-p+1}/(2^n(p-s-1))\);
- \(p=s+1\): \(\log(w/K)\le w/K\), envelope \(w^{n+1}/(2^n K)\);
- \(p<s+1\): \(w^{n+s-p+1}/(2^n(s-p+1))\).

Unknown omitted-profile tails and the reciprocal remainder are contracted
with these nonnegative envelopes. Missing moments are rejected, never
filled with zero.

## Directed denominators

Period length, \(r_{\mathrm{ref}}\) and canonical \(K\) are positive
enclosing intervals. Collapsed upper-only dyadics are rejected as
denominators. Low-bit checks on \(K=7/6\) show that inverting
\(K^{\mathrm{upper}}\) fails to cover the exact reciprocal power, while
the enclosing reciprocal covers it.

## Authentic control, \(K=256\) and \(K=1024\)

Payload
`fc20111d0c72f1ac81d56ee5ea64cdd1bb94a9bfc6e86edfa19d1f8992f0ce33`,
profile identity
`a4c061f297027bcb9934fa2338bb81fb799ca4ec26a025f6bc563384c70cd69d`.
Period \(205/512\), \(\sigma=0.03\), \(r_{\mathrm{ref,min}}=7/5\),
\(\eta_0\le 1.742489018860512\times 10^{-5}<1\). Fourier settings remain
\(p=16\), \(M=65536\), retained index \(8192\). No field or source run.

Uppers below are binary64 displays of exact dyadics. In every listed
kernel the **retained polynomial** dominates; omitted-profile tails and
the reciprocal remainder are many orders smaller. \(K=1024\) is strictly
tighter than \(K=256\) for these nonnegative kernels.

| Kernel | \(K=256\) retained | omitted tail | remainder | total |
|---|---:|---:|---:|---:|
| \(f_{\mathrm{small}}(1,0)\) | \(8.895373565736502\times 10^{-17}\) | \(7.90\times 10^{-36}\) | \(2.67\times 10^{-35}\) | \(8.895373565736502\times 10^{-17}\) |
| \(f_{\mathrm{small}}(5,1)\) | \(3.144648174349598\times 10^{-6}\) | \(3.65\times 10^{-22}\) | \(4.24\times 10^{-24}\) | \(3.144648174349599\times 10^{-6}\) |
| \(f_{\mathrm{shift}}(0)\) | \(9.614077099627725\times 10^{-5}\) | \(4.67\times 10^{-20}\) | \(5.42\times 10^{-22}\) | \(9.614077099627729\times 10^{-5}\) |
| \(f_{\mathrm{shift}}(1)\) | \(4.260291856051376\times 10^{-2}\) | \(3.23\times 10^{-15}\) | \(1.57\times 10^{-19}\) | \(4.260291856051699\times 10^{-2}\) |
| \(f_{\mathrm{taylor}}(0,1,0)\) | \(2.500062951120735\times 10^{-7}\) | \(1.82\times 10^{-22}\) | \(2.12\times 10^{-24}\) | \(2.500062951120737\times 10^{-7}\) |
| \(f_{\mathrm{taylor}}(4,5,1)\) | \(8.255316033087911\times 10^{-5}\) | \(1.57\times 10^{-13}\) | \(2.58\times 10^{-22}\) | \(8.255316048771122\times 10^{-5}\) |

| Kernel | \(K=1024\) retained | omitted tail | remainder | total |
|---|---:|---:|---:|---:|
| \(f_{\mathrm{small}}(1,0)\) | \(1.039330040367922\times 10^{-19}\) | \(3.08\times 10^{-38}\) | \(1.04\times 10^{-37}\) | \(1.039330040367922\times 10^{-19}\) |
| \(f_{\mathrm{small}}(5,1)\) | \(2.568652359828163\times 10^{-7}\) | \(3.65\times 10^{-22}\) | \(4.24\times 10^{-24}\) | \(2.568652359828166\times 10^{-7}\) |
| \(f_{\mathrm{shift}}(0)\) | \(2.303806476183432\times 10^{-6}\) | \(4.67\times 10^{-20}\) | \(5.42\times 10^{-22}\) | \(2.303806476183479\times 10^{-6}\) |
| \(f_{\mathrm{shift}}(1)\) | \(3.044526565441626\times 10^{-3}\) | \(3.23\times 10^{-15}\) | \(1.57\times 10^{-19}\) | \(3.044526565444859\times 10^{-3}\) |
| \(f_{\mathrm{taylor}}(0,1,0)\) | \(1.824450980529068\times 10^{-9}\) | \(4.56\times 10^{-23}\) | \(5.30\times 10^{-25}\) | \(1.824450980529114\times 10^{-9}\) |
| \(f_{\mathrm{taylor}}(4,5,1)\) | \(4.662819314895330\times 10^{-7}\) | \(2.45\times 10^{-15}\) | \(4.04\times 10^{-24}\) | \(4.662819339400346\times 10^{-7}\) |

The largest listed total at \(K=256\) is \(f_{\mathrm{shift}}(1)\), still
dominated by retained modes. High-degree envelopes enlarge the omitted
\(W^p\) tail relative to a constant envelope, but it remains far below the
retained sum. \(f_{\mathrm{taylor}}(0,1,1)\) equals \(f_{\mathrm{shift}}(0)\)
identically.

## What remains OPEN

These sums bound existing radius-control geometry. They are not a UV
certificate and not a physical incoming-gate PASS. Still missing for an
integrated Weyl remainder:

- time integrals of the weighted kernels;
- fourth-order symbol and global-projector coefficients assembled by root;
- the exact Moyal remainder;
- the Bloch trace inequality;
- source-energy quadrature tails (not identified with canonical \(K\)).

Fourteen focused tests pass, including analytic single/cosine modes versus
explicit sums, Neumann remainder plus finite convolution, denominator and
threshold continuity, polynomial envelopes, and the authentic table at
\(K=256\) and \(K=1024\).
