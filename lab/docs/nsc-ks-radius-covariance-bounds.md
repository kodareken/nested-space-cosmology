# Nonlinear radius coefficients and covariance-growth bounds

This record reuses the authenticated p=16 profile table and the unchanged
radius history. It bounds spatial derivatives and Fourier-l1 moments of the
full function `q=1/r_g-1/r_ref` through order9, including the omitted Fourier
tail. It does not truncate q to first order in the history amplitude.

The pointwise radius lower bound and the Fourier inverse bound are kept
distinct. The latter follows from the strictly positive Neumann margin
`r_ref,min - A0(delta r)` and the derivative product rule. All denominator
bounds are enclosed or directed downward; upper norm records are never used
as lower denominators. The physical profile already includes its amplitude.

The first two reciprocal Fourier moments are propagated into the potential
integrals for every group in the original channel inventory. The exact
reference evolution is Fourier diagonal and norm preserving. The resulting
positive comparison matrix controls covariance moments as

\[
\begin{pmatrix}M_0\\M_1\end{pmatrix}_{\Sigma}
\le e^{2K_0}\begin{pmatrix}1&0\\K_1&1\end{pmatrix}
\begin{pmatrix}M_0(D_{\rm up})+\int M_0(R)\\M_1(D_{\rm up})+\int M_1(R)\end{pmatrix}.
\]

The matrix is an error propagator, not a source state. Its zeroth input and
true-defect inputs remain unknown in this record; no zero is assigned to
them. The ell=0 groups have identity comparison matrices because their pure
radius response is zero. No field evolution or source preparation is run.

These bounds identify whether propagation itself amplifies an eventual UV
error budget. They do not yet evaluate that budget: the actual projector
coefficients, the initial auxiliary covariance error and the time-integrated
true Weyl defect still need bounds. The physical local incoming gate remains
OPEN, with the source, action, scales and seam unchanged.
