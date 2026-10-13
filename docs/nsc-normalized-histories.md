# Static-to-causal response and retained history controls

These are two distinct finite comparisons. They retain their original
geometry, state and numerical domains; neither supplies an absolute
renormalized stress or a coupled metric solution.

## Euclidean-frequency and vacuum Kubo match

On the existing smooth axial circle, $R=2,a=1,N=q=1$, the radius channel
$r=r_0/(1+Js)$ gives $H[J]=H_0+JV$ and $V_D=i\beta V$.
The warped-sphere spin connection cancels under the canonical rescaling by
$r$, including its time derivative in this fixed-lapse/radial-metric channel.
For the negative-energy vacuum, Fourier convention $e^{izt}$ and interaction
$+JV$, the retarded response is

$$
\chi_R(z)=\sum_{i-,j+}|V_{ij}|^2
\left[(z-\Delta_{ij})^{-1}-(z+\Delta_{ij})^{-1}\right],
\qquad \Delta_{ij}=E_j-E_i>0,\quad \operatorname{Im}z>0.
$$

Independent full internal-frequency integration of
$\operatorname{tr}(D_{\omega-\nu/2}^{-1}V_D
D_{\omega+\nu/2}^{-1}V_D)$ with measure $d\omega/(2\pi)$ matches
$\chi_R(i\nu)=-2\sum\Delta|V_{ij}|^2/(\Delta^2+\nu^2)$ on finite grids.
The proper-time profile $g(d)=E_1(d^2/\Lambda^2)/2$ has
$g'(d)=-e^{-d^2/\Lambda^2}/d$; its divided-difference Hessian uses $g''$
at coincidence. Its independently integrated zero-frequency functional is

$$
E_\Lambda=\operatorname{Tr}F_\Lambda(H),\qquad
F_\Lambda(E)=\tfrac12[\Lambda e^{-E^2/\Lambda^2}/\sqrt\pi
-|E|\operatorname{erfc}(|E|/\Lambda)].
$$

Periodic-minus-antiperiodic comparison cancels state-independent local
counterterms in the continuum relative response. It does not cancel the
state-dependent operator contact: for $r=r_0e^{Js}$,
$d^2E/dJ^2=\chi_R(0)+\langle H''\rangle$, or
$\operatorname{Tr}[F'_\Lambda(H)H'']$ at finite cutoff.
The representative $\kappa=1$ relative log-radius contact is
`-0.02201411785` before multiplicity $4\kappa$.

Through $\kappa=8$ at proper-time cutoff 4, the relative Euclidean kernels
for $\nu=0,0.6,1.4$ are `0.18577154999`, `0.17252544942`,
`0.13062541226`. The retarded Abel/Richardson differences are
`-1.49e-9`, `-1.41e-9`, `-1.14e-9`. The Abel weight
$e^{-(E_i^2+E_j^2)/K^2}$ is a summation device distinct from the finite
proper-time action. Spatial, angular and frequency refinements are observed
convergence, without rigorous omitted-tower bounds or identification of the
two finite regulators. The result is a matter radius-response control,
not a graviton Hessian or a proof of finite-cutoff Lorentzian causality.

## Targeted retarded history resolution

The separate `14_1` response diagnostic retains the original four signed
energies, coherent three-column source, fixed incoming frame and compact
radius direction. One new 3201-node reference grid exactly nests the
archived 801-node grid; the resampling tolerance is `3e-11`.
It then compares `(points,steps)=(1601,128),(3201,128),(3201,256)`
on PG time `[0,0.3]`, with a 120 CPU-second cap. A failed spatial indicator
stops before the temporal solve. Relative spatial/time limits are
`0.25`, `0.1`; flux, Hermiticity and auxiliary algebra use `3e-11`.

Fields, traces and covariance tangents are compared on nested nodes. The
source covariance, initial/inflow tangents and frame-map derivative remain
unchanged; source phases are not reapplied. The estimated causal trace
interval `(0.05404795482456391,0.24595204517543606)` is a diagnostic,
not a directed bound or acceptance gate. Stationary drift and sampled
projector tangents are reported without subtraction. Replay contracts
saved arrays without new propagation. The declared continuum error bound
remains `null`, continuous $C_0$ matching remains OPEN, and no complete CAR,
physical constraint solution or metric evolution follows from a sampled PASS.
