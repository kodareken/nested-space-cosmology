# Spherical feedback action, first-order chart

This is the checked bulk and boundary primitive for the spherical
Einstein–Maxwell–Weyl chart with lapse and shift still free. It is not a
Cauchy evolution and it is not a regeneration result. \(C_W=0\) is rejected.
There is no pure-Einstein chart. The actual matter-source branch is not
coupled here.

The geometry is the owned \((+,-,-,-)\) reduction

\[
ds^2=N^2dt^2-q^2(dx+\beta dt)^2-r^2d\Omega^2,
\qquad L=N/r,\quad Q=q/r,
\]

\[
h=L^2dt^2-Q^2(dx+\beta dt)^2,\qquad \sqrt{|h|}=LQ.
\]

The accepted two-dimensional identity is

\[
D=\frac{Q_t-(\beta Q)_x}{L},\qquad
\sqrt{|h|}\,R_h=-2D_t+2(\beta D+L_x/Q)_x,
\]

\[
C^2=\frac{(R_h-2)^2}{3r^4}.
\]

No four-dimensional Riemann tensor is rebuilt in this module.

## Bulk coefficients

\[
\alpha=-\frac{4\pi C_W}{3},\qquad
F=-4\pi A r^2+2\alpha\chi,\qquad
Z=-24\pi A,
\]

\[
V=8\pi A r^2-\alpha(4\chi+\chi^2)-2\pi C_F\mathrm{flux}^2.
\]

\(F_r=-8\pi A r\) and \(F_\chi=2\alpha\neq 0\). The value \(\chi=R_h-2\) returns
\(\alpha(R_h-2)^2\) and is the only critical point. The magnetic term is the
owned contraction \(F_{\theta\phi}=q_{\mathrm{mag}}\sin\theta/2\).

The first-order density, after removing
\(\partial_t(-2FD)+\partial_x\big(2F(\beta D+L_x/Q)\big)\), is

\[
\begin{aligned}
\mathcal L=&\,2D(F_r u+F_\chi v)-2F_x L_x/Q\\
&+Z\big[Q/L\,u^2-L/Q\,r_x^2\big]+LQ\,V,
\end{aligned}
\]

with \(u=r_t-\beta r_x\) and \(v=\chi_t-\beta\chi_x\). On a homogeneous slice
the Einstein sector is the owned GHY density

\[
8\pi A\Big(Nq-2r\dot r\dot q/N-q\dot r^2/N\Big).
\]

Flat space and de Sitter give \(R_h=2\) and \(C^2=0\). Schwarzschild gives
\(C^2=48M^2/r^6\).

## Momenta and constraints

\[
\begin{aligned}
p_Q&=2(F_r u+F_\chi v)/L,\\
p_r&=2F_r D+2ZQ u/L,\\
p_\chi&=2F_\chi D,\\
\Pi&=p_r-(F_r/F_\chi)p_\chi.
\end{aligned}
\]

\[
\begin{aligned}
\mathcal C&=\frac{p_Q p_\chi}{2F_\chi}+\frac{\Pi^2}{4ZQ}+\frac{Z r_x^2}{Q}-QV-2\partial_x(F_x/Q),\\
\mathcal D&=p_r r_x+p_\chi\chi_x-Q\partial_x p_Q.
\end{aligned}
\]

For generic jets, including nonconstant \(\beta\),

\[
p_Q Q_t+p_r r_t+p_\chi\chi_t-\mathcal L-(L\mathcal C+\beta\mathcal D)
=\partial_x(\beta Q p_Q+2LF_x/Q).
\]

The same density gives \(\mathrm{EL}_L=-\mathcal C\) and
\(\mathrm{EL}_\beta=-\mathcal D\).

## Explicit Hamilton equations

The generator is \(H=\int(L\mathcal C+\beta\mathcal D)\,dx\). The velocity
equations are

\[
\begin{aligned}
Q_t&=L\frac{p_\chi}{2F_\chi}+\partial_x(\beta Q),\\
r_t&=L\frac{\Pi}{2ZQ}+\beta r_x,\\
\chi_t&=L\frac{p_Q}{2F_\chi}-\frac{F_r}{F_\chi}\frac{L\Pi}{2ZQ}+\beta\chi_x.
\end{aligned}
\]

The momentum equation for \(\chi\) is

\[
\dot p_\chi=LQ\,\partial_\chi V+2F_\chi\partial_x(L_x/Q)+\partial_x(\beta p_\chi).
\]

In the coordinate chart \(F_\chi=-8\pi C_W/3\), the factored right-hand sides
that match these derivatives are

\[
\begin{aligned}
Q_t&=\frac{16\pi C_W\partial_x(\beta Q)-3L p_\chi}{16\pi C_W},\\[4pt]
r_t&=\frac{48\pi A C_W Q\beta r_x+3A L p_\chi r-C_W L p_r}{48\pi A C_W Q},\\[4pt]
\chi_t&=-\frac{3A L p_\chi r^2-16\pi C_W^2 Q\beta\chi_x+3C_W L Q p_Q-C_W L p_r r}{16\pi C_W^2 Q}.
\end{aligned}
\]

\[
\begin{aligned}
\dot p_Q=\frac{1}{96\pi A C_W^2 Q^2}\Big(&
768\pi^2 A^2 C_W^2 L Q^2 r^2-2304\pi^2 A^2 C_W^2 L r_x^2\\
&-1536\pi^2 A^2 C_W^2 L_x r r_x-9A^2 L p_\chi^2 r^2\\
&-192\pi^2 A C_F C_W^2 L\,\mathrm{flux}^2 Q^2
+128\pi^2 A C_W^3 L Q^2\chi^2\\
&+512\pi^2 A C_W^3 L Q^2\chi-512\pi^2 A C_W^3 L_x\chi_x\\
&+96\pi A C_W^2 Q^2\beta\partial_x p_Q+6A C_W L p_\chi p_r r-C_W^2 L p_r^2\Big),
\end{aligned}
\]

\[
\begin{aligned}
\dot p_r=\frac{1}{16\pi C_W^2 Q^2}\Big(&
256\pi^2 A C_W^2 L Q^3 r-768\pi^2 A C_W^2 L Q r_{xx}\\
&+768\pi^2 A C_W^2 L Q_x r_x-768\pi^2 A C_W^2 L_x Q r_x\\
&+256\pi^2 A C_W^2 L_x Q_x r-256\pi^2 A C_W^2 L_{xx} Q r\\
&+3A L Q p_\chi^2 r+16\pi C_W^2 Q^2\partial_x(\beta p_r)-C_W L Q p_\chi p_r\Big),
\end{aligned}
\]

\[
\begin{aligned}
\dot p_\chi=\frac{1}{3Q^2}\Big(&
8\pi C_W L Q^3\chi+16\pi C_W L Q^3+16\pi C_W L_x Q_x\\
&-16\pi C_W L_{xx} Q+3Q^2\partial_x(\beta p_\chi)\Big).
\end{aligned}
\]

## Boundary ledger

The spatial coordinate is periodic. The time interval has two finite caps.
Dirichlet data on those caps are \(Q\), \(r\) and \(\chi\), so the cap form
\(p_Q\delta Q+p_r\delta r+p_\chi\delta\chi\) vanishes. A periodic spatial flux
integrates to zero, so this domain has no spatial corners. Corner terms are
not an additional requirement.

The cap counterterms cancel the bulk divergences:

\[
-4\pi C_E[j^t]+4\pi C_E[j^t]=0,
\qquad
-C_\Box[J^t_R]+C_\Box[J^t_R]=0,
\qquad
-2FD+2FD=0.
\]

The Myers angular charge of the Euler current is \(Q_+=-4\pi j^t\). The
confirmed Euler current, with the full \(\beta\) residual equal to zero, is

\[
j^t=8k A_0-16 v s_x,\qquad
j^x=-8(P+\beta k)A_0+16 v s_t,
\]

where \(k=(q_t-(\beta q)_x)/N\), \(v=(r_t-\beta r_x)/N\), \(s=r_x/q\),
\(P=N_x/q\) and \(A_0=1+v^2-s^2\). The naive current without the mixed
\(16 v s_x\) and \(16 v s_t\) terms is rejected. The section identity used for
that residual is

\[
Nq\,r^2 E_4=Nq\big[-4 R_2 A_0+16(C^2-AB)\big],
\]

with \(A=n(v)-(P/N)s\), \(n(v)=(v_t-\beta v_x)/N\),
\(B=s_x/q-(k/q)v\) and \(C=v_x/q-(k/q)s\).

Einstein is already the GHY completion. Splitting \(F=F_{\mathrm{EH}}+2\alpha\chi\)
with \(F_{\mathrm{EH}}=-4\pi A r^2\), the raw conformal time current is

\[
J^E_t=24\pi A\frac{q r}{N}(r_t-\beta r_x).
\]

The combination \(-J^E_t+2 F_{\mathrm{EH}} D\) equals the homogeneous GHY
primitive \(-8\pi A(r^2\dot q/N+2 q r\dot r/N)\). The only time flux beyond
that GHY completion is the Weyl piece \(4\alpha\chi D\):

\[
2 F_{\mathrm{EH}} D+4\alpha\chi D=2FD.
\]

The Box-\(R\) flux, already containing the angular \(4\pi\), is

\[
J^t_R=4\pi\frac{Q r^2}{L}(R_t-\beta R_x),\qquad
J^x_R=4\pi r^2\Big(-\frac{\beta Q}{L}(R_t-\beta R_x)-\frac{L}{Q}R_x\Big).
\]

Its divergence matches \(\partial_a(\sqrt{-g}\nabla^a R)\) reduced with
\(h^{-1}\). The homogeneous cap is the owned \(4\pi q r^2\dot R/N\).

## Conditional initial slices

These slices solve the constraints. They are not a physical source covariance,
not future imposed pulses, and not the separate source-coupling branch.

Take \(r\), \(Q\), \(\chi\) and \(p_\chi\) constant, \(p_{Q,x}=j/Q\), and

\[
K=\frac{p_Q p_\chi}{2F_\chi}-QV+\rho,\qquad
\Pi=\pm\sqrt{-4ZQK},\qquad
p_r=\frac{F_r}{F_\chi}p_\chi+\Pi.
\]

If \(K>0\) and the periodic mean of \(j\) vanishes, then
\(\mathcal C+\rho=0\) and \(\mathcal D+j=0\).

At the normalization \(A=1/(4\pi)\), \(C_W=3/(4\pi)\), \(C_F=1/\pi\),
\(\mathrm{flux}=1\), \(r=Q=1\), \(\chi=0\) and \(p_\chi=-2\), with
\(\rho=1+\cos x/4\):

- \(j=0\) gives \(p_Q=2\), \(K=2+\cos x/4\) and
  \(\Pi=\pm\sqrt{48+6\cos x}\). The minimum of \(K\) is \(7/4\).
- \(j=\sin x/5\) gives \(p_Q=2-\cos x/5\) and \(K=2+3\cos x/20\).
  The minimum of \(K\) is \(37/20\).

Both signs of \(\Pi\) leave the two constraint residuals at zero.

## Matter

On the canonical half-density \(u=r\sqrt{q}\,\psi\), with that covariance held
fixed, the continuum operator is

\[
H=\sigma_2\tfrac12\{L/Q,P\}+\sigma_1 L\kappa-\tfrac12\{\beta,P\},
\qquad P=-i\partial_x.
\]

A radius variation at fixed \(L,Q\) cancels in this continuum chart. The
existing discrete \(\sqrt{N}\) ordering is not fixed and is not used. The
factor \(4\kappa\) is only the isotropic copy weight of one angular block.
\(\mathrm{Tr}(CH)\) is not added to \(\mathcal C\) or \(\mathcal D\).

## What this does not claim

The six geometric equations have conditional constraint-solving data, not an
evolution. Discrete constraint preservation is not owned. The source-coupling
branch remains a separate owner.

```sh
python scripts/lab.py -m pytest tests/test_nsc_spherical_feedback_action.py -q
```
