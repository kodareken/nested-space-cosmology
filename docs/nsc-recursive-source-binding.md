# The scale-inherited LLL state does not source the black-universe neck

The scale calculation fixes a candidate room ratio, but a scale is not yet a
metric source.  This binding applies the simplest state implied by that ratio
to the actual stored neck, using only the existing massless magnetic
lowest-Landau-level (LLL) stress and PG projection.

For transparent LLL propagation, the inherited Hamiltonian scaling gives

$$
\kappa_c=\Omega\kappa_p,
\qquad
t_u=\frac{|q|\kappa_p^2}{48\pi},
\qquad
t_v=\frac{|q|\kappa_c^2}{48\pi}=\Omega^2t_u.
$$

No incoming thermal reservoir is chosen in this candidate: the deeper room
supplies the incoming state in the common parent energy frame.  The global
clock-transfer geometry has not independently proved this scaling.  Under the
stated inherited-Hamiltonian law, the conserved parent Killing power is

$$
\boxed{P_{\rm parent}=t_u-t_v
=\frac{|q|\kappa_p^2}{48\pi}(1-\Omega^2).}
$$

On the recorded `q=4` scale branch this gives

$$
P_{\rm parent}=-0.02227623165818482.
$$

The sign describes a directed incoming recursive state.  It is not converted
into a cosmological density rate.

## Bind it to the actual neck

Insert the stored black-universe neck data

$$
A=1-\frac{3\pi}{2},\quad A'=6,\quad A''=-3\pi,
\quad r=1,\quad\beta^2=\frac{3\pi}{2}
$$

into the already imported conformal stress and PG projection.  The result is

$$
T^{(4)}_{++}=0.007886485863354125,
\qquad
T^{(4)}_{--}=0.05267065468742516.
$$

Both are positive.  The stored black-universe neck requires both radial null
components to be negative.  This is a sign obstruction, not a resolution or
quadrature uncertainty.

The conclusion is narrow and decisive: the MMP transparent LLL sector cannot
be substituted for the full charged CTP source of the black-universe geometry.
MMP self-sources its own horizonless charged throat; the NSC child continuation
needs the full charged angular and compact covariance, its nonlocal CTP metric
variation, and the actual recursive transmitting domain.  The scale root is
left intact and is not searched again.

The focused reproducer is

```sh
python3 scripts/derive_nsc_recursive_source_binding.py --check
```

It performs no numerical integration and reruns no prior generator.
