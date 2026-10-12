# A finite observer window in the retained parent exterior

## Admission and question

Borrow Schwarzschild radial null propagation and the tortoise coordinate from
[Kuchar, gr-qc/9403003v1](https://arxiv.org/html/gr-qc/9403003v1), Section 2.1,
equations 11–13. Borrow smooth vacuum local geometric uniqueness from
[Sbierski, 1309.7591v2](https://arxiv.org/html/1309.7591v2), Section 2,
Theorems 2.6–2.8. These are imported results, not a new propagation law.

The missing connection is a finite, physically interpreted observer domain
for the [explicit initial slice](nsc-atlas-connection-slice.md). This calculation
decides how long a static observer is guaranteed the same vacuum geometry in
assemblies that preserve its initial data, even if their far exteriors differ.
It unlocks an atlas of locally indistinguishable exterior assemblies.

Fix before calculation: M/ell=0.4, R_obs=4M, R_cut=10M and the exact right
Schwarzschild initial region 2M<R<R_cut, with K=0. This is a chosen retained
region, not a separation bound supplied by a gluing theorem. The gluing must
preserve its closure and a larger open overlap. The observer calculation uses
at most 64 Gauss nodes and a fixed analytic formula; budget 60 CPU seconds.
Stop if the slice lacks that vacuum region, if either radius leaves the
declared domain, or if analytic propagation and independent quadrature disagree.
No larger evolution or fit is admitted here.

## Imported propagation and local consequence

With f=1-2M/R and r*=R+2M log(R/(2M)-1), radial null curves satisfy
t +/- r*=constant. A static observer at R_obs remains in the future domain
of dependence of the retained open initial region for

$$
0\leq t<r_*(R_{cut})-r_*(R_{obs}),\qquad
0\leq\tau_{obs}<\sqrt{f(R_{obs})}
[r_*(R_{cut})-r_*(R_{obs})].
$$

The inward past-null endpoint stays above 2M for every finite t because the
inner optical endpoint is minus infinity. The outward past-null endpoint
reaches R_cut at the displayed limit. Angular motion cannot arrive sooner
than radial null motion. At our fixed radii, the bounds are
t/M=6+2 log(4) and tau_obs/M=(6+2 log(4))/sqrt(2).

Within this open window, different admissible far-parent assemblies that
retain identical initial h and K have isometric local vacuum developments.
This is a finite geometric degeneracy. It does not equate quantum states or
test fields whose initial data have not also been equated. It does not say
a disturbance actually arrives at the endpoint. A clock continues to exist
afterwards, but this certificate no longer determines its environment.
Vacuum local uniqueness is used only on this smooth domain; no global
existence theorem for the distributional shell is inferred.

## Evidence and boundary

The adapter [nsc_atlas_connection_observer.py](../src/recursive_horizons/nsc_atlas_connection_observer.py)
compares the analytic time with direct quadrature of dR/f and checks the
ingoing Kruskal null coordinate independently. The
[tests](../tests/test_nsc_atlas_connection_observer.py) exercise unit rescaling,
invalid radii, certificate exits and incompatibility with the supplied slice.
The successor bundle records exact inputs and residuals in causal-window.json.
No wall-to-parent signal, sibling communication, cross-universe age, dark
component or formation history follows from this calculation. The next
observable must be chosen outside this degeneracy or use separately supplied
field/state information.

Reading boundary: Kuchar Section 2.1 and Sbierski's definitions and theorem
statements were read for this use; neither full proof is re-derived here.
The bounded prior-art comparison identifies an application of standard null
propagation and local uniqueness, with no priority claim.
