# Same-realization conformal local response

On the finalized conformal continuation, the maintained region-0 structure,
resolved surface exchange and a quantitative local response coexist on the
same saved \(T\in[0.2,0.3]\) segment. The original \(T=0\) two-mode
observer is fixed. Its occupation falls while the natural width-2 window
remains the leader. The reducer consumes the generated geometry and state;
it does not integrate a new coupled trajectory. This is a measured finite
result with numerical comparisons, not a certificate for an unknown
geometric path or continuum evolution. Renewal and stress are not claimed.

The consumed episode is `nsc-spherical-conformal-episode-v2`, primary
`nf512_dt_0.0005`. All four source runs reach \(T=0.3\), have a positive
chart and pass their named sampled constraint-forcing assessment. The
successor binder checks the original \(T=0\) Cauchy arrays against v5 and
all eight \(T=0.05\) handoff fields against the predecessor bytes. The
surface assessment is bound to the same source JSON and NPZ. The
[initial feasibility record](nsc-conformal-local-response.md) remains its
short probe and is not rewritten.

## Same source, geometry and observer

The six weights remain \((0.75,0.75,0.5,0.5,0.25,0.25)\). The initial
state is the actual transported \(\Phi(0.2)\), with initial
\(\|C_{AE}\|_F=0.338636575692358\). The observer is the original columns
0 and 1, copied without QR or a phase reset. Exterior means the Hilbert-space
complement of that pair. Both conditional controls and the full reference
use this same preparation and the generated \(Q\) schedule.

The conformal representative generator is
\(H=\sigma_2P+\kappa Q\sigma_1\), with \(L=Q\), \(\beta=0\). Stored
coarse \(Q\) is prolonged to the original 2048-point quadrature and
interpolated linearly between 21 stored frames. The isotropic copy factor
belongs to the episode's geometric source energy and force; it does not
multiply this representative generator or the occupation. Independent AP
Fourier action checks pass at the first and last frames below
\(4.4\times10^{-13}\).

One streamed solve transports twelve columns: the retained and exterior
parts of the same six source columns. Linearity gives the actual amplitude
as their sum. Its covariance contains the actual cross terms. Summing the
two component covariances removes only that cross, and evolving the
retained component alone gives the initial-drive omission. An unsplit
pilot independently checks this superposition. The source preparation is
unchanged by this computational decomposition.

## Local response and own control errors

The two autonomous occupations move from
\((0.10557447719690986,0.10557463752522356)\) to
\((0.03214495758531893,0.03226899862128299)\). The maximum occupation
change is \(0.07342951961159094\); the coherence change is
\(0.0023742079275808343\).

| Comparison on the same segment | Movement or error | Fraction of its effect |
|---|---:|---:|
| Streamed occupation versus conditional full | \(6.3508141397161655\times10^{-6}\) | \(0.00865\%\) |
| Streamed coherence versus conditional full | \(8.478723274307841\times10^{-7}\) | \(0.03571\%\) |
| Conditional full occupation versus stored autonomous frames | \(2.2951084088784768\times10^{-10}\) | \(3.13\times10^{-7}\%\) |
| Memory omission's occupation movement | \(0.016818538289934015\) | — |
| Memory-off output-step halving indicator | \(1.650868287690621\times10^{-6}\) | \(0.00982\%\) of memory effect |
| Memory-off error versus independent triangular exponential reference | \(2.20112157194724\times10^{-6}\) | \(0.01309\%\) of memory effect |
| Initial exterior-drive omission's occupation movement | \(0.048541204982542864\) | — |
| Drive-off error versus independent full retained-initial propagation | \(1.66255239472321\times10^{-5}\) | \(0.03425\%\) of drive effect |
| Actual initial cross omission's occupation movement | \(0.09206143279175863\) | — |
| Cross occupation error versus independent full column split | \(1.2282200189728254\times10^{-5}\) | \(0.01334\%\) of cross effect |

The cross covariance's maximum Frobenius movement is
\(0.13011977207833153\); its independent full-split error is
\(1.709921623148559\times10^{-5}\), or \(0.01314\%\). These are errors
of the respective control, rather than a shared roundoff or headline-error
substitute. Memory's halved-grid comparison is a temporal indicator. The
history vanishes in the memory omission. No synthetic cross is substituted.

The independent memory reference is a separate
[record](../results/development/nsc-conformal-memory-control-v1.json) and
[payload](../results/development/nsc-conformal-memory-control-v1.npz),
preserving the production bytes. With \(P_A=VV^\dagger\),
\(P_E=I-P_A\), it evolves the same six initial columns under

$$
G_{\mathrm{off}}(t)=H(t)-P_EH(t)P_A.
$$

This triangular causal negative control removes the retained-to-exterior
response and keeps the exterior-to-retained drive. Its generator is
non-Hermitian; it is not a physical replacement model. Independent
midpoint exponential action with four substeps gives the memory-off error
above. Two versus four substeps moves its occupation by
\(3.153471234140781\times10^{-11}\), below one percent of the memory
effect. The supplement used \(3.638052\) CPU seconds under a 30-second
cap, stores \(4750\) bytes and no propagator history, and adds two passing
tests of its block ODE and actual correlated-state error. Its owners are
[module](../src/recursive_horizons/nsc_conformal_memory_control.py),
[driver](../scripts/derive_nsc_conformal_memory_control.py), and
[tests](../tests/test_nsc_conformal_memory_control.py). The source geometry
and full production comparison were not rerun.

Doubling the full midpoint substeps moves occupation by
\(1.9253487693049465\times10^{-11}\). A rate-Hermite schedule using the
saved actual lifted \(\dot Q\), which equals the saved \(\dot L\),
moves it by \(2.2907768737479017\times10^{-10}\). Both coherence
movements are also below one percent of the coherence change. The
rate-Hermite comparison is a declared-path indicator, not an enclosure of
unstored geometry.

For geometry refinement, the alternate \(n_f=256\) path is prolonged to
the primary 2048-point grid. The same primary state, observer and
1024-dimensional generator are used, changing only \(Q(t)\). Its
maximum \(Q\) gap is \(4.42\times10^{-13}\) and occupation movement is
\(1.25\times10^{-16}\). The corresponding timestep geometry gap is
\(8.40\times10^{-14}\), with occupation movement
\(7.63\times10^{-17}\). Separately, actual autonomous occupations differ
by \(1.33\times10^{-8}\) across spatial resolution and
\(2.32\times10^{-11}\) across timestep. All reported local numerical
movements remain below one percent of their own occupation or coherence
effect.

## Relationship to the maintained physical window

On this same coordinate-time interval, the time integrals of absolute
normal energy flux across the natural surfaces are:

| Surface | Absolute energy-flow integral | Largest space/time/frame indicator as fraction of that flow |
|---|---:|---:|
| \(x=0\) | \(0.038484363186333706\) | \(0.08449\%\) |
| \(x=2\) | \(0.06757847920532943\) | \(0.02699\%\) |
| \(x=4\) | \(0.060887369296944736\) | \(0.01429\%\) |

These selected integrals reproduce the episode's physical surface ledger,
including its `flux_nodal / dx` normalization. They are energy-flow
quantities; occupation and coherence are separate one-body measurements.
The selected interval carries about 96–98 percent of the corresponding
full-continuation flow. Region 0 remains the shell and probability leader,
with shares at least \(0.6509785261434402\) and
\(0.6558454531429172\). The continuous \(T=0\) normal clock at \(x=1\)
advances from \(0.24120933558087854\) to \(0.361611405711853\); it is
not reset at the local segment's start.

The quantitative connection is that maintained localization and resolved
regional exchange occur on the generated realization whose unresolved
mode complement actively changes the fixed observer. Removing memory,
the initial exterior amplitude or its actual cross correlations produces
a local movement larger than its measured method error. The controls are
conditional field decompositions on that generated geometry; they do not
claim a separately integrated geometry for each omission.

## Budget, data and replay

The two-step shared solve took \(1.897327\) CPU seconds. Before scaling,
the full shared solve was forecast at \(189.733\) seconds and the entire
control plan at \(270.559\) seconds. Production, including independent
controls and refinements, used \(84.442888\) seconds inside the initial
300-second ceiling. The 1800-second confirmation ceiling was not needed.

The shared history is \((21,1022,12)\), \(4120704\) bytes, and its
coherent original-source history has width six. Time-indexed exterior
propagator storage is zero; the avoided dense history would occupy about
351 MB. The 34094-byte JSON/NPZ payload stores local arrays and bindings,
with its own declared 64 MiB cap. It does not duplicate source Cauchy
frames or alter the source trajectory's payload budget.

Owners:

- [Module](../src/recursive_horizons/nsc_conformal_local_response_successor.py)
- [Driver](../scripts/derive_nsc_conformal_local_response_successor.py)
- [Tests](../tests/test_nsc_conformal_local_response_successor.py)
- [Record](../results/development/nsc-conformal-local-response-v2.json)
- [Arrays](../results/development/nsc-conformal-local-response-v2.npz)

```sh
python scripts/lab.py scripts/derive_nsc_conformal_local_response_successor.py --check
python scripts/lab.py -m pytest tests/test_nsc_conformal_local_response_successor.py tests/test_nsc_conformal_local_response.py -q
```

`--check` is read-only. The production entry refuses an existing output;
new scientific runs require separately named successors. Sixteen tests
passed, including independent shared-column algebra and mutations of
readiness, source reset, chart, observer and control data. The conditional
state-dependent Duhamel relation stays in the record: it maps a controlled
Hamiltonian defect to covariance error on the same observer. No unknown
geometric defect or continuum error bound is claimed by this measurement.
