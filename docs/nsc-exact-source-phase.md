# Exact source phases in the incoming field accuracy control

This bounded diagnostic changes only the numerical treatment of the owned
harmonic incident source. The physical source, grid, Dirac generator,
canonical frame and boundary data are unchanged.

The generic prepared adapter freezes its incoming forcing callback at each
time midpoint. For the fixed source law, the forcing phases are already
known exactly. Reuse the augmented phase block from the existing
`evolve_field_jets` owner:

\[
\frac d{d\tau}\begin{pmatrix}\Phi\\Q\end{pmatrix}
=\begin{pmatrix}L_0&B\\0&-i\,\mathrm{diag}(E)\end{pmatrix}
\begin{pmatrix}\Phi\\Q\end{pmatrix},\qquad Q(0)=I.
\]

B is the same already SAT-weighted right-boundary forcing amplitude.
There is no interior forcing correction. The auxiliary phase block is a
numerical representation of the original source phases, not an added
physical degree of freedom or a closed few-frequency bulk surrogate.
The spatial field still has all1602 spin/grid degrees of freedom.

One801-grid reference control samples65 times on the inherited numerical
interval[0,.3], with a30CPU-second ceiling. Since this control's generator
is time independent, one augmented exponential supplies every sample.
The incoming F_z is evaluated from the actual semi-discrete PDE at rho1;
at fixed rho1, partial_tau equals partial_z. The previous cubic derivative
is also evaluated separately to distinguish its contribution.

Compare the resulting raw N,beta matter values with the same original
stationary reference columns. A factor100 reduction relative to the saved
midpoint-forcing drift is the declared numerical discrimination test,
not a physical stationarity tolerance. The full-source and continuum
error bounds remain OPEN, and no nonzero-history accuracy result follows
from a reference control alone.

The result and original arrays are persisted by
`scripts/derive_nsc_exact_source_phase.py --run`; `--check` only authenticates
and replays those arrays. No old field/source producer, new horizon grid,
physical duration, constraint root or metric evolution is introduced.
