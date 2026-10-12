# Mixed geometry enclosure over the whole axial cell

This helper combines the existing coordinate-time jets with certified
spatial derivative norms of the same analytic w and U profiles. It encloses
all mixed geometry derivatives at every axial point in one numerical cell,
for a supplied preparation-time box.

For the normalized coefficient with indices(t,j), the perturbation is

\[
[\chi(s)s]_t\,{w^{(j)}(z)\over j!}
 +[\chi(s)s^3/6]_t\,{U^{(j)}(z)\over j!}.
\]

The time coefficients already include division by t!. The profile derivative
norms, including their original amplitudes, supply symmetric real enclosures
for the spatial factors. The signed background radius and axial jets are
retained. This gives a conservative box for the derivatives of the same
history; it introduces no independent geometric degrees of freedom.

The analytic axial-profile identity must match the supplied norm receipt.
All derivative orders must be present, including U bounds. The numerical
cell must contain the compact profile support. Radius positivity is checked
on the resulting whole-cell enclosure. Normal-window boundaries still
require time subdivision when an interval cannot be classified.

This enclosure replaces a separate spatial-box traversal when bounding the
auxiliary projector. It is not a covariance, source or physical gate
certificate. Tests verify containment of the actual profile's point jets in
the plateau, transitions and exterior, and reject foreign profiles, missing
norms and a cell that cuts through the profile support.
