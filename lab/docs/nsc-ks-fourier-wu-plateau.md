# Fourier-on-I times the owned plateau, same (w,U) class

Chebyshev `n<=32` is a refused discretization of the declared class

\[
\delta r=\chi(s)\bigl[s\,w(z)+s^3U(z)/6\bigr].
\]

This owner keeps those two functions, the locked `I`, axial support and
normal window, and replaces the polynomial basis by

\[
1,\quad \cos(2\pi k\xi),\quad \sin(2\pi k\xi),
\]

with \(\xi\) the affine map of \(I\) onto \([0,1)\), multiplied by the
already-owned axial plateau. On \(I\) the plateau is identically one. In
the collar it kills any computational-envelope period. `n` is 8 or 16
modes per function. `LocalIncomingFamily` stays the Chebyshev regression.

```sh
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_fourier_wu_plateau.py
```
