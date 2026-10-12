# First Newton proposal from the complete retained source

This proposal uses the complete finite-source residual and all16 retarded
history derivatives of the coupled n=8 history. Its phase values use the
192-node rule verified by the directed phase-value certificate. The
32-node edge and tangent are replaced by the corresponding192-node
quantities; source fields and upstream preparation are unchanged.

The existing Newton owner supplies the matrix assembly and geometric
preconditioning. Both rectangular Gauss-Newton and square tau layouts are
recorded for diagnosis. The rectangular layout retains all N and beta solve
nodes and is selected for the first nonlinear trial. The same declared
`w,w',w''` controls at the center are retained.

The output is a proposed coefficient change and a **linear prediction**.
It is not an accepted Newton step or a measured residual at the new g.
Fresh evolution from the same upstream source must evaluate the trial
before it can be accepted. Source, field, finite-tail and between-node
errors remain OPEN, as does the physical local incoming gate.

```sh
python scripts/derive_nsc_ks_coupled_newton_proposal.py --check
```
