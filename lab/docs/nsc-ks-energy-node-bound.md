# Interpolation factor for the actual stored energy nodes

`actual_node_remainder_factor(nodes, interval)` bounds the nodal polynomial
of the actual binary64 nodes. It does not identify those rounded values with
ideal Chebyshev roots. Every supplied binary value enters directed ball
arithmetic exactly, at160 bits by default.

For `n` distinct nodes and a Banach-valued function with a bounded nth
energy derivative, the divided-difference identity and its simplex integral
give

\[
 \sup_E\|f(E)-p_{n-1}(E)\|
 \le {\sup_E|\prod_j(E-E_j)|\over n!}\sup_E\|\partial_E^n f(E)\|.
\]

This reuses the Genocchi–Hermite formula; the norm of its Bochner integral
is bounded by the integral of the norm. A common scalar mean-value point is
unnecessary. See [de Boor, *Divided Differences*, equation52 and the following
interpolation-error identity](https://pages.cs.wisc.edu/~deboor/sat/papers/2/2.pdf).

Scale the declared interval to `E=c+h*x`, with `x` in `[-1,1]`, and form
`prod_j(x-(E_j-c)/h)` in a Chebyshev basis. Multiplication uses `x*T_0=T_1`
and `x*T_k=(T_{k-1}+T_{k+1})/2` for `k>=1`. The sum of absolute coefficient
uppers bounds this polynomial on the whole interval because `|T_k|<=1`.
The recurrence and interval bound are reused from [Trefethen, *Approximation
Theory and Approximation Practice*, equations3.8–3.10](https://people.maths.ox.ac.uk/trefethen/ATAP/first6.pdf).
Multiplication by `h^n/n!` supplies the returned exact dyadic upper.

This procedure accounts for node displacement in the remainder itself.
It applies equally to non-Chebyshev nodes. It does not enclose numerical
barycentric arithmetic, multiplication by the physical source columns, or
errors in the field solves at the nodes. Those remain separate inputs to
the local constraint error budget.

Tests include an exactly sharp quadratic bound, the actual49-node table
used by the degree-48 propagator on `[-320,320]`, and a perturbed table.
The rounded49-node polynomial bound differs from the ideal factor by less
than `1e-10` relatively. No Dirac or source evolution is run by this helper,
and no physical gate result follows from it.
