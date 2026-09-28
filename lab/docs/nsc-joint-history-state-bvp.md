# The joint history/state problem has no selector in the current action owner

The requested problem is

$$
C_g(\tau)=U_g(\tau)C_0U_g(\tau)^\dagger,
\qquad
\frac{\delta\Gamma_{\rm one}^{\rm CTP}[g,C[g]]}
     {\delta g_\Delta^A(\tau)}=0.
$$

The existing records decide both proposed solution routes before a new metric
step is allowed.

## Route A stops on the seed constraint

The serialized state lives on the unit-radius seed.  That exact Cauchy pair has

$$
\mathcal C_H=-0.118944813912720,
\qquad
\mathcal C_M=-0.0135927088490713.
$$

These are initial-data constraints.  Simultaneous evolution cannot start from
that pair without first changing a source or geometry dependency, which this
gate forbids.

## Route B is non-unique under every currently enforced condition

The source-selected endpoint passes its fixed-tensor constraints, and a
declared history produces a unitary covariance map.  But the two authenticated
history controls have identical endpoints and give

$$
\max|U_A-U_B|=1.99966699255,
\qquad
\max|C_A-C_B|=0.192089815976105.
$$

Thus at least two distinct state histories satisfy the current endpoint,
unitarity and CAR conditions.  Their different target covariances precede the
fourth-order/local stress evaluation, so no unique finite $r_\star$ stress can
be assigned.

The existing candidate selectors do not remove this freedom.  The condition
$H_\perp=0$ is shared by both endpoints.  Equal-history CTP normalization and
the energy/work Ward identity hold for every unitary history.  The relational
zero-tadpole law selects $V_{\rm full}$, and scale stationarity selects
$\Omega,\zeta$; neither varies the KS history.  The spherical constraints
restrict a supplied slice but do not supply its missing quantum source.

## Exact missing owner

The selected route is therefore **C: non-uniqueness bound and stop**.  The
missing selector is not a new physical principle.  It is the absent executable
history variation of the already declared action:

$$
\boxed{
\frac{\delta\Gamma_{\rm one}^{\rm CTP}[g,C[g]]}
     {\delta g_\Delta^A(\tau)}=0,
\qquad A\in\{N,\beta,q_{\rm ADM},r\}.
}
$$

The `GeneralKSSameActionHistoryFunctional` must combine the blockwise state
propagator, the general-KS fourth-order reference, the local induced history
once, and boundary/interface variations.  The current causal owner evaluates
a supplied history; it does not select or solve one.

No finite stress, null sign or updated constraint is reported because those
quantities differ with the unselected history.  No control duration, old
tensor or linearized density is promoted to a solution.  Coupled evolution
remains closed.

The focused verifier is

```sh
python3 scripts/derive_nsc_joint_history_state_bvp.py --check
```

It reads the authenticated state and U/history artifacts and performs no
source, metric, MMP or historical generator.
