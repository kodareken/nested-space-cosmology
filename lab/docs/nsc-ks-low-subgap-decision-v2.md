# Low/subgap source decision, v2

This successor corrects the continuation point after inspecting the later
source records. The group-12 order-24 middle bound and the group-32 combined
LOW/middle bound already pass their component tolerances, and V5 has composed
the group-32 `[24,32]` replacement. Those calculations are reused; they are
not tasks to repeat.

The still-missing quantity is the physical reconstruction error of the
stored `rho=1` mode columns on panels below each V5 joined threshold, including
the subgap panels, followed by correlated transport through the final history.
Exact ell-zero radius response does not erase the nonzero baseline, and the
group-13 nested-node movement is not a uniform source-error certificate.

The historical approximately `4.998e-9` aggregate was useful before the
order-24 successors. It is not relabelled as the current remaining low/subgap
bound. The current bound remains `None` until the physical-column enclosures
are supplied.

```sh
python3 scripts/derive_nsc_ks_low_subgap_decision_v2.py --check
PYTHONPATH=src .venv/validation/bin/python -m pytest -q tests/test_nsc_ks_low_subgap_decision_v2.py
```
