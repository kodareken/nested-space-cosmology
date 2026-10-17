"""Evaluate the actual prepared exterior response and its conformal metric test.

One fixed 129-node reference and one 129-node Dirichlet solve. No coupled
evolution, new preparation, field-grid scan or quantum determinant is run.
The response pairing is a boundary diagnostic, not a retarded action.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from recursive_horizons.nsc_atlas_exterior_dtn import (
    kruskal_factor, proper_null_speeds, quiet_outward_curvatures,
    radius_null_derivatives, solve_exterior_dirichlet, weyl_boundary_data,
)
from recursive_horizons.nsc_atlas_field_preparation import (
    preparation_record, radius_from_null, source_profile,
)
from recursive_horizons.nsc_atlas_prepared_response import (
    PreparedFieldBudget, _geometry, _solve,
)


def calculate():
    budget = PreparedFieldBudget(cpu_seconds=120, kernel_values_max=30_000_000)
    geometry = _geometry(129)
    Va, Vb = geometry['Va'], geometry['Vb']
    reference, _, _, dense = _solve(
        129, geometry, preparation_record(), budget, include_dense=True)

    def trace(V):
        coordinate = np.asarray(V, dtype=float)
        row = reference((coordinate-Va)/(Vb-Va))
        return (np.reshape(row['wall_chi_over_epsilon'], coordinate.shape),
                np.reshape(row['wall_dchi_dV_over_epsilon'], coordinate.shape))

    solution = solve_exterior_dirichlet(
        trace, source=lambda V: source_profile(V, Va, Vb), budget=budget)
    cuts = np.unique(np.r_[dense['cuts'], solution.cuts])
    rows = {}
    ward_absolute = 0.
    ward_relative = 0.
    orientation_error = 0.
    for order in (32, 64):
        node, weight = np.polynomial.legendre.leggauss(order)
        total = 0.
        bare_total = 0.
        for lower, upper in zip(cuts[:-1], cuts[1:]):
            half = (Vb-Va)*(upper-lower)/2
            V = Va+(Vb-Va)*(lower+upper)/2+half*node
            row = solution.boundary(V)
            U, prime = geometry['wall'](V), geometry['wall'](V, 1)
            radius = radius_from_null(U, V)
            udot, vdot = proper_null_speeds(prime, radius)
            rU, rV = radius_null_derivatives(U, V, radius)
            beta_plus = .2/radius**2-1.25*radius
            orientation_error = max(orientation_error, float(np.max(np.abs(
                udot*rU-vdot*rV+beta_plus))))
            # Integrate the constant-on-S2 field, not an implicitly normalized Y00.
            measure = 4*np.pi*radius**2*np.sqrt(kruskal_factor(radius)*prime)
            pair = row['phi']*row['conformal_robin']
            total += half*float(weight @ (measure*pair))
            bare_total += half*float(weight @ (measure*row['phi']*row['outward_dphi']))
            # A nonconstant ambient Weyl jet, including its normal derivative.
            omega = .07*U*V+.03*V**2+.02*U
            normal_omega = udot*(.07*V+.02)-vdot*(.07*U+.06*V)
            transformed = weyl_boundary_data(
                row['phi'], row['outward_dphi'], row['outward_K'], omega, normal_omega)
            transformed_pair = np.exp(3*omega)*transformed['phi']*transformed['conformal_robin']
            difference = np.abs(transformed_pair-pair)
            ward_absolute = max(ward_absolute, float(np.max(difference)))
            ward_relative = max(ward_relative, float(np.max(difference/(1+np.abs(pair)))))
        rows[str(order)] = {'robin_pairing_per_epsilon_squared': total,
                            'bare_normal_pairing_per_epsilon_squared': bare_total}
    sample = solution.boundary(np.linspace(Va, Vb, 9))
    record = {
        'schema': 'nsc-atlas-exterior-schur-metric-v1',
        'domain': 'Quiet exterior IV wedge, wall r/ell in [0.45,0.55], monopole.',
        'source': 'Original unit-amplitude C3 source_profile, zero bottom, original prepared wall value and analytic first derivative.',
        'normalization': 'chi=r*phi; tabulated phi=ell*phi_physical/epsilon, Robin=ell^2*B(phi_physical)/epsilon; angular field constant on S2.',
        'outward_curvature_times_ell': quiet_outward_curvatures(),
        'boundary_quadrature': rows,
        'absolute_quadrature_difference': {
            key: abs(rows['32'][key]-rows['64'][key]) for key in rows['64']},
        'nine_wall_samples': [
            {key: float(sample[key][index]) for key in
             ('V', 'phi', 'outward_dphi', 'conformal_robin', 'chi_U')}
            for index in range(9)],
        'metric_checks': {
            'maximum_sampled_weyl_pairing_absolute_residual': ward_absolute,
            'maximum_sampled_weyl_pairing_scaled_residual': ward_relative,
            'maximum_sampled_outward_radial_orientation_residual': orientation_error},
        'classical_added_local_bulk_Lambda_times_Lref_squared': 0,
        'bulk_variation': 'Compact support strictly inside the retained cap; fixed interface, exterior geometry and characteristic data. The classical exterior boundary influence has no support there.',
        'reference_scale': 'Lref=ell; the exact zero is independent of ell. P001 ell=10000 m is a fiducial input, not a derived cosmological ruler.',
        'supplied_cap_Lambda0_times_ell_squared': 3,
        'uncertainty_scope': 'The bulk zero is a support identity for the declared variation. Quadrature differences and Ward residuals are numerical diagnostics, not continuum field or model error bounds.',
        'not_inferred': ['A vanishing boundary stress or full metric response',
                         'A self-adjoint action from the driven retarded pairing',
                         'A quantum covariance, determinant or full renormalized stress',
                         'A new nonzero cosmological Lambda or observational match'],
        'numerical_budget': {'CPU_seconds_limit': 120, 'kernel_values_limit': 30_000_000,
                             'kernel_values_charged': budget.kernel_values},
    }
    budget.check()
    return record, time.process_time()-budget.started


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    record, elapsed = calculate()
    root = Path(__file__).resolve().parents[1]
    paths = ('src/recursive_horizons/nsc_atlas_exterior_dtn.py',
             'src/recursive_horizons/nsc_atlas_prepared_response.py',
             'src/recursive_horizons/nsc_atlas_field_preparation.py',
             'tests/test_nsc_atlas_exterior_dtn.py',
             'scripts/record_nsc_atlas_exterior_schur.py')
    record['source_sha256'] = {
        name: hashlib.sha256((root/name).read_bytes()).hexdigest() for name in paths}
    record['reproduce'] = ('python scripts/run_python.py scripts/record_nsc_atlas_exterior_schur.py '
                           '--output results/summaries/nsc-atlas-exterior-schur-metric-v1.json')
    args.output.write_text(json.dumps(record, indent=2, sort_keys=True, allow_nan=False)+'\n',
                          encoding='utf-8', newline='\n')
    print(json.dumps({'summary': args.output.as_posix(), 'CPU_seconds': elapsed,
                      'kernel_values': record['numerical_budget']['kernel_values_charged'],
                      'classical_added_local_bulk_Lambda_times_Lref_squared': 0,
                      'robin_pairing': record['boundary_quadrature']['64'],
                      'checks': record['metric_checks']}))


if __name__ == '__main__':
    main()
