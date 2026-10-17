"""Save a compact, deterministic summary of the exact Step-7 calculations."""
import argparse
import hashlib
import json
from pathlib import Path
import time

from recursive_horizons.nsc_atlas_effective_vacuum import (
    classical_vacuum_calculation, conformal_anomaly_trace_projection,
    conformal_anomaly_scale_si, constant_curvature_trace_condition,
    conformal_interface_heat_coefficients,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    started = time.process_time()
    products = {}
    for name, calculate in (
        ('classical', classical_vacuum_calculation),
        ('anomaly', conformal_anomaly_trace_projection),
        ('scaled', conformal_anomaly_scale_si),
        ('curvature', constant_curvature_trace_condition),
        ('interface', conformal_interface_heat_coefficients),
    ):
        products[name] = calculate()
        if time.process_time()-started > 60:
            raise RuntimeError('exact vacuum summary exceeded 60 CPU seconds')
    classical, anomaly, scaled, curvature, interface = products.values()
    source_paths = (
        'scripts/record_nsc_atlas_effective_vacuum.py',
        'src/recursive_horizons/nsc_atlas_effective_vacuum.py',
        'tests/test_nsc_atlas_effective_vacuum.py',
    )
    root = Path(__file__).resolve().parents[1]
    record = {
        'schema': 'nsc-atlas-effective-vacuum-v1',
        'domain': 'Classical smooth on-shell bulk; quantum terms on the smooth quiet cap and Schwarzschild exterior; formal quiet-interface polynomial.',
        'classical_added_lambda_Lref_squared': str(classical['induced_lambda_Lref_squared']),
        'quantum_assumption': 'One free real four-dimensional massless conformal scalar, Hadamard bulk state, local covariant renormalization.',
        'quantum_cap_coefficient_of_zeta': str(anomaly['cap_coefficient']),
        'quantum_exterior_coefficient_of_zeta': str(anomaly['exterior_coefficient']),
        'reference_shape': {'h': '1', 'm': '2/5', 'r': '11/20'},
        'fiducial_SI_inputs': {key: str(scaled[key]) for key in
                              ('ell_m', 'G_SI', 'G_standard_uncertainty_SI',
                               'planck_h_SI_exact', 'c_SI_exact')},
        'constant_source': 'https://physics.nist.gov/cuu/pdf/wall_2022.pdf',
        'zeta_at_fiducial_scale': str(scaled['zeta'].evalf(20)),
        'reference_quantum_trace': {
            region: {key: str(value.evalf(20)) for key, value in row.items()}
            for region, row in scaled['regions'].items()},
        'uncertainty_scope': 'Standard input uncertainty in G only; no total model or observational error bound.',
        'constant_curvature_trace_residual': str(curvature['trace_equation_residual']),
        'Lambda0_zero_X_roots': [str(value) for value in curvature['unsourced_X_roots']],
        'curved_root_H2_planck_length2': str(curvature['unsourced_curved_H2_planck_length2']),
        'interface_a2_integrand': str(interface['a2_integrand']),
        'formal_interface_density_times_ell3': str(interface['normalized_density_ell3']),
        'formal_interface_reference': str(interface['reference_normalized_density']),
        'unresolved': ['Full renormalized stress and joint state',
                       'Lorentzian quantum interface construction',
                       'Physical exterior reduction and metric variation',
                       'Cosmological scale selection and observational tests'],
        'zero_residual_checks': {
            name: {key: str(value) for key, value in result['checks'].items()}
            for name, result in products.items() if 'checks' in result},
        'source_sha256': {name: hashlib.sha256((root/name).read_bytes()).hexdigest()
                          for name in source_paths},
        'reproduce': 'python scripts/run_python.py scripts/record_nsc_atlas_effective_vacuum.py --output results/summaries/nsc-atlas-effective-vacuum-v1.json',
    }
    args.output.write_text(json.dumps(record, indent=2, sort_keys=True, allow_nan=False)+'\n',
                           encoding='utf-8', newline='\n')
    print(f'Wrote {args.output}; {time.process_time()-started:.3f} CPU seconds; no evolution.')


if __name__ == '__main__':
    main()
