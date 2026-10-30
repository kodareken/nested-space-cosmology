"""Evaluate the oriented-zero tail, conformal, winding and detector benchmarks."""
import argparse
import hashlib
import json
from pathlib import Path
import time

from recursive_horizons.nsc_atlas_zero_interface import calculation

ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path)
    output=(ROOT/p.parse_args().output).resolve()
    if not output.is_relative_to(ROOT) or output.exists():raise ValueError('fresh repository output required')
    start=time.perf_counter();r=calculation();r['elapsed_seconds']=time.perf_counter()-start
    if r['elapsed_seconds']>90:raise RuntimeError('bounded zero-interface calculation exhausted its budget')
    names=['src/recursive_horizons/nsc_atlas_zero_interface.py','src/recursive_horizons/nsc_atlas_s3_polynomials.py',
           'src/recursive_horizons/nsc_atlas_horizon_response.py','src/recursive_horizons/nsc_atlas_parent_circulation.py',
           'src/recursive_horizons/nsc_atlas_parent_scale_selector.py','src/recursive_horizons/nsc_atlas_inherited_charge.py',
           'scripts/check_nsc_atlas_zero_interface.py']
    r['source_sha256']={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}
    output.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    last=r['tail']['samples'][-1];gauge=r['charge_and_winding']['inherited_KK_connection']
    print(json.dumps(dict(H=r['tail']['H'],carrier_frequency_limit=r['tail']['proper_carrier_frequency_limit'],
        tail_powers=[r['tail']['averaged_excitation_density_exponent_in_a'],r['tail']['spatial_wavelength_exponent_in_a']],
        surviving_state=dict(A=last['amplitude_quadrature'],B=last['phase_quadrature']),
        reciprocal_radiation_normalization=r['conformal_break']['reciprocal_comparison']['radiation_density_times_a_fourth'],
        gauge_gap_relation=dict(ratio=gauge['gauge_coupling_squared_over_kappa_mass_increment'],
            g_squared_over_ggrav=gauge['gauge_coupling_squared_over_g_grav'],
            reference_to_minimum_gauge_ratio=gauge['reference_s_zero_to_minimum_gauge_ratio']),
        expected_internal_charge=r['charge_and_winding']['state_charge']['neutral_state_expected_internal_charge'],
        physical_baryon_photon_eta=None,current_Omega_Lambda_derived=False,
        detector_energy_keV_reference_only=r['detector_comparison']['photon_reference_energy_keV'],
        elapsed=r['elapsed_seconds']),indent=2))


if __name__=='__main__':main()
