"""Connect inherited geometry, actual receiver stress and candidate source gates."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

from recursive_horizons.nsc_atlas_parent_circulation import Parent
from recursive_horizons.nsc_atlas_valve_source_gate import (
    radial_source_gate,geometric_spinor_candidate,mmp_massless_field_gate,
    dirac_clifford_control,neutral_thermal_ec_benchmark,
    axial_pair_correlation_witness,
)

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def calculation():
    start=time.perf_counter();H=Parent.retained().H
    paths=('results/summaries/nsc-atlas-valve-conformal-v1.json',
           'results/summaries/nsc-atlas-valve-restart-v1.json')
    records=[json.loads((ROOT/p).read_text(encoding='utf-8')) for p in paths]
    snapshots=[('old_frame_endpoint',records[0]['levels'][-1]['endpoint']),
               ('new_controlled_cut',records[1]['levels'][-1]['last_actual_cut'])]
    if not records[1]['levels'][-1]['final_cut_admitted']:
        raise ValueError('no admitted restarted cut; do not substitute a coarse diagnostic')
    gates=[]
    for name,s in snapshots:
        gate=radial_source_gate(radius=s['center_radius'],radius_rate=s['center_Rproper'],
            density=s['center_rho'],radial_pressure=s['center_radial_pressure'],
            radial_spacetime_Hessian_over_radius=s['center_spacetime_Hess_ss_over_R_H2']*H*H,
            actual_acceleration=s['center_Rproper_acc'])
        gates.append(dict(name=name,proper_time=s['center_proper_time'],gate=gate))
    spinor=geometric_spinor_candidate()
    sources=('src/recursive_horizons/nsc_atlas_valve_source_gate.py',
             'src/recursive_horizons/nsc_atlas_inherited_charge.py',
             'src/recursive_horizons/nsc_atlas_parent_circulation.py',
             'src/recursive_horizons/nsc_atlas_parent_scale_selector.py',
             'scripts/check_nsc_atlas_valve_source_gate.py')
    return dict(schema='nsc-atlas-valve-source-gate-v1',
        producer_parent_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        source_sha256={name:sha(ROOT/name) for name in sources},
        input_sha256={name:sha(ROOT/name) for name in paths},actual_receiver_gates=gates,
        geometric_spinor_candidate=spinor,Clifford_control=dirac_clifford_control(),
        four_dimensional_minimal_EC_thermal_benchmark=neutral_thermal_ec_benchmark(spinor),
        full_sixteen_dimensional_thermal_contact=dict(
            added_field='Full complex massless16D Dirac, both internal signs and vertices retained',
            benchmark='Free unpolarized nonrotating homogeneous rest-frame mu=0 Fermi state, identical +/- internal occupations, fixed vacuum-Wick subtraction',
            coincident_projectors='P_n,plus(y,y)=P_n,minus(y,y)=d_n/(64 V_gamma) I64',
            coincident_thermal_covariance_difference=0.,fully_vacuum_Wick_ordered_contact=0.,
            full_renormalized_actual_neck_source=None,
            warning='The finite 4D axial benchmark is not a matched full16D torsion truncation. All transformed vertices and KK signs must agree; no generic GH-state or bare composite zero is claimed.'),
        MMP_field_gate=mmp_massless_field_gate(spinor),
        non_Gaussian_pair_correlation_witness=axial_pair_correlation_witness(),
        sources=[dict(reference='Camporesi-Higuchi gr-qc/9505009,3.1 and5.2',url='https://arxiv.org/abs/gr-qc/9505009'),
                 dict(reference='Corrected minimalEC sign,2008.11686v2',url='https://arxiv.org/html/2008.11686v2'),
                 dict(reference='Diakonov-Tumanov-Vladimirov1104.2432v2,64-76',url='https://arxiv.org/html/1104.2432v2'),
                 dict(reference='Maldacena-Milekhin-Popov1807.04726v3,2.3 and5.2',url='https://arxiv.org/html/1807.04726v3')],
        route_status=[
            'Current canonical preparation: no turnaround by the controlled small-radius cut; no all-future exclusion.',
            'Positive vacuum alone has zero direct radial-null stress; changing a cap also requires changed geometry/junction and inherited state.',
            'Existing scalar fields supply no EC spin current. Added neutral free thermal4D fermions give positive NEC; the full balanced16D thermal-Wick contact cancels.',
            'Untwisted positive-curvature KK fermions do not meet the imported MMP massless-field assumption; magnetic topology and return state are also missing.',
            'Non-Gaussian/nonequilibrium inherited stress, or a separately specified geometric transition, remains to be calculated. No negative fluid is inserted.'],
        verdict=dict(regular_expanding_child_proved=False,new_field_selected_by_old_action=False,
            inherited_spinor_population_selected=False,SI_scale_selected=False,
            physical_baryon_photon_eta=None,geometry_derived_mass_and_coupling_relations=True),
        numerical=dict(elapsed_seconds=time.perf_counter()-start,observation_fits=0,raw_arrays_committed=False),
        reproduction='python scripts/run_python.py scripts/check_nsc_atlas_valve_source_gate.py')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'results/summaries/nsc-atlas-valve-source-gate-v1.json')
    args=parser.parse_args();record=calculation();args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(output=str(args.output),gates=record['actual_receiver_gates'],
                         spinor_mass=record['geometric_spinor_candidate']['lowest_mass'],
                         MMP_mass_times_radius=record['MMP_field_gate']['lowest_mass_times_extremal_radius']),indent=2))


if __name__=='__main__':main()
