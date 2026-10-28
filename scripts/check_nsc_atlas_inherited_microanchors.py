"""Fixed SI coupling-ratio lengths versus the still-free parent curvature radius."""
import argparse
import hashlib
import json
import math
from pathlib import Path

from recursive_horizons.nsc_atlas_scale_selection import calculation as scale_calculation

ROOT=Path(__file__).resolve().parents[1]


def calculation():
    scale=scale_calculation()
    kappa=scale['SI_conventions']['kappa_natural_metres_squared']
    ehat,sigmahat,what=.5,.2,1e9
    anchors=dict(charge_over_tension_metres=ehat*math.sqrt(kappa)/sigmahat,
                 inverse_length_frequency_over_charge_metres=what*math.sqrt(kappa)/ehat,
                 sqrt_inverse_length_frequency_over_tension_metres=math.sqrt(what*kappa/sigmahat))
    witnesses=[]
    for unit in (.001,100.):
        e=ehat/(math.sqrt(kappa)*unit)
        sigma=sigmahat/(kappa*unit)
        frequency=what/unit  # canonical inverse-metres frequency, not s^-1
        witnesses.append(dict(arbitrary_model_unit_metres=unit,
            charge_over_tension_metres=e/sigma,
            inverse_length_frequency_over_charge_metres=frequency/e,
            sqrt_inverse_length_frequency_over_tension_metres=math.sqrt(frequency/sigma),
            parent_Lambda_m_inverse2=scale['fixed_parent']['Lambda_model_unit_squared']/unit**2))
    return dict(schema='nsc-atlas-inherited-microanchors-v1',
        fixed_normalized_inputs=dict(charge=ehat,tension=sigmahat,omega=what),
        SI_coupling_ratio_lengths=anchors,
        fixed_dimensionless_microanchor_ratio=what*sigmahat/ehat**2,
        dimensional_completion_witnesses=witnesses,
        inverse_target_hierarchy_not_a_prediction=scale['inverse_target_diagnostic']['Lref_metres']/anchors['charge_over_tension_metres'],
        distinction='The supplied normalized coupling ratios and G,hbar,c fix these microscopic SI combinations. They do not select the parent curvature length L_ref or its hierarchy above those combinations.',
        Lambda_obs_derived=False,
        source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in (
            'scripts/check_nsc_atlas_inherited_microanchors.py',
            'src/recursive_horizons/nsc_atlas_scale_selection.py')})


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();out=(ROOT/args.output).resolve()
    if not out.is_relative_to(ROOT) or out.exists():raise ValueError('new repository output required')
    value=calculation();out.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(value))


if __name__=='__main__':main()
