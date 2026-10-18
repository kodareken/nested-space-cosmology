"""Small reproducible illustration of the written infinite-depth proof."""
import json
import math

import numpy as np

from recursive_horizons.nsc_infinite_inheritance import (
    InheritedTree,geometric_inherited_tree,homogeneous_dilation_source,
    volume_root_observables,
)


def main():
    tree=InheritedTree(np.array([[3.]]),np.array([[1.]]),4.,3)
    exact=2/(3+math.sqrt(6));rows=[]
    for depth in (0,3,6):
        result=tree.response(0.,depth);value=float(result['response'][0,0].real)
        rows.append(dict(depth=depth,response=value,exact_depth_error=exact-value,
            analytic_depth_bound=result['analytic_depth_bound']))
    volume_tree=geometric_inherited_tree(3,.5,spatial_dimension=3)
    volume_mass=.5*(.5**-3+3)
    def calibrated_readout(q):
        result=volume_root_observables(3,q)
        x=.1;upper=1+result['dimensionless_inertial_scale_squared']*x**2
        result['sL']=x
        result['normalized_inverse_response_interval']=[upper-result['dimensionless_s4_remainder_coefficient']*x**4,upper]
        return result
    print(json.dumps(dict(inputs=dict(H=3,B=1,Omega=4,branching=3),
        infinite_static_response=exact,infinite_effective_coefficient=1/exact,
        inherited_local_coefficient_shift=1/exact-3,
        spectral_counting_exponent=tree.counting_exponent,
        second_order_spectral_dimension=2*tree.counting_exponent,
        contraction=tree.contraction,finite_responses=rows,
        spatial_volume_example=dict(branching=3,length_ratio=.5,L=1,A=1,
            static_compliance=1/(volume_tree.gamma*volume_mass),
            compliance_after_doubling_all_spatial_lengths=1/(2*volume_tree.gamma*volume_mass),
            known_force_readout_q_half=calibrated_readout(.5),
            known_force_readout_q_two_fifths=calibrated_readout(.4),
            equal_energy_source=homogeneous_dilation_source(1.,1.)),
        scope='Mathematical inherited-operator example; no metric embedding or cosmological Lambda identification.'),
        indent=2,allow_nan=False))


if __name__=='__main__':
    main()
