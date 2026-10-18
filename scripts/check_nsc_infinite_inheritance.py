"""Small reproducible illustration of the written infinite-depth proof."""
import json
import math

import numpy as np

from recursive_horizons.nsc_infinite_inheritance import InheritedTree


def main():
    tree=InheritedTree(np.array([[3.]]),np.array([[1.]]),4.,3)
    exact=2/(3+math.sqrt(6));rows=[]
    for depth in (0,3,6):
        result=tree.response(0.,depth);value=float(result['response'][0,0].real)
        rows.append(dict(depth=depth,response=value,exact_depth_error=exact-value,
            analytic_depth_bound=result['analytic_depth_bound']))
    print(json.dumps(dict(inputs=dict(H=3,B=1,Omega=4,branching=3),
        infinite_static_response=exact,infinite_effective_coefficient=1/exact,
        inherited_local_coefficient_shift=1/exact-3,
        spectral_counting_exponent=tree.counting_exponent,
        second_order_spectral_dimension=2*tree.counting_exponent,
        contraction=tree.contraction,finite_responses=rows,
        scope='Mathematical inherited-operator example; no metric embedding or cosmological Lambda identification.'),
        indent=2,allow_nan=False))


if __name__=='__main__':
    main()
