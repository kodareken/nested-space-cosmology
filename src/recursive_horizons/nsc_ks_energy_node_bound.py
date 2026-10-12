"""Certified interpolation remainder factor for actual binary energy nodes.

The node polynomial is enclosed in a Chebyshev basis on the full interval.
This bounds displacement of the nodes without pretending they are exact
Chebyshev roots. Numerical interpolation and nodal solve errors are separate.
"""
from math import factorial

import numpy as np
from flint import arb,ctx

from .nsc_ks_ball_trajectory import exact_upper,restored_upper


def actual_node_remainder_factor(nodes,interval,*,bits=160):
    """Upper sup_E |prod(E-x_j)|/n! for the given exact binary nodes.

    For a Banach-valued n-times differentiable function, multiply this by
    sup ||d_E^n f||. Values at the actual nodes are assumed exact here.
    """
    nodes=np.asarray(nodes,float);interval=np.asarray(interval,float)
    if (interval.shape!=(2,) or not np.isfinite(interval).all() or interval[1]<=interval[0]
            or nodes.ndim!=1 or len(nodes)<2 or not np.isfinite(nodes).all()
            or np.any(np.diff(nodes)<=0) or nodes[0]<interval[0] or nodes[-1]>interval[1]):
        raise ValueError('distinct increasing finite nodes inside a positive-length interval required')
    with ctx.workprec(bits):
        lo,hi=map(lambda x:arb(float(x)),interval)
        half=(hi-lo)/2;center=(hi+lo)/2
        coefficients=[arb(1)]
        for node in nodes:
            eta=(arb(float(node))-center)/half
            new=[arb(0) for _ in range(len(coefficients)+1)]
            for k,c in enumerate(coefficients):
                new[k]-=eta*c
                if k==0:new[1]+=c
                else:
                    new[k-1]+=c/2
                    new[k+1]+=c/2
            coefficients=new
        normalized=sum((c.abs_upper() for c in coefficients),arb(0))
        polynomial=half**len(nodes)*normalized
        factor=polynomial/factorial(len(nodes))
        return {
            'nodes':nodes.tolist(),'interval':interval.tolist(),'node_count':len(nodes),
            'degree':len(nodes)-1,'derivative_order':len(nodes),'precision_bits':bits,
            'normalized_node_polynomial_sup_upper':exact_upper(normalized),
            'node_polynomial_sup_upper':exact_upper(polynomial),
            'remainder_factor_upper':exact_upper(factor),
            'exact_binary_nodes_enclosed':True,'ideal_Chebyshev_nodes_assumed':False,
            'node_displacement_included':True,'node_solve_error_bound':None,
            'interpolation_arithmetic_error_bound':None,'physical_local_gate':'OPEN',
        }


def multiply_derivative_bound(factor_record,derivative_upper,*,bits=160):
    """Conditional actual-node remainder; no numerical field error inferred."""
    with ctx.workprec(bits):
        derivative=restored_upper(derivative_upper) if isinstance(derivative_upper,dict) else arb(derivative_upper)
        if not derivative.is_finite() or not derivative>=0:
            raise ValueError('finite nonnegative certified derivative upper required')
        return exact_upper(restored_upper(factor_record['remainder_factor_upper'])*derivative)
