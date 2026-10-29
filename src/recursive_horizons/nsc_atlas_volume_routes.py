"""Evaluate the three proposed selectors in the declared volume actions."""
import math
from .nsc_atlas_parent_scale_selector import VolumeParent,evolve
from .nsc_atlas_volume_family import member
from .nsc_atlas_inherited_matter import exact_initial_geometry

FRONTIER=((12,2),(6,3),(5,4),(4,5),(3,7))


def calculation():
    rows=[]
    for d,j in FRONTIER:
        m=member(d,j)
        minimum=m['stationary_points'][1]
        final,_=evolve(model=VolumeParent(d,j),endpoint=6.)
        rows.append(dict(d=d,j=j,Lambda_Lref0_squared=m['Lambda_Lref0_squared'],
            s_minimum=minimum['s'],vacuum_minimum=minimum['vacuum'],
            de_Sitter_entropy_times_g=24*math.pi**2/minimum['vacuum'],
            endpoint=final,endpoint_distance_to_own_minimum=abs(final['s']-minimum['s']),
            cross_candidate_bounce_to_12_2=None,
            bounce_status='undefined: no common inter-sector action or interface' if (d,j)!=(12,2) else 'same sector; not a decay to itself',
            ends_in_same_discrete_sector=True))
    for row in rows:
        row['entropy_over_reference']=rows[0]['vacuum_minimum']/row['vacuum_minimum']
    model=VolumeParent()
    maximum,minimum=model.vacuum_roots()
    coefficient=24*math.pi**2*(1/minimum['vacuum']-1/maximum['vacuum'])
    tachyon=-maximum['stiffness']/model.K/(maximum['vacuum']/3)
    return dict(schema='nsc-atlas-volume-selection-routes-v1',frontier=rows,
        ranked_routes=[
            dict(route='c_attractor',result='All five fixed-action flows retain d,j and approach their own minimum; none supplies cross-member selection.'),
            dict(route='b_entropy',result='The reference maximizes late de Sitter entropy among the five under common units; entropy ordering supplies neither a transition nor a maximizing inheritance law.'),
            dict(route='a_BT_CDL',result='Cross-member bounce actions are undefined without a common action and interface; the same reference action also has a lower-energy decompactification channel.')],
        same_action_Hawking_Moss=dict(false_vacuum=minimum['vacuum'],barrier=maximum['vacuum'],
            exponent_times_g=coefficient,g_definition='kappa_nat/L_unit(0)^2',
            negative_curvature_over_H_top_squared=tachyon,
            first_nonconstant_scalar_eigenvalue_over_H_top_squared=4-tachyon,
            homogeneous_scalar_negative_modes=1,
            runaway_limit='V(s)->0+ as s->+infinity',
            interpretation='A standard HM saddle in the homogeneous reduced action, not an absolute rate, a dominance proof or a localized higher-dimensional decay calculation.'),
        missing_common_action=['a common state space for d and the ancestor source',
            'inter-sector kinetic terms, interpolating potential or charged interface',
            'inherited initial state and measure; determinants for rates'],
        parameter_free_selection_obtained=False,
        physical_selection_derived=False,
        interpolation_nonidentifiability='Adding lambda*z^4*(1-z)^4 to an interpolation leaves endpoint potential and its first two derivatives fixed but changes its barrier; endpoint data do not fix a bounce.',
        absolute_SI_scale_selected=False)
