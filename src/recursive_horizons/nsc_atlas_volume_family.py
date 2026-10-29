"""Selection audit of the declared parent-volume family, without new dynamics.

Reuse the fixed one-member implementation. General extrema and envelope
identities expose which extra matching coefficients change across members.
The old v1 record and its producer remain reproducible at their source pin.
"""
from fractions import Fraction as F
import math

from scipy.optimize import brentq

from .nsc_atlas_inherited_matter import exact_initial_geometry, outward_decimal_interval
from .nsc_atlas_parent_scale_selector import VolumeParent, matched_coefficients


def initial_stability_polynomial(d, ancestor_steps):
    t = (ancestor_steps+3)**2
    return d*(t-18)-2*t-30


def minimum_stable_dimension(ancestor_steps):
    if not isinstance(ancestor_steps, int) or ancestor_steps < 1:
        raise ValueError('positive integer ancestor depth required')
    t = (ancestor_steps+3)**2
    if t <= 18:
        return None
    return (2*t+30)//(t-18)+1


def member(d, ancestor_steps):
    """Two root brackets valid for every integer d>2; no parameter scan."""
    if not isinstance(d, int) or d <= 2:
        raise ValueError('integer d>2 required for the convex-root proof')
    model = VolumeParent(d=d, ancestor_steps=ancestor_steps)
    c = matched_coefficients(d, ancestor_steps)
    A, B, T = (float(c[k]) for k in ('A','B','T'))
    f = lambda x: d*A-(d+2)*B*x+2*d*T*x**(d/2)
    # F(0)>0, F(1)=-d*E0/2<0, and one positive turning point.
    high = 2.
    for _ in range(12):
        if f(high) > 0:
            break
        high *= 2
    else:
        raise RuntimeError('declared root-bracketing budget exhausted')
    roots = (brentq(f,0.,1.,xtol=5e-15), brentq(f,1.,high,xtol=5e-15))
    rows=[]
    for x in roots:
        s=-math.log(x)/2
        V, ss = model.potential(s), model.potential(s,2)
        rows.append(dict(x=x,s=s,vacuum=V,stiffness=ss,
                         kind='minimum' if ss>0 else 'maximum'))
    minimum=rows[1]
    curvature=float(exact_initial_geometry()['curvature'])
    return dict(d=d,ancestor_depth=ancestor_steps,
        ancestor_flux=str(c['ancestor_flux']),sheet_T=str(c['T']),
        required_A=str(c['A']),required_B=str(c['B']),required_UV_bare=str(c['required_UV_bare']),
        initial_stiffness=str(c['initial_stiffness']),
        initial_homogeneous_stable=c['initial_stiffness']>0,
        stationary_points=rows,Lambda_Lref0_squared=minimum['vacuum']/curvature,
        inherited_discrete_choice=False)


def inverse_matching(A, B):
    """Conditional inference if A and B are independently known at s=0."""
    A,B=F(A),F(B)
    rho,E0=F(9,32),F(3,8)
    T=B+rho-A
    if T<=0 or B<=0:
        raise ValueError('positive sheet and curvature required')
    return dict(T=str(T),d=str(2*B/(rho+T+E0/2)),j=math.sqrt(32*float(T))-3)


def minimum_enclosure(d, j):
    """Exact rational enclosure sufficient to order the five frontier values.

    y=exp(-s)>1 turns stationarity into an integer-power polynomial even
    for odd d. Forty bisections suffice to separate the frontier values.
    """
    c=matched_coefficients(d,j)
    A,B,T=(c[k] for k in ('A','B','T'))
    f=lambda y:d*A-(d+2)*B*y*y+2*d*T*y**d
    lo,hi=F(1),F(2)
    if not f(lo)<0<f(hi):
        raise ValueError('the finite rational bracket does not apply')
    for _ in range(40):
        mid=(lo+hi)/2
        if f(mid)<0:lo=mid
        else:hi=mid
    C=exact_initial_geometry()['curvature']
    lower=(A*lo**d-B*hi**(d+2)+T*lo**(2*d))/C
    upper=(A*hi**d-B*lo**(d+2)+T*hi**(2*d))/C
    return outward_decimal_interval(lower,upper,places=8)


def calculation():
    pairs=((10,2),(12,2),(16,2),(30,2),(12,3),(30,4),(6,3),(5,4),(4,5),(3,7))
    witnesses=[member(d,j) for d,j in pairs]
    frontier_pairs=((12,2),(6,3),(5,4),(4,5),(3,7))
    frontier=[next(r for r in witnesses if (r['d'],r['ancestor_depth'])==pair) for pair in frontier_pairs]
    minimum=min(frontier,key=lambda r:r['Lambda_Lref0_squared'])
    intervals=[minimum_enclosure(r['d'],r['ancestor_depth']) for r in frontier]
    assert all(F(intervals[0][1])<F(row[0]) for row in intervals[1:])
    return dict(schema='nsc-atlas-volume-family-selection-v1',
        reference_source_commit='59ef213',
        inherited_source=dict(rho='9/32',E0='3/8',charge='1/2',parent_flux='3/2',epsilon='1/4'),
        ancestor_meaning='j is the depth of ONE assigned ancestor vacuum sheet, not a sum of j independent sheets.',
        stability_condition='d*((j+3)^2-18)>2*(j+3)^2+30',
        minimum_stable_dimensions=[dict(j=j,d_min=minimum_stable_dimension(j)) for j in range(1,8)],
        witnesses=witnesses,
        stable_integer_frontier=[dict(d=r['d'],j=r['ancestor_depth'],coefficient=r['Lambda_Lref0_squared'],
            exact_outward_enclosure=interval) for r,interval in zip(frontier,intervals)],
        family_energy_minimum=dict(d=minimum['d'],j=minimum['ancestor_depth'],
            coefficient=minimum['Lambda_Lref0_squared'],unique=True,
            proof='Strict coordinate monotonicity reduces all initially stable integer pairs to (12,2),(6,3),(5,4),(4,5),(3,7). Their minimum is (12,2).',
            physical_selection_derived=False,
            limitation='A minimum over rematched effective actions is not a derived transition or ground-state rule of one common parent theory.'),
        rectangle_10_to_30_and_2_to_4=dict(
            lower=witnesses[0]['Lambda_Lref0_squared'],upper=witnesses[5]['Lambda_Lref0_squared'],
            lower_has_stable_initial_point=False,
            lower_with_initial_stability=witnesses[1]['Lambda_Lref0_squared'],
            method='Analytic monotonicity in d and T_j reduces the rectangle extrema to corners; no dynamics or exhaustive parameter scan.'),
        monotonicity=dict(
            variables='u=d*s, x=exp(-2*u/d), S0=rho+T+E0/2',
            dV_dT='exp(-u)*((d/2-1)-(d/2)*x+x^(d/2))>0 at the late minimum',
            dV_dd='exp(-u)*S0/2*(1-(1+2*u/d)*exp(-2*u/d))>0 at the late minimum',
            justification='Convexity of x^(d/2), exp(y)>=1+y, and stationarity eliminate the implicit minimum-location derivative.'),
        conditional_inverse_example=inverse_matching(7,F(15,2)),
        inverse_warning='The reference A=7 and B=15/2 were defined from (12,2). Inverting them does not independently select it.',
        common_action_warning='Changing d changes the parent theory; changing j here also rematches A,B and bulk energy. This audit is a family of conditional actions, not transitions within one completed UV model.',
        conclusion='Inheritance covariance and homogeneous stability do not select d=12 or j=2. 9.7022370093 is the computed late coefficient for that declared member.',
        missing_selection='An independent common parent geometry/action/state must fix d and the ancestor-source coupling or provide a normalized branch-selection rule. No such equation is currently supplied.',
        dimension_selected=False,ancestor_depth_selected=False,SI_scale_selected=False,
        new_field_evolutions=0)
