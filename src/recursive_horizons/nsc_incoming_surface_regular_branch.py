"""Affine normal-data chart from reused cu and constant principal entries."""
from functools import lru_cache
from math import pi as stored_pi

import mpmath as mp
import sympy as sp

from .nsc_incoming_surface_quadratic_lapse import reference_quadratic_group
from .nsc_incoming_source_quadrature_bound import pack,float_enclosure
from .nsc_incoming_vacuum_tail_bound import _precision,_lo,_hi,_range


@lru_cache(maxsize=1)
def affine_identities():
    m,L,Ha,Hr,w,K=sp.symbols('m L H_a H_r w K',real=True)
    a,r,D=sp.symbols('a r D',positive=True)
    M2=m*m+L*L
    Aref=-D*a*L*L/(120*sp.pi*r)*(4*m*m*Hr/M2**3+(Hr-Ha)/M2**2)
    Alocal=-a*r*(K*(Ha-Hr)+(Ha+Hr)/(60*sp.pi))
    dref=-D*a*L*L/(120*sp.pi*r*M2**2);dlocal=a*r*(K+1/(60*sp.pi))
    ar1=sp.diff(Aref.subs(Hr,Hr+w/r),w)
    al1=sp.diff(Alocal.subs(Hr,Hr+w/r),w)
    expected_ref=-D*a*L*L*(L*L+5*m*m)/(120*sp.pi*r*r*M2**3)
    expected_local=a*(K-1/(60*sp.pi))
    values={
        'reference_affine_remainder':Aref.subs(Hr,Hr+w/r)-Aref-w*ar1,
        'local_affine_remainder':Alocal.subs(Hr,Hr+w/r)-Alocal-w*al1,
        'reference_derivative':ar1-expected_ref,'local_derivative':al1-expected_local,
        'reference_reused_d':ar1-dref/r*(1+4*m*m/M2),
        'local_reused_d':al1-dlocal/r+a/(30*sp.pi),
        'reference_independent_C':ar1+dref/r+a*a*reference_quadratic_group(m,L,a,r,D,sp.pi),
        'local_independent_C':al1+dlocal/r+a*a*(-2*K/a),
    }
    A0,A1,d=sp.symbols('A0 A1 d',real=True,nonzero=True)
    A=A0+A1*w;H=Hr+w/r;B=-(2*A+H*d)/a**2;e=-d/a**2
    R0=A0+Hr*d;R1=A1+d/r;det=d*(R0+R1*w)/a**2
    values.update(principal_determinant=A*e-B*d-det,
                  lapse_root=A.subs(w,-A0/A1),principal_root=det.subs(w,-R0/R1),
                  determinant_at_lapse_root=det.subs(w,-A0/A1)-d*d*(Hr-A0/(A1*r))/a**2)
    C=sp.Symbol('C',real=True)
    values['determinant_slope_from_C']=sp.diff(det,w).subs(A1,-a*a*C-d/r)+d*C
    residuals={k:str(sp.factor(v)) for k,v in values.items()}
    if set(residuals.values())!={'0'}:raise ArithmeticError('affine chart identity failed')
    return {'residuals':residuals,'A1_reference':str(expected_ref),'A1_local':str(expected_local),
            'A':'A0+A1*w','B':'-(2*A+(Hr+w/r)*d)/a^2','e':'-d/a^2',
            'determinant':'d*(A0+Hr*d+(A1+d/r)*w)/a^2',
            'roots':{'lapse_pivot':'-A0/A1','principal':'-(A0+Hr*d)/(A1+d/r)'},
            'C_consistency':'A1+d/r=-a^2*C', 'Cu_Cv_producers_called':False}


def _hull(value,stored):
    point=mp.mpf(float(stored))
    return _range(min(_lo(value),point),max(_hi(value),point))


def enclose_branch(lapse,principal,channels,cauchy,*,precision=80):
    if [c['index'] for c in channels]!=list(range(33)):
        raise ValueError('unchanged retained33 inventory required')
    matrix=principal['matrix'];components=matrix['components']
    oldA=lapse['coefficient']['total_coefficient_interval']
    if matrix['reused_cu_interval']!=oldA or not lapse['coefficient']['strictly_positive']:
        raise ValueError('unchanged certified A0 required')
    with _precision(precision):
        pi=_hull(mp.iv.pi,stored_pi)
        a=mp.iv.mpf(components['a_interval']);Hr=mp.iv.mpf(components['H_r_interval'])
        r=_hull(mp.iv.sqrt(2),cauchy['probe_normal_geometry']['intrinsic_r'])
        A0=mp.iv.mpf(oldA);d=mp.iv.mpf(matrix['matrix_intervals'][1][0]);e=mp.iv.mpf(matrix['matrix_intervals'][1][1])
        local_d=mp.iv.mpf(components['local_d_interval'])
        local_A1=local_d/r-a/(30*pi)
        ref_A1=mp.iv.mpf(0);per_group=[]
        rows=components['per_group_reference_d']
        if [row['group'] for row in rows]!=list(range(1,33)):raise ValueError('saved per-group d inventory changed')
        for channel,row in zip(channels[1:],rows):
            m=_hull(channel['compact_level']*mp.iv.pi/2,channel['compact_mass'])
            ell=_hull(mp.iv.sqrt(channel['angular_level']*(channel['angular_level']+4)),channel['angular_eigenvalue'])
            L=ell/r;M2=m*m+L*L
            if _lo(M2)<=0:raise ValueError('gapped reference group required')
            dg=mp.iv.mpf(row['d_interval'])
            value=mp.iv.mpf(0) if channel['angular_eigenvalue']==0 else dg/r*(1+4*m*m/M2)
            ref_A1+=value;per_group.append({'group':channel['index'],'A1_interval':float_enclosure(value)})
        A1=local_A1+ref_A1;R0=A0+Hr*d;R1=A1+d/r
        if not (_lo(A0)>0 and _hi(d)<0 and _hi(A1)<0 and _lo(R0)>0 and _hi(R1)<0):
            raise ArithmeticError('affine root/sign ordering not certified')
        wA=-A0/A1;wD=-R0/R1
        if not (0<_lo(wD) and _hi(wD)<_lo(wA)):
            raise ArithmeticError('distinct positive roots not certified in the expected order')
        A_at_D=A0+A1*wD
        det_at_A=d*d*(Hr+wA/r)/(a*a)
        B_at_A=-(Hr+wA/r)*d/(a*a)
        if _lo(A_at_D)<=0 or _lo(det_at_A)<=0 or _lo(B_at_A)<=0:
            raise ArithmeticError('pivot and principal degeneracy not separated')
        det0=d*R0/(a*a);e_residual=e+d/(a*a)
        importedB=mp.iv.mpf(matrix['matrix_intervals'][0][1]);B0=-(2*A0+Hr*d)/(a*a)
        if not _lo(e_residual)<=0<=_hi(e_residual) or max(_lo(importedB),_lo(B0))>min(_hi(importedB),_hi(B0)):
            raise ArithmeticError('imported d/e/B principal relations inconsistent')
        C_implied=-R1/(a*a)
        values={'A1_local':local_A1,'A1_reference':ref_A1,'A1_total':A1,'R0':R0,'R1':R1,
                'w_lapse_pivot':wA,'w_principal':wD,'A_at_principal_root':A_at_D,
                'determinant_at_lapse_pivot':det_at_A,'B_at_lapse_pivot':B_at_A,
                'determinant_at_zero':det0,'C_implied_for_independent_check':C_implied}
        return {'intervals':{k:float_enclosure(v) for k,v in values.items()},
                'binary_intervals':{k:pack(v) for k,v in values.items()},'A1_reference_groups':per_group,
                'imported_unchanged':{'A0':oldA,'d':matrix['matrix_intervals'][1][0],'e':matrix['matrix_intervals'][1][1]},
                'intrinsic_intervals':{'a':components['a_interval'],'r':float_enclosure(r),'Hr':components['H_r_interval']},
                'strict_root_order':'0 < w_principal < w_lapse_pivot',
                'regular_component_containing_zero':{'exact':'(-infinity,w_principal)',
                    'upper_root_enclosure':float_enclosure(wD),
                    'guaranteed_numeric_subinterval':{'lower':'-infinity','upper_exclusive':float_enclosure(wD)[0]},
                    'signs':{'A':'positive','determinant':'negative','eliminated_K':'negative'}},
                'lapse_pivot_classification':'A=0 but determinant>0 and B>0; this elimination pivot fails, original principal matrix remains invertible',
                'principal_classification':'determinant=0 with A>0 and d<0; rank1 highest-spatial-derivative matrix in this family',
                'root_role':'finite normal-data values delta(r_T), not z or time endpoints',
                'spacetime_singularity_claimed':False,'precision':precision,
                'coefficient_derivation':'new analytic slope only, using stored d components; A0,d,e producers not rerun'}
