#!/usr/bin/env python3
"""One retained initial transfer, independently checked in Bernstein form.

Print only compact evidence. No radial solve, evolution or resolution scan.
Six exact points on each physical polynomial piece certify its identity;
floating agreement at those points is a separate numerical consistency check.
"""
from bisect import bisect_right
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from diagnose_nsc_atlas_boundary_rows import authenticate, Budget, plain
from check_nsc_atlas_cap_background_domain import retained_slice
from recursive_horizons.nsc_atlas_cap_fixed_transfer import transfer_initial_background

ROOT=Path(__file__).resolve().parents[1]
SOURCES=(
    'scripts/check_nsc_atlas_cap_fixed_transfer.py',
    'scripts/check_nsc_atlas_cap_background_domain.py',
    'scripts/diagnose_nsc_atlas_boundary_rows.py',
    'src/recursive_horizons/nsc_atlas_cap_fixed_basis.py',
    'src/recursive_horizons/nsc_atlas_cap_fixed_transfer.py',
    'src/recursive_horizons/nsc_atlas_cap_background_jets.py',
)


class BernsteinSource:
    """Independent exact-binary Hermite data evaluation via de Casteljau."""
    def __init__(self,background,budget):
        data={key:tuple(Q(float(v)) for v in values) for key,values in background.inputs.items()}
        self.labels=data['labels']
        self.quiet=all(data[key][0]==0 for key in ('values','first','second'))
        self.controls=[]
        for i,h in enumerate(b-a for a,b in zip(self.labels,self.labels[1:])):
            y0,d0,s0=(data[key][i] for key in ('values','first','second'))
            y1,d1,s1=(data[key][i+1] for key in ('values','first','second'))
            value=(y0,y0+h*d0/5,y0+2*h*d0/5+h*h*s0/20,
                   y1-2*h*d1/5+h*h*s1/20,y1-h*d1/5,y1)
            first=tuple(5*(b-a)/h for a,b in zip(value,value[1:]))
            second=tuple(4*(b-a)/h for a,b in zip(first,first[1:]))
            self.controls.append((value,first,second))
            budget.check(384)

    def evaluate(self,u,budget):
        u=Q(u)
        if u>self.labels[-1] or (u<self.labels[0] and not self.quiet):
            raise ValueError('Bernstein query exits the exact source domain')
        budget.check(384)
        if u<self.labels[0]: return (Q(0),)*3
        i=min(bisect_right(self.labels,u)-1,len(self.labels)-2)
        x=(u-self.labels[i])/(self.labels[i+1]-self.labels[i])
        result=[]
        for controls in self.controls[i]:
            row=controls
            while len(row)>1:
                row=tuple((1-x)*a+x*b for a,b in zip(row,row[1:]))
            result.append(row[0])
        return tuple(result)


def polynomial_value(poly,x,h,derivative=0):
    if derivative:
        poly=tuple(p*poly[p]/h for p in range(1,len(poly)))
    result=Q(0)
    for coefficient in reversed(poly): result=result*x+coefficient
    return result


def compare_transfer(background,transfer,budget):
    source=BernsteinSource(background,budget)
    basis=transfer['basis'];T=transfer['initial_time'];Z=transfer['wall_coordinate']
    keys=('chi','canonical_density','chi_z','canonical_density_z')
    errors={key:0. for key in keys};peaks=errors.copy()
    equalities=points=0;boundaries={}
    for cell in range(basis.cells):
        low,high=basis.exact_knots[cell:cell+2];h=high-low
        cut=min(Q(1),(Z-low)/h)
        coordinates=[cut*Q(j,7) for j in range(1,7)]
        if cell==0: coordinates.insert(0,Q(0))
        if cell==basis.cells-1: coordinates.append(cut)
        for xi in coordinates:
            z=low+h*xi
            plus=source.evaluate(T+z,budget);minus=source.evaluate(T-z,budget)
            expected=(plus[0]-minus[0],plus[1]-minus[1],
                      plus[1]+minus[1],plus[2]+minus[2])
            polys=transfer['exact_polynomials']
            exact=(polynomial_value(polys['chi'][cell],xi,h),
                   polynomial_value(polys['canonical_density'][cell],xi,h),
                   polynomial_value(polys['chi'][cell],xi,h,1),
                   polynomial_value(polys['canonical_density'][cell],xi,h,1))
            if exact!=expected:
                raise ValueError('independent exact Bernstein transfer identity failed')
            row=basis.values(cell,float(xi),budget=budget);idx=row['indices']
            c=transfer['chi'][idx];d=transfer['canonical_density'][idx]
            numeric=(row['value']@c,row['value']@d,row['first']@c,row['first']@d)
            for key,actual,target in zip(keys,numeric,expected):
                error=abs(float(actual)-float(target))
                if not np.isfinite(error): raise ValueError('nonfinite transfer discrepancy')
                errors[key]=max(errors[key],error);peaks[key]=max(peaks[key],abs(float(target)))
            if z==0 or z==Z:
                boundaries['axis' if z==0 else 'wall']={key:float(v) for key,v in zip(keys,expected)}
            equalities+=4;points+=1;budget.check(256)
    tolerances={key:5e-12*max(1.,peaks[key]) for key in keys}
    if any(errors[key]>tolerances[key] for key in keys):
        raise ValueError('fixed transfer numerical consistency exceeds its declared tolerance')
    return dict(exact_polynomial_identities=4*basis.cells,exact_point_equalities=equalities,evaluation_points=points,
                method='six rational interior points per physical piece plus axis and wall; independent Bernstein Hermite evaluation',
                boundary_targets=boundaries,
                maximum_absolute_error=errors,reference_peaks=peaks,numerical_tolerances=tolerances,
                numerical_check_is_uniform_enclosure=False)


def main():
    budget=Budget(maximum_evaluations=3_000_000,cpu_seconds=30.)
    phase='authenticate'
    try:
        current,_,_,identity=authenticate(ROOT,budget)
        background,selection=retained_slice(current)
        phase='initial transfer'
        transfer=transfer_initial_background(background,selection['T'],selection['Z'],budget=budget)
        phase='independent polynomial check'
        comparison=compare_transfer(background,transfer,budget)
        phase='clipped Gram factor'
        basis=transfer['basis'];gram=basis.gram(transfer['wall_coordinate'],budget=budget)
        budget.check()
    except Exception:
        print(json.dumps(dict(failed_phase=phase,charged_values=budget.evaluations,
                              cpu_seconds=time.process_time()-budget.started),sort_keys=True))
        raise
    selection=dict(selection)
    selection['scope']='retained zero-case initial total fields; recorded alpha is metadata only; no evolved finite state'
    report=dict(schema='nsc-atlas-cap-fixed-transfer-v1',selection=selection,
        input={key:identity[key] for key in ('path','sha256','manifest_sha256','authenticated_input_files')},
        source_bindings={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCES},
        representation=dict(cells=basis.cells,field_dimension=basis.dimension,
                            initial_cut_fraction=0.5,scaled_gram_condition=gram.condition,
                            condition_times_machine_epsilon=gram.condition*np.finfo(float).eps),
        comparison=comparison,coefficient_conversion_bounds=transfer['coefficient_conversion_bounds'],
        numerical_evaluation_enclosure=False,physical_evolution=False,
        execution=dict(charged_values=budget.evaluations,cpu_seconds=time.process_time()-budget.started,
                       transfers=1,gram_factors=1,radial_solves=0,time_steps=0,trajectory_replays=0))
    print(json.dumps(plain(report),indent=2,sort_keys=True,allow_nan=False))


if __name__=='__main__': main()
