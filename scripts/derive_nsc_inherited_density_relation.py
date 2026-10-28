"""An epoch-independent density relation from the frozen inherited-source bound."""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import sympy as sp

from recursive_horizons.nsc_atlas_inherited_matter import (
    exact_initial_geometry,uniform_bound,outward_decimal_interval,
)

ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/'results/summaries/nsc-atlas-inherited-density-relation-v1.json'


def calculation():
    A,H,C,d,h2,R=sp.symbols('A H C d h2 R',positive=True)
    Ochi=d*R/(A**3*H*H);Ov=h2/(H*H);Ok=-C/(A*A*H*H)
    residual=sp.factor((-Ok)**3/(Ochi**2*Ov)-C**3/(d*d*h2*R*R))
    assert residual==0
    g=exact_initial_geometry()
    coefficient=g['curvature']**3/((g['excitation']/3)**2*(g['vacuum']/3))
    rows=[]
    for w in (20,2000,10**9):
        b=uniform_bound(w)
        lower,upper=coefficient/b['upper']**2,coefficient/b['lower']**2
        rows.append(dict(omega=w,
            density_fraction_invariant_interval=outward_decimal_interval(lower,upper,places=18),
            relative_lower_gap=float(1-lower/coefficient),relative_upper_gap=float(upper/coefficient-1)))
    source='results/summaries/nsc-atlas-inherited-matter-forecast-v1.json'
    # The comparison at E/rho_v=1 is descriptive; the interval above controls
    # finite-frequency phase dependence. Omega_chi means the scalar, not a
    # measured total matter density containing additional fields.
    b=float(coefficient)**(1/3)
    return dict(schema='nsc-atlas-inherited-density-relation-v1',
        result='An epoch- and length-unit-independent vacuum/scalar/curvature correlation fixed by the same inherited initial geometry and energy.',
        definitions=dict(Omega_chi='kappa E/(3H^2)',Omega_vacuum='kappa rho_v/(3H^2)',
            Omega_curvature='-C/(A^2 H^2)',comoving_ratio='R=A^3 E/E0'),
        exact_identity='(-Omega_curvature)^3/(Omega_chi^2*Omega_vacuum)=C_Omega/R^2',
        exact_residual=str(residual),coefficient_exact=str(coefficient),coefficient_approximate=float(coefficient),
        frequency_members=rows,
        dust_comparison_energy_equality=dict(excitation_over_vacuum=1,
            Omega_vacuum=1/(2-b),Omega_chi=1/(2-b),Omega_curvature=-b/(2-b),
            role='Approximate decimals for the dust-comparison curve at equal excitation/vacuum, not astronomical values or the finite-frequency unaveraged equality.'),
        fixed_source='Same protected homogeneous scalar/four-form parent; turning-point preparation, from A=1 onward.',
        independent_of=['observer epoch A','Hubble normalization','physical SI length unit'],
        supplied_inputs_remain=['ancestral preparation','flux/mixing sector','initial boundary geometry','frequency member'],
        observations_used=False,
        background_bound=dict(path=source,sha256=hashlib.sha256((ROOT/source).read_bytes()).hexdigest()),
        source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in
            ('src/recursive_horizons/nsc_atlas_inherited_matter.py','scripts/derive_nsc_inherited_density_relation.py')})


if __name__=='__main__':
    if OUTPUT.exists():raise RuntimeError('creation-only derived record already exists')
    r=calculation()
    OUTPUT.write_text(json.dumps(r,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(r,indent=2))
