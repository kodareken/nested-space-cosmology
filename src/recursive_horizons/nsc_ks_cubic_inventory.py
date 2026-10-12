"""Original channel weights for shared cubic coefficients, not a UV certificate."""
import json,math
from pathlib import Path
from flint import arb,ctx
from .nsc_ks_massive_cubic_uv import CAUCHY_RECORD,CUTOFF_BRIDGE_RECORD,INVENTORY_RECORD
from .nsc_ks_cubic_channel_basis import _required_ball

LARGE_CUTOFF_GROUPS={10,11,12,31,32}


def _integer(value,name):
    if type(value) is not int or value<0:raise ValueError('nonnegative integer '+name+' required')
    return value


def verify_inventory(cauchy,bridge,inventory,*,bits=160):
    channels=cauchy['channels']
    for r in channels:_integer(r['index'],'channel index')
    by_id={r['index']:r for r in channels}
    if len(channels)!=33 or set(by_id)!=set(range(33)):raise ValueError('original channel census required')
    # Check even the zero-angular channels before excluding their cubic term.
    with ctx.workprec(bits):
        for r in channels:
            n=_integer(r['compact_level'],'compact level');l=_integer(r['angular_level'],'angular level')
            ell=_required_ball(r['angular_eigenvalue'],'angular label').union(arb(l*(l+4)).sqrt())
            mass=_required_ball(r['compact_mass'],'mass label').union(n*arb.pi()/2)
            if ell.rad()>arb('1e-14') or mass.rad()>arb('1e-14'):
                raise ValueError('original mass/angular label changed')
    expected={(r['index'],s) for r in channels if r['angular_level']>0 for s in (-1,1)}
    if len(expected)!=60:raise ValueError('original angular census required')
    seen=set();rows=[]
    with ctx.workprec(bits):
        for row in bridge['ledger']:
            if type(row['group']) is not int or type(row['positive_angular_sign']) is not int:
                raise ValueError('explicit channel and angular sign required')
            key=(row['group'],row['positive_angular_sign'])
            if key not in expected or key in seen:raise ValueError('missing, duplicate or foreign channel')
            seen.add(key);channel=by_id[key[0]]
            n=_integer(channel['compact_level'],'compact level');l=_integer(channel['angular_level'],'angular level')
            copies=_integer(channel['copy_count'],'copies');deg=_integer(channel['degeneracy'],'degeneracy')
            if not copies or not deg:raise ValueError('positive multiplicity required')
            mu=copies*deg/2
            if isinstance(row['multiplicity_per_signed_family'],bool) or row['multiplicity_per_signed_family']!=mu:
                raise ValueError('signed-family multiplicity changed')
            ell=_required_ball(channel['angular_eigenvalue'],'angular label').union(arb(l*(l+4)).sqrt())
            mass=_required_ball(channel['compact_mass'],'mass label').union(n*arb.pi()/2)
            if ell.rad()>arb('1e-14') or mass.rad()>arb('1e-14'):raise ValueError('original mass/angular label changed')
            if row['absolute_angular']!=channel['angular_eigenvalue']:raise ValueError('angular ledger mismatch')
            weight=mu*channel['angular_eigenvalue']**2
            if not math.isfinite(row['angular_square_weight_per_energy_sign']) or abs(row['angular_square_weight_per_energy_sign']-weight)>4*math.ulp(weight):
                raise ValueError('angular-square weight mismatch')
            cutoff=inventory['quadrature_measure_by_family'][f'{key[0]}_{key[1]}']
            if cutoff!=(320 if key[0] in LARGE_CUTOFF_GROUPS else 160):raise ValueError('original family cutoff changed')
            rows.append({'group':key[0],'angular_sign':key[1],'angular':ell,'mass':mass,
                         'multiplicity':arb(copies)*deg/2,'cutoff':arb(cutoff)})
        if seen!=expected:raise ValueError('missing angular channel')
    return rows


def load_inventory(root,*,bits=160):
    root=Path(root)
    return verify_inventory(*(json.loads((root/p).read_text()) for p in
        (CAUCHY_RECORD,CUTOFF_BRIDGE_RECORD,INVENTORY_RECORD)),bits=bits)


def channel_moments(rows,*,bits=160):
    """Each angular family contributes mu/(2pi); each energy sign stays explicit."""
    with ctx.workprec(bits):
        raw=[arb(0) for _ in range(3)];leading=[arb(0) for _ in range(3)];seen=set()
        for row in rows:
            key=(row['group'],row['angular_sign'])
            if key in seen:raise ValueError('duplicate channel moment')
            seen.add(key)
            l,m,mu,E=(_required_ball(row.get(k),k) for k in ('angular','mass','multiplicity','cutoff'))
            if not l>0 or not m>=0 or not mu>0 or not E>0:raise ValueError('positive channel weights required')
            for i,value in enumerate((l*l,l**4,m*m*l*l)):
                term=mu*value/(2*arb.pi());raw[i]+=term;leading[i]+=term/(2*E*E)
        return {'raw':tuple(raw),'leading':tuple(leading),'angular_families':len(seen),
                'energy_folding_factor':1,'C_M':None,'physical_local_gate':'OPEN'}


def combine_basis(bases,moments,*,bits=160):
    if set(bases)!={-1,1}:raise ValueError('both characteristic signs required')
    with ctx.workprec(bits):
        def endpoints(v):
            v=_required_ball(v,'domain');return v.lower().man_exp(),v.upper().man_exp()
        ref=bases[1]
        for sign,b in bases.items():
            if type(sign) is not int or type(b['sign']) is not int:raise ValueError('explicit characteristic signs required')
            if b['sign']!=sign or b.get('massless_A') is not True:raise ValueError('basis sign or source changed')
            if b['profile_identity']!=ref['profile_identity'] or endpoints(b['rho_up'])!=endpoints(ref['rho_up']) or endpoints(b['z_domain'])!=endpoints(ref['z_domain']):
                raise ValueError('same profile and domain required')
            if b.get('C_M') is not None:raise ValueError('higher remainder not supplied by this method')
        result={}
        for name in ('raw','leading'):
            weights=tuple(_required_ball(v,'moment') for v in moments[name])
            if len(weights)!=3 or any(not w>=0 for w in weights):raise ValueError('three nonnegative moments required')
            N=arb(0);beta=arb(0)
            for b in bases.values():
                N-=sum((w*_required_ball(b.get(k),k) for w,k in zip(weights,('N2','N4','NM'))),arb(0))
                beta+=sum((w*_required_ball(b.get(k),k) for w,k in zip(weights,('J2','J4','JM'))),arb(0))
            result[name]={'N':N,'beta':beta}
        return {**result,'profile_identity':ref['profile_identity'],'z_domain':ref['z_domain'],
                'C_M':None,'complete_UV_tail':None,'physical_residual_modified':False,
                'physical_local_gate':'OPEN'}
