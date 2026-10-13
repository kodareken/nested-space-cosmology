"""Directed one-energy affine-vacuum projector enclosure without state wrapping.

The point column is a numerical approximation, never a replacement for an
archived coherent mode. Each exact frozen Hermitian exponential is enclosed
locally, then recentered; all discarded widths enter an additive error radius.
"""
import mpmath as mp

from .nsc_incoming_vacuum_tail_bound import _precision,_horizon_enclosure,_lo,_hi,_range
from .nsc_incoming_validated_vacuum_propagation import _sinc_minus_one as stable_sinc_minus_one
from .nsc_incoming_source_quadrature_bound import pack,unpack,float_enclosure

ENERGY=0.5562120090641313
STEPS=2048
PRECISION=60
DELTA0='1e-8'


def point(value):return list(mp.mpf(value)._mpf_)
def unpoint(value):return mp.make_mpf(tuple(int(x) for x in value))
def ivpoint(value):
    p=mp.mpf(value)._mpf_
    return mp.iv.make_mpf((p,p))
def cp(value):return [point(value.real),point(value.imag)]
def ucp(value):return mp.mpc(unpoint(value[0]),unpoint(value[1]))
def ivcomplex(value):return mp.iv.mpc(ivpoint(value.real),ivpoint(value.imag))
def cpack(value):return [pack(value.real),pack(value.imag)]
def cunpack(value):return mp.iv.mpc(unpack(value[0]),unpack(value[1]))
def midpoint(value):return (_lo(value)+_hi(value))/2
def upper(value):return ivpoint(_hi(value))
def abs_upper(value):return upper(abs(value))
def vector_norm(values):return mp.iv.sqrt(sum((abs(v)**2 for v in values),mp.iv.mpf(0)))
def describe(value):return {'binary_interval':pack(value),'float_enclosure':float_enclosure(value)}


class SelectedGenerator:
    """Exact real-profile H=[[-D,conj(p)],[p,D]] in y=log(delta)."""
    def __init__(self,channel,config,energy=ENERGY):
        if (channel['index']!=14 or channel['compact_level']!=1 or channel['angular_level']!=1
                or energy!=ENERGY or config['magnetic_flux']!=4):
            raise ValueError('only selected positive group14 low-energy component is owned')
        self.energy=ivpoint(energy)
        self.horizon,self.horizon_derivative,_=_horizon_enclosure(config['horizon_rho'],config['surface_gravity'],mp.iv.dps)
        self.qh=mp.iv.pi/2+mp.iv.atan2(self.horizon,mp.iv.mpf(1))
        self.b0=self.horizon_derivative
        s,c=mp.iv.sin(2*self.qh),mp.iv.cos(2*self.qh)
        self.linear=-3*s-c;self.sinc_coefficient=-3*c+s
        exact_m=mp.iv.pi/2;exact_l=mp.iv.sqrt(5)
        self.mass=_range(min(_lo(exact_m),mp.mpf(channel['compact_mass'])),max(_hi(exact_m),mp.mpf(channel['compact_mass'])))
        self.angular=_range(min(_lo(exact_l),mp.mpf(channel['angular_eigenvalue'])),max(_hi(exact_l),mp.mpf(channel['angular_eigenvalue'])))

    def geometry(self,delta):
        # These functions are real on this real domain. The reused complex
        # sinc enclosure has a harmless imaginary remainder; retain its real
        # projection explicitly rather than giving H complex diagonal entries.
        s1=stable_sinc_minus_one(delta).real
        s2=stable_sinc_minus_one(2*delta).real
        B=(self.b0+self.linear*delta*(1+s1)**2+self.sinc_coefficient*s2).real
        sine=mp.iv.sin(self.qh-delta).real
        if _lo(B)<=0 or _lo(sine)<=0:
            raise ArithmeticError('positive exact-profile geometry not enclosed')
        return B,1/sine

    def components(self,y):
        delta=mp.iv.exp(y);B,r=self.geometry(delta)
        factor=mp.iv.sqrt(delta)/mp.iv.sqrt(B)
        return (self.energy/B,-factor*self.mass*r,factor*self.angular)

    def tail(self,start):
        nominal=mp.iv.mpf(DELTA0);actual=mp.iv.exp(ivpoint(start))
        delta_upper=max(_hi(nominal),_hi(actual))
        domain=_range(mp.mpf(0),delta_upper);B,r=self.geometry(domain)
        A=upper(mp.iv.sqrt((self.mass*r)**2+self.angular**2)/mp.iv.sqrt(B))
        bound=upper(2*A*mp.iv.sqrt(ivpoint(delta_upper)))
        return {'nominal_delta0':pack(nominal),'actual_start_delta':pack(actual),
                'tail_delta_upper':point(delta_upper),'offdiagonal_coefficient_upper':pack(A),
                'affine_projector_tail_error_upper':pack(bound)}

    def binding(self):
        return {name:describe(value) for name,value in (
            ('exact_horizon_rho',self.horizon),('exact_horizon_q',self.qh),
            ('exact_geometric_2kappa',self.b0),('mass_hull',self.mass),('positive_angular_hull',self.angular),
            ('fixed_binary_energy',self.energy))}


def frozen_update(column,D,x,y,step):
    """Outward exp(-i h Hmid) acting on one exact dyadic point column."""
    D,x,y=map(ivpoint,(D,x,y));p=mp.iv.mpc(x,y);pb=mp.iv.mpc(x,-y)
    n=mp.iv.sqrt(D*D+x*x+y*y)
    if _lo(n)<=0:raise ArithmeticError('selected frozen Hamiltonian norm must exclude zero')
    cosine=mp.iv.cos(step*n);sine_over=mp.iv.sin(step*n)/n
    v0,v1=map(ivcomplex,column)
    image=(cosine*v0-mp.iv.j*sine_over*(-D*v0+pb*v1),
           cosine*v1-mp.iv.j*sine_over*(p*v0+D*v1))
    centered=tuple(mp.mpc(midpoint(v.real),midpoint(v.imag)) for v in image)
    rounding=upper(vector_norm([v-ivcomplex(c) for v,c in zip(image,centered)]))
    return centered,image,rounding


def cell_defect(whole,frozen,step):
    # Operator norm of a traceless Hermitian difference. Its three real
    # coordinates are independently enclosed; no commutativity is assumed.
    radii=[abs_upper(value-ivpoint(center)) for value,center in zip(whole,frozen)]
    return upper(step*mp.iv.sqrt(sum((v*v for v in radii),mp.iv.mpf(0))))


def propagate(channel,config,*,progress=None):
    """Exactly one fixed 2048-cell run; no adaptive refinement or normalization."""
    with _precision(PRECISION):
        G=SelectedGenerator(channel,config)
        start=_lo(mp.iv.ln(mp.iv.mpf(DELTA0)))
        endpoint=mp.iv.ln(G.qh-3*mp.iv.pi/4)
        finish=midpoint(endpoint)
        tail=G.tail(start)
        column=(mp.mpc(1),mp.mpc(0));error=mp.iv.mpf(0);rows=[];left=start
        for index in range(STEPS):
            right=finish if index==STEPS-1 else start+(finish-start)*(index+1)/STEPS
            if not left<right:raise ArithmeticError('positive consecutive real steps required')
            step=ivpoint(right)-ivpoint(left)
            whole=G.components(_range(left,right));atmid=G.components(ivpoint((left+right)/2))
            frozen=tuple(midpoint(v) for v in atmid)
            if any(not _lo(box)<=v<=_hi(box) for box,v in zip(atmid,frozen)):
                raise ArithmeticError('frozen exact point left its midpoint enclosure')
            norm=upper(vector_norm([ivcomplex(v) for v in column]));defect=cell_defect(whole,frozen,step)
            after,image,rounding=frozen_update(column,*frozen,step)
            increment=upper(defect*norm+rounding);error=upper(error+increment)
            rows.append({'left':point(left),'right':point(right),'step':pack(step),
                'H_whole':[pack(v) for v in whole],'H_midpoint':[pack(v) for v in atmid],
                'frozen_H':[point(v) for v in frozen],'column_before':[cp(v) for v in column],
                'update_image':[cpack(v) for v in image],'column_after':[cp(v) for v in after],
                'point_norm_upper':pack(norm),'propagator_defect_upper':pack(defect),
                'recenter_error_upper':pack(rounding),'column_error_increment':pack(increment)})
            column=after;left=right
            if progress and (index+1)%256==0:progress(index+1)
        bridge_y=_range(min(finish,_lo(endpoint)),max(finish,_hi(endpoint)))
        endpoint_H=G.components(bridge_y)
        Hnorm=upper(mp.iv.sqrt(sum((abs_upper(v)**2 for v in endpoint_H),mp.iv.mpf(0))))
        distance=abs_upper(endpoint-ivpoint(finish));bridge=upper(2*Hnorm*distance)
        payload={'schema':'NSC-LOW-ENERGY-VACUUM-PROOF-v1','precision':PRECISION,'steps':STEPS,
            'group':14,'energy':ENERGY,'channel':dict(channel),'source_config':dict(config),'geometry':G.binding(),
            'start_y':point(start),'finish_y':point(finish),'physical_endpoint_y':pack(endpoint),'tail':tail,
            'endpoint_bridge':{'range_y':pack(bridge_y),'H_range':[pack(v) for v in endpoint_H],
                'H_norm_upper':pack(Hnorm),'distance_upper':pack(distance),'projector_error_upper':pack(bridge)},
            'rows':rows,'final_column':[cp(v) for v in column],
            'algorithm':'point columns; exact Hermitian frozen matrices; whole-step interval defect; additive unitary recurrence; no wrapping or normalization'}
        payload['result']=replay(payload)
        return payload


def replay(payload):
    """Replay error sums and final predicate only; no geometry/ODE propagation."""
    if payload['precision']!=PRECISION or payload['group']!=14 or payload['energy']!=ENERGY or payload['steps']!=STEPS or len(payload['rows'])!=STEPS:
        raise ValueError('complete selected fixed-size proof required')
    with _precision(payload['precision']):
        previous=unpoint(payload['start_y']);column=(mp.mpc(1),mp.mpc(0))
        error=mp.iv.mpf(0);defect_sum=mp.iv.mpf(0);rounding_sum=mp.iv.mpf(0)
        for row in payload['rows']:
            left,right=unpoint(row['left']),unpoint(row['right'])
            if left!=previous or right<=left:raise ValueError('proof steps have a gap, overlap or reversed endpoint')
            if [cp(v) for v in column]!=row['column_before']:raise ValueError('point-column chain changed')
            step=ivpoint(right)-ivpoint(left)
            if pack(step)!=row['step']:raise ValueError('step interval changed')
            whole=list(map(unpack,row['H_whole']));frozen=list(map(unpoint,row['frozen_H']))
            middle=list(map(unpack,row['H_midpoint']))
            if any(not _lo(box)<=v<=_hi(box) for box,v in zip(middle,frozen)):raise ValueError('frozen point is not in its midpoint range')
            norm=upper(vector_norm([ivcomplex(v) for v in column]));eta=cell_defect(whole,frozen,step)
            after=tuple(map(ucp,row['column_after']));image=tuple(map(cunpack,row['update_image']))
            rounding=upper(vector_norm([v-ivcomplex(c) for v,c in zip(image,after)]))
            increment=upper(eta*norm+rounding)
            for value,key in ((norm,'point_norm_upper'),(eta,'propagator_defect_upper'),(rounding,'recenter_error_upper'),(increment,'column_error_increment')):
                if pack(value)!=row[key]:raise ValueError('directed step bound differs: '+key)
            defect_sum=upper(defect_sum+upper(eta*norm));rounding_sum=upper(rounding_sum+rounding)
            error=upper(error+increment);previous=right;column=after
        if previous!=unpoint(payload['finish_y']) or [cp(v) for v in column]!=payload['final_column']:
            raise ValueError('final endpoint/column changed')
        endpoint=unpack(payload['physical_endpoint_y']);finish=unpoint(payload['finish_y']);bridge_data=payload['endpoint_bridge']
        bridge_range=unpack(bridge_data['range_y'])
        if _lo(bridge_range)>min(finish,_lo(endpoint)) or _hi(bridge_range)<max(finish,_hi(endpoint)):
            raise ValueError('endpoint bridge does not cover the physical endpoint')
        qh=unpack(payload['geometry']['exact_horizon_q']['binary_interval'])
        if pack(mp.iv.ln(qh-3*mp.iv.pi/4))!=pack(endpoint):raise ValueError('physical rho1 endpoint changed')
        H=list(map(unpack,bridge_data['H_range']));Hnorm=upper(mp.iv.sqrt(sum((abs_upper(v)**2 for v in H),mp.iv.mpf(0))))
        distance=abs_upper(endpoint-ivpoint(finish));bridge=upper(2*Hnorm*distance)
        if pack(Hnorm)!=bridge_data['H_norm_upper'] or pack(distance)!=bridge_data['distance_upper'] or pack(bridge)!=bridge_data['projector_error_upper']:
            raise ValueError('endpoint bridge replay differs')
        tail=unpack(payload['tail']['affine_projector_tail_error_upper'])
        nominal=unpack(payload['tail']['nominal_delta0']);actual=unpack(payload['tail']['actual_start_delta'])
        if pack(mp.iv.exp(ivpoint(unpoint(payload['start_y']))))!=pack(actual):raise ValueError('actual start coordinate changed')
        delta=unpoint(payload['tail']['tail_delta_upper']);A=unpack(payload['tail']['offdiagonal_coefficient_upper'])
        if delta<max(_hi(nominal),_hi(actual)) or pack(upper(2*A*mp.iv.sqrt(ivpoint(delta))))!=pack(tail):
            raise ValueError('omitted affine tail does not cover the actual starting coordinate')
        norm=vector_norm([ivcomplex(v) for v in column]);projector_numeric=upper((1+upper(norm))*error)
        total=upper(tail+projector_numeric+bridge)
        P01=(ivcomplex(column[0])*mp.iv.mpc(ivpoint(column[1].real),-ivpoint(column[1].imag))).real
        enclosed=P01+_range(-_hi(total),_hi(total))
        norm_defect=abs_upper(norm*norm-1)
        return {'final_Re_P01':describe(enclosed),'point_outer_product_Re_P01':describe(P01),
            'point_column_norm':describe(norm),'point_column_norm_squared_defect_upper':describe(norm_defect),
            'column_error_upper':describe(error),'summed_unitary_defect_column_error_upper':describe(defect_sum),
            'summed_recenter_error_upper':describe(rounding_sum),
            'numerical_projector_error_upper':describe(projector_numeric),'affine_tail_projector_error_upper':describe(tail),
            'endpoint_bridge_projector_error_upper':describe(bridge),'total_projector_error_upper':describe(total),
            'strict_Re_P01_above_one_quarter':bool(_lo(enclosed)>mp.mpf(1)/4),
            'lossless_step_replay_residual':0.,'column_normalized':False}

from pathlib import Path as _receipt_Path
from .nsc_frozen_sources import selected_source_sha256 as _receipt_source_sha256

# Retained executable control: vacuum_control
# Folded from scripts/derive_nsc_incoming_low_energy_vacuum_enclosure.py at 329ef809; original receipts retain their producer identities.
import argparse as _vacuum_control_argparse
import gzip as _vacuum_control_gzip
from hashlib import sha256 as _vacuum_control_sha256
import json as _vacuum_control_json
from pathlib import Path as _vacuum_control_Path
import signal as _vacuum_control_signal
import sys as _vacuum_control_sys
import time as _vacuum_control_time
import numpy as _vacuum_control_np
from recursive_horizons.nsc_posix import ROOT as _vacuum_control_ROOT
_vacuum_control_sys.path.insert(0, str(_vacuum_control_ROOT / 'src'))
from recursive_horizons.nsc_incoming_low_energy_vacuum_enclosure import ENERGY as _vacuum_control_ENERGY, STEPS as _vacuum_control_STEPS, PRECISION as _vacuum_control_PRECISION, propagate as _vacuum_control_propagate, replay as _vacuum_control_replay
_vacuum_control_OUTPUT = 'results/development/nsc-incoming-low-energy-vacuum-enclosure.json'
_vacuum_control_SOURCES = ('src/recursive_horizons/nsc_incoming_low_energy_vacuum_enclosure.py', 'scripts/derive_nsc_incoming_low_energy_vacuum_enclosure.py', 'tests/test_nsc_incoming_low_energy_vacuum_enclosure.py', 'docs/nsc-incoming-low-energy-vacuum-enclosure.md')
_vacuum_control_INPUTS = ('results/development/nsc-mode-resolved-cauchy-state.json', 'results/development/nsc-compact-matched-restart.json', 'results/development/nsc-pg-ctp-mode-jets.json', 'src/recursive_horizons/nsc_incoming_vacuum_tail_bound.py', 'src/recursive_horizons/nsc_incoming_validated_vacuum_propagation.py', 'src/recursive_horizons/nsc_incoming_source_quadrature_bound.py', 'src/recursive_horizons/nsc_paired_horizon_preparation.py', 'src/recursive_horizons/nsc_pg_massive_modes.py', 'src/recursive_horizons/nsc_unruh_state.py', 'src/recursive_horizons/nsc_lorentzian.py')
_vacuum_control_CPU_CAP = 120.0

def _vacuum_control_digest(path):
    candidate = _receipt_Path(path)
    relative = candidate.relative_to(_vacuum_control_ROOT).as_posix() if candidate.is_absolute() else candidate.as_posix()
    return _receipt_source_sha256(_vacuum_control_ROOT, relative)

def _vacuum_control_signature():
    return {p: _vacuum_control_digest(p) for p in (*_vacuum_control_SOURCES, *_vacuum_control_INPUTS)}

def _vacuum_control_inputs():
    channel = _vacuum_control_json.loads((_vacuum_control_ROOT / _vacuum_control_INPUTS[0]).read_text())['channels'][14]
    config = _vacuum_control_json.loads((_vacuum_control_ROOT / _vacuum_control_INPUTS[1]).read_text())['scattering_provenance']['config']
    reference = _vacuum_control_json.loads((_vacuum_control_ROOT / _vacuum_control_INPUTS[2]).read_text())
    spec = reference['payload']
    if _vacuum_control_digest(spec['path']) != spec['sha256']:
        raise ValueError('selected source-mode artifact changed')
    with _vacuum_control_np.load(_vacuum_control_ROOT / spec['path'], allow_pickle=False) as a:
        if float(a['14_1/energies'][3]) != _vacuum_control_ENERGY or a['14_1/mass'].item() != channel['compact_mass'] or a['14_1/angular'].item() != channel['angular_eigenvalue']:
            raise ValueError('selected energy/channel changed')
    return (channel, config, spec)

def _vacuum_control_make_record(payload, path):
    if payload['signature'] != _vacuum_control_signature():
        raise ValueError('projector proof/input source changed')
    channel, config, spec = _vacuum_control_inputs()
    if payload['channel'] != channel or payload['source_config'] != config or payload['selected_source_artifact'] != spec:
        raise ValueError('physical/source parameter binding changed')
    if payload.get('budget_exceeded'):
        result = None
        status = 'OPEN: fixed CPU budget exhausted before certificate'
    else:
        result = _vacuum_control_replay(payload)
        if result != payload['result']:
            raise ValueError('directed projector error replay differs')
        status = 'PASS: directed physical affine-vacuum Re(P01)>1/4' if result['strict_Re_P01_above_one_quarter'] else 'OPEN: fixed directed enclosure does not establish Re(P01)>1/4'
    return {'schema': 'NSC-INCOMING-LOW-ENERGY-VACUUM-ENCLOSURE-v1', 'accountable_author': 'Douglas Ek', 'status': status, 'group': 14, 'angular_sign': 1, 'energy': _vacuum_control_ENERGY, 'channel': channel, 'source_config': config, 'result': result, 'geometry': payload.get('geometry'), 'initialization': payload.get('tail'), 'endpoint_bridge': payload.get('endpoint_bridge'), 'algorithm': {'steps': _vacuum_control_STEPS, 'precision': _vacuum_control_PRECISION, 'point_column_initial': [1, 0], 'state_wrapping': False, 'column_normalization': False, 'error_recurrence': 'epsilon_next=epsilon+step*sup(norm(H-Hmid))*norm(point_column)+recenter_radius', 'projector_error': 'affine_tail+(1+norm(point_column))*column_error+endpoint_bridge', 'replay': 'directed sum, point chain and coordinate/error accounting; no repeat geometry propagation'}, 'runtime': payload['runtime'], 'selected_source_artifact': spec, 'scope': {'physical_affine_dominant_vacuum_component': True, 'coherent_source_replaced': False, 'coherent_thermal_remainder': 'separate exact owner; not discarded or included in this vacuum projector error', 'archived_modes_or_state_normalized_or_changed': False, 'Jost_phase_changed': False, 'source_occupation_kappa_changed': False, 'group32_bound_applied_to_low_energy': False, 'source_precision_3e_minus_11_claimed': False, 'physical_initial_data_or_metric_evolution': False, 'full_C0_matching_or_extended_stationarity': 'OPEN', 'publication': False}, 'source_hashes': {p: payload['signature'][p] for p in _vacuum_control_SOURCES}, 'input_hashes': {p: payload['signature'][p] for p in _vacuum_control_INPUTS}, 'payload': {'path': path.relative_to(_vacuum_control_ROOT).as_posix(), 'sha256': _vacuum_control_sha256(path.read_bytes()).hexdigest(), 'bytes': path.stat().st_size}, 'reproducer': 'python3 scripts/derive_nsc_incoming_low_energy_vacuum_enclosure.py --check'}

class _vacuum_control_CPUBudget(RuntimeError):
    pass

def _vacuum_control_run():
    before = _vacuum_control_signature()
    channel, config, spec = _vacuum_control_inputs()
    cpu = _vacuum_control_time.process_time()
    wall = _vacuum_control_time.monotonic()
    completed = 0

    def progress(count):
        nonlocal completed
        completed = count
        print(f'{count}/{_vacuum_control_STEPS} directed midpoint cells; CPU {_vacuum_control_time.process_time() - cpu:.2f}s', flush=True)

    def stop(signum, frame):
        raise _vacuum_control_CPUBudget('fixed120CPU-second budget exhausted')
    prior = _vacuum_control_signal.signal(_vacuum_control_signal.SIGPROF, stop)
    _vacuum_control_signal.setitimer(_vacuum_control_signal.ITIMER_PROF, _vacuum_control_CPU_CAP)
    try:
        payload = _vacuum_control_propagate(channel, config, progress=progress)
    except _vacuum_control_CPUBudget:
        payload = {'budget_exceeded': True, 'steps_completed_at_last_checkpoint': completed, 'channel': channel, 'source_config': config}
    finally:
        _vacuum_control_signal.setitimer(_vacuum_control_signal.ITIMER_PROF, 0)
        _vacuum_control_signal.signal(_vacuum_control_signal.SIGPROF, prior)
    if _vacuum_control_signature() != before:
        raise ValueError('proof/input source changed during production')
    payload.update(signature=before, selected_source_artifact=spec, runtime={'CPU_seconds': _vacuum_control_time.process_time() - cpu, 'wall_seconds': _vacuum_control_time.monotonic() - wall, 'CPU_cap_seconds': _vacuum_control_CPU_CAP, 'budget_exceeded': bool(payload.get('budget_exceeded'))})
    raw = _vacuum_control_json.dumps(payload, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()
    compressed = _vacuum_control_gzip.compress(raw, compresslevel=6, mtime=0)
    h = _vacuum_control_sha256(compressed).hexdigest()
    path = _vacuum_control_ROOT / f'results/development/artifacts/nsc-incoming-low-energy-vacuum-enclosure.{h}.json.gz'
    path.write_bytes(compressed)
    record = _vacuum_control_make_record(payload, path)
    (_vacuum_control_ROOT / _vacuum_control_OUTPUT).write_text(_vacuum_control_json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + '\n')
    return record

def _vacuum_control_check():
    old = _vacuum_control_json.loads((_vacuum_control_ROOT / _vacuum_control_OUTPUT).read_text())
    spec = old['payload']
    path = _vacuum_control_ROOT / spec['path']
    if _vacuum_control_digest(spec['path']) != spec['sha256']:
        raise ValueError('directed proof artifact changed')
    payload = _vacuum_control_json.loads(_vacuum_control_gzip.decompress(path.read_bytes()))
    record = _vacuum_control_make_record(payload, path)
    if record != old:
        raise ValueError('directed projector receipt differs')
    return record

from types import SimpleNamespace as _vacuum_control_Namespace
vacuum_control = _vacuum_control_Namespace(
    CPUBudget=_vacuum_control_CPUBudget,
    CPU_CAP=_vacuum_control_CPU_CAP,
    ENERGY=_vacuum_control_ENERGY,
    INPUTS=_vacuum_control_INPUTS,
    OUTPUT=_vacuum_control_OUTPUT,
    PRECISION=_vacuum_control_PRECISION,
    Path=_vacuum_control_Path,
    ROOT=_vacuum_control_ROOT,
    SOURCES=_vacuum_control_SOURCES,
    STEPS=_vacuum_control_STEPS,
    argparse=_vacuum_control_argparse,
    check=_vacuum_control_check,
    digest=_vacuum_control_digest,
    gzip=_vacuum_control_gzip,
    inputs=_vacuum_control_inputs,
    json=_vacuum_control_json,
    make_record=_vacuum_control_make_record,
    np=_vacuum_control_np,
    propagate=_vacuum_control_propagate,
    replay=_vacuum_control_replay,
    run=_vacuum_control_run,
    sha256=_vacuum_control_sha256,
    signal=_vacuum_control_signal,
    signature=_vacuum_control_signature,
    sys=_vacuum_control_sys,
    time=_vacuum_control_time,
)
