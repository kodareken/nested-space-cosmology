#!/usr/bin/env python3
"""Authenticate the new spatial-mode record and reconstruct its field checks.

Routine verification reuses the persisted mode fields. It does not repeat the
Jost solve, the historical seed generators, or any stress/metric calculation.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_pg_massive_modes import MassiveModeSection,physical_sewing
from recursive_horizons.nsc_lorentzian import geometry
from recursive_horizons.nsc_transmitting_dirac_domain import I2,S2
from recursive_horizons.nsc_paired_horizon_preparation import source_covariance

OUTPUT='results/development/nsc-pg-massive-mode-resolution.json'
SOURCES=('src/recursive_horizons/nsc_pg_massive_modes.py','src/recursive_horizons/nsc_massive_jost_modes.py',
         'scripts/derive_nsc_pg_massive_modes.py','scripts/check_nsc_pg_massive_modes.py',
         'tests/test_nsc_pg_massive_modes.py','docs/nsc-pg-massive-mode-resolution.md')
INPUTS=('results/development/nsc-compact-matched-restart.json',
        'results/development/nsc-horizon-paired-pg-map.json',
        'results/development/nsc-massive-signed-preparation.json',
        'results/development/nsc-pg-lll-preparation.json',
        'results/development/nsc-transmitting-ctp-resolvent.json',
        'results/development/nsc-mode-resolved-cauchy-state.json')

def sha(path):return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()

def field_checks(payload):
    if sha(payload['path'])!=payload['sha256'] or (ROOT/payload['path']).stat().st_size!=payload['bytes']:
        raise ValueError('mode payload authentication failed')
    with np.load(ROOT/payload['path'],allow_pickle=False) as a:
        data={k:a[k].copy() for k in a.files}
    for key,value in data.items():
        if key!='metadata_json' and not np.isfinite(value).all():
            raise ValueError('unevaluated or nonfinite mode payload field: '+key)
    meta=json.loads(data['metadata_json'].tobytes())
    for path,digest in meta['source_hashes'].items():
        if sha(path)!=digest:raise ValueError('mode source changed: '+path)
    expected_groups=set(range(1,33));observed_groups=set(data['label'][:,0].astype(int))
    if observed_groups!=expected_groups or len(data['label'])!=3584:
        raise ValueError('retained massive group/node coverage differs')
    config=json.loads((ROOT/'results/development/nsc-compact-matched-restart.json').read_text())['scattering_provenance']['config']
    maxima={k:0. for k in ('interior_current','exterior_current','source_CAR','source_law',
                           'complex_scattering_isometry','spectral_jump_relative')}
    per_group={str(i):dict(maxima) for i in range(1,33)}
    vi=np.array([S2-float(geometry(x)[0])*I2 for x in data['rho_interior']])
    ve=np.array([S2-float(geometry(x)[0])*I2 for x in data['rho_exterior']])
    for i,label in enumerate(data['label']):
        ch,_,_,mass,angular,energy,_=label
        inner,outer,_,projection=physical_sewing(energy,mass,data['reflection'][i].item(),data['transmission'][i].item())
        if not np.array_equal(projection,data['source_projector'][i]):raise ValueError('source rank changed')
        fi=data['interior'][i];fe=data['exterior'][i]
        values={
            'interior_current':float(max(np.linalg.norm(f.conj().T@v@f+inner.conj().T@inner) for f,v in zip(fi,vi))),
            'exterior_current':float(max(np.linalg.norm(f.conj().T@v@f-outer.conj().T@np.diag([1.,-1.])@outer) for f,v in zip(fe,ve))),
            'complex_scattering_isometry':float(np.linalg.norm(data['scattering'][i].conj().T@data['scattering'][i]-projection)),
        }
        eig=np.linalg.eigvalsh(data['source_covariance'][i]);values['source_CAR']=max(0.,-float(eig.min()),float(eig.max())-1.)
        values['source_law']=float(np.linalg.norm(data['source_covariance'][i]-source_covariance(energy,config['surface_gravity'],config['omega'],mass)))
        section=MassiveModeSection(float(energy),float(mass),float(angular),data['rho_interior'],fi,
                                  data['rho_exterior'],fe,data['source_covariance'][i],projection,{},0,
                                  data['interior_fundamental'][i],data['exterior_ingoing'][i])
        values['spectral_jump_relative']=section.spectral_jump_residual()['spectral_jump_relative']
        for key,value in values.items():
            maxima[key]=max(maxima[key],value)
            per_group[str(int(ch))][key]=max(per_group[str(int(ch))][key],value)
    return {'maxima':maxima,'per_group':per_group,'groups':32,'positive_energy_signed_angular_nodes':3584}

def build(evaluation):
    raw=json.loads(Path(evaluation).read_text())
    checks=field_checks(raw['payload'])
    for key,value in checks['maxima'].items():
        tolerance=3e-11 if key in ('source_CAR','source_law','complex_scattering_isometry') else 3e-9
        if value>tolerance:raise ArithmeticError(f'{key}={value} exceeds {tolerance}')
    if raw['maximum_residuals']['outer_radius_doubling_PG_fields']>3e-9:
        raise ArithmeticError('outer boundary refinement remains unresolved')
    record={
        'schema':'NSC-PG-MASSIVE-MODE-RESOLUTION-v1',
        'status':'PASS: massive two-chart mode and whole-line Green construction; spatial covariance integration remains OPEN',
        'source_hashes':{p:sha(p) for p in SOURCES},
        'input_hashes':{**{p:sha(p) for p in INPUTS},**{p:v for p,v in raw['source_hashes'].items() if p not in SOURCES}},
        'locked_inputs':raw['locked_inputs'],'payload':raw['payload'],
        'domain':{
            'operator':'same whole-line self-adjoint PG Dirac H_D on L2(d rho;C2), with the owned angular/compact sectors',
            'neck':'fixed rho=0; T transmission retained',
            'horizon':'independent outgoing exterior/interior coordinates, shared ingoing characteristic',
            'source_rank':'3 for |E|>m; 2 for |E|<m; occupation threshold Omega*m does not set rank',
            'source_covariance':'owned C_H plus inherited incoming occupation; seed covariance is a comparison only',
            'angular_Z':'zero between angular signs; correlations inside C_H retained',
            'negative_energy':'locked paired sigma3 K map, with the opposite angular positive partner',
        },
        'construction':{
            'start_convention':'matched_delta_q',
            'source_generations':'historical and matched seed payloads preserved byte for byte',
            'exterior':'complex Jost amplitude; massive spinor asymptotic expansion; phase evolution for the closed-current branch',
            'incoming_phase':'sqrt(T) at the horizon; complex t_out evaluated in the declared outer Jost coordinate',
            'clock':'exact PG clock F prime=beta/A, F(0)=0; geometry primitives only, no LLL mode reuse',
            'basis_connection':'included in both PG spin-frame conversions',
            'normalization':'whole-line four-support Green jump; dP/dE=Phi_E Phi_E dagger/(2*pi)',
            'sampled_frequency_content':'all existing positive-energy nodes and required angular signs; finite arrays are not a CAR grid',
        },
        'groups':raw['groups'],'measured_residuals':raw['maximum_residuals'],
        'outer_refinements':raw['refinements'],'stored_field_verification':checks,
        'tolerances':{'exact_frame_source_and_scattering_algebra_absolute':3e-11,
                      'spatial_Green_relative':3e-9,'spatial_current_and_outer_fields_absolute':3e-9,
                      'Green_scale':'max(1, Frobenius norm of Phi(x) Phi(y) dagger)',
                      'rationale':'spatial propagation and finite Jost boundary errors are measured separately from exact fiber algebra'},
        'continuum_scope':{
            'no_point_spectrum':'current/L2 argument for the owned whole-line domain, including E=0 and thresholds',
            'Jost_completion':'exact PG clock, Coulomb phase, and radial angular spinor correction; integrable remainder away from thresholds',
            'spectral_identification':'whole-line Green and Stone route stated in the note; finite samples alone do not establish completeness',
            'contact_and_tail':'must be retained explicitly in spatial packet covariance; not supplied by a finite seed-frequency Gram fill',
        },
        'gate':{
            'massive_spatial_mode_map':'PASS',
            'whole_line_Green_support_and_normalization':'PASS',
            'spatial_C_PG_and_full_C1b':'OPEN: packet energy integration and contact/tail error account',
            'C2_transmitting_EndpointBranchJets':'OPEN',
            'C3_Gamma_rest_Weyl':'OPEN',
            'C4_extended_stationarity':'OPEN',
            'PDF_bump':False,'metric_timestep':False,'history_selected':False,
            'stress_computed':False,'scales_refitted':False,'seed_covariance_used_as_input':False,
        },
        'comparison':{'policy':'authenticate every source/input/payload and record field; reconstruct stored-field numerical residuals without repeating Jost or historical generators',
                      'stored_field_atol':3e-13,'stored_field_rtol':3e-7,
                      'fresh_mode_reproduction':'derive_nsc_pg_massive_modes.py produces the mode payload and all measured residuals from the same named owners'},
    }
    (ROOT/OUTPUT).write_text(json.dumps(record,indent=2,sort_keys=True)+'\n')
    return record

def check():
    record=json.loads((ROOT/OUTPUT).read_text())
    for field in ('source_hashes','input_hashes'):
        for path,digest in record[field].items():
            if sha(path)!=digest:raise ValueError('authenticated dependency changed: '+path)
    for flag in ('PDF_bump','metric_timestep','history_selected','stress_computed','scales_refitted','seed_covariance_used_as_input'):
        if record['gate'][flag] is not False:raise ValueError('mode record promoted beyond its scope: '+flag)
    locked=json.loads((ROOT/INPUTS[0]).read_text())['locked_inputs']
    if record['locked_inputs']!=locked:raise ValueError('locked action or scale data changed')
    actual=field_checks(record['payload'])
    def compare(a,b,path=''):
        if isinstance(a,dict):
            if set(a)!=set(b):raise ValueError('record keys differ: '+path)
            for key in a:compare(a[key],b[key],path+'/'+key)
        elif isinstance(a,(int,float)):
            if not np.isclose(a,b,atol=3e-13,rtol=3e-7):raise ValueError(f'residual differs: {path} {a} {b}')
        elif a!=b:raise ValueError('value differs: '+path)
    compare(actual,record['stored_field_verification'])
    return record

def main():
    p=argparse.ArgumentParser();p.add_argument('--build-record');p.add_argument('--check',action='store_true');a=p.parse_args()
    record=build(a.build_record) if a.build_record else check()
    print(json.dumps({'status':record['status'],'residuals':record['stored_field_verification']['maxima'],
                      'full_C1b':record['gate']['spatial_C_PG_and_full_C1b'],'PDF_bumped':False},indent=2))

if __name__=='__main__':main()
