#!/usr/bin/env python3
"""Compute the transmitted packet resolvent and its static metric jets."""
from __future__ import annotations
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_transmitting_resolvent import (
    TransmittingCrossResolvent, quadrature, probes, profile, inverse_response, FIELDS,
)
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
from derive_nsc_transmitting_boundary_binding import compare

OUTPUT='results/development/nsc-transmitting-cross-resolvent.json'
INPUTS={
 'split':'results/development/nsc-common-time-bulk-split.json',
 'state':'results/development/nsc-mode-resolved-cauchy-state.json',
 'T':'results/development/nsc-transmitting-dirac-domain.json',
 'energies':'results/nsc-8-chiral-boundary.json',
 'differential':'results/development/nsc-transmitting-ctp-variation.json',
}
SOURCES=('src/recursive_horizons/nsc_transmitting_resolvent.py','scripts/derive_nsc_transmitting_resolvent.py',
         'tests/test_nsc_transmitting_resolvent.py','docs/nsc-transmitting-cross-resolvent.md')
DEPENDENCIES=('src/recursive_horizons/nsc_lorentzian.py','src/recursive_horizons/nsc_common_time_bulk_split.py',
 'src/recursive_horizons/nsc_transmitting_dirac_domain.py','src/recursive_horizons/nsc_mode_resolved_cauchy_state.py',
 'scripts/derive_nsc_transmitting_boundary_binding.py')


def sha(path):return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
def norm(a):return float(np.linalg.norm(a))
def rel(a,reference):return norm(a)/max(1.,norm(reference))
def pairs(a):
 a=np.asarray(a);return np.stack((a.real,a.imag),axis=-1).tolist()


def calculate():
 records={k:json.loads((ROOT/p).read_text()) for k,p in INPUTS.items()}
 split,state=records['split'],records['state']
 for k in ('A','magnetic_flux','Omega','zeta','V_full'):
  if split['locked_inputs'][k]!=state['locked_inputs'][k]:raise ValueError('locked input changed: '+k)
 old_energies=records['energies']['energy_resolved_maps']
 z=complex(old_energies[0]['energy']['real'],old_energies[0]['energy']['imag'])
 w=complex(old_energies[1]['energy']['real'],old_energies[1]['energy']['imag'])
 if records['energies']['conventions']['ports']!=[1.,0.,-1.]:raise ValueError('inherited ports changed')
 x,weights=quadrature(96)
 gram=sum(weight*probes(rho).conj().T@probes(rho) for rho,weight in zip(x,weights))
 packet_error=float(np.max(abs(gram-np.eye(4))))
 response_arrays={key:[] for key in ('upper_response','lower_response','second_energy_response','upper_metric_jets',
  'inverse_response','inverse_response_metric_jets','upper_fields','lower_fields','upper_field_derivatives','lower_field_derivatives')}
 rows=[];fd_rows=[]
 representatives=[0,max(c['index'] for c in state['channels'] if c['family']==1),max(c['index'] for c in state['channels'] if c['family']==2)]
 for channel in state['channels']:
  owner=TransmittingCrossResolvent(channel['compact_mass'],channel['angular_eigenvalue'])
  upper=owner.solve(z);lower=owner.solve(z.conjugate());other=owner.solve(w,with_jets=False)
  u=np.array([upper.field(rho) for rho in x]);l=np.array([lower.field(rho) for rho in x]);ow=np.array([other.field(rho) for rho in x])
  up=np.array([upper.field_derivative(rho) for rho in x]);lp=np.array([lower.field_derivative(rho) for rho in x])
  true_product=np.einsum('n,nai,naj->ij',weights,l.conj(),ow)
  identity_error=rel(upper.response-other.response-(z-w)*true_product,upper.response-other.response)
  density=np.array([owner.static_metric_kernel(upper,lower,rho) for rho in x])
  s=np.array([profile(rho)[0] for rho in x])
  weak_jets=np.einsum('n,n,naij->aij',weights,s,density)
  k,dk=inverse_response(upper);kl,dkl=inverse_response(lower);kw,_=inverse_response(other)
  weak_dk=np.array([-k@d@k for d in weak_jets])
  im=(upper.response-upper.response.conj().T)/(2j)
  eig=np.linalg.eigvalsh(im)
  v0=owner.metric_parts(0.,np.zeros(4))[0]
  seam_value=upper.field(0.)
  residuals={
   'adjoint_response':rel(lower.response-upper.response.conj().T,upper.response),
   'adjoint_metric_jets':rel(lower.directional_jets-upper.directional_jets.swapaxes(-1,-2).conj(),upper.directional_jets),
   'true_resolvent_identity':identity_error,
   'metric_weak_identity':rel(upper.directional_jets-weak_jets,upper.directional_jets),
   'inverse_metric_weak_identity':rel(dk-weak_dk,dk),
   'inverse_response':norm(k@upper.response-np.eye(4)),
   'cross_generator_reconstruction':norm(upper.response[2:,:2]+np.linalg.solve(k[2:,2:],k[2:,:2])@np.linalg.inv(k[:2,:2])),
   'inverse_response_adjoint':rel(kl-k.conj().T,k),
   'seam_continuity':max(upper.seam_residual,lower.seam_residual,other.seam_residual),
   'reverse_upper_response':norm(upper.response[:2,2:]),
   'reverse_upper_metric_jets':norm(upper.directional_jets[:,:2,2:]),
   'positive_imaginary_violation':max(0.,-float(eig.min())),
   'resolvent_norm_bound_violation':max(0.,float(np.linalg.norm(upper.response,2))-1/z.imag),
  }
  rows.append({'channel':channel,'residuals':residuals,
    'forward_parent_to_child_response':pairs(upper.response[2:,:2]),
    'forward_inverse_response_kernel':pairs(k[2:,:2]),
    'forward_response_norm':norm(upper.response[2:,:2]),
    'forward_kernel_norm':norm(k[2:,:2]),
    'minimum_imaginary_eigenvalue':float(eig.min()),
    'individual_projected_delta_coefficient_norm':norm(1j*v0@seam_value),
    'closed_four_mode_H_test_difference':norm((k+z*np.eye(4))-(kw+w*np.eye(4))),
    'function_evaluations':upper.function_evaluations+lower.function_evaluations+other.function_evaluations})
  for key,value in (('upper_response',upper.response),('lower_response',lower.response),('second_energy_response',other.response),
    ('upper_metric_jets',upper.directional_jets),('inverse_response',k),('inverse_response_metric_jets',dk),
    ('upper_fields',u),('lower_fields',l),('upper_field_derivatives',up),('lower_field_derivatives',lp)):
   response_arrays[key].append(value)
  if channel['index'] in representatives:
   epsilon=2e-5
   fd_owner=TransmittingCrossResolvent(channel['compact_mass'],channel['angular_eigenvalue'],rtol=2e-13,atol=2e-15)
   for field,name in enumerate(FIELDS):
    parameters=np.zeros(4);parameters[field]=epsilon
    plus=fd_owner.solve(z,parameters=parameters,with_jets=False)
    minus=fd_owner.solve(z,parameters=-parameters,with_jets=False)
    dr=(plus.response-minus.response)/(2*epsilon)
    dkernel=(inverse_response(plus)[0]-inverse_response(minus)[0])/(2*epsilon)
    fd_rows.append({'channel':channel['index'],'field':name,'epsilon':epsilon,
      'response_error':rel(dr-upper.directional_jets[field],upper.directional_jets[field]),
      'inverse_response_error':rel(dkernel-dk[field],dk[field])})
 maxima={name:max(row['residuals'][name] for row in rows) for name in rows[0]['residuals']}
 fdmax=max(max(row['response_error'],row['inverse_response_error']) for row in fd_rows)
 if packet_error>3e-11 or max(maxima.values())>3e-11 or fdmax>3e-8:
  raise ArithmeticError(json.dumps({'packet':packet_error,'residuals':maxima,'finite_difference':fdmax}))
 if not all(row['forward_response_norm']>1e-12 and row['individual_projected_delta_coefficient_norm']>1e-12 for row in rows):
  raise ArithmeticError('no resolved forward transmission or projected delta witness')
 payload_arrays={k:np.array(v) for k,v in response_arrays.items()}
 payload_arrays.update(rho_quadrature=x,rho_weights=weights,probe_gram=gram,
                       channel_indices=np.array([c['index'] for c in state['channels']],dtype=np.int32))
 data=deterministic_npz_bytes(payload_arrays);digest=hashlib.sha256(data).hexdigest()
 return {
  'schema':'NSC-TRANSMITTING-CROSS-RESOLVENT-v1',
  'status':'Z2a PASS: domain-correct packet cross-response and inverse-response generator; Z2b static metric jets PASS; time-history/CTP source bridge OPEN',
  'source_hashes':{p:sha(p) for p in SOURCES},'input_hashes':{p:sha(p) for p in (*INPUTS.values(),*DEPENDENCIES)},
  'locked_inputs':split['locked_inputs'],
  'payload':{'path':f'results/development/artifacts/nsc-transmitting-cross-resolvent.{digest}.npz','sha256':digest,'bytes':len(data),
             'description':'packet resolvents, inverse-response generators, four compact static metric jets and uncompressed quadrature fields'},
  'definition':{
   'route':'B: whole-line boundary-resolvent, then normalized packet compression',
   'whole_line_response':'R(z)=(H_D-z)^-1; projectors applied after inversion',
   'packet_response':'R_packet=J^dagger R J, J=[f_parent I2,f_child I2]; order parent,child',
   'probes':'fp=sqrt(30)*rho*(1-rho) on (0,1); fc=fp(-rho) on (-1,0); zero elsewhere',
   'probe_roles':'normalized L2 source/measurement functions in D(H), not a physical initial state or new particle sector',
   'kernel_generator':'K(z)=R_packet(z)^-1; K_cp is an energy-dependent inverse-response cross kernel',
   'ordinary_Hamiltonian_B_claimed':False,
   'not_a_closed_four_mode_resolvent':True,
   'causal_prescription':'Im z>0: forward support toward decreasing rho, zero above source; Im z<0: reverse support, zero below source',
   'convention':'(H-z)^-1=i integral_0^infinity exp(i z t)exp(-iHt)dt in upper half-plane; not directly the -i theta convention without its sign conversion',
   'domain':'actual whole-line transmitted PG Dirac domain; computations restricted to strictly trapped [-1,1], using causal support, not walls',
   'not_evaluated':'numeric response from the asymptotically exterior parent region or through a singular stationary horizon chart',
  },
  'parameters':{'z':[z.real,z.imag],'identity_check_energy':[w.real,w.imag],'rtol':2e-12,'atol':2e-14,
                'finite_difference_rtol':2e-13,'finite_difference_atol':2e-15,'quadrature_points_per_half':96},
  'delta_accounting':{
   'locked_principal_singular_values':[split['Z2_domain_check']['coefficient_smallest_singular_value'],split['Z2_domain_check']['coefficient_largest_singular_value']],
   'individual_terms':'i v0 u(0) delta from P_child and its opposite from P_parent are nonzero',
   'cancellation':'u_parent(0)=u_child(0), including metric sensitivities; they cancel before projection of the full solution',
   'delta_renamed_as_B_or_source':False,
  },
  'metric_derivative':{
   'fields':list(FIELDS),'cut':'fixed rho=0','canonical_probe_functions_held_fixed':True,
   'functional':'delta R_packet=sum_A integral s_A(rho) J_A(rho;z)d rho for smooth compact variations in (-1,1)',
   'kernel':'J_A=i/2[L^dagger V_A U_prime-L_prime^dagger V_A U]-L^dagger M_A U; U=R(z)J, L=R(z*)J',
   'inverse_kernel':'delta K=-K delta R_packet K',
   'numeric_direction':'s=exp(1-1/(1-rho²)) inside |rho|<1, zero outside; no metric history or fitted deformation',
   'static_metric_jets_evaluated':True,'EndpointBranchJets_filled':False,
   'CTP_bridge':'requires the time/spectral reconstruction of the full history response and its state/memory data; a static packet resolvent derivative is not a unitary dU',
  },
  'verification':{'packet_gram_error':packet_error,'maximum_residuals':maxima,'tolerance':3e-11,
                  'finite_difference_rows':fd_rows,'finite_difference_maximum':fdmax,'finite_difference_tolerance':3e-8,
                  'finite_difference_tolerance_reason':'central finite difference at epsilon=2e-5 combines O(epsilon²) truncation with amplified ODE/inversion roundoff',
                  'resolvent_identity_uses_uncompressed_solutions':True,
                  'Z1_projector_and_CAR_prerequisites_reused':split['finite_verification']['maximum_residuals'],
                  'T_current_compatibility_reused':split['T_compatibility']['current_relation_residual']},
  'channels':rows,
  'scope':{'Z2a':'PASS on declared response/probe domain','Z2b':'PASS for compact static metric kernels; OPEN for time-history CTP jets',
           'Z3':'OUT OF SCOPE: fixed cut','Z4':'OPEN','Gamma_rest_boundary_derivative':None,'physical_Weyl_mismatch':None,
           'physical_Vc_or_history_selected':False,'new_physical_terms':False,'parameters_refitted':False,'old_generators_rerun':False,
           'physical_stress_or_nulls_fabricated':False,'metric_timestep_started':False,'optimizer_started':False},
  'comparison':records['T']['comparison'],
 },data


def main():
 p=argparse.ArgumentParser(description=__doc__);g=p.add_mutually_exclusive_group();g.add_argument('--check',action='store_true');g.add_argument('--write',action='store_true');a=p.parse_args()
 record,data=calculate();out=ROOT/OUTPUT;artifact=ROOT/record['payload']['path']
 if a.check:
  compare(json.loads(out.read_text()),record)
  if artifact.read_bytes()!=data:raise AssertionError('resolvent artifact differs')
  print('Transmitting cross-resolvent and static metric jets verified')
 elif a.write:
  if out.exists() or artifact.exists():raise FileExistsError('refusing overwrite')
  artifact.write_bytes(data);out.write_text(json.dumps(record,indent=2,sort_keys=True,allow_nan=False)+'\n');print(record['status'])
 else:print(json.dumps(record,indent=2,sort_keys=True,allow_nan=False))


if __name__=='__main__':main()
