import json
from recursive_horizons.nsc_atlas_causal_observers import causal_observers_record, centre_emission_to_wall
r=causal_observers_record()
print(json.dumps(r['named_events'],indent=2))
print(json.dumps(r['independent_checks'],indent=2))
print('domain',r['trajectory']['t_flat_over_ell_domain'])
w=r['named_events']['W_early']
print('centre_to_early',centre_emission_to_wall(.5,w['t_flat_over_ell']))
print('nodes',r['admission']['quadrature_evaluations'])
