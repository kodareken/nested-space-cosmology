import json
from recursive_horizons.nsc_atlas_connection_slice import connection_slice_record
r=connection_slice_record()
print(json.dumps(r['wall'], indent=2))
print(json.dumps(r['independent_checks'], indent=2))
print('transition_length', r['transition']['proper_length_over_ell'])
print('throat_l', r['bridge']['throat_l_over_ell'])
print('outer_l', r['bridge']['samples']['l_over_ell'][-1])
