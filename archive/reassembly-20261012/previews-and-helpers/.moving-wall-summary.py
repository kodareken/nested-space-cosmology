from recursive_horizons.nsc_atlas_moving_boundaries import moving_boundary_record
import json
r = moving_boundary_record()
print(json.dumps(r['independent_checks'], indent=2))
print(r['vacuum_wall_control']['domain'])
print({k: [v[0], v[-1]] for k, v in r['vacuum_wall_control']['samples'].items()})
