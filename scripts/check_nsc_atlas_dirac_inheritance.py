"""Print the analytic quiet-cap inheritance scales; no grids or data files."""
import json
import math

from recursive_horizons.nsc_atlas_dirac_preparation import (
    cap_dirac_detector_shape,cap_inheritance_scales,
)
from recursive_horizons.nsc_atlas_nested_spacetime import expanding_cap_wall_clock
from recursive_horizons.nsc_atlas_dirac_trace import inherited_late_cap_vacuum


def main():
    print(json.dumps(dict(
        quiet_cap=dict(H_ell=1.,M_parent_over_ell=.4,k_wall_ell=2.,compact_cap_mass=0.),
        horizon_wall_clock=expanding_cap_wall_clock(
            radius=1.,H=1.,beta_minus=.95,Rdot=.95),
        after_one_expansion_efold=dict(two_point_difference_weight=math.exp(-3),
            stress_difference_weight=math.exp(-4),leading_stress_error_weight=math.exp(-5)),
        unit_gap_spectrum_shape=float(cap_dirac_detector_shape(1.)),
        conditional_physical_scale=cap_inheritance_scales(10000.),
        conditional_late_vacuum_source=inherited_late_cap_vacuum(ell_m=10000.),
        scope=('Full massless 4D Dirac local late limit for regular inherited states after '
            'cap horizon isolation. The length is supplied; no finite parent covariance, '
            'source-dependent Psi/Phi or generated classical cap Lambda is calculated.')),
        indent=2,allow_nan=False))


if __name__=='__main__':
    main()
