"""An incomplete or unconverged calculation cannot become endpoint evidence."""
import copy

import pytest

from analyze_nsc_atlas_fixed_endpoint import compare, GAPS
from recursive_horizons.nsc_atlas_endpoint_comparison import ENDPOINT_SOURCE_X


def records():
    out=[]
    for i,(radius,gap) in enumerate(zip([.54,.545,.546],[1.,.5,.25])):
        out.append(dict(schema='nsc-atlas-fixed-endpoint-v1',level=i,
            status='completed_numerical_endpoint',source_commit=str(i)*40,
            radial_accuracy=dict(rtol=2.5e-14,atol=1e-15,adjoint_order=32,gradient_absolute_tolerance=2e-8),
            execution=dict(last_accepted_x=ENDPOINT_SOURCE_X,raw_history_saved=False,final_state_saved=False),
            geometry=dict(x_left=.4,x_right=ENDPOINT_SOURCE_X,reference_radius_endpoint=.55,
                          cap_cells=[2,4,8][i],exterior_cells=[32,64,128][i]),
            input=dict(sha256='a'*64,alpha=.0001,source_rows=257),
            preparation=dict(input_sha256='b'*64,original_input_sha256='c'*64),
            endpoint=dict(einstein_radius=radius,jordan_radius=radius,c=0.,frame_formula_needs_check=False),
            samples=[dict(x=.4+(ENDPOINT_SOURCE_X-.4)*j/8) for j in range(9)],
            maximum_absolute={key:gap for key in GAPS.values()},
            transport_comparison={key:gap for key in ('mass_absolute_gap','scalar_absolute_gap',
                'max_abs_running_mass_gap','max_abs_running_scalar_gap')}))
    return out


def test_three_completed_runs_produce_a_scaled_resolution_envelope():
    data=records();before=copy.deepcopy(data)
    result=compare(data[::-1])
    assert data==before
    assert result['requested_mass_and_flux_contraction'] is True
    assert result['physical_observable']['value_metres']==pytest.approx(5460.)
    assert result['physical_observable']['resolution_half_width_metres']==pytest.approx(60.)
    assert result['physical_evolutions']==0
    assert result['independent_observational_test'] is False


@pytest.mark.parametrize('mutation',[
    lambda r:r.update(status='stopped'),
    lambda r:r['execution'].update(last_accepted_x=.45),
    lambda r:r['input'].update(alpha=.0002),
    lambda r:r['preparation'].update(input_sha256='d'*64),
    lambda r:r['endpoint'].update(jordan_radius=.6),
    lambda r:r['endpoint'].update(frame_formula_needs_check=True),
    lambda r:r['samples'].pop(),
    lambda r:r['radial_accuracy'].update(rtol=2e-13),
])
def test_wrong_or_incomplete_evidence_is_rejected(mutation):
    data=records();mutation(data[2])
    with pytest.raises(ValueError):compare(data)


def test_nonshrinking_gaps_are_reported_without_a_success_claim():
    data=records()
    data[2]['maximum_absolute'][GAPS['wall_flux']]=.7
    result=compare(data)
    assert result['requested_mass_and_flux_contraction'] is False
    assert result['prediction_status']=='unconverged_numerical_estimate'


def test_three_receipts_must_share_the_confirmed_panel_step_rule():
    data=records()
    data[0]['radial_accuracy']['maximum_nonzero_panel_step_fraction']=.25
    with pytest.raises(ValueError,match='settings differ'):
        compare(data)
    for row in data:
        row['radial_accuracy']['maximum_nonzero_panel_step_fraction']=.25
    assert compare(data)['radial_accuracy']['maximum_nonzero_panel_step_fraction']==.25
    data[2]['radial_accuracy']['maximum_nonzero_panel_step_fraction']=.5
    with pytest.raises(ValueError,match='confirmed radial accuracy'):
        compare(data)


def test_panel_polynomial_rule_must_be_shared():
    data=records()
    data[0]['radial_accuracy']['background_panel_sides']=True
    with pytest.raises(ValueError,match='settings differ'):compare(data)
    for row in data:row['radial_accuracy']['background_panel_sides']=True
    assert compare(data)['radial_accuracy']['background_panel_sides'] is True
    data[1]['radial_accuracy']['background_panel_sides']=False
    with pytest.raises(ValueError,match='confirmed radial accuracy'):compare(data)


def test_known_final_serializer_failure_preserves_endpoint_and_marks_missing_ledger():
    data=records();coarse=data[0]
    coarse.update(status='stopped',source_commit='ae2675980de38f547efa6b4fb1afbd2d12eee62e',
        failure=dict(phase='endpoint summary',
            reason="dict() got multiple values for keyword argument 'state_corrected'"),
        schedule=dict(accepted_step_budget=1024),integrated_defect={})
    coarse['execution']['accepted_steps']=1024
    coarse['samples'][-1].update(coarse['endpoint'])
    del coarse['transport_comparison']
    before=copy.deepcopy(data)
    result=compare(data)
    assert data==before
    assert result['all_endpoints_completed'] is True
    assert result['summary_coverage'][0]['original_receipt_status']=='stopped'
    assert result['accumulated_transport_gaps']['mass_endpoint']['values']==[None,.5,.25]
    assert result['accumulated_transport_gaps']['mass_endpoint']['strictly_decreasing'] is None
    assert result['requested_mass_and_flux_contraction'] is True
    coarse['execution']['accepted_steps']=1023
    with pytest.raises(ValueError,match='complete'):compare(data)
    coarse['execution']['accepted_steps']=1024
    coarse['failure']['reason']='a physical constraint failed'
    with pytest.raises(ValueError,match='complete'):compare(data)
