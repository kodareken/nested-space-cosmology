from check_nsc_atlas_matched_trace import calculate


def test_matching_the_strong_characteristic_stays_matched_through_adapter():
    row=calculate()
    assert row['execution']['time_steps']==0
    assert abs(row['direct']['incoming_discrepancy'])<2e-12
    assert abs(row['hybrid']['incoming_discrepancy'])<2e-12
    assert abs(row['hybrid']['scalar_channel_identity_defect'])<2e-12
    assert row['physical_input']['alpha']>0
