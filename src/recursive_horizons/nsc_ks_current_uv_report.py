"""UV v3 report separated from the immutable v2 proof owner."""
from .nsc_ks_finite_history_uv_remainder import (
    CAUCHY_RECORD,
    CUTOFF_BRIDGE_CONTROL,
    FIRST_NONCANCELLING_PAIRED_ORDER,
    HISTORY_RECORD,
    INVENTORY_RECORD,
    LocalIncomingFamily,
    MISSING_L0_A4,
    MISSING_UPSTREAM_TAIL,
    MappingProxyType,
    PILOT_RECORD,
    PILOT_SCHEMA,
    PREPARATION_OWNER,
    PROVISIONAL_UV_ALLOCATION,
    Path,
    QUANTITY_NAMES,
    REMAINDER_ORDER,
    RESTART_RECORD,
    STATE_LAW_OWNER,
    SUBTRACTION_OWNER,
    SUCCESSOR_SCHEMA,
    V2_RECORD,
    V3_RECORD,
    V4_RECORD,
    _as_map,
    _digest,
    _quantity,
    _quantity_binding_slice,
    _require_missing,
    _root,
    authenticate_changed_history_uv_bindings,
    changed_history_uv_bindings,
    conditional_m4_remainder_majorant,
    current_history_commutator_integrals,
    fermi_thermal_tail_bound,
    first_noncancelling_paired_coefficient,
    frozen_coefficient_records,
    json,
    manufactured_defect_majorant_control,
    negative_pair_and_omission_controls,
    np,
    remainder_recurrence_identities,
    smallest_representative_pair,
    thermal_tail_direct_window_check,
    truncated_envelope_L0_defect,
)

from .nsc_ks_current_uv_transport import (
    C4_SCOPE,
    NEXT_PRIMITIVE as TRANSPORT_NEXT_PRIMITIVE,
    current_history_characteristic_transport_report,
    validate_characteristic_transports,
)

V3_SUCCESSOR_SCHEMA = 'NSC-KS-CURRENT-UV-REMAINDER-v3'

V2_REMAINDER_RECORD = 'results/development/nsc-ks-current-uv-remainder-v2.json'

V3_SUCCESSOR_STATUS = (
    'OPEN: current-history characteristic transports of major A2, A3 and n4 '
    'owned under the switched state law; history-minus-reference e^{-3} C4 '
    'defined; massless history-minus-reference Im(major A2) is exactly 0; '
    'Re(major A2) difference and q have proved interval majorants; numerical '
    'C4, C_M and the vacuum tail remain missing pending minor A3 / L0 A2 '
    'second jets, production L0 A4 H2 integrals and the same-column upstream '
    'H2 remainder')

def changed_history_uv_v3_quantities(root=None, bindings=None, transport=None):
    """v3 quantities: same names as v2, with proved transport fields filled in."""
    root = Path(root) if root is not None else _root()
    bindings = (
        changed_history_uv_bindings(root) if bindings is None
        else authenticate_changed_history_uv_bindings(bindings, root))
    transport = (
        current_history_characteristic_transport_report(root)
        if transport is None else transport)
    validate_characteristic_transports(transport, root)
    slice_ = _quantity_binding_slice(bindings)
    identities = remainder_recurrence_identities()
    first = first_noncancelling_paired_coefficient()
    pair = smallest_representative_pair(root)
    family = LocalIncomingFamily(np.array(
        json.loads((root / HISTORY_RECORD).read_text())['history']['coefficients']))
    commutators = current_history_commutator_integrals(
        family, pair['absolute_angular'])
    majorant = conditional_m4_remainder_majorant(
        None, None, commutators['Bz_integral'], commutators['Bzz_integral'])
    imag = transport['imag_major_A2']
    majorants = transport['majorants']
    c4 = transport['c4']
    formulas = dict(first['formulas'])
    formulas['C4'] = c4['formulas']['C4']
    quantities = {
        'A2': _quantity(
            'A2',
            formula=(
                'Im(T_s major A2)=-m ell/4 (s a^2 r_rho + r_z)/r^2; '
                'Re(T_s major A2)=[-a_rho m^2 a r^3 + ell^2 (-a_rho a r + '
                'r_rho a^2 + s r_z + 2 q r) + 2 m^2 q r^3]/(4 r^3). '
                'History-minus-reference. Not a local metric jet'),
            derivation=(
                'nsc_ks_current_uv_transport characteristic_transport_identities '
                'and history_minus_reference_imag_major_A2'),
            numerical_value=None,
            status='OPEN',
            missing_primitive=(
                'full-envelope major A2, including the unowned upstream '
                'imaginary datum; the massless difference is closed'),
            bindings=slice_,
            imag_major_transport_rhs_this_massless_pair=imag[
                'transport_rhs_this_massless_pair'],
            imag_major_transport_rhs_vanishes_by_m_0_identity=pair['mass'] == 0.0,
            imag_major_A2_numerical_value=None,
            upstream_imag_major_A2_owned=False,
            real_major_A2_numerical_value=None,
            minor_A2_numerical_value=None,
            history_minus_reference_imag_major_A2=imag['numerical_value'],
            history_minus_reference_imag_major_A2_proved_exact=imag['proved_exact'],
            history_minus_reference_real_major_A2_abs_upper=majorants[
                'history_minus_reference_real_major_A2_abs_upper']['binary64_upper'],
            q_abs_upper=majorants['q_abs_upper']['binary64_upper'],
            delta_q_abs_upper=majorants['delta_q_abs_upper']['binary64_upper'],
            A3_locally_algebraic_in_metric_jets=False,
            coefficient_scope='history-minus-reference',
        ),
        'A3': _quantity(
            'A3',
            formula=(
                'minor A3 = a^2/(2 i) Pi_{-s} L0 A2 on transported major A2; '
                'massless Re(T_s major A3)=ell^2 h/(2 r^2); Im(T_s major A3) '
                'retains Re(major A2). Not a local metric jet'),
            derivation=(
                'nsc_ks_current_uv_transport Re_Ts_major_A3_massless_is_ell2_h_over_2r2'),
            numerical_value=None,
            status='OPEN',
            missing_primitive=TRANSPORT_NEXT_PRIMITIVE,
            bindings=slice_,
            depends_on_major_A2=True,
            A3_locally_algebraic_in_metric_jets=False,
            massless_Re_Ts_depends_on_unowned_h=True,
        ),
        'n4': _quantity(
            'n4',
            formula=formulas['n4_transport'],
            derivation=(
                'first_noncancelling continuity_e_minus4_retains_minor_A3; '
                'massless |minor A2|^2 is independent of the unowned h'),
            numerical_value=None,
            status='OPEN',
            missing_primitive='minor A3 and the n4 characteristic integral',
            bindings=slice_,
            retains_minor_A3=True,
            n4_eliminated_without_minor_A3=False,
        ),
        'C4': _quantity(
            'C4',
            formula=formulas['C4'],
            derivation=(
                'production history-minus-reference e^{-3} contraction on I; '
                'massless J_3 drops the unowned imag major A2. The E^{-3} first '
                'allowed order is not this coefficient'),
            numerical_value=None,
            status='OPEN',
            missing_primitive=c4['missing_primitive'],
            bindings=slice_,
            first_noncancelling_paired_inverse_energy_order=first[
                'first_noncancelling_paired_inverse_energy_order'],
            complete_e_minus3_coefficient_proven_nonzero=False,
            C4_or_C_M_bounded=False,
            history_minus_reference_coefficient_formed=True,
            full_envelope_coefficient_formed=False,
            definition_scope=C4_SCOPE,
            massless_J3_independent_of_upstream_imag_major_A2=True,
        ),
        'C_M': _quantity(
            'C_M',
            formula=(
                '||X-X^{(4)}||_H2 <= C_M e^{-4} after propagate_h2 of the '
                'upstream H2 remainder, current-history L0 A4 H2 integrals, '
                'and owned B_z, B_zz'),
            derivation=(
                'conditional_m4_remainder_majorant on the declared slab; '
                'manufactured H2 fixtures are not production L0 A4 or upstream '
                'H2 inputs'),
            numerical_value=None,
            status='OPEN',
            missing_primitive=tuple(majorant['missing_inputs']),
            bindings=slice_,
            C_M_h2=None,
            Bz_integral=commutators['Bz_integral'],
            Bzz_integral=commutators['Bzz_integral'],
            certified_L0_A4=False,
            manufactured_fixture_used_as_production_constant=False,
            applicability=(
                'current-history slab of 0b0e4ced with owned radius commutators; '
                'C_M remains unbounded until production L0 A4 and the '
                'same-column upstream H2 remainder exist'),
        ),
        'current_history_L0_A4_H2_integrals': _quantity(
            'current_history_L0_A4_H2_integrals',
            formula='H2 integrals of L0 A4 on the declared current-history slab',
            derivation=(
                'truncated_envelope_L0_defect: formal defect is e^{-4} L0 A4; '
                'production A4 is absent. Manufactured periodic/small-matrix '
                'controls are not this integral'),
            numerical_value=None,
            status='OPEN',
            missing_primitive=MISSING_L0_A4,
            bindings=slice_,
            M=REMAINDER_ORDER,
            production_L0_A4_control=False,
        ),
        'upstream_higher_order_H2_remainder': _quantity(
            'upstream_higher_order_H2_remainder',
            formula=(
                'H2 bound on e^4 (X-X^{(4)}) at rho_up for inverse-energy '
                'orders greater than 4, same mathematical vacuum column'),
            derivation=(
                'conditional_m4_remainder_majorant initial_h2; the frozen '
                'background tail certificate is not a proved bridge'),
            numerical_value=None,
            status='OPEN',
            missing_primitive=MISSING_UPSTREAM_TAIL,
            bindings=slice_,
            frozen_background_tail_copied=False,
            same_column_as_current_history_L0_A4=True,
        ),
        'transported_majors': _quantity(
            'transported_majors',
            formula=(
                'characteristic integrals of major A2 (real and imaginary) and '
                'the subsequent major-A3 transport; only the local '
                'Im(T_s major A2) right-hand side vanishes for m=0'),
            derivation=(
                'nsc_ks_current_uv_transport owns the T_s equations; the '
                'massless difference of Im(major A2) is exact 0'),
            numerical_value=None,
            status='OPEN',
            missing_primitive=TRANSPORT_NEXT_PRIMITIVE,
            bindings=slice_,
            imag_major_A2_this_massless_pair=None,
            imag_major_A2_characteristic_integral_this_massless_pair=None,
            imag_major_A2_transport_rhs_this_massless_pair=imag[
                'transport_rhs_this_massless_pair'],
            real_major_A2=None,
            major_A3=None,
            characteristic_integrals_owned=False,
            characteristic_transport_equations_owned=True,
            history_minus_reference_imag_major_A2=imag['numerical_value'],
            imag_major_A2_vanishes_for_this_massless_pair=pair['mass'] == 0.0,
        ),
        'vacuum_N_beta_tail': _quantity(
            'vacuum_N_beta_tail',
            formula=(
                'occupied vacuum-envelope N,beta tail above the archived '
                '160/320 splits after C4 and C_M; not the Fermi occupation'),
            derivation=(
                'C4 and C_M remain missing, so the vacuum tail stays None; '
                'thermal occupation is bounded separately'),
            numerical_value=None,
            status='OPEN',
            missing_primitive=(
                'numerical C4 and C_M after minor A3, L0 A4 H2 and the '
                'same-column upstream H2 remainder; the E^{-2} cancellation '
                'is not this tail'),
            bindings=slice_,
            N=None,
            beta=None,
            vacuum_integrated_tail_160=None,
            vacuum_integrated_tail_320=None,
            thermal_occupation_bounded_separately=True,
        ),
    }
    if tuple(quantities) != QUANTITY_NAMES:
        raise ArithmeticError('changed-history UV quantity names drifted')
    return MappingProxyType(quantities)

def validate_changed_history_uv_v3(report, root=None):
    """v3 validator: proved transport zeros/majorants, still no invented C4/C_M."""
    report = _as_map(report, 'v3 successor report')
    if report.get('schema') in (PILOT_SCHEMA, SUCCESSOR_SCHEMA):
        raise ValueError('historical v1/v2 remainder schemas are not the v3 successor')
    if report.get('schema') != V3_SUCCESSOR_SCHEMA:
        raise ValueError('unexpected changed-history UV v3 schema')
    status = report.get('status', V3_SUCCESSOR_STATUS)
    if not isinstance(status, str) or not status.startswith('OPEN'):
        raise ValueError('OPEN calculation may not be reported as PASS or NON_EXISTENCE')
    if 'NON-EXISTENCE' in status or status.startswith('PASS'):
        raise ValueError('OPEN calculation may not be reported as PASS or NON_EXISTENCE')
    if report.get('physical_EXISTENCE_certificate') or report.get(
            'physical_NONEXISTENCE_certificate'):
        raise ValueError('OPEN calculation may not be reported as PASS or NON_EXISTENCE')
    if report.get('physical_local_gate') == 'PASS':
        raise ValueError('OPEN calculation may not be reported as PASS or NON_EXISTENCE')
    if report.get('certificate_from_leading_e_minus2_cancellation'):
        raise ValueError('E^{-2} cancellation is not a tail bound')
    if report.get('C_M_or_physical_gate_from_leading_cancellation'):
        raise ValueError('leading cancellation does not supply C_M')
    if report.get('first_noncancelling_paired_order') != FIRST_NONCANCELLING_PAIRED_ORDER:
        raise ValueError('E^{-3} first allowed paired order was not preserved')
    if report.get('leading_paired_e_minus2_cancels_for_equal_mu') is False:
        raise ValueError('E^{-2} equal-mu cancellation was not preserved')
    if report.get('field_or_source_evolutions') not in (0, None):
        raise ValueError('successor control may not launch field or source evolutions')
    if report.get('v1_pilot_rewritten') or report.get('v2_record_rewritten'):
        raise ValueError('historical v1/v2 remainder bytes are immutable')
    bindings = _as_map(report.get('bindings'), 'bindings')
    authenticate_changed_history_uv_bindings(bindings, root)
    transport = _as_map(report.get('transport'), 'transport')
    validate_characteristic_transports(transport, root)
    quantities = _as_map(report.get('quantities'), 'quantities')
    if set(quantities) != set(QUANTITY_NAMES):
        raise ValueError('changed-history UV quantity names drifted')
    slice_ = _quantity_binding_slice(bindings)
    for name in QUANTITY_NAMES:
        item = _as_map(quantities[name], name)
        if item.get('name') != name:
            raise ValueError('quantity name mismatch: ' + name)
        if item.get('status') != 'OPEN':
            raise ValueError('quantity ' + name + ' must remain OPEN')
        if not item.get('leading_e_minus2_cancellation_is_not_this_value'):
            raise ValueError('E^{-2} cancellation is not ' + name)
        if not item.get('e_minus3_first_allowed_order_is_not_a_tail_bound'):
            raise ValueError('E^{-3} first allowed order is not a tail bound')
        if item.get('background_tail_certificate_transferred'):
            raise ValueError('background tail certificate must not transfer onto ' + name)
        if item.get('replaced_by_local_metric_jet'):
            raise ValueError(name + ' may not be replaced by a local metric jet')
        if item.get('fitted_decay'):
            raise ValueError(name + ' may not use a fitted decay')
        bound = item.get('bindings')
        if _as_map(bound, name + ' bindings') != dict(slice_):
            raise ValueError(
                name + ' is not bound to state law/profile/preparation/inventory/subtraction')
        _require_missing(item.get('numerical_value'), name)
    vacuum = quantities['vacuum_N_beta_tail']
    for key in ('N', 'beta', 'vacuum_integrated_tail_160', 'vacuum_integrated_tail_320'):
        _require_missing(vacuum.get(key), 'vacuum ' + key)
    _require_missing(quantities['C4']['numerical_value'], 'C4')
    _require_missing(quantities['C_M']['numerical_value'], 'C_M')
    _require_missing(quantities['C_M'].get('C_M_h2'), 'C_M_h2')
    _require_missing(
        quantities['A2'].get('real_major_A2_numerical_value'), 'real major A2')
    _require_missing(
        quantities['A2'].get('imag_major_A2_numerical_value'), 'imag major A2')
    _require_missing(quantities['transported_majors'].get('real_major_A2'), 'real major A2')
    _require_missing(quantities['transported_majors'].get('major_A3'), 'major A3')
    _require_missing(quantities['transported_majors'].get(
        'imag_major_A2_this_massless_pair'), 'imag major A2')
    _require_missing(quantities['transported_majors'].get(
        'imag_major_A2_characteristic_integral_this_massless_pair'),
        'imag major A2 characteristic value')
    if quantities['A2'].get('imag_major_transport_rhs_this_massless_pair') not in (0.0, 0):
        raise ValueError('massless Im(T_s major A2) identity was not preserved')
    if quantities['A2'].get('history_minus_reference_imag_major_A2') not in (0.0, 0):
        raise ValueError('massless history-minus-reference Im(major A2) must be exact 0')
    if not quantities['A2'].get('history_minus_reference_imag_major_A2_proved_exact'):
        raise ValueError('massless history-minus-reference Im(major A2) must be proved exact')
    if quantities['A2'].get('coefficient_scope') != 'history-minus-reference':
        raise ValueError('A2 has the wrong coefficient scope')
    if not quantities['C4'].get('history_minus_reference_coefficient_formed'):
        raise ValueError('production history-minus-reference C4 definition was not formed')
    if quantities['C4'].get('full_envelope_coefficient_formed'):
        raise ValueError('C4 may not use the generic full-envelope coefficient')
    if quantities['C4'].get('definition_scope') != C4_SCOPE:
        raise ValueError('C4 has the wrong coefficient scope')
    if quantities['A3'].get('A3_locally_algebraic_in_metric_jets'):
        raise ValueError('A3 may not be replaced by a local metric jet')
    if quantities['upstream_higher_order_H2_remainder'].get('frozen_background_tail_copied'):
        raise ValueError('frozen-background tail certificate copied without a proved bridge')
    if quantities['C_M'].get('manufactured_fixture_used_as_production_constant'):
        raise ValueError('manufactured fixture cannot supply current-history C_M')
    if quantities['current_history_L0_A4_H2_integrals'].get('production_L0_A4_control'):
        raise ValueError('manufactured fixture cannot supply current-history L0 A4')
    identities_report = _as_map(report.get('identities'), 'identity summary')
    if identities_report.get('dummy_symbol_L0_A_j'):
        raise ValueError('v3 transport may not use dummy L0 A_j symbols')
    if identities_report.get('production_A3_A4_constructed'):
        raise ValueError('formal M=4 telescope may not claim production A3/A4')
    manufactured = _as_map(report.get('manufactured'), 'manufactured control')
    if (manufactured.get('production_L0_A4_control')
            or manufactured.get('certificate_use')
            or manufactured.get('current_history_C_M') is not None):
        raise ValueError('manufactured fixture cannot supply current-history C_M')
    negative = _as_map(report.get('negative_controls'), 'negative controls')
    if (negative.get('production_C4_coefficient_control')
            or negative.get('certificate_use')):
        raise ValueError('manufactured parity fixture cannot supply C4')
    thermal = _as_map(report.get('thermal'), 'thermal occupation')
    if thermal.get('vacuum_envelope_used') or thermal.get(
            'background_tail_certificate_transferred'):
        raise ValueError('thermal occupation cannot become a vacuum-tail bound')
    for row in thermal.get('cutoffs', ()):
        row = _as_map(row, 'thermal cutoff')
        _require_missing(row.get('N'), 'thermal N float')
        _require_missing(row.get('beta'), 'thermal beta float')
        if (row.get('N_decimal_upper') in (None, '0.0')
                or row.get('beta_decimal_upper') in (None, '0.0')):
            raise ValueError('thermal occupation upper must remain directed and nonzero')
    windows = report.get('direct_windows')
    if windows is not None:
        windows = _as_map(windows, 'direct_windows')
        _require_missing(windows.get('vacuum_integrated_tail_N_beta'), 'vacuum tail')
        if windows.get('background_tail_certificate_transferred'):
            raise ValueError('background tail certificate must not transfer')
    _require_missing(report.get('numerical_C4'), 'C4')
    _require_missing(report.get('numerical_C_M'), 'C_M')
    _require_missing(report.get('vacuum_integrated_tail_N_beta'), 'vacuum tail')
    proved = _as_map(report.get('proved'), 'proved values')
    if proved.get('history_minus_reference_imag_major_A2_massless') not in (0.0, 0):
        raise ValueError('proved massless Im(major A2) difference must be exact 0')
    for name in (
            'q_abs_upper', 'delta_q_abs_upper',
            'history_minus_reference_real_major_A2_abs_upper'):
        value = proved.get(name)
        if isinstance(value, (bool, np.bool_)) or value is None or float(value) <= 0:
            raise ValueError('proved positive majorant required: ' + name)
    if report.get('next_primitive') != TRANSPORT_NEXT_PRIMITIVE:
        raise ValueError('v3 successor must name the exact next primitive')
    return True

def changed_history_uv_v3_report(root=None):
    """v3 successor: own transports, form history-minus-reference C4, keep C_M None."""
    root = Path(root) if root is not None else _root()
    bindings = changed_history_uv_bindings(root)
    transport = current_history_characteristic_transport_report(root)
    quantities = changed_history_uv_v3_quantities(root, bindings, transport)
    first = first_noncancelling_paired_coefficient()
    identities = remainder_recurrence_identities()
    transport_identities = transport['identities']
    defect = truncated_envelope_L0_defect()
    manufactured = manufactured_defect_majorant_control()
    negative = negative_pair_and_omission_controls()
    pair = smallest_representative_pair(root)
    thermal = fermi_thermal_tail_bound(
        pair['mass'], pair['absolute_angular'],
        pair['multiplicity_per_signed_family'],
        cutoffs=pair['archived_splits'])
    windows = thermal_tail_direct_window_check(
        pair['mass'], pair['absolute_angular'],
        pair['multiplicity_per_signed_family'],
        cutoffs=pair['archived_splits'])
    family = LocalIncomingFamily(np.array(
        json.loads((root / HISTORY_RECORD).read_text())['history']['coefficients']))
    commutators = current_history_commutator_integrals(
        family, pair['absolute_angular'])
    c_m = conditional_m4_remainder_majorant(
        None, None, commutators['Bz_integral'], commutators['Bzz_integral'])
    if manufactured['production_L0_A4_control'] or c_m['numerical_C_M'] is not None:
        raise ArithmeticError('manufactured fixture leaked into production C_M')
    proved = {
        'history_minus_reference_imag_major_A2_massless': quantities['A2'][
            'history_minus_reference_imag_major_A2'],
        'q_abs_upper': quantities['A2']['q_abs_upper'],
        'delta_q_abs_upper': quantities['A2']['delta_q_abs_upper'],
        'history_minus_reference_real_major_A2_abs_upper': quantities['A2'][
            'history_minus_reference_real_major_A2_abs_upper'],
        'history_minus_reference_C4_definition_formed': True,
        'numerical_C4': None,
        'numerical_C_M': None,
        'vacuum_integrated_tail_N_beta': None,
    }
    report = MappingProxyType({
        'schema': V3_SUCCESSOR_SCHEMA,
        'status': V3_SUCCESSOR_STATUS,
        'profile_identity': bindings['profile']['identity'],
        'pair': dict(pair),
        'bindings': bindings,
        'transport': transport,
        'quantities': quantities,
        'proved': proved,
        'next_primitive': TRANSPORT_NEXT_PRIMITIVE,
        'first_noncancelling_paired_order': first[
            'first_noncancelling_paired_inverse_energy_order'],
        'leading_paired_e_minus2_cancels_for_equal_mu': first[
            'leading_paired_e_minus2_cancels_for_equal_mu'],
        'e_minus3_even_in_ell': first['e_minus3_even_in_ell'],
        'e_minus3_first_allowed_order_is_not_a_tail_bound': True,
        'leading_e_minus2_cancellation_is_not_a_tail_bound': True,
        'identities': MappingProxyType({
            'expansion_order': identities['expansion_order'],
            'defect': identities['defect'],
            'A3_locally_algebraic_in_metric_jets': False,
            'major_A2_replaced_by_local_metric_jet': False,
            'dummy_symbol_L0_A_j': transport_identities['dummy_symbol_L0_A_j'],
            'formal_symbolic_telescope': identities['formal_symbolic_telescope'],
            'production_A2_transport_constructed': transport_identities[
                'production_A2_transport_constructed'],
            'production_A3_A4_constructed': False,
            'telescope': identities['telescope'],
        }),
        'defect': MappingProxyType({
            'M': defect['M'],
            'defect': defect['defect'],
            'dummy_symbol_L0_A_j': False,
            'formal_symbolic_telescope': defect['formal_symbolic_telescope'],
            'production_A3_A4_constructed': False,
        }),
        'commutators': dict(commutators),
        'manufactured': manufactured,
        'negative_controls': negative,
        'thermal': {
            'cutoffs': [dict(row) for row in thermal['cutoffs']],
            'source_law': thermal['source_law'],
            'vacuum_envelope_used': thermal['vacuum_envelope_used'],
            'background_tail_certificate_transferred': False,
        },
        'direct_windows': windows,
        'numerical_C4': None,
        'numerical_C_M': None,
        'vacuum_integrated_tail_N_beta': None,
        'field_or_source_evolutions': 0,
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
        'physical_local_gate': 'OPEN',
        'certificate_from_leading_e_minus2_cancellation': False,
        'C_M_or_physical_gate_from_leading_cancellation': False,
        'new_action_term': False,
        'Gamma_rest_assigned': False,
        'existing_source_ad759_evaluation_changed': False,
        'v2_v4_records_rewritten': False,
        'v1_pilot_schema': PILOT_SCHEMA,
        'v1_pilot_record': PILOT_RECORD,
        'v1_pilot_sha256': _digest(root, PILOT_RECORD),
        'v1_pilot_rewritten': False,
        'v2_record': V2_REMAINDER_RECORD,
        'v2_record_sha256': _digest(root, V2_REMAINDER_RECORD),
        'v2_record_rewritten': False,
        'provisional_UV_allocation': PROVISIONAL_UV_ALLOCATION,
        'budget_implication': (
            'Gate 2 remains OPEN. Proved massless history-minus-reference '
            'Im(major A2)=0 and interval majorants of q and Re(major A2) do '
            'not close the UV allocation; C4, C_M and the vacuum tail are '
            'still missing'),
    })
    validate_changed_history_uv_v3(report, root)
    return report
