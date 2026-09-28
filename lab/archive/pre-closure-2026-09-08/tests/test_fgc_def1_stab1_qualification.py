from __future__ import annotations

from dataclasses import replace
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from recursive_horizons.fgc.constraint_system import (  # noqa: E402
    ACTIVE_SPHERICAL_GAUGE_COMPONENTS,
)
from recursive_horizons.fgc.def1_geometry_error import INPUT_NAMES  # noqa: E402
from recursive_horizons.fgc.def1_stab1 import (  # noqa: E402
    DEF1_BOOLEAN_NAMES,
    ERROR_BUDGET_COMPONENTS,
    PREMISE_STATUS_CONDITIONAL,
    SOURCE_IMP1_ADMISSION_DEBIT,
)
from recursive_horizons.fgc.def1_stab1_qualification import (  # noqa: E402
    ADM_METRIC_TWO_JET_ORDER,
    ADM_SOURCE_PHYSICAL_ARGUMENT_ORDER,
    BASE_METRIC_TWO_JET_ORDER,
    CONVERSION_RED1_MHG2,
    Def1Stab1QualificationError,
    FO1_GROUP_TO_JET_SLOT,
    IMP1_18_CHANNEL_ORDER,
    INVENTORY_ADM_METRIC_TWO_JET,
    INVENTORY_ADM_SOURCE_PHYSICAL,
    INVENTORY_ADM_SOURCE_STATE,
    INVENTORY_FO1_BASE_PHYSICAL,
    INVENTORY_FO1_BASE_STATE,
    INVENTORY_FO1_BASE_TO_ADM,
    PROVIDER_ROUTES,
    RED1_TO_MHG2_EQUATION_PAIRS,
    SLOT_EXACT_INVERSE,
    bind_qualification_provenance,
    convert_base_to_adm_geometry_slots,
    convert_fo1_to_geometry_slots,
    convert_imp1_channels_to_q,
    execute_qualification_coverage_matrix,
    fo1_adm_geometry_slot_map,
    geometry_slot_zero_justification,
    imp1_channel_geometry_assignments,
    invert_base_metric_two_jets,
    invert_fo1_physical_arguments_to_adm,
    invert_ordered_base_metric_two_jets,
    live_imp1_18_channel_order,
    map_equation_order_to_mhg2,
    owner_file_sha256,
    phi_chi_geometry_status,
    red1_mhg2_equation_identity_map,
)
from recursive_horizons.fgc.modified_harmonic import (  # noqa: E402
    MHG_EQUATION_ORDER,
    MHG_GAUGE_CONSTRAINT_ORDER,
)
from recursive_horizons.fgc.modified_harmonic_constraints import (  # noqa: E402
    PHYSICAL_PROJECTION_ORDER,
)
from recursive_horizons.fgc.modified_harmonic_first_order import (  # noqa: E402
    FO1_PHYSICAL_ARGUMENT_ORDER,
    FO1_PHYSICAL_EQUATION_ORDER,
    FO1_REDUCTION_CONSTRAINT_ORDER,
    FO1_STATE_ORDER,
)
from recursive_horizons.fgc.modified_harmonic_implicit import (  # noqa: E402
    IMPLICIT_EQUATION_ORDER,
)
from recursive_horizons.fgc.modified_harmonic_reference import (  # noqa: E402
    MHG2_FULL_EQUATION_ORDER,
)
from recursive_horizons.fgc.regular_center import REGULAR_EQUATION_ORDER  # noqa: E402
from recursive_horizons.fgc.sgb1_ctl1_source import SOURCE_EQUATION_ORDER  # noqa: E402
from recursive_horizons.fgc.spherical_reduction import (  # noqa: E402
    INDEPENDENT_EQUATION_ORDER,
    Jet2,
    state_from_generalized_adm_pg_fixture,
)


def _ones(names: tuple[str, ...]) -> dict[str, Q]:
    return {name: Q(index + 1, 8) for index, name in enumerate(names)}


def _nonzero_shift_adm() -> tuple[Jet2, Jet2, Jet2, Jet2]:
    return (
        Jet2(2, Q(1, 3), Q(-1, 5), Q(1, 7), Q(1, 11), Q(-1, 13)),
        Jet2(Q(1, 2), Q(1, 4), Q(1, 6), Q(-1, 8), Q(1, 9), Q(1, 10)),
        Jet2(3, Q(1, 2), Q(-1, 3), Q(1, 5), Q(-1, 6), Q(1, 7)),
        Jet2(4, 0, 1, 0, 0, 0),
    )


def _jet_mapping(jet: Jet2) -> dict[str, Q]:
    return {
        part: getattr(jet, part)
        for part in ("value", "dt", "dr", "dtt", "dtr", "drr")
    }


def _base_from_adm(
    alpha: Jet2, shift: Jet2, lam: Jet2, radius: Jet2
):
    return state_from_generalized_adm_pg_fixture(
        {
            "state": {
                "alpha": _jet_mapping(alpha),
                "shift": _jet_mapping(shift),
                "lambda": _jet_mapping(lam),
                "areal_radius": _jet_mapping(radius),
                "phi": _jet_mapping(Jet2.constant(0)),
                "chi": _jet_mapping(Jet2.constant(0)),
            },
            "branch": "GR-0",
        }
    )


def _fo1_physical_from_base(state) -> dict[str, Q]:
    values: dict[str, Q] = {}
    for group, slot in FO1_GROUP_TO_JET_SLOT:
        for field in (
            "h_tt",
            "h_tr",
            "h_rr",
            "areal_radius",
            "phi",
            "chi",
        ):
            values[f"{group}.{field}"] = getattr(getattr(state, field), slot)
    return values


class EquationOrderIdentityTests(unittest.TestCase):
    def test_identity_is_exact_six_row_owner_order(self) -> None:
        identity = red1_mhg2_equation_identity_map()
        self.assertEqual(identity.pairs, RED1_TO_MHG2_EQUATION_PAIRS)
        self.assertEqual(identity.red1_order, INDEPENDENT_EQUATION_ORDER)
        self.assertEqual(identity.mhg2_order, MHG2_FULL_EQUATION_ORDER)
        self.assertEqual(len(identity.pairs), 6)
        self.assertEqual(
            map_equation_order_to_mhg2(INDEPENDENT_EQUATION_ORDER, owner_id="RED1"),
            MHG2_FULL_EQUATION_ORDER,
        )
        for owner_id, order in (
            ("MHG1", MHG_EQUATION_ORDER),
            ("MHG2-REF1", MHG2_FULL_EQUATION_ORDER),
            ("MHG3-IMP1", IMPLICIT_EQUATION_ORDER),
            ("FO1-physical", FO1_PHYSICAL_EQUATION_ORDER),
            ("SGB1-source", SOURCE_EQUATION_ORDER),
            ("SRC1-NL1", MHG2_FULL_EQUATION_ORDER),
            ("CON4-unredefined-residual-names", INDEPENDENT_EQUATION_ORDER),
        ):
            self.assertEqual(
                map_equation_order_to_mhg2(order, owner_id=owner_id),
                MHG2_FULL_EQUATION_ORDER,
            )
        self.assertFalse(identity.claims_col1_values)
        self.assertFalse(identity.global_pde_error_certified)

    def test_omissions_aliases_and_reordering_refuse(self) -> None:
        swapped = (
            INDEPENDENT_EQUATION_ORDER[1],
            INDEPENDENT_EQUATION_ORDER[0],
        ) + INDEPENDENT_EQUATION_ORDER[2:]
        with self.assertRaisesRegex(Def1Stab1QualificationError, "reordering"):
            map_equation_order_to_mhg2(swapped, owner_id="RED1")
        with self.assertRaisesRegex(Def1Stab1QualificationError, "aliases"):
            map_equation_order_to_mhg2(
                {name: index for index, name in enumerate(INDEPENDENT_EQUATION_ORDER)},
                owner_id="RED1",
            )
        with self.assertRaisesRegex(Def1Stab1QualificationError, "omission"):
            map_equation_order_to_mhg2(
                INDEPENDENT_EQUATION_ORDER[:-1], owner_id="RED1"
            )
        with self.assertRaisesRegex(Def1Stab1QualificationError, "aliases"):
            map_equation_order_to_mhg2(MHG2_FULL_EQUATION_ORDER, owner_id="RED1")
        with self.assertRaisesRegex(Def1Stab1QualificationError, "aliases"):
            map_equation_order_to_mhg2(
                INDEPENDENT_EQUATION_ORDER, owner_id="MHG2-REF1"
            )

    def test_non_six_row_owners_are_not_mhg2_aliases(self) -> None:
        for owner_id, order in (
            ("CON4-physical-projection", PHYSICAL_PROJECTION_ORDER),
            ("MHG1-gauge-constraint", MHG_GAUGE_CONSTRAINT_ORDER),
            ("CON4-gauge-vector", ACTIVE_SPHERICAL_GAUGE_COMPONENTS),
            ("FO1-reduction-constraint", FO1_REDUCTION_CONSTRAINT_ORDER),
            ("REG1-regular-center", REGULAR_EQUATION_ORDER),
        ):
            with self.subTest(owner_id=owner_id):
                with self.assertRaisesRegex(
                    Def1Stab1QualificationError, "cannot be mapped"
                ):
                    map_equation_order_to_mhg2(order, owner_id=owner_id)
        with self.assertRaisesRegex(Def1Stab1QualificationError, "unknown"):
            map_equation_order_to_mhg2(
                INDEPENDENT_EQUATION_ORDER, owner_id="COL1"
            )

    def test_replace_cannot_promote_or_reorder(self) -> None:
        identity = red1_mhg2_equation_identity_map()
        with self.assertRaisesRegex(Def1Stab1QualificationError, "COL1|PDE"):
            replace(identity, claims_col1_values=True)
        with self.assertRaisesRegex(Def1Stab1QualificationError, "COL1|PDE"):
            replace(identity, global_pde_error_certified=True)
        with self.assertRaisesRegex(Def1Stab1QualificationError, "live owners"):
            replace(identity, pairs=tuple(reversed(identity.pairs)))


class Fo1AdmGeometrySlotTests(unittest.TestCase):
    def test_contract_classifies_every_geometry_slot(self) -> None:
        contract = fo1_adm_geometry_slot_map()
        self.assertEqual(contract.jet_slot_identity, FO1_GROUP_TO_JET_SLOT)
        self.assertFalse(contract.claims_col1_values)
        self.assertFalse(contract.global_pde_error_certified)
        for _name, assignments in contract.inventories:
            self.assertEqual(
                tuple(item.geometry_slot for item in assignments), INPUT_NAMES
            )
        for slot in INPUT_NAMES:
            self.assertIsNone(geometry_slot_zero_justification(slot))
        self.assertIn("not_a_geometry_slot", phi_chi_geometry_status())

    def test_fo1_base_fills_only_areal_radius_identity(self) -> None:
        contract = fo1_adm_geometry_slot_map()
        base_slots = dict(contract.inventories)[INVENTORY_FO1_BASE_PHYSICAL]
        statuses = {item.geometry_slot: item.status for item in base_slots}
        self.assertEqual(statuses["alpha.value"], SLOT_EXACT_INVERSE)
        self.assertEqual(statuses["shift.dtt"], SLOT_EXACT_INVERSE)
        self.assertEqual(statuses["lambda.dr"], SLOT_EXACT_INVERSE)
        converted = convert_fo1_to_geometry_slots(
            _ones(FO1_PHYSICAL_ARGUMENT_ORDER),
            inventory=INVENTORY_FO1_BASE_PHYSICAL,
            require_complete=False,
        )
        supplied = dict(converted.supplied)
        self.assertEqual(
            set(supplied),
            {f"R.{jet}" for jet in ("value", "dt", "dr", "dtt", "dtr", "drr")},
        )
        self.assertIn("alpha.value", converted.refused)
        self.assertIn("shift.dt", converted.refused)
        self.assertIn("lambda.dtt", converted.refused)
        self.assertIn("k.t", converted.refused)
        self.assertFalse(converted.complete)
        self.assertFalse(converted.used_implicit_zero)
        with self.assertRaisesRegex(Def1Stab1QualificationError, "refused="):
            convert_fo1_to_geometry_slots(
                _ones(FO1_PHYSICAL_ARGUMENT_ORDER),
                inventory=INVENTORY_FO1_BASE_PHYSICAL,
            )

    def test_fo1_first_order_state_refuses_second_jets(self) -> None:
        converted = convert_fo1_to_geometry_slots(
            _ones(FO1_STATE_ORDER),
            inventory=INVENTORY_FO1_BASE_STATE,
            require_complete=False,
        )
        self.assertEqual(
            set(dict(converted.supplied)), {"R.value", "R.dt", "R.dr"}
        )
        self.assertIn("R.dtt", converted.refused)
        self.assertIn("R.dtr", converted.refused)
        self.assertIn("R.drr", converted.refused)

    def test_adm_source_physical_plus_tangent_is_complete(self) -> None:
        converted = convert_fo1_to_geometry_slots(
            _ones(ADM_SOURCE_PHYSICAL_ARGUMENT_ORDER),
            inventory=INVENTORY_ADM_SOURCE_PHYSICAL,
            tangent=(("k.t", Q(3)), ("k.r", Q(4))),
        )
        self.assertTrue(converted.complete)
        self.assertEqual(tuple(name for name, _ in converted.supplied), INPUT_NAMES)
        self.assertEqual(converted.supplied[INPUT_NAMES.index("k.t")][1], Q(3))
        self.assertFalse(converted.used_implicit_zero)
        with self.assertRaisesRegex(Def1Stab1QualificationError, "k.t"):
            convert_fo1_to_geometry_slots(
                _ones(ADM_SOURCE_PHYSICAL_ARGUMENT_ORDER),
                inventory=INVENTORY_ADM_SOURCE_PHYSICAL,
            )

    def test_adm_two_jet_identity_and_alias_refusal(self) -> None:
        converted = convert_fo1_to_geometry_slots(
            _ones(ADM_METRIC_TWO_JET_ORDER),
            inventory=INVENTORY_ADM_METRIC_TWO_JET,
            tangent={"k.t": Q(1), "k.r": Q(2)},
        )
        self.assertTrue(converted.complete)
        with self.assertRaisesRegex(Def1Stab1QualificationError, "extra|missing"):
            convert_fo1_to_geometry_slots(
                {"areal_radius.value": Q(1)},
                inventory=INVENTORY_ADM_METRIC_TWO_JET,
                tangent={"k.t": 1, "k.r": 1},
                require_complete=False,
            )
        with self.assertRaisesRegex(Def1Stab1QualificationError, "extra|missing"):
            convert_fo1_to_geometry_slots(
                _ones(ADM_SOURCE_PHYSICAL_ARGUMENT_ORDER) | {"u.v": Q(1)},
                inventory=INVENTORY_ADM_SOURCE_PHYSICAL,
                require_complete=False,
            )
        reordered = tuple(reversed(_ones(ADM_METRIC_TWO_JET_ORDER).items()))
        with self.assertRaisesRegex(Def1Stab1QualificationError, "reordering"):
            convert_fo1_to_geometry_slots(
                reordered,
                inventory=INVENTORY_ADM_METRIC_TWO_JET,
                require_complete=False,
            )

    def test_first_order_adm_state_does_not_zero_second_jets(self) -> None:
        converted = convert_fo1_to_geometry_slots(
            _ones(tuple(
                f"{group}.{field}"
                for group in ("u", "p", "q")
                for field in ("alpha", "shift", "lambda", "R", "phi", "chi")
            )),
            inventory=INVENTORY_ADM_SOURCE_STATE,
            require_complete=False,
        )
        supplied = set(dict(converted.supplied))
        self.assertIn("alpha.value", supplied)
        self.assertIn("shift.dr", supplied)
        self.assertNotIn("alpha.dtt", supplied)
        self.assertIn("alpha.dtt", converted.refused)
        self.assertFalse(converted.used_implicit_zero)


class BaseToAdmTwoJetInverseTests(unittest.TestCase):
    def test_nonzero_shift_roundtrip_and_complete_26_input(self) -> None:
        alpha, shift, lam, radius = _nonzero_shift_adm()
        self.assertNotEqual(shift.value, 0)
        base = _base_from_adm(alpha, shift, lam, radius)
        inverse = invert_base_metric_two_jets(
            h_tt=base.h_tt,
            h_tr=base.h_tr,
            h_rr=base.h_rr,
            areal_radius=base.areal_radius,
            lapse_root=alpha.value,
            radial_scale_root=lam.value,
        )
        self.assertEqual(inverse.alpha, alpha)
        self.assertEqual(inverse.shift, shift)
        self.assertEqual(inverse.lambda_jet, lam)
        self.assertEqual(inverse.areal_radius, radius)
        self.assertTrue(inverse.roundtrip_holds)
        self.assertFalse(inverse.claims_col1_values)
        self.assertFalse(inverse.global_pde_error_certified)
        self.assertFalse(inverse.used_implicit_zero)
        packed = convert_base_to_adm_geometry_slots(
            inverse, tangent=(("k.t", Q(5, 2)), ("k.r", Q(-1, 4)))
        )
        self.assertEqual(packed.inventory, INVENTORY_FO1_BASE_TO_ADM)
        self.assertTrue(packed.complete)
        self.assertEqual(tuple(name for name, _ in packed.supplied), INPUT_NAMES)
        values = dict(packed.supplied)
        self.assertEqual(values["alpha.value"], alpha.value)
        self.assertEqual(values["shift.dtr"], shift.dtr)
        self.assertEqual(values["lambda.drr"], lam.drr)
        self.assertEqual(values["R.dr"], radius.dr)
        self.assertEqual(values["k.t"], Q(5, 2))
        self.assertEqual(values["k.r"], Q(-1, 4))
        self.assertNotIn("phi.value", values)
        self.assertNotIn("chi.value", values)
        self.assertFalse(packed.used_implicit_zero)
        statuses = {item.geometry_slot: item.status for item in packed.assignments}
        self.assertEqual(statuses["alpha.dtt"], SLOT_EXACT_INVERSE)
        self.assertEqual(statuses["R.value"], "documented_name_identity")

    def test_fo1_physical_arguments_invert_with_supplied_roots(self) -> None:
        alpha, shift, lam, radius = _nonzero_shift_adm()
        base = _base_from_adm(alpha, shift, lam, radius)
        inverse = invert_fo1_physical_arguments_to_adm(
            _fo1_physical_from_base(base),
            lapse_root=alpha.value,
            radial_scale_root=lam.value,
        )
        self.assertEqual(inverse.alpha, alpha)
        self.assertEqual(inverse.shift, shift)
        self.assertEqual(inverse.lambda_jet, lam)

    def test_wrong_root_sign_singular_incomplete_and_reorder_refuse(self) -> None:
        alpha, shift, lam, radius = _nonzero_shift_adm()
        base = _base_from_adm(alpha, shift, lam, radius)
        with self.assertRaisesRegex(Def1Stab1QualificationError, "square exactly"):
            invert_base_metric_two_jets(
                h_tt=base.h_tt,
                h_tr=base.h_tr,
                h_rr=base.h_rr,
                areal_radius=base.areal_radius,
                lapse_root=alpha.value,
                radial_scale_root=lam.value + 1,
            )
        with self.assertRaisesRegex(Def1Stab1QualificationError, "positive branch"):
            invert_base_metric_two_jets(
                h_tt=base.h_tt,
                h_tr=base.h_tr,
                h_rr=base.h_rr,
                areal_radius=base.areal_radius,
                lapse_root=alpha.value,
                radial_scale_root=-lam.value,
            )
        with self.assertRaisesRegex(Def1Stab1QualificationError, "positive branch"):
            invert_base_metric_two_jets(
                h_tt=base.h_tt,
                h_tr=base.h_tr,
                h_rr=base.h_rr,
                areal_radius=base.areal_radius,
                lapse_root=-alpha.value,
                radial_scale_root=lam.value,
            )
        with self.assertRaisesRegex(Def1Stab1QualificationError, "strictly positive"):
            invert_base_metric_two_jets(
                h_tt=base.h_tt,
                h_tr=base.h_tr,
                h_rr=Jet2(0, 1, 1, 1, 1, 1),
                areal_radius=base.areal_radius,
                lapse_root=1,
                radial_scale_root=1,
            )
        with self.assertRaisesRegex(Def1Stab1QualificationError, "Lorentzian"):
            invert_base_metric_two_jets(
                h_tt=Jet2(1),
                h_tr=Jet2(0),
                h_rr=Jet2(1),
                areal_radius=Jet2(2),
                lapse_root=1,
                radial_scale_root=1,
            )
        incomplete = {
            "value": base.h_rr.value,
            "dt": base.h_rr.dt,
            "dr": base.h_rr.dr,
            "dtt": base.h_rr.dtt,
            "dtr": base.h_rr.dtr,
        }
        with self.assertRaisesRegex(Def1Stab1QualificationError, "complete two-jet"):
            invert_base_metric_two_jets(
                h_tt=base.h_tt,
                h_tr=base.h_tr,
                h_rr=incomplete,
                areal_radius=base.areal_radius,
                lapse_root=alpha.value,
                radial_scale_root=lam.value,
            )
        ordered = (
            ("h_rr", base.h_rr),
            ("h_tt", base.h_tt),
            ("h_tr", base.h_tr),
            ("areal_radius", base.areal_radius),
        )
        with self.assertRaisesRegex(Def1Stab1QualificationError, "reordering"):
            invert_ordered_base_metric_two_jets(
                ordered,
                lapse_root=alpha.value,
                radial_scale_root=lam.value,
            )
        self.assertEqual(
            BASE_METRIC_TWO_JET_ORDER,
            ("h_tt", "h_tr", "h_rr", "areal_radius"),
        )
        invert_ordered_base_metric_two_jets(
            (
                ("h_tt", base.h_tt),
                ("h_tr", base.h_tr),
                ("h_rr", base.h_rr),
                ("areal_radius", base.areal_radius),
            ),
            lapse_root=alpha.value,
            radial_scale_root=lam.value,
        )

    def test_replace_cannot_drop_roundtrip_or_claim_col1(self) -> None:
        alpha, shift, lam, radius = _nonzero_shift_adm()
        base = _base_from_adm(alpha, shift, lam, radius)
        inverse = invert_base_metric_two_jets(
            h_tt=base.h_tt,
            h_tr=base.h_tr,
            h_rr=base.h_rr,
            areal_radius=base.areal_radius,
            lapse_root=alpha.value,
            radial_scale_root=lam.value,
        )
        with self.assertRaisesRegex(Def1Stab1QualificationError, "COL1|PDE"):
            replace(inverse, claims_col1_values=True)
        with self.assertRaisesRegex(Def1Stab1QualificationError, "COL1|PDE"):
            replace(inverse, global_pde_error_certified=True)
        with self.assertRaisesRegex(Def1Stab1QualificationError, "roundtrip|zeros"):
            replace(inverse, roundtrip_holds=False)
        with self.assertRaisesRegex(Def1Stab1QualificationError, "live Jet2 inverse"):
            replace(inverse, shift=Jet2.constant(0))
        packed = convert_base_to_adm_geometry_slots(
            inverse, tangent={"k.t": 1, "k.r": 1}
        )
        self.assertTrue(packed.complete)
        with self.assertRaisesRegex(Def1Stab1QualificationError, "missing"):
            convert_base_to_adm_geometry_slots(inverse, tangent={"k.t": 1})


class Imp1ToQConversionTests(unittest.TestCase):
    def test_exact_lipschitz_sum_stays_conditional(self) -> None:
        channels = live_imp1_18_channel_order()
        self.assertEqual(channels, IMP1_18_CHANNEL_ORDER)
        debits = tuple((name, Q(index + 1, 16)) for index, name in enumerate(channels))
        lipschitz = tuple((name, Q(2, 5)) for name in channels)
        converted = convert_imp1_channels_to_q(
            channel_debits=debits,
            lipschitz=lipschitz,
            context="synthetic IMP1-to-Q control",
            unit="complete-Q units",
        )
        expected = sum((Q(index + 1, 16) * Q(2, 5) for index in range(18)), Q(0))
        self.assertEqual(converted.additive_q, expected)
        self.assertEqual(converted.status, PREMISE_STATUS_CONDITIONAL)
        self.assertEqual(converted.source, SOURCE_IMP1_ADMISSION_DEBIT)
        self.assertFalse(converted.global_pde_error_certified)
        self.assertFalse(converted.def1_error_map_passed)
        self.assertFalse(converted.used_measured_q)
        self.assertFalse(converted.imp1_admission_debit_treated_as_global_pde_error)
        self.assertFalse(converted.claims_col1_values)

    def test_universal_factor_aliases_and_reordering_refuse(self) -> None:
        channels = IMP1_18_CHANNEL_ORDER
        debits = [Q(1, 8)] * 18
        with self.assertRaisesRegex(Def1Stab1QualificationError, "universal"):
            convert_imp1_channels_to_q(
                channel_debits=debits,
                lipschitz=1,
                context="synthetic IMP1-to-Q control",
                unit="complete-Q units",
            )
        with self.assertRaisesRegex(Def1Stab1QualificationError, "reordering"):
            convert_imp1_channels_to_q(
                channel_debits=tuple(reversed(list(zip(channels, debits)))),
                lipschitz=debits,
                context="synthetic IMP1-to-Q control",
                unit="complete-Q units",
            )
        with self.assertRaisesRegex(Def1Stab1QualificationError, "extra|missing"):
            convert_imp1_channels_to_q(
                channel_debits={"u:shift": Q(1)},
                lipschitz=debits,
                context="synthetic IMP1-to-Q control",
                unit="complete-Q units",
            )
        with self.assertRaisesRegex(
            Def1Stab1QualificationError, "exactly 18|missing"
        ):
            convert_imp1_channels_to_q(
                channel_debits=debits[:-1],
                lipschitz=debits,
                context="synthetic IMP1-to-Q control",
                unit="complete-Q units",
            )

    def test_replace_cannot_certify_pde_or_claim_col1(self) -> None:
        converted = convert_imp1_channels_to_q(
            channel_debits=[Q(0)] * 18,
            lipschitz=[Q(0)] * 18,
            context="synthetic IMP1-to-Q control",
            unit="complete-Q units",
        )
        self.assertEqual(converted.additive_q, 0)
        self.assertEqual(converted.status, PREMISE_STATUS_CONDITIONAL)
        with self.assertRaisesRegex(Def1Stab1QualificationError, "global PDE"):
            replace(converted, global_pde_error_certified=True)
        with self.assertRaisesRegex(Def1Stab1QualificationError, "global PDE"):
            replace(converted, def1_error_map_passed=True)
        with self.assertRaisesRegex(Def1Stab1QualificationError, "global PDE"):
            replace(converted, claims_col1_values=True)
        assignments = imp1_channel_geometry_assignments()
        self.assertEqual(tuple(item.channel for item in assignments), IMP1_18_CHANNEL_ORDER)
        geometry = {
            item.channel: item.geometry_slot
            for item in assignments
            if item.geometry_slot is not None
        }
        self.assertEqual(geometry["u:v"], "shift.value")
        self.assertEqual(geometry["q:R"], "R.dr")
        self.assertIsNone(
            next(item.geometry_slot for item in assignments if item.channel == "u:phi")
        )


class ProvenanceTests(unittest.TestCase):
    def test_binds_owner_bytes_without_col1_claim(self) -> None:
        record = bind_qualification_provenance(
            artifact_id="FGC-1-HYP1-RED1",
            owner_path="src/recursive_horizons/fgc/spherical_reduction.py",
            conversion_name=CONVERSION_RED1_MHG2,
        )
        self.assertEqual(
            record.owner_sha256,
            owner_file_sha256("src/recursive_horizons/fgc/spherical_reduction.py"),
        )
        self.assertTrue(record.authenticates_owner_bytes)
        self.assertFalse(record.claims_col1_values)
        self.assertFalse(record.authenticates_future_col1)
        with self.assertRaisesRegex(Def1Stab1QualificationError, "COL1"):
            replace(record, claims_col1_values=True)
        with self.assertRaisesRegex(Def1Stab1QualificationError, "COL1"):
            replace(record, authenticates_future_col1=True)
        with self.assertRaisesRegex(Def1Stab1QualificationError, "live owner"):
            replace(record, owner_sha256="0" * 64)
        with self.assertRaisesRegex(Def1Stab1QualificationError, "conversion_identity"):
            replace(record, conversion_identity="0" * 64)
        with self.assertRaisesRegex(Def1Stab1QualificationError, "escaped|relative"):
            bind_qualification_provenance(
                artifact_id="FGC-1-HYP1-RED1",
                owner_path="../spherical_reduction.py",
                conversion_name=CONVERSION_RED1_MHG2,
            )


class CoverageMatrixTests(unittest.TestCase):
    def test_every_provider_route_has_positive_and_injected_failure_controls(self) -> None:
        matrix = execute_qualification_coverage_matrix()
        self.assertEqual(matrix.components, ERROR_BUDGET_COMPONENTS)
        self.assertEqual(len(ERROR_BUDGET_COMPONENTS), 14)
        self.assertEqual(
            tuple((item.component, item.route_id) for item in matrix.routes),
            tuple((item.component, item.route_id) for item in PROVIDER_ROUTES),
        )
        self.assertEqual(
            {item.component for item in matrix.routes},
            set(ERROR_BUDGET_COMPONENTS),
        )
        for item in matrix.routes:
            self.assertTrue(item.positive_control_passed)
            self.assertTrue(item.injected_failure_refused)
            self.assertFalse(item.trajectory_values_evaluated)
            self.assertFalse(item.def1_booleans_evaluated)
            self.assertEqual(item.def1_boolean_names, DEF1_BOOLEAN_NAMES)
        self.assertFalse(matrix.trajectory_values_evaluated)
        self.assertFalse(matrix.def1_booleans_evaluated)
        self.assertFalse(matrix.def1_error_map_passed)
        self.assertFalse(matrix.global_pde_error_certified)
        self.assertFalse(matrix.claims_col1_values)
        self.assertFalse(matrix.used_measured_q)
        self.assertEqual(matrix.provenance.conversion_name, "def1_14_component_route_coverage")

    def test_replace_cannot_evaluate_booleans_or_pass_the_gate(self) -> None:
        matrix = execute_qualification_coverage_matrix()
        with self.assertRaisesRegex(Def1Stab1QualificationError, "promote|DEF1"):
            replace(matrix, def1_error_map_passed=True)
        with self.assertRaisesRegex(Def1Stab1QualificationError, "promote|DEF1"):
            replace(matrix, def1_booleans_evaluated=True)
        with self.assertRaisesRegex(Def1Stab1QualificationError, "promote|DEF1"):
            replace(matrix, trajectory_values_evaluated=True)
        with self.assertRaisesRegex(Def1Stab1QualificationError, "promote|DEF1"):
            replace(matrix, global_pde_error_certified=True)
        broken = replace(matrix.routes[0], positive_control_passed=False)
        with self.assertRaisesRegex(Def1Stab1QualificationError, "lacks positive"):
            replace(matrix, routes=(broken,) + matrix.routes[1:])


if __name__ == "__main__":
    unittest.main()
