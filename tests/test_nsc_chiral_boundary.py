"""Small independent probes; the runner separately reproduces the full record."""
from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.nsc_chiral_boundary import (
    ALPHA, BETA, DEFAULT_RESULT, ENERGIES, GAMMA5, I4, PORTS, SCHEMA,
    algebra_and_domain, clifford_basis, comm, continuum_data, decompose_link,
    eliminate_center, frame_control, green_from_fundamental, joined_green,
    mass_control, norm, regulator_checks, response_row, transfer_green,
)
from recursive_horizons.nsc_boundary import compare


class DomainAndAlgebraTests(unittest.TestCase):
    def test_paired_chirality_and_current_domain(self):
        result = algebra_and_domain()
        for key in ("gamma5_from_Clifford_orientation_residual", "gamma5_square_residual",
                    "kinetic_commutator_residual", "angular_commutator_residual",
                    "compatible_wall_commutator_residual", "compatible_wall_spatial_flux_residual",
                    "mass_anticommutator_residual", "identity_transmission_oriented_flux_residual"):
            self.assertEqual(result[key], 0.)
        self.assertGreater(result["legacy_wall_chirality_leakage_norm"], 1.)
        self.assertGreater(result["same_wall_PG_throat_flux_norm"], 1.)
        endpoints = {r["location"]: r for r in result["oriented_PG_characteristics"]}
        self.assertEqual(endpoints["child_outer_cut"]["incoming_full_spinor_channels"], 0)
        self.assertEqual(endpoints["parent_outer_cut"]["incoming_full_spinor_channels"], 2)
        self.assertEqual(endpoints["parent_throat"]["outgoing_full_spinor_channels"], 4)
        self.assertEqual(endpoints["child_throat"]["incoming_full_spinor_channels"], 4)

    def test_clifford_basis_orthogonal_and_scalar_benchmark(self):
        matrices = list(clifford_basis().values())
        gram = np.array([[np.trace(a.conj().T@b)/4 for b in matrices] for a in matrices])
        np.testing.assert_allclose(gram, np.eye(16), atol=1e-15)
        decomposition = decompose_link(.37*BETA)
        self.assertAlmostEqual(decomposition["action_clifford_coefficients"]["scalar"], .37)
        self.assertEqual(decomposition["chiral_even_link_norm"], 0.)
        control = mass_control(.37)
        self.assertEqual(control["Pi_sheet_commutator_norm"], 0.)
        self.assertLess(control["dispersion_square_residual"], 1e-14)

    def test_schur_symmetry_theorem_with_noncommuting_blocks(self):
        # Independent finite Hermitian example; chirality-invariant blocks
        # need not mutually commute. Distinct rows ensure ordering is tested.
        from recursive_horizons.nsc_chiral_boundary import ANGULAR
        rr = 2.1*I4+.3*ALPHA
        ee = 1.7*I4+.2*ANGULAR
        re = .8*ALPHA+.4j*ANGULAR
        er = re.conj().T
        schur = rr-re@np.linalg.solve(ee, er)
        self.assertLess(norm(comm(schur, GAMMA5)), 1e-14)
        self.assertGreater(norm(schur-(rr-np.linalg.solve(ee, re@er))), .01)


class EvaluatedMapTests(unittest.TestCase):
    def test_direct_joined_transfer_and_elimination(self):
        direct, residual = joined_green(40, ENERGIES[0])
        transfer = transfer_green(40, ENERGIES[0])
        self.assertLess(residual, 1e-12)
        np.testing.assert_allclose(direct, transfer, rtol=1e-9, atol=1e-11)
        two, _ = joined_green(40, ENERGIES[0], points=(1., -1.))
        np.testing.assert_allclose(eliminate_center(direct), np.linalg.inv(two), atol=1e-12)
        # Correct source orientation also enforces the resolvent sign.
        imaginary_part = (direct-direct.conj().T)/(2j)
        self.assertGreater(np.min(np.linalg.eigvalsh(imaginary_part)), 0.)

    def test_continuum_resolution_and_exact_chirality(self):
        continuum = green_from_fundamental(continuum_data(ENERGIES[0]), PORTS)
        target = eliminate_center(continuum)
        errors = []
        for intervals in (40, 80):
            green, _ = joined_green(intervals, ENERGIES[0])
            response = eliminate_center(green)
            errors.append(norm(response-target)/norm(target))
            self.assertLess(norm(comm(response, np.kron(np.eye(2), GAMMA5))), 1e-12)
        self.assertLess(errors[-1], .001)
        self.assertGreater(errors[0]/errors[1], 3.8)
        self.assertLess(errors[0]/errors[1], 4.2)

    def test_wall_and_inserted_mass_controls_detect_breaking(self):
        base = response_row(ENERGIES[0], 40)
        wall = response_row(ENERGIES[0], 40, wall="legacy_radial")
        mass = response_row(ENERGIES[0], 40, mass=.35)
        self.assertLess(base["link"]["chiral_odd_link_norm"], 1e-12)
        self.assertGreater(wall["link"]["chiral_odd_link_norm"], .01)
        self.assertGreater(mass["link"]["chiral_odd_link_norm"], .01)
        self.assertGreater(base["Pi_sheet_commutator_relative_norm"], .01)
        self.assertLess(base["link"]["action_scalar_pseudoscalar_tensor_coefficient_norm"], 1e-12)

    def test_frame_transport_exposes_false_mass_inference(self):
        green, _ = joined_green(40, ENERGIES[0])
        link = eliminate_center(green)[:4, 4:]
        result = frame_control(link)
        self.assertLess(result["misleading_raw_same_gamma_anticommutator_relative_norm"], 1e-12)
        self.assertLess(result["correct_intertwining_relative_norm"], 1e-12)
        self.assertAlmostEqual(result["common_frame_anticommutator_relative_norm"], 2.)
        self.assertEqual(result["restored_common_frame_relative_error"], 0.)
        self.assertTrue(result["passive_frame_transport_only"])
        self.assertFalse(result["same_Q_plus_outer_walls_are_reflection_invariant"])

    def test_regulator_current_and_complex_green_identity(self):
        result = regulator_checks()
        for key in ("real_energy_cell_current_residual", "complex_energy_cell_Green_identity_residual",
                    "cell_chiral_commutator_residual", "resolvent_Schwarz_reflection_relative_error"):
            self.assertLess(result[key], 1e-12)
        self.assertGreater(result["upper_half_plane_green_imaginary_min_eigenvalue"], 0.)


class RecordTests(unittest.TestCase):
    def test_every_field_comparison_rejects_mutations(self):
        reference = {"map": [{"real": .37, "imag": -.1}], "domain": "transmission", "gate": True}
        for changed in ({**reference, "domain": "wall"}, {**reference, "gate": False},
                        {**reference, "map": [{"real": .38, "imag": -.1}]},
                        {**reference, "extra": 1}):
            with self.assertRaises(RuntimeError):
                compare(reference, changed)

    def test_output_is_exclusive(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"existing.json"
            path.write_bytes(b"untouched")
            completed = subprocess.run([sys.executable, str(ROOT/"scripts/check_nsc_chiral_boundary.py"),
                                        "--output", str(path)], capture_output=True, text=True)
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("refusing to overwrite", completed.stderr)
            self.assertEqual(path.read_bytes(), b"untouched")

    def test_record_has_scoped_positive_gates(self):
        if not DEFAULT_RESULT.exists():
            self.skipTest("result not generated yet")
        record = json.loads(DEFAULT_RESULT.read_bytes())
        self.assertEqual(record["schema"], SCHEMA)
        self.assertTrue(all(record["gate"].values()))
        self.assertFalse(any(record["nonclaims"].values()))
        self.assertEqual(len(record["energy_resolved_maps"]), 3)


if __name__ == "__main__":
    unittest.main()
