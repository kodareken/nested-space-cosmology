from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import unittest


ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = ROOT / "results" / "manifest.json"


def load_script(name: str):
    path = ROOT / "scripts" / name
    spec = importlib.util.spec_from_file_location(name.removesuffix(".py"), path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class PublicationManifestTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
        cls.steps = cls.manifest["steps"]

    HISTORICAL_OUTPUT_SHA256 = {
        "results/nsc-1-bps-self-equal-energy.json": "2dbc0ee428cfbd0205c41af6f7641723d30f304379742c2280117513d14f29e0",
        "results/nsc-1-exact-black-universe-defocusing.json": "285a3c657ca70547ec3ad92b776f122572725f2ac9772861ec03e0392ab65fe2",
        "results/nsc-1-gravitating-bps-observation-link.json": "04054b5cc6bed8fe49c75922ad39e9306f49d3d158368fbc94885ec8616f9ec3",
        "results/nsc-1-s-one-constant-dictionary.json": "0dbb3d47aa9ef74b6fb898da5441d673a254ecf7ca5c7ccbde4f0c756a6ced1c",
        "results/nsc-1-s-one-nested-pair-closure.json": "494995c310b6b44e82f5724b8d815d81804a6fb4ded4a6e267776e1f1851c654",
        "results/nsc-1-s-one-particle-wave-identity.json": "6dffd874d93d5bcce867af845d205adeb7f7455d7c9020d5c456522c200e31c0",
        "results/nsc-1-s-one-reciprocal-closure.json": "9f52811ad1b85af1d93a20cdff5f8ee00e59bb98114afd7c6870218654b18b59",
        "results/nsc-1-s-one-recursive-outside-kernel.json": "bf1063999658e9d9a43aacc9baa6b221d539466bf09eae164765460f5de20dae",
        "results/nsc-1-unified-gradient-boundary-dev3.json": "f8259c836993c928df1afeba79b7dc9c83dedb3c98f8e24fc926015f32096e93",
        "results/nsc-1-s-one-invariant-closure-bind.json": "77edf049951f12c525c3cde2dbf29ccb9d7d6e3a537204ecc7d445c3ad2ea5ac",
        "results/nsc-1-s-one-relevant-direction-count.json": "9afeabf93a804204b97e5e96a87281608cc87d01ff4fad6e726bb1d7fd42e1d7",
        "results/nsc-1-s-one-outside-black-bind.json": "31948f9c7a1229e314e62931f0d152856e604fa5e967f78baf13d6653033656b",
        "results/nsc-1-s-one-shadow-matter-bind.json": "b6789d6eb7e8609a1592e3e5bdd581e817eee062fc9b5e08a2c9ddfe4b656fcd",
        "results/nsc-1-s-one-shadow-lensing-observation.json": "520cada50d6c1fd7f2b8b5fe176ce5f2b4d43f6f1dd54240507d367e96f840c3",
        "results/nsc-1-s-one-resolution-channel-closure.json": "2d517ae3a957f93706becc33da8844d2d009a8f8f34050b4c0fc2f3fd2fd9406",
        "results/nsc-1-s-one-transition-link-coefficient.json": "798c486efb5f4fde5da16616f9e228693ce131aa137652e4976e1e09909c93e0",
        "results/nsc-1-s-one-warped-resolution-bind.json": "404edddbffa0273f1215ff8225a0f328929dc50b4c0d964ec87851c84a882dbf",
        "results/nsc-1-s-one-spectral-room-bind.json": "ffda536a4a090ae9000f09c7775a2424d1a1ef76de6818e226746317e64766f8",
        "results/nsc-1-s-one-spectral-boundary-constraint.json": "802af6cfdd4f6e084f9fd20513e67d34ea67dc90b639b80aa455e8eb85001a4f",
        "results/nsc-1-s-one-two-sheet-black-geometry.json": "90a9d8b3e8a33edbc5e87b48fff71f91f0c113cff600dc635632264d15aa2778",
        "results/nsc-1-s-one-one-operator-bind.json": "a8a2a22a41f05e1d317f8da49d396d6a15cbb65a161553b76c000dad1250b8e4",
        "results/nsc-1-s-one-two-sheet-5d-geometry.json": "08adb805bf70d37de22e8dadf53886bd0c127f00aa9d052450c6ba49dacc9563",
        "results/nsc-1-s-one-two-sheet-5d-null-source.json": "ee3e09df3d74e78397a061a523359682d7a491d6f37a75eb425ab75bf8cfad57",
        "results/nsc-1-s-one-gauss-bonnet-null-closure.json": "79ad1f29b9be355fe6858119bb00aef8a413f78c8343818ffe7c6d70b2df08ab",
        "results/nsc-1-s-one-gauss-bonnet-component-closure.json": "76ba23b812d2a703b86c625f1408a291e63cbee65d0e85ae57800d445e02b9ea",
        "results/nsc-1-s-one-dark-black-invariant.json": "1e20a7ca0259790359eabdc697cc6025adc69581035b5bf361ae2e6a2f5f333d",
        "results/nsc-1-s-one-gauss-bonnet-criticality.json": "fa69cd5c600dc68d8c0c00811fc21dca60a78a3903d9c4de9a953445125cff0d",
        "results/nsc-1-s-one-gauss-bonnet-kinetic-rank.json": "4558792032532344faa682929aa7d35eb1b08370f4fba62459f52c0be639a1b3",
        "results/nsc-1-s-one-induced-interface-rank.json": "f37cdfd56975c40b5a02679273a7c3c946d51574061ca6ac6d3d3dfabd2d60d9",
        "results/nsc-1-s-one-interface-characteristic-scan.json": "967af0ee446f5cd3c8021062009e4fffde4d1a86fbef645832358759b9bc518a",
        "results/nsc-1-s-one-spectral-resolution-wall-bind.json": "36bb0bb6f990c13a7724cb1b51c80723277228c3cba98dfb70e7cca89bfb6435",
        "results/nsc-1-s-one-spectral-profile-closure.json": "94d6134079caff00d9a8a662c75b836ce48b489324392aaf2c4c3096c094a3be",
        "results/nsc-1-s-one-finite-momentum-graviton.json": "6851dc74c16cab50fc11dd2ee0800498d7e8d48482868972943bac8c415cb11b",
        "results/nsc-1-s-one-flat-spectral-poles.json": "0bd0d95d507263d9632e6a0b47b253ed713de29703fa3a0019351500e53db47c",
        "results/nsc-1-s-one-background-adjusted-graviton.json": "629ad9e11671ecef8719dcfefb3ec9a98022c7f374495f7401451723cd91ad4c",
        "results/nsc-1-s-one-reflection-positivity-obstruction.json": "8dd1ba9c8629504deb963e0ce64b586f64b59fc4ce09f1e0324c45658b4c40bb",
        "results/nsc-1-s-one-fermionic-geometry-bind.json": "8c567f1070c181369995ad76bc6049f417381bcb0de986c887645ebe698ce4d3",
        "results/nsc-1-s-one-invariant-partition-bind.json": "59bf475d3bef3d1e42cfb0beee62ae7cfe8fbb827b8ede25dc3acc5c942cd996",
        "results/nsc-1-s-one-local-two-sheet-anomaly.json": "d3d3bdfe2a5421ef7a522ac32c9fbd0117bca4f395cdb5b1be4288ba26a77ba7",
        "results/nsc-1-s-one-doubled-flrw-closure.json": "6d22ec7d1b876aea1af327861414252b7629370899412bb51326a9643e55f3d0",
        "results/nsc-1-s-one-full-exponential-relative-kernel.json": "1b9cd7106bd0377ed2a80f8ffdb67ddf0228b49a8a29765a7602e9dc1d56aae8",
        "results/nsc-1-s-one-gap-transition-closure.json": "f13159e86bf024961332f59e61a568dda16c5060a6f84b32096f69fc157cb738",
        "results/nsc-1-s-one-child-orientation.json": "82e3b32d49c42f4ba29bc5ddd8df55d67b41617efa72661742641b5279abbb62",
        "results/nsc-1-s-one-boundary-retarded-pole.json": "d45114887d2507ed2a7280446194f0c264867d362e740dc29b3a20661ae82e7a",
        "results/nsc-1-s-one-relative-reflection-test.json": "d12756930ed3ce225d439d9f294955b36c39f14c3f2e66e14049ae20bfe83ade",
        "results/nsc-1-s-one-fermionic-relative-observable.json": "6a5fff0ff2d781ee92b71ca3960be244c2a3a05874a2c7442250c27ffe32c318",
        "results/nsc-1-s-one-child-scale-correction.json": "19c50a14a9902b372a1868c147dd1506570b40a1ac3bd1f1200d03cfe5b8763b",
        "results/nsc-2-zeta1-foliation.json": "39bd1c160c8b72a740f8f26b62c8c202b6f90d357949737981a4e35e1a9d30a1",
        "results/nsc-2-zeta1-self-adjoint-domain.json": "77167b67ef6b9f2260a76a6fcc5782eeb6f4e61507a83992222acf1980825263",
        "results/nsc-2-zeta1-lowest-mode.json": "37e19c2b611cc6ff31f0f866e2ac9a07e4dd0c8d9a21f21c6b1a08b79d7238e0",
        "results/nsc-2-zeta1-angular-tower.json": "1246bee679d33bdefa6e18efddc62bc7f59439e805ea1485847d86211a65e0ad",
        "results/nsc-2-zeta1-regulated-determinant.json": "5dc40ccbc5f5e27cc40dc5c6be574b3713d0bbc91f5a97c74bcb4b10ad27ff10",
        "results/nsc-2-zeta1-y-boundary-sensitivity.json": "970a9e26a8f70cde46a67cb3a48f97255d2ff0b246a7818e9e06601f6078b01e",
        "results/nsc-2-zeta1-orbifold-parity.json": "79871e46cde9ba3e7f6b05df2b20074db0ca07cfbf84cb6b274bf7e743d8c747",
        "results/nsc-2-zeta1-warped-y.json": "243dbdc8d6085815c588ced0280ad05b9b2e7437ce59e4962e560b3945fb90e5",
        "results/nsc-2-zeta1-anomaly-decomposition.json": "052750d49e969a1d82b683f6405a1568aff5f56352183eb9da948a6a2ab0b750",
        "results/nsc-2-zeta1-anomaly-owner-correction.json": "02a40da4985f572fe3c535b6b31329b25da017b85e9ce4fcb7870edd6b13c65c",
        "results/nsc-2-zeta1-recursion-map.json": "4b2c7c3cea91bc980a3de3750e89d32c1b3b105adbdaf7844be1ab70c9565774",
    }

    def test_exact_frontier_and_count(self) -> None:
        self.assertEqual(59, self.manifest["result_count"])
        self.assertEqual(59, len(self.steps))
        self.assertEqual(58, self.manifest["historical_result_count"])
        self.assertEqual(
            "NSC-2-ZETA1-RECURSION-MAP", self.manifest["frontier_artifact_id"]
        )
        self.assertEqual(
            "NSC-2-ZETA1-UNIT-CLOSURE-CHECK",
            self.manifest["follow_up_artifact_id"],
        )
        current = [
            step["artifact_id"]
            for step in self.steps
            if step["category"] == "current_frontier"
        ]
        self.assertEqual(["NSC-2-ZETA1-RECURSION-MAP"], current)

    def test_historical_v010_bytes_are_preserved(self) -> None:
        historical = self.steps[:58]
        self.assertEqual(58, len(historical))
        self.assertEqual(
            "NSC-2-ZETA1-RECURSION-MAP", historical[-1]["artifact_id"]
        )
        observed = {step["output"]: step["output_sha256"] for step in historical}
        self.assertEqual(self.HISTORICAL_OUTPUT_SHA256, observed)
        for step in historical:
            path = ROOT / step["output"]
            self.assertEqual(
                step["output_sha256"],
                hashlib.sha256(path.read_bytes()).hexdigest(),
            )
            generator = ROOT / step["generator"]
            self.assertEqual(
                step["generator_sha256"],
                hashlib.sha256(generator.read_bytes()).hexdigest(),
            )
            self.assertNotIn("auxiliary_inputs", step)

    def test_v010_generator_bytes_match_tagged_snapshot(self) -> None:
        original = json.loads(
            subprocess.check_output(
                ["git", "show", "v0.1.0:results/manifest.json"],
                cwd=ROOT,
            )
        )
        original_steps = {
            step["output"]: step for step in original["steps"]
        }
        self.assertEqual(58, len(original_steps))
        for step in self.steps[:58]:
            previous = original_steps[step["output"]]
            self.assertEqual(previous["output_sha256"], step["output_sha256"])
            self.assertEqual(previous["generator_sha256"], step["generator_sha256"])
            self.assertEqual(previous["generator"], step["generator"])

    def test_unit_closure_follow_up_is_all_fields(self) -> None:
        follow_up = self.steps[-1]
        self.assertEqual("NSC-2-ZETA1-UNIT-CLOSURE-CHECK", follow_up["artifact_id"])
        self.assertEqual("diagnostic_nonpass", follow_up["category"])
        self.assertEqual("all_fields", follow_up["comparison_policy"]["kind"])
        self.assertEqual([], follow_up["dependencies"])
        self.assertTrue(follow_up["auxiliary_inputs"])
        record = json.loads((ROOT / follow_up["output"]).read_text(encoding="utf-8"))
        self.assertIs(
            False, record["nonclaims"]["physical_zeta_or_particle_mass_derived"]
        )
        self.assertIs(
            False,
            record["finite_family_sign_theorem"][
                "continuum_or_general_recursive_theorem"
            ],
        )
        self.assertEqual(
            "82f506316c8b73f47f692eb6573c0a8c5156c98944704fccd3cb46f3e2fd0928",
            follow_up["output_sha256"],
        )

    def test_manifest_is_topological(self) -> None:
        complete: set[str] = set()
        for step in self.steps:
            self.assertLessEqual(set(step["dependencies"]), complete)
            complete.add(step["output"])

    def test_declared_hashes_match_checkout(self) -> None:
        for step in self.steps:
            for path_key, hash_key in (
                ("output", "output_sha256"),
                ("generator", "generator_sha256"),
            ):
                path = ROOT / step[path_key]
                self.assertTrue(path.is_file(), path)
                self.assertEqual(
                    step[hash_key], hashlib.sha256(path.read_bytes()).hexdigest()
                )

    def test_every_result_is_terminal_and_identified(self) -> None:
        for step in self.steps:
            result = json.loads((ROOT / step["output"]).read_text(encoding="utf-8"))
            self.assertIs(True, result["terminal"])
            self.assertEqual(step["artifact_id"], result["artifact_id"])

    def test_determinant_scale_roots_are_not_frontier(self) -> None:
        superseded = {
            step["artifact_id"]
            for step in self.steps
            if step["category"] == "superseded_candidate"
        }
        self.assertTrue(
            {
                "NSC-2-ZETA1-LOWEST-MODE",
                "NSC-2-ZETA1-REGULATED-DETERMINANT",
                "NSC-2-ZETA1-ORBIFOLD-PARITY",
                "NSC-2-ZETA1-WARPED-Y",
            }
            <= superseded
        )

    def test_known_black_universe_is_imported(self) -> None:
        categories = {step["artifact_id"]: step["category"] for step in self.steps}
        self.assertEqual(
            "imported_benchmark_reproduction",
            categories["NSC-1-EXACT-BLACK-UNIVERSE-DEFOCUSING"],
        )

    def test_frontier_preserves_open_scale_claims(self) -> None:
        frontier = json.loads(
            (ROOT / self.manifest["frontier_output"]).read_text(encoding="utf-8")
        )
        self.assertIs(False, frontier["gate"]["mode_resolved_tail_solved"])
        self.assertIs(False, frontier["gate"]["zeta_derived"])
        self.assertIs(False, frontier["nonclaims"]["Omega_value_selected"])
        self.assertIs(False, frontier["nonclaims"]["physical_zeta_promoted"])

    def test_manifest_builder_round_trip(self) -> None:
        builder = load_script("build_result_manifest.py")
        self.assertEqual(self.manifest, builder.build())

    def test_numeric_comparison_policy_propagates_through_dynamic_steps(self) -> None:
        policies = {
            step["artifact_id"]: step["comparison_policy"]["kind"]
            for step in self.steps
        }
        self.assertEqual(
            "portable_numeric",
            policies["NSC-1-S-ONE-FULL-EXPONENTIAL-RELATIVE-KERNEL"],
        )
        self.assertEqual(
            "portable_numeric",
            policies["NSC-1-S-ONE-GAP-TRANSITION-CLOSURE"],
        )
        self.assertEqual(
            "portable_numeric",
            policies["NSC-2-ZETA1-ANOMALY-DECOMPOSITION"],
        )
        self.assertEqual(
            "all_fields",
            policies["NSC-2-ZETA1-UNIT-CLOSURE-CHECK"],
        )

    def test_reproduction_validator_accepts_frozen_checkout(self) -> None:
        reproducer = load_script("reproduce_public_results.py")
        reproducer.validate_checkout(self.manifest)

    def test_all_authenticated_field_variants_are_validated_and_normalized(self) -> None:
        reproducer = load_script("reproduce_public_results.py")
        for name in (
            "nsc-1-s-one-spectral-room-bind.json",
            "nsc-1-s-one-shadow-lensing-observation.json",
            "nsc-1-s-one-resolution-channel-closure.json",
            "nsc-1-s-one-flat-spectral-poles.json",
        ):
            value = json.loads((ROOT / "results" / name).read_text(encoding="utf-8"))
            reproducer.validate_authenticated_inputs(ROOT, value)
            normalized = reproducer.normalize_dynamic_hashes(value)
            entries = list(reproducer.authenticated_entries(normalized))
            self.assertTrue(entries)
            self.assertTrue(
                all(
                    entry["sha256"] == "<validated-generated-result>"
                    for entry in entries
                    if entry["path"].startswith("results/")
                )
            )

    def test_all_fields_comparison_checks_undeclared_numerics(self) -> None:
        reproducer = load_script("reproduce_public_results.py")
        expected = {"classification": "same", "scan": 5.0}
        actual = {"classification": "same", "scan": 9.0}
        with self.assertRaises(reproducer.ReproductionError):
            reproducer.compare_portable(
                expected,
                actual,
                path="",
                relative_tolerance=1e-8,
                absolute_tolerance=1e-8,
                compare_numbers=True,
            )

    def test_portable_comparison_ignores_only_undeclared_numerics(self) -> None:
        reproducer = load_script("reproduce_public_results.py")
        expected = {
            "classification": "same",
            "gate": {"passes": True},
            "scan": 5.0,
            "argmin_index": 2,
        }
        actual = {
            "classification": "same",
            "gate": {"passes": True},
            "scan": 9.0,
            "argmin_index": 3,
        }
        reproducer.compare_portable(
            expected,
            actual,
            path="",
            relative_tolerance=1e-8,
            absolute_tolerance=1e-10,
            compare_numbers=False,
        )
        with self.assertRaises(reproducer.ReproductionError):
            reproducer.compare_portable(
                expected,
                {
                    "classification": "same",
                    "gate": {"passes": False},
                    "scan": 5.0,
                    "argmin_index": 2,
                },
                path="",
                relative_tolerance=1e-8,
                absolute_tolerance=1e-10,
                compare_numbers=False,
            )

    def test_headline_pointer_is_resolved(self) -> None:
        reproducer = load_script("reproduce_public_results.py")
        value = {"outer": {"rows": [{"zeta": 4.75}]}}
        self.assertEqual(4.75, reproducer.resolve_pointer(value, "/outer/rows/0/zeta"))

    def test_headline_observables_are_granular_numeric_paths(self) -> None:
        for step in self.steps:
            for observable in step["headline_observables"]:
                self.assertEqual(
                    {"pointer", "relative_tolerance", "absolute_tolerance"},
                    set(observable),
                )
                value = json.loads(
                    (ROOT / step["output"]).read_text(encoding="utf-8")
                )
                pointer_value = load_script("reproduce_public_results.py").resolve_pointer(
                    value, observable["pointer"]
                )
                self.assertIsInstance(pointer_value, (int, float))


if __name__ == "__main__":
    unittest.main()
