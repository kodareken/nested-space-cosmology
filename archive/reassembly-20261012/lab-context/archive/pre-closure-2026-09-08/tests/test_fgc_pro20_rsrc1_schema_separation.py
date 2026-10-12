"""RSRC1 schema/store separation from the sealed PRO20 namespace."""

from __future__ import annotations

import ast
from copy import deepcopy
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from recursive_horizons.fgc.evolution import pro20_ev1_protocol as closed_protocol
from recursive_horizons.fgc.evolution import pro20_ev1_store as closed_store
from recursive_horizons.fgc.evolution import pro20_rsrc1_protocol as protocol
from recursive_horizons.fgc.evolution import pro20_rsrc1_member as member_owner
from recursive_horizons.fgc.evolution import pro20_rsrc1_seed as seed_owner
from recursive_horizons.fgc.evolution import pro20_rsrc1_store as store
from tests.test_fgc_pro20_rsrc1_protocol import CAMPAIGN_ID, _receipt, _seed_live
from tests.test_fgc_pro20_rsrc1_seed import _bundle, _request


ROOT = Path(__file__).resolve().parents[1]


class _NormalizeSibling(ast.NodeTransformer):
    TEXT = (
        (
            "runs/fgc-2-sf1/pro20-rsrc1-event1/event",
            "runs/fgc-2-sf1/pro20-event1/calibration",
        ),
        ("FGC-1-PRO20-EV1-RSRC1", "FGC-1-PRO20-EV1"),
        ("pro20-rsrc1-event1", "pro20-event1"),
        ("production_rsrc1_event1", "production_event1"),
        ("PRO20 RSRC1", "PRO20-EV1"),
        ("pro20_rsrc1_protocol", "pro20_ev1_protocol"),
        ("PRO20RSRC1", "PRO20EV1"),
    )

    @classmethod
    def text(cls, value: str) -> str:
        for source, target in cls.TEXT:
            value = value.replace(source, target)
        return value

    def visit_Name(self, node: ast.Name):
        return ast.copy_location(ast.Name(self.text(node.id), node.ctx), node)

    def visit_Attribute(self, node: ast.Attribute):
        self.generic_visit(node)
        node.attr = self.text(node.attr)
        return node

    def visit_ClassDef(self, node: ast.ClassDef):
        self.generic_visit(node)
        node.name = self.text(node.name)
        return node

    def visit_FunctionDef(self, node: ast.FunctionDef):
        self.generic_visit(node)
        node.name = self.text(node.name)
        return node

    def visit_ImportFrom(self, node: ast.ImportFrom):
        self.generic_visit(node)
        if node.module is not None:
            node.module = self.text(node.module)
        return node

    def visit_alias(self, node: ast.alias):
        return ast.alias(self.text(node.name), self.text(node.asname) if node.asname else None)

    def visit_Constant(self, node: ast.Constant):
        if isinstance(node.value, str):
            return ast.copy_location(ast.Constant(self.text(node.value)), node)
        return node

    def visit_Assign(self, node: ast.Assign):
        self.generic_visit(node)
        if any(
            isinstance(target, ast.Name) and target.id == "ARTIFACT_ID"
            for target in node.targets
        ):
            node.value = ast.Constant("FGC-1-PRO20-EV1-FRZ1")
        if any(
            isinstance(target, ast.Name) and target.id == "PRODUCTION_STORE_LEAF"
            for target in node.targets
        ):
            node.value = ast.Constant("calibration")
        return node

    def visit_AnnAssign(self, node: ast.AnnAssign):
        self.generic_visit(node)
        if isinstance(node.target, ast.Name) and node.target.id == "ARTIFACT_ID":
            node.value = ast.Constant("FGC-1-PRO20-EV1-FRZ1")
        return node


def _normalized_tree(path: Path) -> str:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    if (
        tree.body
        and isinstance(tree.body[0], ast.Expr)
        and isinstance(tree.body[0].value, ast.Constant)
        and isinstance(tree.body[0].value.value, str)
    ):
        tree.body.pop(0)
    tree = _NormalizeSibling().visit(tree)
    ast.fix_missing_locations(tree)
    return ast.dump(tree, include_attributes=False)


class RSRC1SchemaSeparationTests(unittest.TestCase):
    def test_protocol_and_store_are_mechanical_siblings_except_new_identity(self) -> None:
        for closed, current in (
            (
                ROOT / "src/recursive_horizons/fgc/evolution/pro20_ev1_protocol.py",
                ROOT / "src/recursive_horizons/fgc/evolution/pro20_rsrc1_protocol.py",
            ),
            (
                ROOT / "src/recursive_horizons/fgc/evolution/pro20_ev1_store.py",
                ROOT / "src/recursive_horizons/fgc/evolution/pro20_rsrc1_store.py",
            ),
        ):
            with self.subTest(current=current.name):
                self.assertEqual(_normalized_tree(current), _normalized_tree(closed))

    def test_new_schema_identity_is_exact_and_disjoint(self) -> None:
        self.assertEqual(protocol.PROTOCOL_ARTIFACT_ID, "FGC-2-SF1-PROTO19")
        self.assertEqual(
            protocol.FREEZE_ARTIFACT_ID,
            "FGC-1-PRO20-EV1-RSRC1-FRZ1",
        )
        self.assertEqual(protocol.ARTIFACT_ID, "FGC-1-PRO20-EV1-RSRC1")
        self.assertEqual(store.ARTIFACT_ID, "FGC-1-PRO20-EV1-RSRC1")
        self.assertEqual(
            protocol.PRODUCTION_NAMESPACE,
            "runs/fgc-2-sf1/pro20-rsrc1-event1/event",
        )
        self.assertEqual(protocol.STORE_KIND_PRODUCTION, "production_rsrc1_event1")
        self.assertEqual(
            store.PRODUCTION_CONTAINER,
            "runs/fgc-2-sf1/pro20-rsrc1-event1",
        )
        self.assertNotEqual(protocol.PRODUCTION_NAMESPACE, closed_protocol.PRODUCTION_NAMESPACE)
        self.assertNotEqual(store.PRODUCTION_CONTAINER, closed_store.PRODUCTION_CONTAINER)

    def test_old_and_new_root_schemas_refuse_each_other(self) -> None:
        identities = {
            "campaign_id": CAMPAIGN_ID,
            "authority_sha256": "1" * 64,
            "implementation_sha256": "2" * 64,
            "config_sha256": "3" * 64,
            "source_sha256": "4" * 64,
            "origin_sha256": "5" * 64,
            "environment_sha256": "6" * 64,
        }
        new_root = protocol.build_root_manifest(protocol.build_authority_receipt(**identities))
        old_root = closed_protocol.build_root_manifest(
            closed_protocol.build_authority_receipt(**identities)
        )
        with self.assertRaises(closed_protocol.PRO20EV1ProtocolError):
            closed_protocol.validate_root_manifest(deepcopy(new_root))
        with self.assertRaises(protocol.PRO20RSRC1ProtocolError):
            protocol.validate_root_manifest(deepcopy(old_root))

    def test_new_seed_publish_leaves_an_existing_old_container_byte_untouched(self) -> None:
        _fixture, seed, _checkpoints = _seed_live()
        with TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            parent = root / "runs" / "fgc-2-sf1"
            old = parent / "pro20-event1"
            old.mkdir(parents=True)
            sentinel = old / "immutable-sentinel"
            sentinel.write_bytes(b"closed-pro20-evidence\n")
            before = sentinel.read_bytes()
            published = store.PRO20RSRC1CampaignStore.publish_seed_store(
                root,
                receipt=_receipt(campaign_id=CAMPAIGN_ID),
                bundles=seed.bundles,
            )
            view = published.authenticate()
            self.assertEqual(sentinel.read_bytes(), before)
            self.assertTrue(view.store_root.is_relative_to(parent / "pro20-rsrc1-event1"))
            self.assertFalse(view.terminal)
            self.assertFalse(view.first_event_complete)

    def test_store_consumes_the_complete_rsrc1_seed_cohort_without_translation(self) -> None:
        outcomes = {}
        for index, key in enumerate(protocol.HLT17_MEMBER_KEYS, start=1):
            request = _request(key)
            raw = seed_owner.encode_seed_response(
                request,
                bundle=_bundle(key),
                construction_sha256=f"{index:x}" * 64,
                peak_rss_bytes=1,
                wall_seconds=1.0,
            )
            outcomes[key] = seed_owner.decode_seed_response(request, raw)
        cohort = seed_owner.bind_seed_cohort(outcomes)
        with TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            (root / "runs" / "fgc-2-sf1").mkdir(parents=True)
            published = store.PRO20RSRC1CampaignStore.publish_seed_store(
                root,
                receipt=_receipt(campaign_id=member_owner.CAMPAIGN_ID),
                bundles=cohort.bundles,
            )
            view = published.authenticate()
            self.assertEqual(
                tuple(view.bundles),
                protocol.HLT17_MEMBER_KEYS,
            )
            self.assertEqual(
                view.root_manifest["campaign_id"],
                member_owner.CAMPAIGN_ID,
            )

    def test_new_store_import_graph_never_imports_the_closed_store_or_protocol(self) -> None:
        for relative in (
            "src/recursive_horizons/fgc/evolution/pro20_rsrc1_protocol.py",
            "src/recursive_horizons/fgc/evolution/pro20_rsrc1_store.py",
        ):
            tree = ast.parse((ROOT / relative).read_text(encoding="utf-8"))
            imports: list[str] = []
            for node in ast.walk(tree):
                if isinstance(node, ast.ImportFrom):
                    imports.append(node.module or "")
                elif isinstance(node, ast.Import):
                    imports.extend(alias.name for alias in node.names)
            joined = " ".join(imports)
            self.assertNotIn("pro20_ev1_protocol", joined)
            self.assertNotIn("pro20_ev1_store", joined)


if __name__ == "__main__":
    unittest.main()
