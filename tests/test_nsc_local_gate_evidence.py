"""Transitive integrity and local-gate arithmetic are distinct obligations."""
from hashlib import sha256
from fractions import Fraction
import json
import math
from pathlib import Path
import subprocess

import pytest
from recursive_horizons.nsc_local_gate_evidence import (
    ERROR_COMPONENTS, build_closure, enclosure_arithmetic, verify_closure)


def node(root, name, value, dependencies=()):
    raw = json.dumps(value).encode() if isinstance(value, dict) else value
    (root / name).write_bytes(raw)
    return {"id": name, "path": name, "sha256": sha256(raw).hexdigest(),
            "bytes": len(raw), "dependencies": list(dependencies)}


def chain(root):
    source = node(root, "source.py", b"saved implementation\n")
    payload = node(root, "payload.npz", b"immutable numerical payload")
    child = node(root, "child.json", {
        "source_hashes": {"source.py": source["sha256"]},
        "payload": {"path": "payload.npz", "sha256": payload["sha256"]},
    }, ["source.py", "payload.npz"])
    parent = node(root, "parent.json", {
        "input_hashes": {"child.json": child["sha256"]},
    }, ["child.json"])
    return {"schema": "NSC-LOCAL-GATE-EVIDENCE-v1",
            "roots": ["parent.json"], "nodes": [source, payload, child, parent]}


def test_entire_declared_chain_verified_without_claiming_physics(tmp_path):
    result = verify_closure(tmp_path, chain(tmp_path))
    assert len(result["verified_nodes"]) == 4
    assert result["integrity_status"] == "PASS"
    assert result["physical_claim_verified"] is False


@pytest.mark.parametrize("name", ["source.py", "payload.npz", "child.json"])
def test_changed_descendant_fails_even_when_top_record_unchanged(tmp_path, name):
    manifest = chain(tmp_path)
    (tmp_path / name).write_bytes(b"changed")
    with pytest.raises(ValueError, match="evidence bytes changed"):
        verify_closure(tmp_path, manifest)


def test_missing_and_undeclared_descendants_fail(tmp_path):
    manifest = chain(tmp_path)
    manifest["nodes"][2]["dependencies"] = ["source.py"]
    with pytest.raises(ValueError, match="undeclared transitive"):
        verify_closure(tmp_path, manifest)
    manifest = chain(tmp_path)
    manifest["nodes"] = manifest["nodes"][1:]
    with pytest.raises(ValueError, match="missing dependency"):
        verify_closure(tmp_path, manifest)


def test_cycles_and_unreachable_nodes_fail(tmp_path):
    manifest = chain(tmp_path)
    manifest["nodes"][0]["dependencies"] = ["parent.json"]
    with pytest.raises(ValueError, match="cyclic"):
        verify_closure(tmp_path, manifest)
    manifest = chain(tmp_path)
    manifest["nodes"].append(node(tmp_path, "extra.py", b"unused"))
    with pytest.raises(ValueError, match="unreachable"):
        verify_closure(tmp_path, manifest)


def test_pinned_historical_bytes_do_not_use_modified_worktree(tmp_path):
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    source = node(tmp_path, "source.py", b"version one")
    subprocess.run(["git", "-C", str(tmp_path), "add", "source.py"], check=True)
    subprocess.run(["git", "-C", str(tmp_path), "-c", "user.name=Test",
                    "-c", "user.email=test@example.invalid", "commit", "-qm", "fixture"], check=True)
    source["git_commit"] = subprocess.check_output(
        ["git", "-C", str(tmp_path), "rev-parse", "HEAD"], text=True).strip()
    (tmp_path / "source.py").write_bytes(b"version two")
    manifest = {"schema": "NSC-LOCAL-GATE-EVIDENCE-v1",
                "roots": ["source.py"], "nodes": [source]}
    assert verify_closure(tmp_path, manifest)["integrity_status"] == "PASS"
    assert (tmp_path / "source.py").read_bytes() == b"version two"


def test_missing_bound_and_single_failing_component_block_criterion():
    bounds = {name: [1e-13, 1e-13] for name in ERROR_COMPONENTS}
    result = enclosure_arithmetic([1e-12, 1e-12], bounds, [1.12, 1.18])
    assert result["criterion_satisfied"] is True
    assert result["physical_claim_verified"] is False
    bounds["between_node"] = None
    result = enclosure_arithmetic([0, 0], bounds, [1.12, 1.18])
    assert result["continuous_upper"] is None
    assert result["criterion_satisfied"] is False
    bounds["between_node"] = [1e-13, 4e-11]
    result = enclosure_arithmetic([0, 0], bounds, [1.12, 1.18])
    assert result["within_tolerance"] == [True, False]
    assert result["criterion_satisfied"] is False


@pytest.mark.parametrize("bad", [-1, float("nan"), float("inf"), True])
def test_invalid_bound_rejected(bad):
    bounds = {name: [0, 0] for name in ERROR_COMPONENTS}
    bounds["field_space_time"] = [bad, 0]
    with pytest.raises(ValueError):
        enclosure_arithmetic([0, 0], bounds, [0, 1])


def test_pointwise_or_reversed_interval_rejected():
    bounds = {name: [0, 0] for name in ERROR_COMPONENTS}
    for interval in ([1, 1], [1, 0]):
        with pytest.raises(ValueError, match="positive length"):
            enclosure_arithmetic([0, 0], bounds, interval)


def test_tolerance_uses_exact_decimal_and_sum_is_outward():
    bounds = {name: [0, 0] for name in ERROR_COMPONENTS}
    at_float = 3e-11
    result = enclosure_arithmetic([at_float, 0], bounds, [0, 1])
    assert result["within_tolerance"][0] == (
        Fraction.from_float(at_float) <= Fraction(3, 10**11))
    bounds["field_space_time"] = [2**-100, 0]
    result = enclosure_arithmetic([1e-11, 0], bounds, [0, 1])
    exact = Fraction.from_float(1e-11) + Fraction(1, 2**100)
    assert Fraction.from_float(result["continuous_upper"][0]) >= exact
    assert result["continuous_upper"][0] == math.nextafter(1e-11, math.inf)


def test_discovered_closure_includes_nested_payload_and_source_dependencies(tmp_path):
    chain(tmp_path)
    manifest = build_closure(tmp_path, ["parent.json"])
    assert {n["path"] for n in manifest["nodes"]} == {
        "parent.json", "child.json", "source.py", "payload.npz"}
    assert verify_closure(tmp_path, manifest)["integrity_status"] == "PASS"
    (tmp_path / "payload.npz").write_bytes(b"changed payload")
    with pytest.raises(ValueError, match="no matching current"):
        build_closure(tmp_path, ["parent.json"])


def test_builder_preserves_two_pinned_versions_and_refuses_implicit_git_lookup(tmp_path):
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    old = node(tmp_path, "source.py", b"old method")
    subprocess.run(["git", "-C", str(tmp_path), "add", "source.py"], check=True)
    subprocess.run(["git", "-C", str(tmp_path), "-c", "user.name=Test",
                    "-c", "user.email=test@example.invalid", "commit", "-qm", "old method"], check=True)
    commit = subprocess.check_output(
        ["git", "-C", str(tmp_path), "rev-parse", "HEAD"], text=True).strip()
    new = node(tmp_path, "source.py", b"new method")
    node(tmp_path, "old.json", {"source_hashes": {"source.py": old["sha256"]}})
    node(tmp_path, "new.json", {"source_hashes": {"source.py": new["sha256"]}})
    with pytest.raises(ValueError, match="no matching current"):
        build_closure(tmp_path, ["old.json", "new.json"])
    manifest = build_closure(tmp_path, ["old.json", "new.json"], historical_commits=[commit])
    versions = [n for n in manifest["nodes"] if n["path"] == "source.py"]
    assert len(versions) == 2
    assert {n.get("git_commit") for n in versions} == {None, commit}
    assert verify_closure(tmp_path, manifest)["integrity_status"] == "PASS"
    assert (tmp_path / "source.py").read_bytes() == b"new method"


@pytest.mark.parametrize("bad", ["../outside.json", "/tmp/outside.json", "a/../b.json"])
def test_builder_refuses_noncanonical_paths(tmp_path, bad):
    with pytest.raises(ValueError, match="repository-relative"):
        build_closure(tmp_path, [bad])


def test_builder_rejects_duplicate_keys_and_escaping_symlinks(tmp_path):
    (tmp_path / "duplicate.json").write_text('{"value":1,"value":2}')
    with pytest.raises(ValueError, match="duplicate"):
        build_closure(tmp_path, ["duplicate.json"])
    (tmp_path / "link.json").symlink_to(tmp_path.parent / "outside.json")
    with pytest.raises(ValueError, match="escaping"):
        build_closure(tmp_path, ["link.json"])


def test_historical_descriptor_lists_remain_bound(tmp_path):
    source = node(tmp_path, "old.py", b"historical source")
    node(tmp_path, "old.json", {"source_hashes": [
        {"path": "old.py", "sha256": source["sha256"]}]})
    manifest = build_closure(tmp_path, ["old.json"])
    assert len(manifest["nodes"]) == 2
    assert verify_closure(tmp_path, manifest)["integrity_status"] == "PASS"
    (tmp_path/"old.py").write_bytes(b"changed")
    with pytest.raises(ValueError, match="evidence bytes changed"):
        verify_closure(tmp_path, manifest)
    node(tmp_path, "bad.json", {"input_hashes": ["unbound path"]})
    with pytest.raises(ValueError, match="path/sha256 descriptors"):
        build_closure(tmp_path, ["bad.json"])
