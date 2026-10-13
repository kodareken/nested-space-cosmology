"""Fresh owners and explicit saved-source selection remain distinct."""
import importlib.util
from pathlib import Path

import pytest

from recursive_horizons.nsc_ks_evaluation_binding import implementation_hashes

ROOT = Path(__file__).resolve().parents[1]
ORIGINAL = "329ef809e466a86ed4a290ee32c19485f85add9d"


def driver(stem):
    spec = importlib.util.spec_from_file_location(stem, ROOT / "scripts" / ("derive_nsc_" + stem + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("stem", ["subgap_upstream_covariance", "vacuum_source_correction",
                                 "vacuum_source_remainder", "vacuum_spinor_remainder"])
def test_fresh_owners_hash_actual_retained_files(stem):
    module = driver(stem)
    actual = implementation_hashes(ROOT, owners=module.OWNERS)
    assert "docs/nsc-ks-source-operator-majorant-v7.md" in actual
    assert module.ORIGINAL_OWNERS[-1] not in actual
    assert all((ROOT / name).is_file() for name in actual)
    old = implementation_hashes(ROOT, owners=module.ORIGINAL_OWNERS, source_context=ORIGINAL)
    assert module.ORIGINAL_OWNERS[-1] in old
    assert old[module.ORIGINAL_OWNERS[0]] != actual[module.ORIGINAL_OWNERS[0]]


def test_rank_fresh_record_has_current_owner_and_saved_replay_retains_old_identity():
    module = driver("upstream_compensation_rank")
    fresh = module.build_record()
    old = module.build_record(source_context=ORIGINAL)
    assert "docs/nsc-preparation-domain-decision.md" in fresh["input_hashes"]
    assert "docs/nsc-upstream-metric-compensation.md" not in fresh["input_hashes"]
    assert "docs/nsc-upstream-metric-compensation.md" in old["input_hashes"]
    for key in ("bounds", "identities", "selected_input", "scope"):
        assert fresh[key] == old[key]


def test_original_parent_selection_does_not_replace_missing_numerical_inputs(tmp_path, monkeypatch):
    module = driver("subgap_upstream_covariance")
    source = ROOT / module.SOURCE
    target = tmp_path / module.SOURCE
    target.parent.mkdir(parents=True)
    target.write_bytes(source.read_bytes())
    monkeypatch.setattr(module, "ROOT", tmp_path)
    # Isolate the source authenticator; input/payload checks still read only
    # the explicit current root, even when the original source is selected.
    monkeypatch.setattr(module, "source_sha256",
                        lambda root, name, context: implementation_hashes(
                            ROOT, owners=(name,), source_context=context)[name])
    with pytest.raises(FileNotFoundError):
        module.authenticated_parent(module.SOURCE, source_context=ORIGINAL)
