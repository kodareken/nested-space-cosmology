"""A retrospective prediction stays bound to its saved values and limitations."""
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("prediction_check", ROOT/"scripts/check_prediction_register.py")
checker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)


def old_register_view(root):
    """Drop only an optional finite-amplitude extension and clear its SI flag."""
    path = root/checker.REGISTER
    data = json.loads(path.read_text(encoding="utf-8"))
    entry = data["entries"][0]
    changed = "finite_amplitude" in entry or len(data["entries"]) != 1
    data["entries"] = [entry]
    entry.pop("finite_amplitude", None)
    if entry["observational_test"].get("SI_scale_selected") is not False:
        entry["observational_test"]["SI_scale_selected"] = False
        changed = True
    if changed:
        path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n",
                        encoding="utf-8", newline="\n")


@pytest.fixture
def register_root(tmp_path):
    for name in (checker.REGISTER, checker.DOCUMENT, *checker.EVIDENCE.values()):
        path = tmp_path/name
        path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT/name, path)
    old_register_view(tmp_path)
    return tmp_path


def edit(root, mutate):
    path = root/checker.REGISTER
    data = json.loads(path.read_text(encoding="utf-8"))
    mutate(data["entries"][0])
    path.write_text(json.dumps(data, allow_nan=False), encoding="utf-8")


def test_saved_projection_requires_no_subprocess_or_scientific_replay(register_root, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("prediction check must not run a generator")
    monkeypatch.setattr(subprocess, "run", forbidden)
    before = {str(path):path.read_bytes() for path in register_root.rglob("*") if path.is_file()}
    assert checker.check(register_root) == dict(predictions=1, evidence_files=7,
        coefficient_projections=4, exact_identity_projections=2, scientific_replays=0)
    assert before == {str(path):path.read_bytes() for path in register_root.rglob("*") if path.is_file()}


def test_p003_binds_the_observable_and_preserves_synthetic_scope(tmp_path):
    entry=json.loads((ROOT/checker.REGISTER).read_text(encoding='utf-8'))['entries'][2]
    path=entry['evidence']['path'];target=tmp_path/path
    target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/path,target)
    checker.check_p003(tmp_path,entry)
    bad=copy.deepcopy(entry);bad['prediction']['null_weights'][0]+=.01
    with pytest.raises(ValueError):checker.check_p003(tmp_path,bad)
    bad=copy.deepcopy(entry);bad['observational_validation']=True
    with pytest.raises(ValueError):checker.check_p003(tmp_path,bad)


@pytest.mark.parametrize("mutation", [
    lambda e:e["coefficients"]["rho"].update(value=-5.0),
    lambda e:e["coefficients"]["eta"].update(field_refinement_maximum_absolute_difference=0.),
    lambda e:e["coefficients"]["rho"].update(value=True),
    lambda e:e["prediction_point"].update(pointer="/solutions/0/wall_samples/8"),
    lambda e:e["identities"]["delta_mean_K"].update(value="1"),
    lambda e:e["parameters"]["m"].update(role="theory"),
    lambda e:e["registration"].update(preregistered=True),
    lambda e:e["uncertainty"].update(rigorous_error_enclosure=True),
    lambda e:e["observational_test"].update(status="confirmed"),
    lambda e:e["observational_test"].update(SI_scale_selected=True),
    lambda e:e.update(open_limits=[]),
])
def test_values_provenance_and_scope_cannot_be_promoted(register_root, mutation):
    edit(register_root, mutation)
    with pytest.raises(ValueError):
        checker.check(register_root)


def test_changed_evidence_bytes_fail_before_projection(register_root):
    path = register_root/checker.EVIDENCE["wall_field"]
    path.write_bytes(path.read_bytes()+b"\n")
    with pytest.raises(ValueError, match="evidence bytes changed: wall_field"):
        checker.check(register_root)


def test_rehashed_record_must_still_match_original_manifest(register_root):
    path = register_root/checker.EVIDENCE["cap_wall"]
    data = json.loads(path.read_bytes())
    data["exact_algebra"]["zero_residuals"]["linear_mean_curvature_variation_vanishes"] = "1"
    raw = json.dumps(data).encode(); path.write_bytes(raw)
    edit(register_root, lambda e:e["evidence"]["cap_wall"].update(
        bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest()))
    with pytest.raises(ValueError, match="original manifest"):
        checker.check(register_root)


@pytest.mark.parametrize("name", ["../outside.json", "C:/outside.json", "docs\\model.json"])
def test_evidence_paths_cannot_escape_or_be_retargeted(register_root, name):
    edit(register_root, lambda e:e["evidence"]["model"].update(path=name))
    with pytest.raises(ValueError, match="model path differs"):
        checker.check(register_root)
    with pytest.raises(ValueError, match="path"):
        checker.local_file(register_root, name)


def test_missing_evidence_and_duplicate_prediction_are_not_accepted(register_root):
    path = register_root/checker.REGISTER
    data = json.loads(path.read_bytes())
    data["entries"].append(copy.deepcopy(data["entries"][0]))
    path.write_text(json.dumps(data), encoding="utf-8")
    with pytest.raises(ValueError, match="one currently registered"):
        checker.check(register_root)
    data["entries"].pop(); path.write_text(json.dumps(data), encoding="utf-8")
    (register_root/checker.EVIDENCE["model"]).unlink()
    assert checker.main(["--root", str(register_root)]) == 1


def test_historical_check_strips_an_extension_already_present_in_the_source(register_root, tmp_path):
    source = tmp_path/"source-with-extension"
    for name in (checker.REGISTER, checker.DOCUMENT, *checker.EVIDENCE.values()):
        destination = source/name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(register_root/name, destination)
    path = source/checker.REGISTER
    data = json.loads(path.read_text(encoding="utf-8"))
    original_coefficients = copy.deepcopy(data["entries"][0]["coefficients"])
    data["entries"][0]["finite_amplitude"] = {"schema": "already-installed"}
    data["entries"][0]["observational_test"]["SI_scale_selected"] = True
    path.write_text(json.dumps(data), encoding="utf-8")
    view = tmp_path/"historical-view"
    for name in (checker.REGISTER, checker.DOCUMENT, *checker.EVIDENCE.values()):
        destination = view/name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source/name, destination)
    old_register_view(view)
    saved = json.loads((view/checker.REGISTER).read_text(encoding="utf-8"))["entries"][0]
    assert "finite_amplitude" not in saved
    assert saved["observational_test"]["SI_scale_selected"] is False
    assert saved["observational_test"]["calibration"] == "none"
    assert saved["observational_test"]["physical_source_selected"] is False
    assert saved["coefficients"] == original_coefficients
    assert checker.check(view) == dict(predictions=1, evidence_files=7,
        coefficient_projections=4, exact_identity_projections=2, scientific_replays=0)


def test_duplicate_json_keys_and_nonfinite_coefficients_are_rejected():
    for raw in ('{"value":1,"value":2}', '{"value":NaN}', '{"value":Infinity}'):
        with pytest.raises(ValueError):
            checker.read_json(raw)


@pytest.mark.parametrize("mutation", [
    lambda e:e.update(observational_validation=True),
    lambda e:e.update(Lambda_Lref_squared=[1.125,.5,.126]),
    lambda e:e.update(inputs={}),
    lambda e:e.update(limits="Observed Lambda explained"),
])
def test_charged_result_keeps_inputs_and_conditional_scope(mutation):
    entry=copy.deepcopy(json.loads((ROOT/checker.REGISTER).read_text(encoding="utf-8"))["entries"][1])
    checker.check_p002(ROOT,entry)
    mutation(entry)
    with pytest.raises(ValueError):
        checker.check_p002(ROOT,entry)
