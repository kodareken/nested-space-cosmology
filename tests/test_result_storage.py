"""Generic compact-result and legacy-byte storage controls; no science runs."""
import copy
import json
from pathlib import Path

import pytest

import check_result_storage as storage
import summarize_result as summary


def fixture(tmp_path):
    raw=tmp_path/".build/results/raw.json"
    raw.parent.mkdir(parents=True)
    producer=tmp_path/"scripts/producer.py"
    producer.parent.mkdir(); producer.write_bytes(b"# synthetic producer\n")
    dependency=tmp_path/"inputs/model.json"
    dependency.parent.mkdir(); dependency.write_bytes(b'{"model":"fixed"}\n')
    payload=dict(schema="synthetic-science-v1",runtime={"python":"fixed"},
        execution={"evaluated_points":42},input={"sha256":"a"*64,"path":"inputs/model.json"},
        audit_source_bindings={"scripts/producer.py":summary.sha256(producer.read_bytes())},
        scientific=dict(values=list(range(40)),small=[1,2,3],
            rows=[dict(name="important label",samples=list(range(20))),dict(name="retained caveat",samples=[4,5])],
            status="conditional; no remainder bound"))
    raw.write_bytes(summary.encoded(payload))
    return raw,payload,producer,dependency


def create(tmp_path,**kwargs):
    raw,payload,producer,dependency=fixture(tmp_path)
    output=tmp_path/"results/summaries/one.json"
    row=summary.create_summary(raw,output,root=tmp_path,
        argv=["python","scripts/producer.py","--output",".build/results/raw.json"],
        producer_paths=[producer],input_paths=[dependency],**kwargs)
    return raw,payload,producer,dependency,output,row


def test_arrays_collapsed_without_losing_structure_strings_or_provenance(tmp_path):
    raw,payload,producer,dependency,output,row=create(tmp_path)
    values=row["data"]["scientific"]["values"]
    assert values==dict(compacted_numeric_array=True,length=40,scalar_count=40,
        sha256=summary.sha256(summary.encoded(list(range(40)))),minimum=0,maximum=39)
    assert row["data"]["scientific"]["small"]==[1,2,3]
    assert row["data"]["scientific"]["rows"][0]["name"]=="important label"
    assert row["data"]["scientific"]["rows"][1]==payload["scientific"]["rows"][1]
    assert row["reported_provenance"]["runtime"]==payload["runtime"]
    assert row["reported_provenance"]["input"]==payload["input"]
    assert row["raw_payload"]==dict(path=".build/results/raw.json",bytes=raw.stat().st_size,
        sha256=summary.sha256(raw.read_bytes()))
    assert row["regeneration"]["argv"]==["python","scripts/producer.py","--output",".build/results/raw.json"]
    assert not row["authentication"]["scientific_validation_executed"]
    assert summary.check_summary(output,root=tmp_path)["raw_available"]
    before=output.read_bytes()
    with pytest.raises(FileExistsError):
        summary.create_summary(raw,output,root=tmp_path,argv=["never"],producer_paths=[producer])
    assert output.read_bytes()==before


def test_json_pointer_selection_retains_provenance_and_nested_numeric_hash(tmp_path):
    raw,payload,producer,dependency,output,row=create(tmp_path,selection="/scientific/rows/0")
    assert row["data"]["name"]=="important label"
    assert "runtime" in row["reported_provenance"]
    assert summary.select({"a/b":{"~name":3}},"/a~1b/~0name")==3
    nested=[[i for i in range(8)] for _ in range(3)]
    compact=summary.compact(nested)
    assert compact["length"]==3 and compact["scalar_count"]==24
    assert compact["sha256"]==summary.sha256(summary.encoded(nested))
    for pointer in ("missing","/no", "/scientific/rows/01", "/a~2b"):
        with pytest.raises(ValueError): summary.select(payload,pointer)


@pytest.mark.parametrize("kind",["raw","producer","dependency","data"])
def test_byte_or_summary_tampering_fails(tmp_path,kind):
    raw,payload,producer,dependency,output,row=create(tmp_path)
    if kind=="raw": raw.write_bytes(raw.read_bytes()+b" ")
    elif kind=="producer": producer.write_bytes(b"changed")
    elif kind=="dependency": dependency.write_bytes(b"changed")
    else:
        row["data"]["scientific"]["values"]["maximum"]=41
        output.write_bytes(summary.encoded(row))
    with pytest.raises(ValueError): summary.check_summary(output,root=tmp_path)


def test_absent_ignored_raw_reported_without_regeneration(tmp_path,monkeypatch):
    raw,payload,producer,dependency,output,row=create(tmp_path)
    raw.unlink()  # synthetic disposable fixture only
    monkeypatch.setattr(storage.subprocess,"run",lambda *a,**k:pytest.fail("no command may be executed"))
    receipt=summary.check_summary(output,root=tmp_path)
    assert not receipt["raw_available"] and receipt["scientific_replays"]==0
    assert output.exists()


def test_path_finite_json_and_summary_size_guards(tmp_path,monkeypatch):
    raw,payload,producer,dependency=fixture(tmp_path)
    for output in (tmp_path.parent/"outside.json","../outside.json",".build/summary.json"):
        with pytest.raises(ValueError):
            summary.create_summary(raw,output,root=tmp_path,argv=["python"],producer_paths=[producer])
    for body in (b'{"x":NaN}',b'{"x":1,"x":2}',b'{"x":1e999}'):
        with pytest.raises(ValueError): summary.load_json(body)
    monkeypatch.setattr(summary,"MAX_SUMMARY_BYTES",20)
    with pytest.raises(ValueError,match="exceeds"):
        summary.create_summary(raw,"results/summary.json",root=tmp_path,argv=["python"],producer_paths=[producer])
    assert not (tmp_path/"results/summary.json").exists()


def legacy(tmp_path,monkeypatch):
    folder=tmp_path/"results"; folder.mkdir()
    old=folder/"old.json"; old.write_bytes(b"["+b"0,"*(storage.LIMIT//2)+b"0]")
    small=folder/"small.json"; small.write_bytes(b"{}")
    untracked=folder/"not-tracked.json"; untracked.write_bytes(old.read_bytes())
    monkeypatch.setattr(storage,"tracked_results",lambda root:["results/old.json","results/small.json"])
    baseline=storage.create_baseline(root=tmp_path)
    assert set(baseline["entries"])=={"results/old.json"}
    return old,small,untracked,baseline


def test_new_large_payload_blocked_and_only_tracked_legacy_baselined(tmp_path,monkeypatch):
    old,small,untracked,baseline=legacy(tmp_path,monkeypatch)
    with pytest.raises(ValueError,match="new oversized") as error: storage.check(tmp_path)
    assert error.value.violations==[dict(path="results/not-tracked.json",bytes=untracked.stat().st_size,reason="new oversized result payload")]
    untracked.unlink()
    receipt=storage.check(tmp_path)
    assert receipt["unchanged_legacy_payloads"]==1
    with pytest.raises(FileExistsError): storage.create_baseline(root=tmp_path)


@pytest.mark.parametrize("new_body",[b"{}",None])
def test_changed_legacy_blocked_even_if_smaller(tmp_path,monkeypatch,new_body):
    old,small,untracked,baseline=legacy(tmp_path,monkeypatch)
    untracked.unlink()
    old.write_bytes(new_body if new_body is not None else old.read_bytes().replace(b"0",b"1",1))
    with pytest.raises(ValueError,match="changed legacy"): storage.check(tmp_path)


def test_ignored_build_payloads_and_baseline_path_guard(tmp_path,monkeypatch):
    old,small,untracked,baseline=legacy(tmp_path,monkeypatch)
    untracked.unlink()
    scratch=tmp_path/".build/results/raw.json"; scratch.parent.mkdir(parents=True)
    scratch.write_bytes(old.read_bytes())
    assert storage.check(tmp_path)["checked_json"]==2
    path=tmp_path/storage.BASELINE
    corrupted=copy.deepcopy(baseline)
    corrupted["entries"]["../escape.json"]=corrupted["entries"].pop("results/old.json")
    path.write_bytes(summary.encoded(corrupted))
    with pytest.raises(ValueError,match="traversal"): storage.check(tmp_path)


@pytest.mark.parametrize("name", ["results/.build/escaped.json", "results/escaped.JSON"])
def test_nested_scratch_and_uppercase_json_do_not_evade_storage_guard(tmp_path,monkeypatch,name):
    old,small,untracked,baseline=legacy(tmp_path,monkeypatch)
    untracked.unlink()
    target=tmp_path/name; target.parent.mkdir(parents=True,exist_ok=True)
    target.write_bytes(old.read_bytes())
    with pytest.raises(ValueError,match="new oversized") as error:
        storage.check(tmp_path)
    assert error.value.violations[0]["path"]==name


@pytest.mark.parametrize('stem',[
    'nsc-atlas-cap-mass-direction-v1','nsc-atlas-cap-mass-balance-v1',
    'nsc-atlas-cap-transport-balance-v1',
])
def test_tracked_summary_authenticates_when_ignored_raw_is_absent(monkeypatch,stem):
    path=summary.ROOT/'results/summaries'/str(stem+'.json')
    before=path.read_bytes();saved=summary.load_json(before)
    raw=summary.ROOT/saved['raw_payload']['path']
    # Model a clean checkout without moving or deleting actual local evidence.
    is_file,exists=Path.is_file,Path.exists
    monkeypatch.setattr(Path,'is_file',lambda self:False if self==raw else is_file(self))
    monkeypatch.setattr(Path,'exists',lambda self:False if self==raw else exists(self))
    receipt=summary.check_summary(path)
    assert not receipt['raw_available'] and receipt['scientific_replays']==0
    assert receipt['producer_context']=='reviewed-successors'
    assert path.read_bytes()==before
