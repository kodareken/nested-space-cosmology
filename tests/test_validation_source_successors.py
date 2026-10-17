"""Source-successor authentication controls; no scientific reconstruction."""
import hashlib
import json
from pathlib import Path
import subprocess

import pytest

import validation_source_successors as tool


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def write(root, name, raw):
    path = root/name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(raw)


@pytest.fixture
def history(tmp_path, monkeypatch):
    old = b"# original validator\n"
    new = b"# reviewed validation-only successor\n"
    name = "scripts/validator.py"
    write(tmp_path, name, old)
    for item in tool.RESOLVER_SOURCES:
        write(tmp_path, item, ("# fixture " + item).encode())
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True, capture_output=True)
    subprocess.run(["git", "-C", str(tmp_path), "add", "."], check=True, capture_output=True)
    subprocess.run(["git", "-C", str(tmp_path), "-c", "user.name=Control",
                    "-c", "user.email=control@example.invalid", "commit", "-qm", "original"],
                   check=True, capture_output=True)
    commit = subprocess.check_output(["git", "-C", str(tmp_path), "rev-parse", "HEAD"]).decode().strip()
    monkeypatch.setattr(tool, "ORIGIN_COMMIT", commit)
    write(tmp_path, name, new)
    registry = dict(schema=tool.SCHEMA, origin_commit=commit,
        resolver_bindings={p:sha((tmp_path/p).read_bytes()) for p in tool.RESOLVER_SOURCES},
        successors={name:dict(original_sha256=sha(old), successor_sha256=sha(new),
            classification="validation-only", reason="Reviewed comparison policy only.")})
    save_registry(tmp_path, registry)
    return tmp_path, name, old, new, registry


def save_registry(root, registry):
    write(root, tool.REGISTRY, json.dumps(registry).encode())


def test_exact_files_need_no_registry_or_git(tmp_path):
    raw = b'{"retained_data":true}\n'
    write(tmp_path, "results/record.json", raw)
    assert tool.verify_saved_source("results/record.json", sha(raw), tmp_path) == raw
    assert tool.verify_saved_bindings({"results/record.json":sha(raw)}, tmp_path,
                                      expected_paths=["results/record.json"]) is True


def test_registered_successor_returns_original_receipt_bytes(history):
    root, name, old, new, registry = history
    assert len(old) != len(new)
    assert tool.verify_saved_source(name, sha(old), root) == old
    assert (root/name).read_bytes() == new
    assert tool.verify_saved_source(name, sha(new), root) == new
    tool.verify_saved_bindings({name:sha(old)}, root, expected_paths=[name])


def test_only_successful_origin_reads_cached_current_and_registry_still_checked(history,monkeypatch):
    root,name,old,new,registry=history
    original=tool.subprocess.run
    calls=[]
    def counted(*args,**kwargs):
        calls.append(args[0])
        return original(*args,**kwargs)
    monkeypatch.setattr(tool.subprocess,"run",counted)
    assert tool.verify_saved_source(name,sha(old),root)==old
    assert tool.verify_saved_source(name,sha(old),root)==old
    assert len(calls)==1
    write(root,name,b"# unreviewed current drift\n")
    with pytest.raises(ValueError,match="unregistered"):
        tool.verify_saved_source(name,sha(old),root)
    write(root,name,new)
    registry["successors"][name]["successor_sha256"]=sha(b"wrong successor")
    save_registry(root,registry)
    with pytest.raises(ValueError,match="unregistered"):
        tool.verify_saved_source(name,sha(old),root)
    registry["successors"][name]["successor_sha256"]=sha(new)
    save_registry(root,registry)
    registry["successors"][name]["original_sha256"]=sha(b"different origin")
    save_registry(root,registry)
    with pytest.raises(ValueError,match="pinned origin"):
        tool.verify_saved_source(name,sha(b"different origin"),root)
    assert len(calls)==2
    registry["successors"][name]["original_sha256"]=sha(old)
    save_registry(root,registry)
    write(root,tool.RESOLVER_SOURCES[0],b"# current resolver drift\n")
    with pytest.raises(ValueError,match="resolver drift"):
        tool.verify_saved_source(name,sha(old),root)
    assert len(calls)==2


@pytest.mark.parametrize("failure",("git_error","wrong_bytes"))
def test_failed_origin_reads_are_not_cached(history,monkeypatch,failure):
    root,name,old,_,_=history
    original=tool.subprocess.run
    calls=[]
    def first_fails(*args,**kwargs):
        calls.append(args[0])
        if len(calls)==1:
            return subprocess.CompletedProcess(args[0],1 if failure=="git_error" else 0,
                                               stdout=b"unverified",stderr=b"fixture")
        return original(*args,**kwargs)
    monkeypatch.setattr(tool.subprocess,"run",first_fails)
    with pytest.raises(ValueError,match="pinned origin"):
        tool.verify_saved_source(name,sha(old),root)
    assert tool.verify_saved_source(name,sha(old),root)==old
    assert tool.verify_saved_source(name,sha(old),root)==old
    assert len(calls)==2


def test_unexpected_current_drift_and_unlisted_source_are_rejected(history):
    root, name, old, _, registry = history
    write(root, name, b"# unreviewed edit\n")
    with pytest.raises(ValueError, match="unregistered"):
        tool.verify_saved_source(name, sha(old), root)
    write(root, "scripts/other.py", b"changed")
    with pytest.raises(ValueError, match="unregistered"):
        tool.verify_saved_source("scripts/other.py", sha(b"old"), root)


@pytest.mark.parametrize("name", ["results/record.json", "docs/scientific-owner.md",
                                   "src/recursive_horizons/scientific_solver.py"])
def test_data_owner_and_unreviewed_package_cannot_use_registry(history, name):
    root, _, old, new, registry = history
    write(root, name, new)
    registry["successors"][name] = dict(original_sha256=sha(old), successor_sha256=sha(new),
        classification="validation-only", reason="Not an allowed source.")
    save_registry(root, registry)
    with pytest.raises(ValueError, match="saved data/source"):
        tool.verify_saved_source(name, sha(old), root)


def test_origin_git_bytes_are_authenticated_not_merely_declared(history):
    root, name, _, _, registry = history
    invented = sha(b"invented original")
    registry["successors"][name]["original_sha256"] = invented
    save_registry(root, registry)
    with pytest.raises(ValueError, match="pinned origin"):
        tool.verify_saved_source(name, invented, root)


def test_resolver_drift_inventory_and_reason_are_checked(history):
    root, name, old, _, registry = history
    with pytest.raises(ValueError, match="inventory"):
        tool.verify_saved_bindings({name:sha(old)}, root, expected_paths=[name,"tests/missing.py"])
    registry["successors"][name]["reason"] = ""
    save_registry(root, registry)
    with pytest.raises(ValueError, match="entry"):
        tool.verify_saved_source(name, sha(old), root)
    registry["successors"][name]["reason"] = "reviewed"
    save_registry(root, registry)
    write(root, tool.RESOLVER_SOURCES[0], b"unreviewed resolver")
    with pytest.raises(ValueError, match="resolver drift"):
        tool.verify_saved_source(name, sha(old), root)


@pytest.mark.parametrize("name", ["../outside.py", "scripts/../outside.py", "/absolute.py",
                                   "C:/absolute.py", "scripts\\validator.py", ".git/config"])
def test_unsafe_paths_rejected(tmp_path, name):
    with pytest.raises(ValueError, match="safe repository-relative"):
        tool.verify_saved_source(name, "a"*64, tmp_path)


def test_only_explicit_package_path_is_eligible():
    assert tool._source_path("src/recursive_horizons/nsc_atlas_characteristic_trace.py")
    assert not tool._source_path("src/recursive_horizons/nsc_atlas_characteristic_initialization.py")
    assert tool._source_path("src/recursive_horizons/nsc_atlas_refined_cache.py")
    assert not tool._source_path("src/recursive_horizons/nsc_atlas_nonlinear_source.py")
    assert tool._source_path("src/recursive_horizons/nsc_atlas_cap_response.py")
    assert tool._source_path("src/recursive_horizons/nsc_atlas_connection_slice.py")
    assert not tool._source_path("src/recursive_horizons/nsc_atlas_connection_observer.py")
    assert not tool._source_path("scripts/data.json")
    for name in ("nsc_atlas_cap_field", "nsc_atlas_background_profile", "nsc_atlas_prepared_response",
                 "nsc_atlas_characteristic_background"):
        assert tool._source_path("src/recursive_horizons/"+name+".py")


def test_implementation_successor_keeps_original_calculation_bytes(history):
    root,_,old,new,registry=history
    name="src/recursive_horizons/nsc_atlas_cap_field.py"
    write(root,name,old)
    subprocess.run(["git","-C",str(root),"add",name],check=True,capture_output=True)
    subprocess.run(["git","-C",str(root),"-c","user.name=Control",
                    "-c","user.email=control@example.invalid","commit","--amend","--no-edit","-q"],
                   check=True,capture_output=True)
    commit=subprocess.check_output(["git","-C",str(root),"rev-parse","HEAD"]).decode().strip()
    # The test fixture's origin now contains this original package source.
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(tool,"ORIGIN_COMMIT",commit)
        registry['origin_commit']=commit
        registry['successors'][name]=dict(original_sha256=sha(old),successor_sha256=sha(new),
            classification='implementation-successor',reason='Changed calculation, original evidence remains at its source.')
        write(root,name,new);save_registry(root,registry)
        assert tool.verify_saved_source(name,sha(old),root)==old
        assert tool.verify_saved_source(name,sha(new),root)==new
        write(root,name,b'unreviewed further implementation')
        with pytest.raises(ValueError,match='unregistered'):
            tool.verify_saved_source(name,sha(old),root)
        write(root,name,new)
        registry['successors']['scripts/validator.py']['classification']='implementation-successor'
        save_registry(root,registry)
        with pytest.raises(ValueError,match='entry'):
            tool.verify_saved_source(name,sha(old),root)


def test_active_cap_replay_policy_keeps_counts_structure_and_scientific_changes_strict():
    import diagnose_nsc_atlas_boundary_rows as boundary
    import diagnose_nsc_atlas_cap_consistency as consistency
    import diagnose_nsc_atlas_cap_mass_direction as direction
    for checker in (boundary.same_recomputed,consistency.same_recomputed,direction.same_recomputed):
        expected={"residual":[-.661250447],"calls":20,"quiet":True,"scope":"fixed"}
        assert checker(dict(expected,residual=[-.661250447+1e-14]),expected)
        assert not checker(dict(expected,residual=[-.661250447+1e-6]),expected)
        assert not checker(dict(expected,calls=20.0),expected)
        assert not checker(dict(expected,quiet=1),expected)
        assert not checker(dict(expected,scope="changed"),expected)


def test_later_source_origin_is_explicit_immutable_and_authenticated(history):
    root,_,_,_,registry=history
    name='src/recursive_horizons/nsc_atlas_cap_mass_derivatives.py'
    old=b'# later original derivative\n';new=b'# reviewed degree-seven derivative\n'
    write(root,name,old)
    subprocess.run(['git','-C',str(root),'add',name],check=True,capture_output=True)
    subprocess.run(['git','-C',str(root),'-c','user.name=Control','-c','user.email=control@example.invalid',
                    'commit','-qm','later source'],check=True,capture_output=True)
    commit=subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD']).decode().strip()
    registry['successors'][name]=dict(original_sha256=sha(old),successor_sha256=sha(new),
        classification='implementation-successor',reason='Later numerical representation.',original_commit=commit)
    write(root,name,new);save_registry(root,registry)
    assert tool.verify_saved_source(name,sha(old),root)==old
    registry['successors'][name].pop('original_commit');save_registry(root,registry)
    with pytest.raises(ValueError,match='pinned origin'):
        tool.verify_saved_source(name,sha(old),root)
    for invalid in ('HEAD',commit[:7],None,123):
        registry['successors'][name]['original_commit']=invalid;save_registry(root,registry)
        with pytest.raises(ValueError,match='entry'): tool.verify_saved_source(name,sha(old),root)
    registry['successors'][name]['original_commit']=commit
    registry['successors'][name]['original_sha256']=sha(b'invented original');save_registry(root,registry)
    with pytest.raises(ValueError,match='pinned origin'):
        tool.verify_saved_source(name,sha(b'invented original'),root)


@pytest.mark.parametrize('raw_available',[False,True])
def test_compact_summary_preserves_reviewed_original_with_or_without_raw(history,raw_available):
    import summarize_result as summary
    root,name,old,new,registry=history
    raw=root/'.build/results/raw.json';raw.parent.mkdir(parents=True)
    raw.write_bytes(summary.encoded({'samples':list(range(32))}))
    write(root,name,old)
    output=root/'results/summaries/control.json'
    summary.create_summary(raw,output,argv=['never-executed'],producer_paths=[name],root=root)
    saved=output.read_bytes();write(root,name,new)
    if not raw_available: raw.unlink()  # disposable fixture only
    receipt=summary.check_summary(output,root=root)
    assert receipt['raw_available']==raw_available and receipt['scientific_replays']==0
    assert receipt['producer_context']=='reviewed-successors'
    assert receipt['original_producer_paths']==[name] and output.read_bytes()==saved
    write(root,name,b'unreviewed change')
    with pytest.raises(ValueError,match='unregistered'): summary.check_summary(output,root=root)
    write(root,name,new)
    changed=json.loads(saved);changed['producer_bindings'][name]['bytes']+=1
    output.write_bytes(summary.encoded(changed))
    with pytest.raises(ValueError,match='byte count'): summary.check_summary(output,root=root)


def test_connection_slice_accepts_only_the_reviewed_hash_pair():
    root = Path(__file__).resolve().parents[1]
    name = "src/recursive_horizons/nsc_atlas_connection_slice.py"
    original = "ef7870830057258c762d6679168e03c34bcd24ab300e29f803aae9ea1828b832"
    current = (root / name).read_bytes()
    authenticated = tool.verify_saved_source(name, original, root)
    assert sha(authenticated) == original and authenticated != current
    assert tool.verify_saved_source(name, sha(current), root) == current
    with pytest.raises(ValueError, match="unregistered"):
        tool.verify_saved_source(name, sha(b"unreviewed connection slice"), root)
