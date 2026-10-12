"""Integrity closure and componentwise arithmetic for local-gate evidence.

Integrity is necessary but is not a mathematical proof. These helpers never
assign a physical EXISTENCE or NON-EXISTENCE verdict from unverified numbers.
Proof owners must establish coverage and the supplied continuous remainders.
"""
from hashlib import sha256
from fractions import Fraction
import json
import math
from pathlib import Path, PurePosixPath
import re
import subprocess


ERROR_COMPONENTS = (
    "field_space_time", "changed_history_UV_tail", "baseline_low_subgap",
    "upstream", "energy_interpolation", "covered_regions", "phase_value",
    "between_node", "arithmetic",
)
TOLERANCE = 3e-11
EXACT_TOLERANCE = Fraction(3, 10**11)
_HASH = re.compile(r"[0-9a-f]{64}\Z")
_COMMIT = re.compile(r"[0-9a-f]{40}\Z")


def _path(value):
    if not isinstance(value, str) or not value:
        raise ValueError("nonempty repository-relative evidence path required")
    path = PurePosixPath(value)
    if path.is_absolute() or ".." in path.parts or str(path) != value or "\\" in value:
        raise ValueError("evidence path must be canonical and repository-relative")
    return value


def _reference(path, expected):
    path = _path(path)
    if not isinstance(expected, str) or not _HASH.fullmatch(expected):
        raise ValueError("SHA-256 required for evidence: " + path)
    return path, expected


def declared_references(value):
    """All hash maps and path/hash descriptors, including nested payloads."""
    references = set()

    def visit(item):
        if isinstance(item, dict):
            if "path" in item and "sha256" in item:
                references.add(_reference(item["path"], item["sha256"]))
            for name, child in item.items():
                if name in ("source_hashes", "input_hashes"):
                    if isinstance(child, dict):
                        references.update(_reference(path, expected) for path, expected in child.items())
                    elif isinstance(child, list):
                        # Original NSC operator records use explicit descriptors
                        # rather than a path-to-hash map. Preserve those bytes.
                        if any(not isinstance(v, dict) or "path" not in v or "sha256" not in v
                               for v in child):
                            raise ValueError("evidence hash list requires path/sha256 descriptors")
                        visit(child)
                    else:
                        raise ValueError("evidence hashes must be a map or descriptor list")
                else:
                    visit(child)
        elif isinstance(item, list):
            for child in item:
                visit(child)

    visit(value)
    return references


def build_closure(root, roots, *, historical_commits=()):
    """Expand declared references into a reproducible evidence manifest.

    Roots are repository-relative paths. Current bytes must match every
    declared hash. A mismatch may use ONLY one of the explicitly supplied
    full Git commits, whose identity is saved on the resulting node. This
    does not silently repair a changed payload or search arbitrary history.
    Each (path, hash) is a distinct node, so two historical implementations
    can coexist. The files and Git checkout are never modified.

    This discovers integrity dependencies, not mathematical obligations.
    Unknown record semantics must still be checked by the proof owner.
    """
    root = Path(root).resolve()
    commits = tuple(historical_commits)
    if any(not isinstance(c, str) or not _COMMIT.fullmatch(c) for c in commits):
        raise ValueError("historical evidence requires full pinned commits")
    if not isinstance(roots, (list, tuple)) or not roots:
        raise ValueError("nonempty evidence root paths required")
    nodes, active = {}, set()
    current = {}

    def read_current(name):
        if name not in current:
            path = root / _path(name)
            if path.is_symlink() or not path.resolve().is_relative_to(root):
                raise ValueError("escaping evidence path: " + name)
            if not path.is_file():
                current[name] = None
            else:
                checksum = sha256()
                size = 0
                with path.open("rb") as stream:
                    while block := stream.read(1024 * 1024):
                        checksum.update(block)
                        size += len(block)
                # Only JSON contents are needed for dependency discovery;
                # do not retain large numerical payloads in memory.
                raw = path.read_bytes() if name.endswith(".json") else None
                if raw is not None and sha256(raw).hexdigest() != checksum.hexdigest():
                    raise ValueError("evidence changed while building closure: " + name)
                current[name] = (checksum.hexdigest(), size, raw)
        return current[name]

    def obtain(name, expected):
        found = read_current(name)
        if found is not None and found[0] == expected:
            return found[1], found[2], None
        for commit in commits:
            completed = subprocess.run(
                ["git", "-C", str(root), "show", f"{commit}:{name}"],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=30)
            if completed.returncode == 0 and sha256(completed.stdout).hexdigest() == expected:
                raw = completed.stdout
                return len(raw), raw if name.endswith(".json") else None, commit
        raise ValueError("no matching current or explicitly pinned evidence: "
                         + name + " (sha256 " + expected + ")")

    def visit(name, expected):
        name, expected = _reference(name, expected)
        key = (name, expected)
        identity = sha256((name + "\0" + expected).encode()).hexdigest()
        if key in active:
            raise ValueError("cyclic evidence dependencies")
        if key in nodes:
            return identity
        active.add(key)
        size, raw, commit = obtain(name, expected)
        dependencies = []
        if raw is not None:
            def unique(pairs):
                result = {}
                for k, v in pairs:
                    if k in result:
                        raise ValueError("duplicate evidence JSON key")
                    result[k] = v
                return result
            value = json.loads(raw, object_pairs_hook=unique)
            try:
                references = declared_references(value)
            except ValueError as error:
                raise ValueError("invalid references in " + name + ": " + str(error)) from error
            dependencies = [visit(path, checksum) for path, checksum in sorted(references)]
        node = {"id": identity, "path": name, "sha256": expected,
                "bytes": size, "dependencies": dependencies}
        if commit is not None:
            node["git_commit"] = commit
        nodes[key] = node
        active.remove(key)
        return identity

    root_ids = []
    for name in roots:
        name = _path(name)
        found = read_current(name)
        if found is None:
            raise ValueError("missing current evidence root: " + name)
        identity = visit(name, found[0])
        if identity not in root_ids:
            root_ids.append(identity)
    return {"schema": "NSC-LOCAL-GATE-EVIDENCE-v1", "roots": root_ids,
            "nodes": [nodes[key] for key in sorted(nodes)]}


def verify_closure(root, manifest):
    """Verify a pinned DAG, including dependencies declared inside JSONs.

    Nodes have id, path, sha256, bytes, dependencies, and optional git_commit.
    Distinct historical versions of a path can coexist under distinct node ids.
    Pinned Git sources are read directly; no checkout or historical rewrite occurs.
    """
    root = Path(root).resolve()
    if manifest.get("schema") != "NSC-LOCAL-GATE-EVIDENCE-v1":
        raise ValueError("explicit local-gate evidence schema required")
    nodes = {}
    for node in manifest["nodes"]:
        identity = node["id"]
        if not isinstance(identity, str) or not identity or identity in nodes:
            raise ValueError("distinct nonempty evidence node ids required")
        _reference(node["path"], node["sha256"])
        if type(node["bytes"]) is not int or node["bytes"] < 0:
            raise ValueError("evidence byte count must be nonnegative")
        deps = node["dependencies"]
        if not isinstance(deps, list) or any(not isinstance(v, str) for v in deps) or len(deps) != len(set(deps)):
            raise ValueError("distinct dependency node ids required")
        nodes[identity] = node
    roots = manifest["roots"]
    if not isinstance(roots, list) or not roots or any(v not in nodes for v in roots):
        raise ValueError("existing evidence roots required")
    visiting, verified = set(), {}

    def read(node):
        if node.get("git_commit") is not None:
            commit = node["git_commit"]
            if not isinstance(commit, str) or not _COMMIT.fullmatch(commit):
                raise ValueError("historical evidence requires a full pinned commit")
            raw = subprocess.check_output(
                ["git", "-C", str(root), "show", f"{commit}:{node['path']}"],
                stderr=subprocess.PIPE, timeout=30)
        else:
            path = root / node["path"]
            if path.is_symlink() or not path.resolve().is_relative_to(root) or not path.is_file():
                raise ValueError("missing or escaping evidence file: " + node["path"])
            raw = path.read_bytes()
        if len(raw) != node["bytes"] or sha256(raw).hexdigest() != node["sha256"]:
            raise ValueError("evidence bytes changed: " + node["path"])
        return raw

    def visit(identity):
        if identity in verified:
            return
        if identity not in nodes:
            raise ValueError("missing dependency node: " + identity)
        if identity in visiting:
            raise ValueError("cyclic evidence dependencies")
        visiting.add(identity)
        node = nodes[identity]
        raw = read(node)
        for dep in node["dependencies"]:
            visit(dep)
        if node["path"].endswith(".json"):
            def unique(pairs):
                value = {}
                for key, content in pairs:
                    if key in value:
                        raise ValueError("duplicate evidence JSON key")
                    value[key] = content
                return value
            value = json.loads(raw, object_pairs_hook=unique)
            actual = declared_references(value)
            covered = {(nodes[dep]["path"], nodes[dep]["sha256"]) for dep in node["dependencies"]}
            if not actual.issubset(covered):
                raise ValueError("undeclared transitive evidence: " + repr(sorted(actual-covered)[0]))
        visiting.remove(identity)
        verified[identity] = {"path": node["path"], "sha256": node["sha256"],
                              "bytes": node["bytes"], "git_commit": node.get("git_commit")}

    for identity in roots:
        visit(identity)
    if set(verified) != set(nodes):
        raise ValueError("manifest contains unreachable evidence nodes")
    return {
        "integrity_status": "PASS", "verified_nodes": verified,
        "physical_claim_verified": False,
        "scope": "declared dependency integrity; mathematical proof replay remains required",
    }


def enclosure_arithmetic(nodal_maxima, component_bounds, interval):
    """Evaluate the local sup criterion after proof owners bound each term.

    Every bound is an (N,beta) pair. The between_node component must enclose
    the full residual over I, not merely geometric interpolation or samples.
    """
    def pair(value, name, *, nonnegative=True):
        if not isinstance(value, (list, tuple)) or len(value) != 2:
            raise ValueError(name + " must be an (N,beta) pair")
        if any(isinstance(v, bool) or not isinstance(v, (int, float)) for v in value):
            raise ValueError(name + " must contain real numbers")
        numbers = tuple(float(v) for v in value)
        if any(not math.isfinite(v) or (nonnegative and v < 0) for v in numbers):
            raise ValueError(name + " must be finite and nonnegative")
        return numbers

    endpoints = pair(interval, "interval", nonnegative=False)
    if endpoints[0] >= endpoints[1]:
        raise ValueError("connected interval of positive length required")
    maxima = pair(nodal_maxima, "nodal maxima")
    if set(component_bounds) != set(ERROR_COMPONENTS):
        raise ValueError("all declared error components must be present")
    missing = [name for name in ERROR_COMPONENTS if component_bounds[name] is None]
    known = {name: pair(value, name) for name, value in component_bounds.items() if value is not None}
    if missing:
        upper = None
        within = [False, False]
    else:
        # Sum the binary input bounds exactly, then round outward. Compare
        # with the exact decimal tolerance, not its rounded binary float.
        upper = []
        for i in range(2):
            exact = sum((Fraction.from_float(value) for value in
                         [maxima[i], *(known[name][i] for name in ERROR_COMPONENTS)]),
                        Fraction())
            rounded = float(exact)
            if Fraction.from_float(rounded) < exact:
                rounded = math.nextafter(rounded, math.inf)
            upper.append(rounded)
        within = [Fraction.from_float(value) <= EXACT_TOLERANCE for value in upper]
    return {
        "constraint_order": ["N", "beta"], "interval": list(endpoints),
        "tolerance": TOLERANCE, "nodal_maxima": list(maxima),
        "missing_components": missing, "continuous_upper": upper,
        "within_tolerance": within, "criterion_satisfied": all(within),
        "physical_claim_verified": False,
        "scope": "componentwise arithmetic conditional on verified enclosures and source coverage",
    }
