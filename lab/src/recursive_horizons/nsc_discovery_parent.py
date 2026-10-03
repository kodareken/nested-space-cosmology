"""Global parent initial data on the unchanged leading Einstein/Dirac action.

The sealed collar at w=0.0013 is the interior seed. This module builds one
periodic slice: band-limited geometry, a rank-2 source, and both retained
momenta. It does not step a trajectory, add a force, or write a campaign.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import hashlib
import io
import json
import os
import resource
from pathlib import Path
import subprocess
import time

import numpy as np
from scipy.integrate import solve_ivp
from threadpoolctl import threadpool_limits

from . import nsc_discovery_backend as backend
from . import nsc_discovery_extent as extent
from . import nsc_discovery_leading_einstein as leading
from . import nsc_discovery_parent_traction as collar
from . import nsc_discovery_parent_step_control as step_control
from . import nsc_spherical_coupling as coupling
from . import nsc_spherical_galerkin_coupling as galerkin
from . import provenance


SCHEMA = "NSC-DISCOVERY-PARENT-v1"
LAB = Path(__file__).resolve().parents[2]
ROOT = LAB.parent
OUTPUT = LAB / "results/development/nsc-discovery-parent-v1"
SEALED_COMMIT = "b7aa5f5e0f87607e16c22bd8b5955993ce42608d"
SEALED_PAYLOAD = "e7d306105830662ff83f8bbbc4bcb363a49a15a9734a8d6ab0e05b55e669a5b8"
SEED_WEIGHT = 0.0013
LENGTH = 8.0
CENTER = LENGTH / 2.0
CHILD_RADIUS = 0.5
COLLAR_RADIUS = 1.0
PARENT_RADIUS = 3.0
TRANSITION_END = 1.2
ENVELOPE_ZERO = 2.8
EXTERIOR_Q = 2.0
ANNULUS = (1.2, 3.0)
POPULATIONS = (0, 1, 2)
SIGNS = (1, -1)
IMBALANCE_CAP = 0.5
K_MARGIN = 0.05
MAX_ACCEPTED = 12
MAX_BACKTRACKS = 8
FALLBACK_STEPS = 8
CPU_LIMIT = 30.0
CHUNK_LIMIT = 64 * 1024 ** 2
CHILD_INTERVAL = (CENTER - CHILD_RADIUS, CENTER + CHILD_RADIUS)
PARENT_INTERVAL = (CENTER - PARENT_RADIUS, CENTER + PARENT_RADIUS)
CLOCK_LOCATIONS = (1.0, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0, 5.2, 7.0)
OWNERS = (
    "lab/src/recursive_horizons/nsc_discovery_parent.py",
    "lab/scripts/derive_nsc_discovery_parent.py",
    "lab/tests/test_nsc_discovery_parent.py",
    "lab/docs/nsc-discovery-parent.md",
    "lab/src/recursive_horizons/nsc_discovery_leading_einstein.py",
    "lab/src/recursive_horizons/nsc_discovery_dynamic_preparation.py",
    "lab/src/recursive_horizons/nsc_discovery_parent_traction.py",
    "lab/src/recursive_horizons/nsc_discovery_parent_step_control.py",
    "lab/src/recursive_horizons/nsc_discovery_leading_step_control.py",
    "lab/src/recursive_horizons/nsc_discovery_backend.py",
    "lab/src/recursive_horizons/nsc_discovery_extent.py",
    "lab/src/recursive_horizons/nsc_discovery_tidal.py",
    "lab/src/recursive_horizons/nsc_spherical_coupling.py",
    "lab/src/recursive_horizons/nsc_spherical_galerkin_coupling.py",
    "lab/src/recursive_horizons/provenance.py",
)

_CACHE = {}


@dataclass(frozen=True)
class ParentPair:
    """Parent holder. W may be the identity; every geometry mode stays active."""

    grid: galerkin.GalerkinGrid
    geometry_map: np.ndarray
    weights: np.ndarray
    reference_columns: np.ndarray
    source_columns: np.ndarray
    source_metadata: dict
    geometry_metadata: dict
    child_interval: tuple = CHILD_INTERVAL
    parent_interval: tuple = PARENT_INTERVAL
    clock_locations: tuple = CLOCK_LOCATIONS

    @property
    def source_phi0(self):return self.source_columns[:self.grid.nf]

    @property
    def source_phi1(self):return self.source_columns[self.grid.nf:]


def _cpu_time():
    """Parent plus completed authentication subprocess CPU; no scientific pool."""
    children=resource.getrusage(resource.RUSAGE_CHILDREN)
    return time.process_time()+children.ru_utime+children.ru_stime


def _sha_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _sha_array(array):
    values = np.ascontiguousarray(array)
    digest = hashlib.sha256(str(values.dtype).encode())
    digest.update(str(values.shape).encode())
    digest.update(values.tobytes())
    return digest.hexdigest()


def source_hashes():
    return dict(leading.source_hashes(),**{path: _sha_file(ROOT / path) for path in OWNERS})


def _plain(value):
    if isinstance(value, dict):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if isinstance(value, np.ndarray):
        return _plain(value.tolist())
    if isinstance(value, (np.floating, float)):
        number = float(value)
        if not np.isfinite(number):
            raise ValueError("nonfinite report value")
        return number
    if isinstance(value, (np.integer,)) or (isinstance(value, int) and not isinstance(value, bool)):
        return int(value)
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if value is None or isinstance(value, str):
        return value
    raise TypeError("report value is not JSON data: " + type(value).__name__)


def _maximum(values):
    return float(np.max(np.abs(values))) if np.size(values) else 0.0


def _summary(values):
    array = np.asarray(values, float)
    return {"min": float(np.min(array)), "max": float(np.max(array)), "max_abs": _maximum(array)}


def signed_distance(x, length=LENGTH):
    """Signed distance from L/2, in (-L/2, L/2]. Centre is 4 when L is 8."""
    return (np.asarray(x, float) % length) - length / 2.0


def _profile(profile):
    selected = {
        "transition_end": TRANSITION_END,
        "envelope_zero": ENVELOPE_ZERO,
        "exterior_Q": EXTERIOR_Q,
        "k_margin": K_MARGIN,
        "exterior_radius": None,
    }
    if profile is not None:
        if not isinstance(profile, dict):
            raise ValueError("profile must be a mapping or None")
        unknown = set(profile) - set(selected)
        if unknown:
            raise ValueError("unknown parent profile fields: " + ",".join(sorted(unknown)))
        selected.update(profile)
    end = float(selected["transition_end"])
    zero = float(selected["envelope_zero"])
    if not COLLAR_RADIUS < end < 0.5 * np.pi:
        raise ValueError("geometry transition must end after the collar and before pi/2")
    if not end < zero < PARENT_RADIUS:
        raise ValueError("field envelope must stay 1 through the transition and reach 0 before 3")
    if float(selected["exterior_Q"]) <= 0 or float(selected["k_margin"]) < 0:
        raise ValueError("exterior lapse and k margin must be positive")
    if selected["exterior_radius"] is not None and float(selected["exterior_radius"]) <= 0:
        raise ValueError("exterior radius must be positive")
    selected["transition_end"] = end
    selected["envelope_zero"] = zero
    selected["exterior_Q"] = float(selected["exterior_Q"])
    selected["k_margin"] = float(selected["k_margin"])
    return selected


def _require_case(nf, population, sign, cpu_limit):
    if int(nf) != nf or int(nf) < 32 or int(nf) % 2:
        raise ValueError("parent preparation needs an even nf of at least 32")
    if int(population) not in POPULATIONS:
        raise ValueError("population must be 0, 1 or 2")
    if int(sign) not in SIGNS:
        raise ValueError("momentum sign must be +1 or -1")
    if not 0 < float(cpu_limit) <= CPU_LIMIT:
        raise ValueError("parent CPU limit must be in (0, 30]")
    return int(nf), int(population), int(sign), float(cpu_limit)


def authenticate_collar():
    """Bind the sealed fine collar payload, its producer commit and input closure."""
    directory = collar.OUTPUT
    record_path = directory / "measurement.json"
    payload_path = directory / "measurement.npz"
    record = json.loads(record_path.read_text())
    payload_hash = _sha_file(payload_path)
    if record.get("schema") != collar.SCHEMA or record.get("producing_commit") != SEALED_COMMIT:
        raise ValueError("collar measurement is not the sealed b7aa5f5 record")
    if payload_hash != record["payload_sha256"] or payload_hash != SEALED_PAYLOAD:
        raise ValueError("sealed collar payload hash mismatch")
    with np.load(payload_path, allow_pickle=False) as saved:
        arrays = {key: saved[key].copy() for key in saved.files}
    for key, digest in record["array_sha256"].items():
        if hashlib.sha256(np.ascontiguousarray(arrays[key]).tobytes()).hexdigest() != digest:
            raise ValueError("sealed collar array binding failed: " + key)
    for path, digest in record["producers"].items():
        provenance.resolve_pinned_source_bytes(ROOT, path, digest, commit=SEALED_COMMIT)
    for path, digest in record["input_hashes"].items():
        if _sha_file(path) != digest:
            raise ValueError("sealed collar input closure mismatch")
    prediction_path = directory / "prediction.json"
    prepare_path = directory / "prepare.json"
    if _sha_file(prediction_path) != record["prediction_json_sha256"]:
        raise ValueError("sealed collar prediction binding failed")
    prediction = json.loads(prediction_path.read_text())
    if _sha_file(prepare_path) != prediction["prepare_json_sha256"]:
        raise ValueError("sealed collar preparation binding failed")
    if prediction.get("locked_before_held_nonlinear") is not True:
        raise ValueError("sealed collar prediction was not locked before measurement")
    held = record["tighter_measurement"]
    if float(held["w"]) != SEED_WEIGHT or float(record["measurement"]["w"]) != SEED_WEIGHT:
        raise ValueError("sealed collar seed weight is not 0.0013")
    sigma = arrays["fine_x"]
    state = arrays["fine_state"]
    child = int(np.argmin(np.abs(sigma - CHILD_RADIUS)))
    if abs(float(sigma[child]) - CHILD_RADIUS) > 1e-12:
        raise ValueError("sealed fine collar grid misses the child cut")
    eigenvalue = SEED_WEIGHT * 2.0 * float(state[child, 6])
    stored = float(held["readouts"]["child"]["finite_collar_covariance_eigenvalue"])
    if abs(eigenvalue - stored) > 1e-12:
        raise ValueError("sealed fine child covariance does not match its state integral")
    constants, input_hashes, binding = collar.load_inputs()
    if abs(constants["g"] - 8.0 * np.pi * constants["A"]) > 1e-12:
        raise ValueError("sealed collar coupling g is not 8 pi A")
    return {
        "record": record,
        "arrays": arrays,
        "constants": constants,
        "input_hashes": {str(Path(path).resolve().relative_to(ROOT)): digest for path, digest in input_hashes.items()},
        "input_binding": binding,
        "measurement_sha256": _sha_file(record_path),
        "payload_sha256": payload_hash,
        "prediction_sha256": _sha_file(prediction_path),
        "prepare_sha256": _sha_file(prepare_path),
        "sealed_child_covariance": stored,
        "sealed_collar_covariance": float(held["readouts"]["parent"]["finite_collar_covariance_eigenvalue"]),
        "center_radius": float(state[0, 0]),
        "center_proper_frequency": float(held["center_proper_frequency"]),
    }


def _reference_extension(auth, transition_end, deadline):
    constants = auth["constants"]
    sigma = np.asarray(auth["arrays"]["fine_x"], float)
    state = np.asarray(auth["arrays"]["fine_state"], float)
    if abs(float(sigma[-1]) - COLLAR_RADIUS) > 1e-12:
        raise ValueError("sealed fine collar must end at the collar radius")
    started = _cpu_time()

    def equation(x, y):
        if _cpu_time() > deadline:
            raise RuntimeError("parent CPU budget exhausted during collar reference extension")
        return collar.rhs(x, y, SEED_WEIGHT, constants)

    solved = solve_ivp(equation, (float(sigma[-1]), float(transition_end)), state[-1],
                       method="DOP853", dense_output=True, rtol=1e-10, atol=1e-12)
    if not solved.success:
        raise RuntimeError("collar static reference did not reach the blend point")
    probes = np.linspace(float(sigma[-1]), float(transition_end), 9)
    rows = solved.sol(probes).T
    gap = max(abs(collar.first_integral(row, SEED_WEIGHT, constants)) for row in rows)
    return solved, {"first_integral_max": float(gap), "CPU_seconds": _cpu_time() - started,
                    "extension": "sealed endpoint to the blend point; not a new held measurement"}


def _values_at(sigma, auth, extension, transition_end):
    """Collar tables on sigma>=0. Geometry uses the reference; spinor uses blended Q outside the core."""
    sigma = np.asarray(sigma, float)
    sealed_x = auth["arrays"]["fine_x"]
    sealed = auth["arrays"]["fine_state"]
    r = np.empty(sigma.shape, float)
    h = np.empty(sigma.shape, float)
    u = np.empty(sigma.shape, float)
    v = np.empty(sigma.shape, float)
    inner = sigma <= COLLAR_RADIUS + 1e-15
    r[inner] = np.interp(sigma[inner], sealed_x, sealed[:, 0])
    h[inner] = np.interp(sigma[inner], sealed_x, sealed[:, 2])
    u[inner] = np.interp(sigma[inner], sealed_x, sealed[:, 4])
    v[inner] = np.interp(sigma[inner], sealed_x, sealed[:, 5])
    outer = ~inner
    if np.any(outer):
        rows = extension.sol(np.clip(sigma[outer], COLLAR_RADIUS, transition_end))
        r[outer] = rows[0]
        h[outer] = rows[2]
        u[outer] = rows[4]
        v[outer] = rows[5]
    return r, np.exp(h), u, v


def _blend_geometry(sigma, reference_r, reference_Q, exterior_r, exterior_Q, transition_end):
    sigma = np.asarray(sigma, float)
    blend = np.zeros(sigma.shape, float)
    mid = (sigma > COLLAR_RADIUS) & (sigma < transition_end)
    blend[mid] = coupling.smooth_step((sigma[mid] - COLLAR_RADIUS) / (transition_end - COLLAR_RADIUS))
    blend[sigma >= transition_end] = 1.0
    return (1.0 - blend) * reference_r + blend * exterior_r, (1.0 - blend) * reference_Q + blend * exterior_Q, blend


def field_envelope(sigma, transition_end=TRANSITION_END, zero_at=ENVELOPE_ZERO):
    """One through the blend, then a flat decay that is zero before the parent edge."""
    sigma = np.abs(np.asarray(sigma, float))
    out = np.ones(sigma.shape, float)
    mid = (sigma > transition_end) & (sigma < zero_at)
    out[mid] = 1.0 - coupling.smooth_step((sigma[mid] - transition_end) / (zero_at - transition_end))
    out[sigma >= zero_at] = 0.0
    return out


def annulus_packet(sigma, inner=ANNULUS[0], outer=ANNULUS[1]):
    """Even C∞ bump on the parent annulus. Both spinor components use this amplitude."""
    sigma = np.abs(np.asarray(sigma, float))
    out = np.zeros(sigma.shape, float)
    mid = (sigma > inner) & (sigma < outer)
    t = (sigma[mid] - inner) / (outer - inner)
    out[mid] = coupling.smooth_step(t) * coupling.smooth_step(1.0 - t)
    return out


def _dirac_extension(auth, extension, profile, exterior_r, deadline):
    constants = auth["constants"]
    end = profile["transition_end"]
    zero = profile["envelope_zero"]
    sealed = auth["arrays"]["fine_state"]
    y0 = sealed[-1, 4:6].astype(float)

    def lapse(sigma):
        point = np.array([float(sigma)])
        reference_r, reference_Q, _u, _v = _values_at(point, auth, extension, end)
        _radius, blended, _weight = _blend_geometry(
            point, reference_r, reference_Q, exterior_r, profile["exterior_Q"], end)
        return float(blended[0])

    def equation(sigma, y):
        if _cpu_time() > deadline:
            raise RuntimeError("parent CPU budget exhausted during source continuation")
        dummy = np.array([1.0, 0.0, np.log(lapse(sigma)), 0.0, y[0], y[1], 0.0, 0.0, 0.0])
        return collar.rhs(sigma, dummy, SEED_WEIGHT, constants)[4:6]

    solved = solve_ivp(equation, (COLLAR_RADIUS, zero), y0, method="DOP853", dense_output=True,
                       rtol=1e-10, atol=1e-12, max_step=0.02)
    if not solved.success:
        raise RuntimeError("source Dirac continuation on the blended lapse failed")
    probes = np.linspace(COLLAR_RADIUS + 1e-4, min(end, zero) - 1e-4, 5)
    gaps = []
    for sigma in probes:
        step = 1e-6
        derivative = (solved.sol(sigma + step) - solved.sol(sigma - step)) / (2 * step)
        image = equation(float(sigma), solved.sol(sigma))
        gaps.append(_maximum(derivative - image))
    return solved, {"rhs_gap_max": float(max(gaps)), "follows_blended_Q_outside_core": True,
                    "core_uses_sealed_spinor": True}


def parent_continuum(profile=None, *, cpu_limit=CPU_LIMIT):
    """Authenticated collar blend. No momentum solve and no held-weight remeasurement."""
    selected = _profile(profile)
    key = json.dumps(selected, sort_keys=True)
    cached = _CACHE.get(("continuum", key))
    if cached is not None:
        return cached
    deadline = _cpu_time() + float(cpu_limit)
    auth = authenticate_collar()
    exterior_r = (float(np.sqrt(auth["constants"]["mag"] / auth["constants"]["g"]))
                  if selected["exterior_radius"] is None else float(selected["exterior_radius"]))
    extension, extension_report = _reference_extension(auth, selected["transition_end"], deadline)
    dirac, dirac_report = _dirac_extension(auth, extension, selected, exterior_r, deadline)
    sigma = np.unique(np.concatenate([
        auth["arrays"]["fine_x"], np.linspace(0.0, PARENT_RADIUS, 481),
        [selected["transition_end"], selected["envelope_zero"], ANNULUS[0], ANNULUS[1]]]))
    reference_r, reference_Q, sealed_u, sealed_v = _values_at(sigma, auth, extension, selected["transition_end"])
    radius, lapse, blend = _blend_geometry(
        sigma, reference_r, reference_Q, exterior_r, selected["exterior_Q"], selected["transition_end"])
    envelope = field_envelope(sigma, selected["transition_end"], selected["envelope_zero"])
    packet = annulus_packet(sigma)
    core = sigma <= COLLAR_RADIUS + 1e-15
    spinor_u = np.array(sealed_u, copy=True)
    spinor_v = np.array(sealed_v, copy=True)
    outside = ~core
    if np.any(outside):
        continued = dirac.sol(np.clip(sigma[outside], COLLAR_RADIUS, selected["envelope_zero"]))
        spinor_u[outside] = continued[0]
        spinor_v[outside] = continued[1]
    spinor_u *= envelope
    spinor_v *= envelope
    interior = auth["arrays"]["fine_x"] <= COLLAR_RADIUS + 1e-12
    sealed_r = np.interp(auth["arrays"]["fine_x"][interior], sigma, radius)
    sealed_Q = np.interp(auth["arrays"]["fine_x"][interior], sigma, lapse)
    tables = {
        "sigma": sigma, "radius": radius, "Q": lapse, "blend": blend, "envelope": envelope, "packet": packet,
        "spinor_u": spinor_u, "spinor_v": spinor_v, "exterior_radius": exterior_r,
        "profile": selected, "constants": auth["constants"], "authentication": {
            "producing_commit": SEALED_COMMIT, "payload_sha256": auth["payload_sha256"],
            "measurement_sha256": auth["measurement_sha256"], "prediction_sha256": auth["prediction_sha256"],
            "prepare_sha256": auth["prepare_sha256"], "input_hashes": auth["input_hashes"],
            "input_binding": auth["input_binding"], "sealed_child_covariance": auth["sealed_child_covariance"],
            "sealed_collar_covariance": auth["sealed_collar_covariance"],
            "center_radius": auth["center_radius"], "center_proper_frequency": auth["center_proper_frequency"],
        },
        "collar_radius_gap": _maximum(sealed_r - auth["arrays"]["fine_state"][interior, 0]),
        "collar_Q_gap": _maximum(sealed_Q - np.exp(auth["arrays"]["fine_state"][interior, 2])),
        "reference_extension": extension_report, "dirac_extension": dirac_report,
        "held_nonlinear_remeasured": False,
    }
    _CACHE[("continuum", key)] = tables
    return tables


def _sample_circle(grid, table, name):
    sigma = np.abs(signed_distance(grid.xi_q, grid.length))
    return np.interp(sigma, table["sigma"], table[name])


def _standing_spinor(grid, table):
    """Swap the sealed components across the centre so n and S stay even."""
    sigma = np.abs(signed_distance(grid.xi_q, grid.length))
    positive_u = np.interp(sigma, table["sigma"], table["spinor_u"])
    positive_v = np.interp(sigma, table["sigma"], table["spinor_v"])
    left = signed_distance(grid.xi_q, grid.length) < 0
    u = np.where(left, positive_v, positive_u)
    v = np.where(left, positive_u, positive_v)
    return u, v


def _project_spinor(grid, u, v):
    half0 = np.sqrt(grid.dx_q) * np.asarray(u, float)
    half1 = np.sqrt(grid.dx_q) * np.asarray(v, float)
    phi0 = grid.U_f.conj().T @ half0
    phi1 = grid.U_f.conj().T @ half1
    norm2 = float(np.vdot(phi0, phi0).real + np.vdot(phi1, phi1).real)
    if norm2 <= 0 or not np.isfinite(norm2):
        raise ValueError("AP projection of a parent source column vanished")
    return phi0, phi1, norm2


def _region_mass(grid, phi0, phi1, mask):
    density = np.abs(grid.U_f @ phi0) ** 2 + np.abs(grid.U_f @ phi1) ** 2
    return float(np.sum(density[mask]))


def _weighted_eigenvalues(gram, weights):
    scale = np.sqrt(np.asarray(weights, float))
    weighted = scale[:, None] * gram * scale[None, :]
    hermitian = 0.5 * (weighted + weighted.conj().T)
    return np.linalg.eigvalsh(hermitian).real


def population_weights(probabilities, raw_norms, gram, population, *, seed_weight=SEED_WEIGHT):
    """Actual regional imbalance at fixed trace T=2*w0*ZC; no weight reversal."""
    regional=np.asarray(probabilities,float);norms=np.asarray(raw_norms,float)
    if regional.shape!=(2,2) or norms.shape!=(2,) or not np.isfinite(regional).all() or np.any(regional<0) or np.any(norms<=0):
        raise ValueError("rank-two source needs actual child and annulus contents")
    if population not in POPULATIONS:raise ValueError("population must be 0, 1 or 2")
    if seed_weight!=SEED_WEIGHT:raise ValueError("the sealed source normalization weight stays fixed")
    child,annulus=regional;difference=child-annulus;total=child+annulus
    if np.any(total<=0):raise ValueError("regional source column has no observed content")
    endpoints=difference/total
    if not np.min(endpoints)<0<np.max(endpoints):
        raise ValueError("preparation obstruction: actual regional column endpoints do not straddle zero")
    delta_star=min(IMBALANCE_CAP,float(np.max(endpoints)),float(-np.min(endpoints)))
    target=(delta_star,0.,-delta_star)[population];trace=2*seed_weight*norms[0]
    fractions=np.linalg.solve(np.vstack((np.ones(2),difference-target*total)),np.array([1.,0.]))
    weights=trace*fractions
    if np.min(weights)<-1e-14:raise ValueError("regional imbalance requires negative source weight")
    weights=np.maximum(weights,0.)
    eigenvalues=_weighted_eigenvalues(gram,weights)
    if np.max(eigenvalues)>1+1e-10 or np.min(eigenvalues)<-1e-12:raise ValueError("weighted nonorthogonal source left CAR interval")
    achieved=float(np.dot(weights,difference)/np.dot(weights,total))
    return weights,{"trace":float(trace),"imbalance":achieved,"delta_star":delta_star,"alpha":delta_star,
       "imbalance_cap":IMBALANCE_CAP,"endpoint_imbalances":endpoints.tolist(),"target_imbalance":target,
       "population":population,"population_name":("child_heavy","balanced","parent_heavy")[population],"rank":2,
       "child_weight":float(weights[0]),"base_child_normalization_weight":float(seed_weight*norms[0]),
       "base_c_equals_wZ":True,"population_intervention_changes_base_weights":True,
       "trace_rule":"2*w0*projected_raw_child_norm_squared","lowdin_applied":False,
       "uniform_weights_substituted":False,"reversed_probabilities_substituted":False}


def _parity_basis(ng):
    if int(ng) % 2 == 0:
        raise ValueError("even/odd parent momenta need the odd geometry count")
    even = [np.eye(ng)[0]]
    odd = []
    for index in range(1, (ng - 1) // 2 + 1):
        right = np.eye(ng)[index]
        left = np.eye(ng)[ng - index]
        even.append(right + left)
        odd.append(right - left)
    even_basis = np.linalg.qr(np.column_stack(even))[0]
    odd_basis = np.linalg.qr(np.column_stack(odd))[0]
    return even_basis, odd_basis


def _band_geometry(grid, table):
    radius = galerkin.pull_geometry(grid, _sample_circle(grid, table, "radius"))
    lapse = galerkin.pull_geometry(grid, _sample_circle(grid, table, "Q"))
    even, _odd = _parity_basis(grid.ng)
    radius = even @ (even.T @ radius)
    lapse = even @ (even.T @ lapse)
    if min(float(np.min(grid.A_g @ radius)), float(np.min(grid.A_g @ lapse))) <= 0:
        raise ValueError("band-limited parent geometry left the positive chart")
    return radius, lapse


def _source_columns(grid, table):
    standing_u, standing_v = _standing_spinor(grid, table)
    packet = _sample_circle(grid, table, "packet")
    packet_gap = _maximum(packet - _sample_circle(grid, table, "packet"))
    del packet_gap
    child_raw = _project_spinor(grid, standing_u, standing_v)
    packet_raw = _project_spinor(grid, packet, packet)
    raw = (child_raw, packet_raw)
    units = []
    for phi0, phi1, norm2 in raw:
        scale = np.sqrt(norm2)
        units.append((phi0 / scale, phi1 / scale, norm2))
    child_window=np.abs(signed_distance(grid.xi_q,grid.length))<=CHILD_RADIUS
    probe0=grid.U_f.conj().T@(np.sqrt(grid.dx_q)*standing_u*child_window)
    probe1=grid.U_f.conj().T@(np.sqrt(grid.dx_q)*standing_v*child_window)
    probes=np.zeros((2*grid.nf,2),complex);probes[:grid.nf,0]=probe0;probes[grid.nf:,1]=probe1
    reference=np.linalg.qr(probes)[0]
    source = np.vstack((
        np.column_stack((units[0][0], units[1][0])),
        np.column_stack((units[0][1], units[1][1])),
    ))
    sigma = np.abs(signed_distance(grid.xi_q, grid.length))
    child = sigma <= CHILD_RADIUS
    annulus = (sigma >= ANNULUS[0]) & (sigma <= ANNULUS[1])
    annulus_masses=np.array([_region_mass(grid,item[0],item[1],annulus) for item in units])
    child_masses = np.array([
        _region_mass(grid, units[0][0], units[0][1], child),
        _region_mass(grid, units[1][0], units[1][1], child),
    ])
    probabilities=np.vstack((child_masses,annulus_masses))
    raw_child = _region_mass(grid, raw[0][0], raw[0][1], child)
    return {
        "reference_columns": reference, "source_columns": source,
        "phi0": np.column_stack((units[0][0], units[1][0])),
        "phi1": np.column_stack((units[0][1], units[1][1])),
        "raw_norms": np.array([raw[0][2], raw[1][2]]),
        "probabilities": probabilities, "unit_child_mass": child_masses,
        "raw_child_mass": raw_child,
        "unit_annulus_mass":annulus_masses,"observer_Gram_gap":_maximum(reference.conj().T@reference-np.eye(2)),
        "observer_source_distinct":True,
        "packet_components_equal_before_projection": True,
        "packet_preprojection_difference": 0.0,
    }


def _zero_constraint(pair, radius, lapse, phi0, phi1):
    zeros = np.zeros(pair.grid.ng)
    nodal = leading.State(lapse, radius, zeros, zeros, phi0, phi1)
    state = leading.encode(pair, nodal)
    fine = leading.fine_state(pair, state)
    system = leading.active_system(pair.grid, fine)
    source = coupling.source_from_columns(system, fine)
    hamilton, shift = leading.constraint_arrays(pair.grid, fine, system, source)
    return fine, system, source, hamilton, shift


def _primitive(values, length, origin):
    count = values.shape[0]
    coefficients = np.fft.fft(values) / count
    modes = np.fft.fftfreq(count) * count
    mean = float(coefficients[0].real)
    nyquist = float(coefficients[count // 2].real) if count % 2 == 0 else 0.0
    coefficients = coefficients.copy()
    coefficients[0] = 0
    if count % 2 == 0:
        coefficients[count // 2] = 0
    omega = 2 * np.pi * modes / length
    integrated = np.zeros_like(coefficients)
    nonzero = modes != 0
    integrated[nonzero] = coefficients[nonzero] / (1j * omega[nonzero])

    def evaluate(coordinate):
        coordinate = np.mod(np.asarray(coordinate, float), length)
        phase = np.exp(2j * np.pi * np.outer(coordinate, modes) / length)
        return np.real(phase @ integrated)

    baseline = float(np.asarray(evaluate(origin)).reshape(-1)[0])
    return lambda coordinate: evaluate(coordinate) - baseline, {
        "integrand_mean_removed": mean, "nyquist_cosine_removed": nyquist,
        "removed_modes_reported": True,
    }


def _offset_seed(pair, radius, lapse, phi0, phi1, margin):
    fine, system, source, hamilton, _shift = _zero_constraint(pair, radius, lapse, phi0, phi1)
    coupling_g = 8.0 * np.pi * float(system.A)
    slope = system.derivative @ fine.r
    integrand = 4.0 * coupling_g * slope * hamilton / (fine.r ** 2 * fine.Q)
    primitive, removal = _primitive(integrand, pair.grid.length, CENTER)
    jay = primitive(pair.grid.xi_q)
    reproduced = system.derivative @ jay
    removal["integrand_reproduction_gap"] = _maximum(reproduced - (integrand - removal["integrand_mean_removed"]))
    kmin = max(0.0, -float(np.min(jay)))
    scale = max(kmin, _maximum(jay), 1e-8)
    return {
        "J": jay, "kmin": kmin, "margin_scale": scale, "margin_fraction": float(margin),
        "suggested_k": kmin + max(float(margin) * scale, 1e-8), "C_s": hamilton, "fine": fine,
        "system": system, "source": source, "g": coupling_g, "removal": removal,
    }


def common_k(nf=64, profile=None, *, cpu_limit=CPU_LIMIT):
    """One offset above every population's kmin. Both momentum signs reuse it."""
    nf, _population, _sign, cpu_limit = _require_case(nf, 0, 1, cpu_limit)
    selected = _profile(profile)
    key = ("common", nf, json.dumps(selected, sort_keys=True))
    if key in _CACHE:
        return _CACHE[key]
    started = _cpu_time()
    deadline = started + cpu_limit
    table = parent_continuum(selected, cpu_limit=cpu_limit)
    with threadpool_limits(limits=1), backend.fft_thread_limit(1):
        grid = backend.make_fft_grid(galerkin.build_grid(nf, gauge="conformal"))
        radius, lapse = _band_geometry(grid, table)
        columns = _source_columns(grid, table)
        gram = columns["phi0"].conj().T @ columns["phi0"] + columns["phi1"].conj().T @ columns["phi1"]
        populations = []
        kmins = []
        for population in POPULATIONS:
            if _cpu_time() > deadline:
                raise RuntimeError("parent CPU budget exhausted while comparing population offsets")
            weights, meta = population_weights(columns["probabilities"], columns["raw_norms"], gram, population)
            grid.fine = replace(grid.fine, occupations=np.array(weights, dtype=float, copy=True))
            pair = _holder(grid, weights, columns, meta, table, population)
            seed = _offset_seed(pair, radius, lapse, columns["phi0"], columns["phi1"], selected["k_margin"])
            kmins.append(seed["kmin"])
            populations.append({"population": population, "weights": weights.tolist(), **meta,
                                "kmin": seed["kmin"], "J_abs_max": _maximum(seed["J"]),
                                "integrand_mean_removed": seed["removal"]["integrand_mean_removed"]})
        kmin = float(max(kmins))
        scale = max([kmin, 1e-8] + [item["J_abs_max"] for item in populations])
        suggested = max(item["kmin"] + max(selected["k_margin"] * max(item["kmin"], item["J_abs_max"], 1e-8), 1e-8)
                        for item in populations)
    report = {"k": float(suggested), "kmin": kmin, "kmin_by_population": [float(value) for value in kmins],
              "margin": float(suggested - kmin), "margin_fraction": selected["k_margin"],
              "signs": list(SIGNS), "populations": _plain(populations), "trace": populations[0]["trace"],
              "alpha": populations[0]["alpha"], "sign_independent": True, "static_not_forced": True,
              "CPU_seconds": _cpu_time() - started}
    _CACHE[key] = report
    return report


def _holder(grid, weights, columns, meta, table, population):
    identity = np.eye(grid.ng)
    geometry = {
        "map": "identity", "all_geometry_degrees_retained": True, "ambient_field_retained": True,
        "centre": CENTER, "child_radius": CHILD_RADIUS, "collar_radius": COLLAR_RADIUS,
        "parent_radius": PARENT_RADIUS, "transition_end": table["profile"]["transition_end"],
        "transition_before_pi_over_2": True, "exterior_Q": table["profile"]["exterior_Q"],
        "exterior_radius": table["exterior_radius"], "collar_radius_gap": table["collar_radius_gap"],
        "collar_Q_gap": table["collar_Q_gap"],
    }
    source = {
        "layout": "collar-standing-plus-even-annulus-packet", "rank": 2, "orthogonalized": False,
        "lowdin_applied": False, "packet_components_equal_before_projection": True,
        "packet_support": list(ANNULUS), "envelope_one_through": table["profile"]["transition_end"],
        "envelope_zero_at": table["profile"]["envelope_zero"],
        "source_plateau_is_maintained_wall": False, "source_plateau_enters_timestep": True,
        "probabilities": columns["probabilities"].tolist(), "raw_norms": columns["raw_norms"].tolist(),
        "raw_child_mass": columns["raw_child_mass"], "unit_child_mass": columns["unit_child_mass"].tolist(),
        "unit_annulus_mass":columns['unit_annulus_mass'].tolist(),
        "observer_Gram_gap":columns['observer_Gram_gap'],"observer_source_distinct":True,
        **meta,
    }
    weights = np.array(weights, float, copy=True)
    return ParentPair(
        grid, _readonly(identity), _readonly(weights), _readonly(columns["reference_columns"]),
        _readonly(columns["source_columns"]), source, geometry)


def _readonly(value):
    array = np.array(value, copy=True)
    array.setflags(write=False)
    return array


class _MomentumProblem:
    """Retained even h=QC, odd shift, and the seed anchor. Q, r, Phi and c stay fixed."""

    def __init__(self, pair, radius, lapse, phi0, phi1, sign, anchor, coupling_g):
        self.pair = pair
        self.grid = pair.grid
        self.radius = np.array(radius, float, copy=True)
        self.lapse = np.array(lapse, float, copy=True)
        self.phi0 = phi0
        self.phi1 = phi1
        self.sign = int(sign)
        self.anchor = float(anchor)
        self.g = float(coupling_g)
        self.even, self.odd = _parity_basis(self.grid.ng)
        self.count = self.even.shape[1]
        zeros = np.zeros(self.grid.ng)
        nodal = leading.State(self.lapse, self.radius, zeros, zeros, phi0, phi1)
        fine = leading.fine_state(pair, leading.encode(pair, nodal))
        self.system = leading.active_system(self.grid, fine)
        self.source = coupling.source_from_columns(self.system, fine)
        self.fine_r = fine.r
        self.fine_Q = fine.Q
        a, zeta, _mag = leading.coefficients(self.system)
        self.a = a
        self.zeta = zeta
        self.b = a * fine.r
        self.rx = self.system.derivative @ fine.r
        self.synthesis = self.grid.A_g @ np.eye(self.grid.ng)
        self.geometry_to_fine = self.system.derivative @ self.synthesis
        scales,self.frequency,_kg,_kf=step_control.physical_scales(self.grid,fine,self.system)
        self.momentum_scales=scales[:,2:4]
        lift=self.synthesis@self.even
        gp=self.grid.dx_q*(lift.T/scales[:,2]**2)@lift
        gr=self.grid.dx_q*(lift.T/scales[:,3]**2)@lift
        self.metric=np.block([[gp,np.zeros_like(gp)],[np.zeros_like(gr),gr]])
        self.whitener=np.linalg.solve(np.linalg.cholesky(self.metric).T,np.eye(2*self.count))
        zero_C,_D=leading.constraint_arrays(self.grid,fine,self.system,self.source)
        source_h=fine.Q*self.source['force_L']/self.grid.dx_q
        self.zero_gravity_h=fine.Q*zero_C-source_h
        self.source_h=source_h

    def nodal(self, theta):
        alpha = theta[:self.count]
        beta = theta[self.count:]
        return self.even @ alpha, self.even @ beta

    def fine_momenta(self, theta):
        momentum, radial = self.nodal(theta)
        return self.grid.A_g @ momentum, self.grid.A_g @ radial

    def residual(self, theta):
        momentum, radial = self.fine_momenta(theta)
        fine = leading.State(self.fine_Q, self.fine_r, momentum, radial, self.grid.U_f @ self.phi0, self.grid.U_f @ self.phi1)
        hamilton, shift = leading.constraint_arrays(self.grid, fine, self.system, self.source)
        density = self.fine_Q * hamilton
        retained_h = galerkin.pull_geometry(self.grid, density)
        retained_d = galerkin.pull_geometry(self.grid, shift)
        mean = float(np.mean(momentum ** 2 / self.fine_r ** 3))
        return np.concatenate((self.even.T @ retained_h, self.odd.T @ retained_d, np.array([mean - self.anchor])))

    def jacobian(self, theta):
        """Analytic momentum derivative of the finite constraints. D(F) stays inside the residual."""
        momentum, radial = self.fine_momenta(theta)
        dh_dp = self.fine_Q * (radial / (2 * self.b) - self.zeta * self.fine_Q * momentum / (2 * self.b ** 2))
        dh_dv = self.fine_Q * (momentum / (2 * self.b))
        weight = self.grid.weight
        synthesis = self.synthesis
        analysis = weight * synthesis.T
        mass_p = (analysis * dh_dp) @ synthesis
        mass_v = (analysis * dh_dv) @ synthesis
        shift_v = (analysis * self.rx) @ synthesis
        shift_p = analysis @ (-self.fine_Q[:, None] * self.geometry_to_fine)
        anchor_p = (2 * momentum / (self.fine_r ** 3 * self.grid.nq)) @ synthesis
        blocks = (
            (self.even.T @ mass_p @ self.even, self.even.T @ mass_v @ self.even),
            (self.odd.T @ shift_p @ self.even, self.odd.T @ shift_v @ self.even),
            (anchor_p @ self.even, np.zeros(self.count)),
        )
        return np.block([[blocks[0][0], blocks[0][1]], [blocks[1][0], blocks[1][1]],
                         [blocks[2][0], blocks[2][1]]])

    def sign_floor(self, theta):
        momentum, _radial = self.fine_momenta(theta)
        scale = max(1e-12, float(np.max(np.abs(momentum))))
        return float(np.min(self.sign * momentum)), scale

    def row_scales(self,theta):
        p,v=self.fine_momenta(theta);Q=self.fine_Q;r=self.fine_r
        hscale=max(_maximum(self.zero_gravity_h),_maximum(self.source_h),
                   _maximum(Q*p*v/(2*self.b)),_maximum(self.zeta*Q*Q*p*p/(4*self.b*self.b)),
                   _maximum(abs(self.a)*Q*Q*r*r),1e-15)
        dscale=max(_maximum(v*self.rx),_maximum(Q*(self.system.derivative@p)),
                   _maximum(self.source['force_beta']/self.grid.dx_q),
                   _maximum(abs(self.a)*r*r*Q*self.frequency),1e-15)
        jscale=max(self.anchor,float(np.mean(self.momentum_scales[:,0]**2/r**3)),1e-15)
        return np.r_[np.full(self.count,hscale),np.full(self.odd.shape[1],dscale),jscale]


def _svd_step(jacobian, residual):
    """Minimum-norm step in the orthonormal parity metric, which is the uniform nodal metric."""
    factor = np.linalg.svd(jacobian, full_matrices=False)
    left, singular, right = factor
    cutoff = 1e-10 * float(singular[0]) if singular.size else 1.0
    inverse = np.zeros_like(singular);active=singular>max(cutoff,1e-12);inverse[active]=1/singular[active]
    step = right.T @ (inverse * (left.T @ (-residual)))
    return step, int(np.count_nonzero(singular > max(cutoff, 1e-12))), singular


def _seed_theta(problem, seed, sign, *, anchor_override=None):
    safe = np.maximum(seed["J"] + seed["suggested_k"], 0.0)
    momentum = sign * np.sqrt(np.maximum(seed["fine"].r ** 3 * (seed["suggested_k"] + seed["J"]), 0.0))
    if float(np.min(safe)) <= 0:
        raise ValueError("offset does not keep the continuum square root positive")
    radial = (3 * seed["fine"].Q * momentum / (2 * seed["fine"].r)
              + 2 * seed["g"] * seed["fine"].r * seed["C_s"] / momentum)
    nodal_p = problem.even @ (problem.even.T @ galerkin.pull_geometry(problem.grid, momentum))
    nodal_v = problem.even @ (problem.even.T @ galerkin.pull_geometry(problem.grid, radial))
    theta = np.concatenate((problem.even.T @ nodal_p, problem.even.T @ nodal_v))
    prolonged, _radial_fine = problem.fine_momenta(theta)
    problem.anchor = float(np.mean(prolonged ** 2 / problem.fine_r ** 3)) if anchor_override is None else float(anchor_override)
    return theta


def correct_momenta(problem, theta, deadline, *, max_accepted=MAX_ACCEPTED):
    """Sign-preserving damped correction. A 1e-8 residual is not a veto."""
    accepted = 0
    backtracks = 0
    history = []
    converged = False
    blocker = None
    best = np.array(theta, copy=True)
    rows=problem.row_scales(theta)
    best_norm = float(np.linalg.norm(problem.residual(best)/rows))
    floor_seed, scale_seed = problem.sign_floor(best)
    while accepted < max_accepted:
        if _cpu_time() > deadline:
            blocker = "preparation CPU cap"
            break
        residual = problem.residual(theta)
        jacobian = problem.jacobian(theta)
        dimensionless, rank, _singular = _svd_step((jacobian/rows[:,None])@problem.whitener,residual/rows)
        step=problem.whitener@dimensionless
        correction = float(np.linalg.norm(dimensionless))
        floor = max(1e-11, 64*np.finfo(float).eps*max(1.,np.linalg.norm(np.linalg.solve(problem.whitener,theta))))
        history.append({"accepted": accepted, "residual": float(np.linalg.norm(residual)),
                        "scaled_residual":float(np.linalg.norm(residual/rows)),
                        "correction": correction, "numerical_rank": rank, "floor": floor,
                        "sign_margin":problem.sign_floor(theta)[0]})
        if correction <= floor:
            converged = True
            blocker = None
            break
        improved = False
        for halving in range(MAX_BACKTRACKS):
            if _cpu_time() > deadline:
                blocker = "preparation CPU cap"
                break
            backtracks += 1
            trial = theta + step * 2.0 ** (-halving)
            trial_floor, trial_scale = problem.sign_floor(trial)
            if trial_floor < min(floor_seed, 0.0) - 1e-8 * max(scale_seed, trial_scale):
                continue
            if not np.isfinite(trial).all():
                continue
            trial_residual = problem.residual(trial)
            trial_norm = float(np.linalg.norm(trial_residual/rows))
            if trial_norm < float(np.linalg.norm(residual/rows)):
                theta = trial
                accepted += 1
                improved = True
                if trial_norm < best_norm:
                    best = np.array(trial, copy=True)
                    best_norm = trial_norm
                break
        if blocker == "preparation CPU cap":
            break
        if not improved:
            blocker = "momentum Newton stalled"
            break
    else:
        residual = problem.residual(theta)
        if max_accepted and correction <= floor:
            converged = True
        else:
            blocker = blocker or "finite momentum iteration cap"
    if float(np.linalg.norm(problem.residual(theta)/rows)) > best_norm:
        theta = best
    return theta, {"converged": converged, "blocker": blocker, "accepted": accepted,
                   "backtracks": backtracks, "history": history, "best_residual": best_norm,
                   "universal_1e-8_veto": False, "max_accepted": MAX_ACCEPTED,
                   "max_backtracks": MAX_BACKTRACKS,"physical_metric":"dxq Bfine.T diag(SP^-2,SR^-2) Bfine; Cholesky-whitened SVD",
                   "row_scales":rows.tolist(),"actual_residual_max":_maximum(problem.residual(theta)),
                   "stationary_correction_is_not_continuous_constraint_certificate":True}


def even_geometry_update(grid, base_radius, base_lapse, theta):
    """Fallback chart: even mean and first cosine of Q and r. Four parameters."""
    if np.asarray(theta).shape != (4,):
        raise ValueError("geometry fallback has four even parameters")
    cosine = np.cos(2 * np.pi * (grid.xi_g - CENTER) / grid.length)
    radius = base_radius + float(theta[2]) + float(theta[3]) * cosine
    lapse = base_lapse + float(theta[0]) + float(theta[1]) * cosine
    even, _odd = _parity_basis(grid.ng)
    radius = even @ (even.T @ radius)
    lapse = even @ (even.T @ lapse)
    if min(float(np.min(grid.A_g @ radius)), float(np.min(grid.A_g @ lapse))) <= 0:
        raise ValueError("fallback geometry left the positive chart")
    return radius, lapse


def _fallback(pair, radius, lapse, phi0, phi1, sign, coupling_g, k, margin, deadline, anchor):
    """One four-parameter even correction. Phi and the weights stay fixed; the source is rebuilt."""
    theta = np.zeros(4)
    best = (np.array(radius, copy=True), np.array(lapse, copy=True), None, None)
    best_norm = None
    steps = 0
    for _step in range(FALLBACK_STEPS):
        if _cpu_time() > deadline:
            break
        current_r, current_Q = even_geometry_update(pair.grid, radius, lapse, theta)
        seed = _offset_seed(pair, current_r, current_Q, phi0, phi1, margin)
        seed["suggested_k"] = float(k)
        if float(np.min(seed["J"] + k)) <= 0:
            break
        problem = _MomentumProblem(pair, current_r, current_Q, phi0, phi1, sign, 0.0, coupling_g)
        seeded = _seed_theta(problem, seed, sign,anchor_override=anchor)
        residual = problem.residual(seeded)
        rows=problem.row_scales(seeded);norm = float(np.linalg.norm(residual/rows))
        if best_norm is None or norm < best_norm:
            best = (current_r, current_Q, seeded, float(problem.anchor))
            best_norm = norm
        columns = []
        base = residual/rows
        for index in range(4):
            if _cpu_time() > deadline:
                break
            step = np.zeros(4)
            step[index] = 1e-5
            try:
                bumped_r, bumped_Q = even_geometry_update(pair.grid, radius, lapse, theta + step)
            except ValueError:
                columns.append(np.zeros_like(base))
                continue
            bumped = _offset_seed(pair, bumped_r, bumped_Q, phi0, phi1, margin)
            bumped["suggested_k"] = float(k)
            if float(np.min(bumped["J"] + k)) <= 0:
                columns.append(np.zeros_like(base))
                continue
            probe = _MomentumProblem(pair, bumped_r, bumped_Q, phi0, phi1, sign, problem.anchor, coupling_g)
            columns.append((probe.residual(_seed_theta(probe, bumped, sign,anchor_override=anchor))/rows - base) / 1e-5)
        if len(columns) < 4:
            break
        jacobian = np.column_stack(columns)
        modes=np.column_stack((np.ones(pair.grid.ng),np.cos(2*np.pi*(pair.grid.xi_g-CENTER)/pair.grid.length)))
        lift=pair.grid.A_g@modes
        Gq=pair.grid.dx_q*(lift.T/problem.fine_Q**2)@lift
        Gr=pair.grid.dx_q*(lift.T/problem.fine_r**2)@lift
        G=np.block([[Gq,np.zeros_like(Gq)],[np.zeros_like(Gr),Gr]])
        whitening=np.linalg.solve(np.linalg.cholesky(G).T,np.eye(4))
        step,_rank,_singular=_svd_step(jacobian@whitening,base);correction=whitening@step
        improved = False
        for halving in range(MAX_BACKTRACKS):
            trial = theta + correction * 2.0 ** (-halving)
            try:
                trial_r, trial_Q = even_geometry_update(pair.grid, radius, lapse, trial)
                trial_seed = _offset_seed(pair, trial_r, trial_Q, phi0, phi1, margin)
            except ValueError:
                continue
            trial_seed["suggested_k"] = float(k)
            if float(np.min(trial_seed["J"] + k)) <= 0:
                continue
            trial_problem = _MomentumProblem(pair, trial_r, trial_Q, phi0, phi1, sign, 0.0, coupling_g)
            trial_norm = float(np.linalg.norm(trial_problem.residual(_seed_theta(trial_problem, trial_seed, sign,anchor_override=anchor))/rows))
            if trial_norm < norm:
                theta = trial
                steps += 1
                improved = True
                break
        if not improved:
            break
    return best, {"steps": steps, "max_steps": FALLBACK_STEPS, "parameters": 4,
                  "phi_and_weights_fixed": True, "source_recomputed": True, "best_seed_residual": best_norm,
                  "original_kinetic_anchor_fixed":anchor,"physical_geometry_metric":True}


def _diagnose(pair, state, problem):
    fine = leading.fine_state(pair, state)
    system = leading.active_system(pair.grid, fine)
    source = coupling.source_from_columns(system, fine)
    hamilton, shift = leading.constraint_arrays(pair.grid, fine, system, source)
    density = fine.Q * hamilton
    retained_h = galerkin.pull_geometry(pair.grid, density)
    retained_d = galerkin.pull_geometry(pair.grid, shift)
    tail_h = density - pair.grid.A_g @ retained_h
    tail_d = shift - pair.grid.A_g @ retained_d
    even_h = problem.even.T @ retained_h
    odd_h = problem.odd.T @ retained_h
    even_d = problem.even.T @ retained_d
    odd_d = problem.odd.T @ retained_d
    if _maximum(problem.even @ even_h + problem.odd @ odd_h - retained_h) > 1e-8:
        raise ValueError("parity split discarded part of the retained lapse constraint")
    if _maximum(problem.even @ even_d + problem.odd @ odd_d - retained_d) > 1e-8:
        raise ValueError("parity split discarded part of the retained shift constraint")
    current = source["force_beta"] / pair.grid.dx_q
    rate, bundle = leading.rates(pair, state, return_bundle=True)
    curvature, velocity, _acceleration = leading.metric_jets(pair, state, rate, bundle)
    clocks = leading.clock_rates(pair, state)
    coordinate_angular=velocity.r/fine.r
    expansion = velocity.r/(fine.r**2*fine.Q)
    radial_normal=(velocity.r/fine.r+velocity.Q/fine.Q)/(fine.r*fine.Q)
    centre = extent.real_periodic_values(pair.grid, expansion, (CENTER,))[0]
    restriction = step_control.step_restriction(pair, state, 1.0)
    energies = leading.energy(pair, state)
    leading_constraints = leading.constraints(pair, state)
    return {
        "full_h_max": _maximum(density), "full_D_max": _maximum(shift),
        "retained_h_max": _maximum(retained_h), "retained_D_max": _maximum(retained_d),
        "tail_h_max": _maximum(tail_h), "tail_D_max": _maximum(tail_d),
        "retained_even_h_max": _maximum(even_h), "retained_odd_D_max": _maximum(odd_d),
        "excluded_odd_h_max": _maximum(odd_h), "excluded_even_D_max": _maximum(even_d),
        "excluded_parity_retained_in_report": True, "current_max": _maximum(current),
        "current_mean": float(np.mean(current)), "current_deleted": False,
        "momentum_includes_current": True, "anchor": problem.anchor,
        "anchor_residual": float(np.mean((fine.p_Q ** 2) / fine.r ** 3) - problem.anchor),
        "leading_constraints": leading_constraints, "energy": energies,
        "expansion": {**_summary(expansion), "centre": float(centre),"scope":"normal angular H=projected rdot/(r^2 Q)"},
        "coordinate_angular_rate":{**_summary(coordinate_angular),"scope":"projected rdot/r; coordinate T"},
        "normal_radial_rate":{**_summary(radial_normal),"scope":"(projected rdot/r+projected Qdot/Q)/(r Q)"},
        "rates": {name: _summary(getattr(velocity,name)) for name in ("Q", "r", "p_Q", "p_r")},
        "rates_scope":"actual decoded and prolonged projected rates; coordinate T",
        "raw_unprojected_rates":{name: _summary(bundle["unprojected_rates"][index]) for index,name in enumerate(("Q","r","p_Q","p_r"))},
        "raw_unprojected_rates_used_for_acceptance":False,
        "tides": {name: _summary(curvature["tides"][name]) for name in ("R_0101", "R_0202", "R4", "owned_W")},
        "curvature": {name: _summary(curvature["base"][name]) for name in ("Ricci2", "K")},
        "proper_clocks": [{"x": float(place), "N": float(value)} for place, value in zip(pair.clock_locations, clocks)],
        "velocity_r": _summary(velocity.r), "timestep": restriction[1],
        "source_plateau_is_maintained_wall": False,
        "CAR_eigenvalues": _weighted_eigenvalues(
            state.phi0.conj().T @ state.phi0 + state.phi1.conj().T @ state.phi1, pair.weights).tolist(),
    }


def _state_from_theta(pair, problem, theta):
    momentum, radial = problem.nodal(theta)
    nodal = leading.State(problem.lapse, problem.radius, momentum, radial, problem.phi0, problem.phi1)
    return leading.encode(pair, nodal)


def prepare_parent(nf=64, population=0, sign=+1, k_override=None, profile=None, cpu_limit=30):
    """Return one parent pair, canonical leading state and diagnostic report. No trajectory."""
    nf, population, sign, cpu_limit = _require_case(nf, population, sign, cpu_limit)
    selected = _profile(profile)
    override = None if k_override is None else float(k_override)
    key = (nf, population, sign, override, json.dumps(selected, sort_keys=True))
    if key in _CACHE and _CACHE[key][0] <= cpu_limit:
        pair, state, report = _CACHE[key][1]
        return pair, state.copy(), _plain(report)
    started = _cpu_time()
    hard_deadline=started+cpu_limit
    deadline=hard_deadline-min(.5,cpu_limit*.1)
    with threadpool_limits(limits=1), backend.fft_thread_limit(1):
        table = parent_continuum(selected, cpu_limit=max(1e-8,deadline-_cpu_time()))
        family = common_k(nf, selected, cpu_limit=max(1e-8,deadline-_cpu_time()))
        kmin = float(family["kmin"])
        k = family["k"] if override is None else override
        if not k > kmin:
            raise ValueError("k must be strictly above the common kmin")
        grid = backend.make_fft_grid(galerkin.build_grid(nf, gauge="conformal"))
        radius, lapse = _band_geometry(grid, table)
        columns = _source_columns(grid, table)
        gram = columns["phi0"].conj().T @ columns["phi0"] + columns["phi1"].conj().T @ columns["phi1"]
        weights, meta = population_weights(columns["probabilities"], columns["raw_norms"], gram, population)
        from dataclasses import replace as replace_dataclass
        grid.fine = replace_dataclass(grid.fine, occupations=np.array(weights, dtype=float, copy=True))
        pair = _holder(grid, weights, columns, meta, table, population)
        seed = _offset_seed(pair, radius, lapse, columns["phi0"], columns["phi1"], selected["k_margin"])
        seed["suggested_k"] = float(k)
        if float(np.min(seed["J"] + k)) <= 0:
            raise ValueError("common offset does not keep this population's square root positive")
        problem = _MomentumProblem(pair, radius, lapse, columns["phi0"], columns["phi1"], sign, 0.0, seed["g"])
        theta = _seed_theta(problem, seed, sign)
        theta, correction = correct_momenta(problem, theta, deadline,max_accepted=MAX_ACCEPTED-4)
        fallback = {"steps": 0, "used": False, "max_steps": FALLBACK_STEPS}
        if not correction["converged"] and _cpu_time() < deadline:
            (radius, lapse, seeded, fallback_anchor), fallback = _fallback(
                pair, radius, lapse, columns["phi0"], columns["phi1"], sign, seed["g"], k,
                selected["k_margin"], deadline,problem.anchor)
            fallback["used"] = True
            problem = _MomentumProblem(pair, radius, lapse, columns["phi0"], columns["phi1"], sign, fallback_anchor or problem.anchor, seed["g"])
            theta = seeded if seeded is not None else theta
            theta, second = correct_momenta(problem, theta, deadline,max_accepted=MAX_ACCEPTED-correction['accepted'])
            second['accepted_before_fallback']=correction['accepted']
            second['accepted_total']=second['accepted']+correction['accepted']
            second['pre_fallback_correction']=correction
            correction = second
            correction["fallback_preceded"] = True
        state = _state_from_theta(pair, problem, theta)
        diagnosis = _diagnose(pair, state, problem) if _cpu_time() < hard_deadline else None
    child_target = SEED_WEIGHT * float(columns["raw_child_mass"])
    child_column = float(weights[0] * columns["unit_child_mass"][0])
    child_region = float(np.dot(weights, columns["unit_child_mass"]))
    report = {
        "schema": SCHEMA, "evolved": False, "trajectory": False, "new_force": False,
        "nf": nf, "population": population, "sign": sign, "k": float(k), "kmin": kmin,
        "k_offset": float(k - kmin), "k_common": float(family["k"]),
        "kmin_by_population": family["kmin_by_population"], "declared_margin_fraction": selected["k_margin"],
        "static_not_forced": True, "universal_1e-8_veto": False,
        "weights": weights.tolist(), "trace": meta["trace"], "imbalance": meta["imbalance"],
        "alpha": meta["alpha"], "imbalance_cap": IMBALANCE_CAP, "weight_metadata": meta,
        "probabilities": columns["probabilities"].tolist(), "raw_norms": columns["raw_norms"].tolist(),
        "gram_offdiag_max": _maximum(gram - np.diag(np.diag(gram))),
        "column_norms": np.linalg.norm(np.vstack((columns["phi0"], columns["phi1"])), axis=0).tolist(),
        "child_column_covariance": child_column, "child_region_covariance": child_region,
        "child_covariance_target": child_target,
        "sealed_child_covariance": table["authentication"]["sealed_child_covariance"],
        "base_child_normalization_identity_gap": float(SEED_WEIGHT*columns["raw_norms"][0]*columns["unit_child_mass"][0]-child_target),
        "population_changes_child_covariance": True,
        "child_column_preserved": bool(abs(child_column-child_target)<=1e-12),
        "occupations_match_weights": True, "geometry_map": "identity",
        "intervals": {"child": list(CHILD_INTERVAL), "collar_radius": COLLAR_RADIUS,
                      "parent": list(PARENT_INTERVAL), "annulus": list(ANNULUS)},
        "clock_locations": list(CLOCK_LOCATIONS), "correction": correction, "fallback": fallback,
        "observer_Gram_gap":columns['observer_Gram_gap'],"observer_is_source_frame":False,
        "source_field_sha256":{"phi0":_sha_array(state.phi0),"phi1":_sha_array(state.phi1)},
        "seed_removal": seed["removal"], "continuum_kmin": seed["kmin"],
        "collar_radius_gap": table["collar_radius_gap"], "collar_Q_gap": table["collar_Q_gap"],
        "reference_extension": table["reference_extension"], "dirac_extension": table["dirac_extension"],
        "authentication": table["authentication"], "centre_proper_frequency_sealed": table["authentication"]["center_proper_frequency"],
        "diagnosis": diagnosis, "CPU_seconds": _cpu_time() - started, "CPU_limit_seconds": cpu_limit,
        "held_nonlinear_remeasured": False, "old_build_pair_used": False,
        "converged": bool(correction["converged"]), "blocker": correction["blocker"],
        "status":"NUMERICAL_STATIONARY_CANDIDATE" if correction['converged'] else "OPEN_PREPARATION",
        "continuous_constraint_certificate":False,"global_static_parent_claimed":False,
        "unresolved": [
            "finite SBP residual and band-limit gaps are indicators, not a continuum certificate",
            "physical vacuum identification and strong-curvature EFT validity remain open",
            "no trajectory or autonomous renewal is claimed",
        ],
    }
    if diagnosis is None:
        report["blocker"] = report["blocker"] or "CPU cap before diagnostics"
        report["converged"] = False
    plain = _plain(report)
    _CACHE[key] = (cpu_limit, (pair, state, plain))
    return pair, state.copy(), plain


def linearization_error(pair, state, seed=1):
    """One directional finite-difference check of the analytic momentum Jacobian."""
    zeros = np.zeros(pair.grid.ng)
    nodal = leading.decode(pair, state)
    problem = _MomentumProblem(pair, nodal.Q, nodal.r, nodal.phi0, nodal.phi1, 1, 0.0, 1.0)
    theta = np.concatenate((problem.even.T @ nodal.p_Q, problem.even.T @ nodal.p_r))
    problem.anchor = float(np.mean(problem.fine_momenta(theta)[0] ** 2 / problem.fine_r ** 3))
    jacobian = problem.jacobian(theta)
    generator = np.random.default_rng(seed)
    direction = generator.normal(size=theta.size)
    direction /= np.linalg.norm(direction)
    step = 1e-6
    difference = (problem.residual(theta + step * direction) - problem.residual(theta - step * direction)) / (2 * step)
    return _maximum(jacobian @ direction - difference)


def preview():
    auth=authenticate_collar();selected=_profile(None)
    return {"schema": SCHEMA, "mode": "preview", "evolved": False, "trajectory": False, "new_force": False,
            "api": "prepare_parent(nf=64, population=0, sign=+1, k_override=None, profile=None, cpu_limit=30)",
            "populations": list(POPULATIONS), "signs": list(SIGNS), "imbalance_cap": IMBALANCE_CAP,
            "centre": CENTER, "intervals": {"child": list(CHILD_INTERVAL), "parent": list(PARENT_INTERVAL),
                                            "annulus": list(ANNULUS), "collar_radius": COLLAR_RADIUS},
            "clock_locations": list(CLOCK_LOCATIONS), "sealed_commit": SEALED_COMMIT,
            "sealed_payload_sha256": auth["payload_sha256"], "held_w": SEED_WEIGHT,
            "transition_end": selected["transition_end"], "envelope_zero": selected["envelope_zero"],
            "exterior_radius":float(np.sqrt(auth['constants']['mag']/auth['constants']['g'])),"exterior_Q":selected['exterior_Q'],
            "constructor_called":False,"spatial_ODE_solved":False,
            "source_plateau_is_maintained_wall": False, "held_nonlinear_remeasured": False,
            "production": "root freezes producers before --prepare; this preview writes nothing"}


def producing_commit(pins):
    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
        for path, digest in pins.items():
            provenance.resolve_pinned_source_bytes(ROOT, path, digest, commit=commit)
        return commit
    except (OSError, subprocess.SubprocessError, RuntimeError):
        return None


def _write(directory, report, arrays):
    directory = Path(directory)
    record_path = directory / "parent.json"
    payload_path = directory / "parent.npz"
    if record_path.exists() or payload_path.exists():
        raise FileExistsError("parent record already exists")
    commit = report.get("producing_commit")
    if not commit:
        raise ValueError("record creation requires a frozen non-None producing commit")
    if report["producers"] != source_hashes():
        raise ValueError("record creation requires the current producer bytes")
    for path, digest in report["producers"].items():
        provenance.resolve_pinned_source_bytes(ROOT, path, digest, commit=commit)
    for path,digest in report['input_hashes'].items():
        if _sha_file(ROOT/path)!=digest:raise ValueError('parent input changed before immutable creation: '+path)
    stream = io.BytesIO()
    np.savez_compressed(stream, **arrays)
    payload = stream.getvalue()
    bound = dict(report, payload_sha256=hashlib.sha256(payload).hexdigest(),
                 array_sha256={key: hashlib.sha256(np.ascontiguousarray(value).tobytes()).hexdigest()
                               for key, value in arrays.items()})
    text = (json.dumps(_plain(bound), sort_keys=True, indent=2, allow_nan=False) + "\n").encode()
    if len(payload) + len(text) > CHUNK_LIMIT:
        raise RuntimeError("parent record exceeds 64 MiB")
    directory.mkdir(parents=True, exist_ok=True)
    for path, blob in ((payload_path, payload), (record_path, text)):
        descriptor = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o444)
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(blob)
    return bound


def prepare(directory, *, nf=64, population=0, sign=+1, k_override=None, profile=None,
            cpu_limit=30, producer_commit=None):
    """Explicit frozen record. Root calls this after the producer commit is immutable."""
    pins = source_hashes()
    commit = producer_commit or producing_commit(pins)
    if not commit:
        raise ValueError("explicit parent preparation requires a frozen producer commit")
    pair, state, report = prepare_parent(nf, population, sign, k_override, profile, cpu_limit)
    arrays = {
        "Q": np.array(state.Q, copy=True), "r": np.array(state.r, copy=True),
        "pi_Q": np.array(state.p_Q, copy=True), "pi_r": np.array(state.p_r, copy=True),
        "phi0": np.array(state.phi0, copy=True), "phi1": np.array(state.phi1, copy=True),
        "W": np.array(pair.geometry_map, copy=True), "weights": np.array(pair.weights, copy=True),
        "source_columns": np.array(pair.source_columns, copy=True),
        "reference_columns": np.array(pair.reference_columns, copy=True),
    }
    report = dict(report, mode="prepared_parent", producers=pins, producing_commit=commit,
                  evidence_written=True, input_hashes=report["authentication"]["input_hashes"])
    report["input_hashes"] = dict(report["input_hashes"], **{
        "lab/results/development/nsc-discovery-parent-traction-v1/measurement.json": report["authentication"]["measurement_sha256"],
        "lab/results/development/nsc-discovery-parent-traction-v1/measurement.npz": report["authentication"]["payload_sha256"],
    })
    for name in ('prepare.json','prepare.npz','prediction.json','prediction.npz'):
        path=collar.OUTPUT/name;report['input_hashes'][str(path.relative_to(ROOT))]=_sha_file(path)
    return _write(directory, report, arrays)


def pilot(directory, *, population=0, sign=1, cpu_limit=30., producer_commit=None):
    """Root-only explicit two-resolution feasibility, one aggregate CPU allowance."""
    _require_case(64,population,sign,cpu_limit);directory=Path(directory)
    ledger=directory/'pilot.json'
    if ledger.exists() or any((directory/f'nf{nf}'/'parent.json').exists() or (directory/f'nf{nf}'/'parent.npz').exists() for nf in (64,128)):
        raise FileExistsError('immutable parent pilot prefix already exists')
    started=_cpu_time();pins=source_hashes();commit=producer_commit or producing_commit(pins)
    if not commit:raise ValueError('parent pilot requires a frozen producer commit')
    for path,digest in pins.items():provenance.resolve_pinned_source_bytes(ROOT,path,digest,commit=commit)
    remaining=lambda:float(cpu_limit)-(_cpu_time()-started)
    families={}
    for nf in (64,128):
        if remaining()<=0:raise RuntimeError('aggregate pilot budget exhausted during seed comparison')
        families[nf]=common_k(nf,cpu_limit=remaining())
    common=max(families[nf]['k'] for nf in families);records=[]
    for index,nf in enumerate((64,128)):
        allowance=(remaining()-.1)/(2-index)
        if allowance<=0:break
        record=prepare(directory/f'nf{nf}',nf=nf,population=population,sign=sign,k_override=common,
                       cpu_limit=allowance,producer_commit=commit)
        records.append({'nf':nf,'directory':f'nf{nf}','parent_json_sha256':_sha_file(directory/f'nf{nf}'/'parent.json'),
            'k':record['k'],'kmin':record['kmin'],'kinetic_anchor':None if record['diagnosis'] is None else record['diagnosis']['anchor'],
            'converged':record['converged'],'blocker':record['blocker'],'CPU_seconds':record['CPU_seconds']})
    elapsed=_cpu_time()-started
    result={'schema':SCHEMA+'-PILOT','mode':'bounded_preparation_pilot','producing_commit':commit,'producers':pins,
            'population':population,'sign':sign,'common_k':common,'same_declared_offset':True,
            'seed_families':{str(nf):families[nf] for nf in families},'records':records,
            'aggregate_CPU_seconds':elapsed,'CPU_limit_seconds':float(cpu_limit),'complete':len(records)==2,
            'budget_exceeded':elapsed>cpu_limit,'trajectory':False,'evolved':False,
            'status':'NUMERICAL_PREPARATIONS_RETAINED','global_parent_PASS':False,
            'refinement_scope':'same declared k; actual projected kinetic anchors and source/projection readouts are reported separately'}
    if pins!=source_hashes():raise ValueError('parent pilot producer changed')
    directory.mkdir(parents=True,exist_ok=True)
    data=(json.dumps(_plain(result),sort_keys=True,indent=2,allow_nan=False)+'\n').encode()
    if len(data)>CHUNK_LIMIT:raise RuntimeError('pilot ledger exceeds 64 MiB')
    with os.fdopen(os.open(ledger,os.O_CREAT|os.O_EXCL|os.O_WRONLY,0o444),'wb') as handle:handle.write(data)
    return result


def _load(directory):
    directory = Path(directory)
    record = json.loads((directory / "parent.json").read_text())
    payload = directory / "parent.npz"
    if record.get("schema") != SCHEMA or _sha_file(payload) != record["payload_sha256"]:
        raise ValueError("parent record binding failed")
    with np.load(payload, allow_pickle=False) as saved:
        arrays = {key: saved[key].copy() for key in saved.files}
    for key, digest in record["array_sha256"].items():
        if hashlib.sha256(np.ascontiguousarray(arrays[key]).tobytes()).hexdigest() != digest:
            raise ValueError("parent array binding failed: " + key)
    for path,digest in record['input_hashes'].items():
        if _sha_file(ROOT/path)!=digest:raise ValueError('parent frozen input drift: '+path)
    if not record.get('producing_commit'):raise ValueError('parent evidence lacks frozen producer')
    for path,digest in record['producers'].items():
        provenance.resolve_pinned_source_bytes(ROOT,path,digest,commit=record['producing_commit'])
    return record, arrays


def check(directory):
    """Read-only replay of stored hashes and the saved state's constraints. Writes nothing."""
    directory=Path(directory)
    if (directory/'pilot.json').exists():
        ledger=json.loads((directory/'pilot.json').read_text())
        for path,digest in ledger['producers'].items():provenance.resolve_pinned_source_bytes(ROOT,path,digest,commit=ledger['producing_commit'])
        children=[]
        for row in ledger['records']:
            if _sha_file(directory/row['directory']/'parent.json')!=row['parent_json_sha256']:raise ValueError('pilot child-record binding failed')
            child,arrays=_load(directory/row['directory'])
            if child['k']!=ledger['common_k']:raise ValueError('pilot common offset differs')
            children.append(check(directory/row['directory']))
        return {'schema':SCHEMA+'-PILOT','mode':'read_only_check','ok':True,'evolved':False,'bytes_written':0,
                'common_k':ledger['common_k'],'children':children,'aggregate_CPU_seconds':ledger['aggregate_CPU_seconds']}
    record, arrays = _load(directory)
    report = {"schema": SCHEMA, "mode": "read_only_check", "ok": True, "evolved": False, "bytes_written": 0}
    if record["producers"] != source_hashes():
        report.update(mode="historical_authentication_only", numerical_replay=False,
                      reason="current producers differ from the frozen parent record")
        return report
    weights = arrays["weights"]
    gram = arrays["phi0"].conj().T @ arrays["phi0"] + arrays["phi1"].conj().T @ arrays["phi1"]
    eigenvalues = _weighted_eigenvalues(gram, weights)
    observer_gap=_maximum(arrays['reference_columns'].conj().T@arrays['reference_columns']-np.eye(2))
    if observer_gap>1e-10:raise ValueError('parent local observer is not orthonormal')
    diagnosis=record.get('diagnosis')
    if diagnosis is not None and _maximum(eigenvalues - np.asarray(diagnosis["CAR_eigenvalues"], float)) > 1e-10:
        raise ValueError("stored CAR eigenvalues do not match the saved columns")
    grid = backend.make_fft_grid(galerkin.build_grid(int(record["nf"]), gauge="conformal"))
    from dataclasses import replace as replace_dataclass
    grid.fine = replace_dataclass(grid.fine, occupations=np.array(weights, float, copy=True))
    pair = ParentPair(grid, arrays["W"], weights, arrays["reference_columns"], arrays["source_columns"],
                      record["weight_metadata"], {"map": "identity"}, tuple(record["intervals"]["child"]),
                      tuple(record["intervals"]["parent"]), tuple(record["clock_locations"]))
    state = leading.State(arrays["Q"], arrays["r"], arrays["pi_Q"], arrays["pi_r"], arrays["phi0"], arrays["phi1"])
    if not np.array_equal(np.vstack((state.phi0,state.phi1)),arrays['source_columns']):raise ValueError('parent source columns differ from saved initial fields')
    measured = leading.constraints(pair, state)
    if diagnosis is not None:
        stored = diagnosis["leading_constraints"]
        for key in ("raw_C_max", "D_max", "source_current_max"):
            if abs(measured[key] - stored[key]) > 1e-8 * max(1.0, abs(stored[key])):
                raise ValueError("stored constraint summary does not match a read-only replay")
    else:report['finite_checkpoint_diagnostics_recomputed']=True
    report.update(numerical_replay=True, converged=record["converged"], blocker=record["blocker"])
    return report
