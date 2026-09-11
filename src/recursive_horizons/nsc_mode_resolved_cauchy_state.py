"""Mode-resolved Gaussian Cauchy-state serialization for the NSC source.

The accepted source kernels remain byte-immutable.  A return-frame adapter
captures their actual final canonical vectors during one authorized base run,
then stores physical 2x2 covariances and the existing reference allocation.
The public loader is ordinary, blockwise, pickle-free code; extraction
instrumentation is not used by consumers.
"""
from __future__ import annotations

from dataclasses import dataclass
from io import BytesIO
import hashlib
from math import pi, sqrt
from pathlib import Path
import sys
import zipfile

import numpy as np
from scipy.special import expit

from .nsc_angular_stress import (
    curvature_squared_tensor,
    physical_source,
    profile_jets,
    static_cylinder_density,
)
from .nsc_charged_ctp_neck import (
    ChargedCTPNeckConfig,
    _massive_angular_source,
)
from .nsc_compact_ctp_neck import (
    CompactCTPConfig,
    compact_ctp_source,
    compact_mode_source,
)
from .nsc_horizon_source import conformal_stress, pg_source
from .nsc_unruh_state import frequency_grid


FAMILY_LLL = 0
FAMILY_ANGULAR = 1
FAMILY_COMPACT = 2
ARRAY_NAMES = (
    "channel_index",
    "frequency",
    "quadrature_weight",
    "covariance_seed",
    "covariance_reference_seed",
    "hx_seed",
    "hy_seed",
    "hz_seed",
    "reference_energy",
    "reference_parallel",
    "transmission",
    "occupation_outgoing",
    "occupation_incoming",
)


def _copy_local(value):
    return np.array(value, copy=True) if isinstance(value, np.ndarray) else value


def _capture_one_return(function, call, names):
    """Capture selected locals without changing the hash-pinned v1 kernel."""
    captured = {}
    code = function.__code__

    def profile(frame, event, _argument):
        if event == "return" and frame.f_code is code:
            captured.update({
                name: _copy_local(frame.f_locals[name])
                for name in names
            })

    previous = sys.getprofile()
    sys.setprofile(profile)
    try:
        result = call()
    finally:
        sys.setprofile(previous)
    if set(captured) != set(names):
        raise RuntimeError("immutable charged kernel did not expose expected state locals")
    return result, captured


def _capture_many_returns(function, call, names):
    rows = []
    code = function.__code__

    def profile(frame, event, _argument):
        if event == "return" and frame.f_code is code:
            rows.append({
                name: _copy_local(frame.f_locals[name])
                for name in names
            })

    previous = sys.getprofile()
    sys.setprofile(profile)
    try:
        result = call()
    finally:
        sys.setprofile(previous)
    if not rows or any(set(row) != set(names) for row in rows):
        raise RuntimeError("immutable compact kernel did not expose expected state locals")
    return result, rows


def _outer_covariance(upper, lower, incoming_upper, incoming_lower):
    occupied = np.stack([upper, lower], axis=-1)
    incoming = np.stack([incoming_upper, incoming_lower], axis=-1)
    return (
        np.einsum("...i,...j->...ij", occupied, occupied.conj())
        + np.einsum("...i,...j->...ij", incoming, incoming.conj())
    )


def extract_charged_angular(config: ChargedCTPNeckConfig):
    names = (
        "frequencies", "weights", "masses_1d", "degeneracies",
        "canonical_positive", "canonical_negative",
        "canonical_incoming_positive", "canonical_incoming_negative",
        "energy2", "energy4", "pressure2", "pressure4",
        "momentum", "transmission", "outgoing_occupation",
        "incoming_occupation", "jets",
    )
    result, values = _capture_one_return(
        _massive_angular_source,
        lambda: _massive_angular_source(config),
        names,
    )
    covariance = _outer_covariance(
        values["canonical_positive"], values["canonical_negative"],
        values["canonical_incoming_positive"],
        values["canonical_incoming_negative"],
    )
    channels = []
    for index, (mass, degeneracy) in enumerate(zip(
        values["masses_1d"], values["degeneracies"]
    )):
        # The immutable reducer stores one longitudinal row and broadcasts it
        # over every angular mass.
        momentum = values["momentum"][0]
        channels.append({
            "family": FAMILY_ANGULAR,
            "compact_level": 0,
            "angular_level": index+1,
            "angular_eigenvalue": float(mass),
            "compact_mass": 0.0,
            "degeneracy": int(round(float(degeneracy))),
            "copy_count": 1,
            "reference_scheme": "massive-angular E2+E4/P2+P4",
            "frequency": values["frequencies"],
            "weight": values["weights"],
            "covariance": covariance[index],
            "reference_covariance": np.zeros_like(covariance[index]),
            "hx": np.full_like(momentum, float(mass)),
            "hy": np.zeros_like(momentum),
            "hz": momentum,
            "reference_energy": (
                values["energy2"][index]/mass
                + values["energy4"][index]/mass**3
            ),
            "reference_parallel": (
                values["pressure2"][index]/mass
                + values["pressure4"][index]/mass**3
            ),
            "transmission": values["transmission"][index],
            "outgoing": values["outgoing_occupation"],
            "incoming": values["incoming_occupation"],
            "factor": float(degeneracy/(4*pi*pi*sqrt(values["jets"][0]))),
        })
    return result, channels, values["jets"]


def extract_positive_compact(config: CompactCTPConfig):
    names = (
        "compact_level", "angular_level", "compact_mass", "angular",
        "degeneracy", "compact_copies", "frequencies", "weights",
        "upper", "lower", "incoming_upper", "incoming_lower",
        "differences", "hx", "hy", "hz", "transmission",
        "outgoing_occupation", "incoming_occupation", "factor",
    )
    result, captured = _capture_many_returns(
        compact_mode_source,
        lambda: compact_ctp_source(config),
        names,
    )
    channels = []
    for values in captured:
        covariance = _outer_covariance(
            values["upper"], values["lower"],
            values["incoming_upper"], values["incoming_lower"],
        )
        reference = covariance-values["differences"]
        channels.append({
            "family": FAMILY_COMPACT,
            "compact_level": int(values["compact_level"]),
            "angular_level": int(values["angular_level"]),
            "angular_eigenvalue": float(values["angular"]),
            "compact_mass": float(values["compact_mass"]),
            "degeneracy": int(values["degeneracy"]),
            "copy_count": int(values["compact_copies"]),
            "reference_scheme": "positive-compact fourth-order superadiabatic Bloch projector",
            "frequency": values["frequencies"],
            "weight": values["weights"],
            "covariance": covariance,
            "reference_covariance": reference,
            "hx": values["hx"],
            "hy": values["hy"],
            "hz": values["hz"],
            "reference_energy": np.zeros_like(values["frequencies"]),
            "reference_parallel": np.zeros_like(values["frequencies"]),
            "transmission": values["transmission"],
            "outgoing": values["outgoing_occupation"],
            "incoming": values["incoming_occupation"],
            "factor": float(values["factor"]),
        })
    channels.sort(key=lambda row: (row["compact_level"], row["angular_level"]))
    return result, channels


def extract_lll(config: ChargedCTPNeckConfig):
    frequencies, weights = frequency_grid(
        config.frequency_edges, config.points_per_frequency_interval
    )
    outgoing = expit(-2*pi*frequencies/config.surface_gravity)
    incoming = expit(
        -2*pi*frequencies/(config.omega*config.surface_gravity)
    )
    covariance = np.zeros((len(frequencies), 2, 2), complex)
    covariance[:, 0, 0] = 1.0-outgoing
    covariance[:, 1, 1] = incoming
    reference = np.zeros_like(covariance)
    reference[:, 0, 0] = 1.0
    return {
        "family": FAMILY_LLL,
        "compact_level": 0,
        "angular_level": 0,
        "angular_eigenvalue": 0.0,
        "compact_mass": 0.0,
        "degeneracy": abs(config.magnetic_flux),
        "copy_count": 1,
        "reference_scheme": "transparent LLL filled-sea covariance plus analytic conformal local term",
        "frequency": frequencies,
        "weight": weights,
        "covariance": covariance,
        "reference_covariance": reference,
        "hx": np.zeros_like(frequencies),
        "hy": np.zeros_like(frequencies),
        "hz": -frequencies/sqrt(3*pi/2-1),
        "reference_energy": np.zeros_like(frequencies),
        "reference_parallel": np.zeros_like(frequencies),
        "transmission": np.ones_like(frequencies),
        "outgoing": outgoing,
        "incoming": incoming,
        "factor": abs(config.magnetic_flux)/pi,
    }


def pack_channels(channels):
    channels = sorted(
        channels,
        key=lambda row: (row["compact_level"], row["angular_level"]),
    )
    offsets = [0]
    arrays = {name: [] for name in ARRAY_NAMES}
    metadata = []
    for channel_index, row in enumerate(channels):
        count = len(row["frequency"])
        offsets.append(offsets[-1]+count)
        arrays["channel_index"].append(np.full(count, channel_index, dtype="<i4"))
        arrays["frequency"].append(np.asarray(row["frequency"], dtype="<f8"))
        arrays["quadrature_weight"].append(np.asarray(row["weight"], dtype="<f8"))
        arrays["covariance_seed"].append(np.asarray(row["covariance"], dtype="<c16"))
        arrays["covariance_reference_seed"].append(
            np.asarray(row["reference_covariance"], dtype="<c16")
        )
        for target, source in (
            ("hx_seed", "hx"), ("hy_seed", "hy"), ("hz_seed", "hz"),
            ("reference_energy", "reference_energy"),
            ("reference_parallel", "reference_parallel"),
            ("transmission", "transmission"),
            ("occupation_outgoing", "outgoing"),
            ("occupation_incoming", "incoming"),
        ):
            arrays[target].append(np.asarray(row[source], dtype="<f8"))
        metadata.append({
            "index": channel_index,
            "family": int(row["family"]),
            "compact_level": int(row["compact_level"]),
            "angular_level": int(row["angular_level"]),
            "angular_eigenvalue": float(row["angular_eigenvalue"]),
            "compact_mass": float(row["compact_mass"]),
            "degeneracy": int(row["degeneracy"]),
            "copy_count": int(row["copy_count"]),
            "factor": float(row["factor"]),
            "reference_scheme": row["reference_scheme"],
            "sample_offset": offsets[-2],
            "sample_count": count,
        })
    packed = {
        name: np.concatenate(values, axis=0)
        for name, values in arrays.items()
    }
    packed["sample_offsets"] = np.asarray(offsets, dtype="<i8")
    return packed, metadata


def _child_tensor_from_two_d(metric_a, metric_beta, uu, uv, vv):
    area = 4*pi
    pg = pg_source(metric_a, metric_beta, 1.0, uu, uv, vv)
    coordinate = np.array([
        [pg["T_tau_tau_2D"], pg["T_tau_rho_2D"]],
        [pg["T_tau_rho_2D"], pg["T_rho_rho_2D"]],
    ])/area
    pg_coframe = np.array([[1.0, 0.0], [metric_beta, 1.0]])
    inverse_pg = np.linalg.inv(pg_coframe)
    pg_frame = inverse_pg.T@coordinate@inverse_pg
    interior_f = -metric_a
    child_from_pg = np.array([
        [metric_beta/sqrt(interior_f), -1.0/sqrt(interior_f)],
        [-1.0/sqrt(interior_f), metric_beta/sqrt(interior_f)],
    ])
    inverse_child = np.linalg.inv(child_from_pg)
    return inverse_child.T@pg_frame@inverse_child


def reconstruct_seed_tensor(arrays, channels, compact_local):
    sigma1 = np.array([[0, 1], [1, 0]], complex)
    sigma2 = np.array([[0, -1j], [1j, 0]], complex)
    sigma3 = np.diag([1, -1]).astype(complex)
    angular_energy, angular_parallel, angular_power = [], [], 0.0
    compact_total = dict(rho=0.0, T01=0.0, p_parallel=0.0,
                         p_sphere=0.0, parent_Killing_power=0.0)
    lll_channel = None
    for channel in channels:
        start = channel["sample_offset"]
        end = start+channel["sample_count"]
        sl = slice(start, end)
        c = arrays["covariance_seed"][sl]
        w = arrays["quadrature_weight"][sl]
        f = arrays["frequency"][sl]
        if channel["family"] == FAMILY_LLL:
            lll_channel = (channel, sl)
        elif channel["family"] == FAMILY_ANGULAR:
            mass = channel["angular_eigenvalue"]
            momentum = arrays["hz_seed"][sl]
            energy = np.sqrt(mass*mass+momentum*momentum)
            cosine = mass/np.sqrt(2*energy*(energy-momentum))
            sine = np.sqrt((energy-momentum)/(2*energy))
            rotation = np.empty((len(f), 2, 2))
            rotation[:, 0, 0] = cosine
            rotation[:, 0, 1] = sine
            rotation[:, 1, 0] = -sine
            rotation[:, 1, 1] = cosine
            ce = np.einsum("fai,fij,fbj->fab", rotation, c, rotation)
            occupation = ce[:, 0, 0].real
            tangent = -2*ce[:, 0, 1].real
            trace = np.trace(c, axis1=1, axis2=2).real
            excess = 2*energy*occupation+energy*(1-trace)
            pressure = (
                momentum*momentum/energy**2*excess
                + momentum*mass/energy*tangent
            )
            factor = channel["factor"]
            angular_energy.append(factor*np.dot(
                w, excess-arrays["reference_energy"][sl]
            ))
            angular_parallel.append(factor*np.dot(
                w, pressure-arrays["reference_parallel"][sl]
            ))
            angular_power += channel["degeneracy"]/pi*np.dot(
                w*f, arrays["transmission"][sl]
                *(arrays["occupation_outgoing"][sl]
                  -arrays["occupation_incoming"][sl])
            )
        else:
            difference = c-arrays["covariance_reference_seed"][sl]
            h = (
                arrays["hx_seed"][sl, None, None]*sigma1
                + arrays["hy_seed"][sl, None, None]*sigma2
                + arrays["hz_seed"][sl, None, None]*sigma3
            )
            energy_density = np.einsum("fij,fji->f", difference, h).real
            trace3 = np.einsum("fij,ji->f", difference, sigma3).real
            trace2 = np.einsum("fij,ji->f", difference, sigma2).real
            factor = channel["factor"]
            compact_total["rho"] += factor*np.dot(w, energy_density)
            compact_total["p_parallel"] += factor*np.dot(
                w, arrays["hz_seed"][sl]*trace3
            )
            compact_total["p_sphere"] += factor*np.dot(
                w, channel["angular_eigenvalue"]*trace2/2
            )
            power = channel["copy_count"]*channel["degeneracy"]/pi*np.dot(
                w*f, arrays["transmission"][sl]
                *(arrays["occupation_outgoing"][sl]
                  -arrays["occupation_incoming"][sl])
            )
            compact_total["parent_Killing_power"] += power

    if lll_channel is None:
        raise ArithmeticError("LLL channel absent")
    lll, sl = lll_channel
    w = arrays["quadrature_weight"][sl]
    f = arrays["frequency"][sl]
    outgoing = arrays["occupation_outgoing"][sl]
    incoming = arrays["occupation_incoming"][sl]
    t_u = lll["factor"]*np.dot(w*f, outgoing)
    t_v = lll["factor"]*np.dot(w*f, incoming)
    metric_a = 1-3*pi/2
    metric_beta = sqrt(1-metric_a)
    geometry = conformal_stress(
        metric_a, 6.0, -3*pi, 0.0, 0.0,
        central_charge=lll["degeneracy"],
    )
    lll_tensor = _child_tensor_from_two_d(
        metric_a, metric_beta,
        geometry[0]+t_u, geometry[1], geometry[2]+t_v,
    )

    jets = profile_jets(pi/2)
    h_energy, h_parallel = curvature_squared_tensor(jets)
    harmonic = float(__import__("mpmath").euler)/2
    bar_rho = (
        static_cylinder_density()
        + harmonic*h_energy/(480*pi*pi)
        + float(sum(angular_energy))
    )
    bar_parallel = (
        -static_cylinder_density()
        - harmonic*h_parallel/(480*pi*pi)
        + float(sum(angular_parallel))
    )
    w0, w1, w2, w3, w4 = jets
    scalar2 = -w2
    box_scalar = -w0*w4-w1*w3
    trace = (
        -(scalar2-2)**2/60-11*scalar2/90-box_scalar/30
    )/(16*pi*pi)
    angular_source = physical_source({
        "q": pi/2, "rho_coordinate": 0.0, "sphere_radius": 1.0,
        "W_jets": jets.tolist(), "bar_rho": bar_rho,
        "bar_p_parallel": bar_parallel, "bar_trace": trace,
        "bar_p_sphere": (bar_rho-bar_parallel-trace)/2,
    })
    interior_f = 3*pi/2-1
    angular_t01 = -angular_power/(4*pi*interior_f)
    compact_total["T01"] = -compact_total["parent_Killing_power"]/(4*pi*interior_f)
    total = {
        "rho": angular_source["rho"]+lll_tensor[0, 0]+compact_local["rho"]+compact_total["rho"],
        "T01": angular_t01+lll_tensor[0, 1]+compact_total["T01"],
        "p_parallel": angular_source["p_parallel"]+lll_tensor[1, 1]+compact_local["p_parallel"]+compact_total["p_parallel"],
        "p_sphere": angular_source["p_sphere"]+compact_local["p_sphere"]+compact_total["p_sphere"],
        "parent_Killing_power": angular_power-4*pi*interior_f*lll_tensor[0, 1]+compact_total["parent_Killing_power"],
    }
    details = {
        "lll_quadrature": {"t_u": float(t_u), "t_v": float(t_v)},
        "charged_angular": {
            "rho": angular_source["rho"], "T01": angular_t01,
            "p_parallel": angular_source["p_parallel"],
            "p_sphere": angular_source["p_sphere"],
            "parent_Killing_power": angular_power,
        },
        "positive_compact": compact_total,
    }
    return {key: float(value) for key, value in total.items()}, details


def deterministic_npz_bytes(arrays):
    output = BytesIO()
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_STORED) as archive:
        for name in sorted(arrays):
            value = np.ascontiguousarray(arrays[name])
            if value.dtype.hasobject:
                raise TypeError("object arrays are forbidden")
            stream = BytesIO()
            np.lib.format.write_array(stream, value, allow_pickle=False)
            info = zipfile.ZipInfo(name+".npy", date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_STORED
            info.create_system = 3
            info.external_attr = 0o100644 << 16
            archive.writestr(info, stream.getvalue())
    return output.getvalue()


def array_digest(name, value):
    value = np.ascontiguousarray(value)
    header = f"{name}\0{value.dtype.str}\0{value.shape}\0C\0".encode()
    return hashlib.sha256(header+value.tobytes(order="C")).hexdigest()


@dataclass(frozen=True)
class ModeResolvedCauchyState:
    arrays: dict[str, np.ndarray]
    channels: tuple[dict, ...]

    @classmethod
    def load(cls, path: Path, channels):
        with np.load(path, allow_pickle=False) as payload:
            names = set(payload.files)
            expected = set(ARRAY_NAMES) | {"sample_offsets"}
            if names != expected:
                raise ValueError("payload array names differ from schema")
            arrays = {name: np.array(payload[name], copy=True) for name in payload.files}
        state = cls(arrays=arrays, channels=tuple(channels))
        state.validate()
        return state

    def validate(self, tolerance=2e-10):
        offsets = self.arrays["sample_offsets"]
        if offsets.dtype.kind != "i" or offsets[0] != 0 or np.any(np.diff(offsets) <= 0):
            raise ValueError("invalid channel offsets")
        samples = int(offsets[-1])
        if len(self.channels) != len(offsets)-1:
            raise ValueError("channel metadata and offsets differ")
        if self.arrays["covariance_seed"].shape != (samples, 2, 2):
            raise ValueError("covariance shape differs from schema")
        for name, value in self.arrays.items():
            if value.dtype.hasobject or not np.isfinite(value).all():
                raise ValueError(f"unsafe or nonfinite array: {name}")
            if name != "sample_offsets" and value.shape[0] != samples:
                raise ValueError(f"sample count differs: {name}")
        covariance = self.arrays["covariance_seed"]
        hermitian = float(np.max(abs(covariance-covariance.swapaxes(1, 2).conj())))
        eigenvalues = np.linalg.eigvalsh(covariance)
        if hermitian > tolerance or eigenvalues.min() < -tolerance or eigenvalues.max() > 1+tolerance:
            raise ValueError("physical covariance violates Hermiticity/CAR")
        return {
            "samples": samples,
            "channels": len(self.channels),
            "hermiticity_residual": hermitian,
            "minimum_eigenvalue": float(eigenvalues.min()),
            "maximum_eigenvalue": float(eigenvalues.max()),
        }
