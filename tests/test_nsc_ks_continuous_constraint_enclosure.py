"""257/256 continuous enclosure: coverage, hashes and sample-not-proof."""
from copy import deepcopy
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import json
import math
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from recursive_horizons.nsc_ks_continuous_constraint_enclosure import (
    CELL_COUNT,
    PROFILE_IDENTITY,
    SCHEMA,
    STATE_LAW,
    VERIFICATION_NODE_COUNT,
    EnclosureBinding,
    enclose_continuous_constraints,
    manufactured_linear_enclosures,
    production_between_node_component,
    require_verification_nodes,
    serialize_fraction,
    validate_continuous_enclosure_record,
)
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily
import derive_nsc_ks_continuous_constraint_enclosure as D


def digest(raw):
    return sha256(raw).hexdigest()


def linear_grid():
    nodes = tuple(Fraction(index, 256) for index in range(257))
    manufactured = manufactured_linear_enclosures(nodes)
    return nodes, manufactured["nodal_abs"], manufactured["derivative_uppers"]


def test_linear_enclosure_is_exact_on_257_nodes():
    nodes, values, derivatives = linear_grid()
    result = enclose_continuous_constraints(
        nodes, values, derivatives, interval=(0, 1))
    assert result["node_count"] == VERIFICATION_NODE_COUNT
    assert result["cell_count"] == CELL_COUNT
    assert result["N"]["continuous_maximum_upper"] == serialize_fraction(1)
    assert result["beta"]["continuous_maximum_upper"] == serialize_fraction(1)
    assert result["N"]["between_node_remainder_upper"] == serialize_fraction(0)
    assert result["samples_are_not_proof"]
    assert not result["dense_sampling_used_as_bound"]


def test_missing_cell_and_reordered_nodes_are_rejected():
    nodes, values, derivatives = linear_grid()
    with pytest.raises(ValueError, match="256"):
        enclose_continuous_constraints(
            nodes, values, derivatives[:-1], interval=(0, 1))
    missing = np.array(derivatives, dtype=object, copy=True)
    missing[10, 0] = None
    with pytest.raises(ValueError, match="missing"):
        enclose_continuous_constraints(nodes, values, missing, interval=(0, 1))
    reordered = list(nodes)
    reordered[3], reordered[4] = reordered[4], reordered[3]
    with pytest.raises(ValueError, match="strictly increasing"):
        require_verification_nodes(reordered, (0, 1))
    with pytest.raises(ValueError, match="cover I"):
        require_verification_nodes(nodes, (0, 2))


def test_inconsistent_derivative_nan_and_negative_are_rejected():
    nodes, values, derivatives = linear_grid()
    too_small = np.array(derivatives, dtype=object, copy=True)
    too_small[0, 0] = Fraction(1, 10 ** 9)
    with pytest.raises(ValueError, match="endpoint movement"):
        enclose_continuous_constraints(
            nodes, values, too_small, interval=(0, 1))
    with pytest.raises(ValueError, match="finite"):
        enclose_continuous_constraints(
            nodes, np.array(values, dtype=object), np.full((256, 2), float("nan")),
            interval=(0, 1))
    with pytest.raises(ValueError, match="nonnegative"):
        enclose_continuous_constraints(
            nodes, values, np.full((256, 2), -1), interval=(0, 1))


def test_samples_are_not_a_bound():
    nodes, values, derivatives = linear_grid()
    samples = np.zeros((257, 2))
    with pytest.raises(ValueError, match="samples are not a continuous bound"):
        enclose_continuous_constraints(
            nodes, values, derivatives, interval=(0, 1), samples=samples,
            samples_used_as_bound=True)
    result = enclose_continuous_constraints(
        nodes, values, derivatives, interval=(0, 1), samples=samples)
    assert result["samples_are_not_proof"]
    null = production_between_node_component(samples=samples)
    assert null["bound"] is None
    assert null["certificate_use"] is False
    assert "outward nodal" in null["missing_inputs"][0]


def test_profile_state_law_and_dependency_hash_mismatch():
    family = LocalIncomingFamily(np.zeros((2, 8)))
    identity = profile_identity(family, include_normal_window=True)
    dep = "docs/active-code-map.md"
    dep_hash = sha256((ROOT / dep).read_bytes()).hexdigest()
    binding = EnclosureBinding(
        identity, STATE_LAW, digest(b"source"), family.interval,
        ((dep, dep_hash),))
    with pytest.raises(ValueError, match="state law"):
        EnclosureBinding(
            identity, "frozen C0", digest(b"source"), family.interval,
            ((dep, dep_hash),))
    nodes = family.collocation_nodes(257)
    manufactured = manufactured_linear_enclosures(nodes)
    wrong_profile = EnclosureBinding(
        PROFILE_IDENTITY, STATE_LAW, digest(b"source"), family.interval,
        ((dep, dep_hash),))
    with pytest.raises(ValueError, match="profile identity mismatch"):
        enclose_continuous_constraints(
            nodes, manufactured["nodal_abs"], manufactured["derivative_uppers"],
            interval=family.interval, family=family, binding=wrong_profile)
    invented = EnclosureBinding(
        identity, STATE_LAW, digest(b"source"), family.interval,
        ((dep, digest(b"dep")),))
    with pytest.raises(ValueError, match="does not match file bytes"):
        enclose_continuous_constraints(
            nodes, manufactured["nodal_abs"], manufactured["derivative_uppers"],
            interval=family.interval, family=family, binding=invented)
    wrong_source = EnclosureBinding(
        identity, STATE_LAW, digest(b"other-source"), family.interval,
        ((dep, dep_hash),))
    enclosed = enclose_continuous_constraints(
        nodes, manufactured["nodal_abs"], manufactured["derivative_uppers"],
        interval=family.interval, family=family, binding=binding)
    assert enclosed["node_count"] == 257
    assert wrong_source.source_identity != binding.source_identity


def test_outward_float_rounds_up_and_false_pass_is_rejected():
    from recursive_horizons.nsc_ks_continuous_constraint_enclosure import (
        binary_rational,
        outward_float,
    )
    value = Fraction(1, 10)
    assert binary_rational(outward_float(value)) >= value
    third = Fraction(1, 3)
    assert binary_rational(outward_float(third)) >= third
    fifteenth = Fraction(1, 15)
    exact = binary_rational(outward_float(fifteenth))
    assert exact >= fifteenth
    record = D.calculate()
    assert record["schema"] == SCHEMA
    assert record["production"]["between_node"] is None
    assert record["production"]["arithmetic"] is None
    assert record["diagnostic"]["certificate_use"] is False
    D.check(record)
    forged = deepcopy(record)
    forged["physical_EXISTENCE_certificate"] = True
    with pytest.raises(ValueError, match="EXISTENCE"):
        validate_continuous_enclosure_record(forged)
    forged_pass = deepcopy(record)
    forged_pass["status"] = "PASS: invented"
    with pytest.raises(ValueError, match="OPEN"):
        validate_continuous_enclosure_record(forged_pass)
    forged_diag = deepcopy(record)
    forged_diag["diagnostic"]["certificate_use"] = True
    with pytest.raises(ValueError, match="diagnostic"):
        validate_continuous_enclosure_record(forged_diag)


def test_recorder_binds_current_profile_and_lists_missing_inputs():
    record = D.calculate()
    assert record["profile_identity"] == PROFILE_IDENTITY
    assert record["state_law"] == STATE_LAW
    assert record["verification_node_count"] == 257
    assert record["cell_count"] == 256
    assert len(record["nodes_hex"]) == 257
    assert record["missing_scientific_inputs"]
    assert "between_node" in record["error_components_not_filled"]
    assert record["phase_value_arithmetic_excluded"]
    assert math.isfinite(int(record["diagnostic"]["continuous"]["N"]
                             ["outward_binary64"]["continuous_maximum_upper"]))
    history = json.loads((ROOT / "results/development/nsc-ks-gate-history-lm-broyden.json").read_text())
    family = LocalIncomingFamily(np.asarray(history["history"]["coefficients"], float))
    assert profile_identity(family, include_normal_window=True) == PROFILE_IDENTITY
