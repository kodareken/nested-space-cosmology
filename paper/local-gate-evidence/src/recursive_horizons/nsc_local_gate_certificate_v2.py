"""Decision arithmetic for the source-fixed local incoming-gate certificate v2."""
from fractions import Fraction


SCHEMA = "NSC-LOCAL-INCOMING-GATE-CERTIFICATE-v2"
CONSTRAINTS = ("N", "beta")
ERROR_COMPONENTS = (
    "field_space_time",
    "changed_history_UV_tail",
    "baseline_low_subgap",
    "upstream",
    "energy_interpolation",
    "covered_regions",
    "phase_value",
    "between_node",
    "arithmetic",
)
REQUIRED_MUTATIONS = (
    "drop-delta-C",
    "drop-coherence",
    "double-weights",
    "k=-E",
    "drop-minus-E-sector",
)
TOLERANCE = Fraction(3, 10**11)


def exact_nonnegative(value, name):
    if value is None:
        return None
    if isinstance(value, bool):
        raise ValueError(name + " must be a nonnegative numerical bound")
    if isinstance(value, dict) and set(value) == {"mantissa", "exponent"}:
        exponent = int(value["exponent"])
        value = Fraction(int(value["mantissa"]) * (1 << max(exponent, 0)),
                         1 << max(-exponent, 0))
    else:
        value = Fraction(str(value)) if isinstance(value, float) else Fraction(value)
    if value < 0:
        raise ValueError(name + " must be nonnegative")
    return value


def _pair(value, name):
    if value is None or len(value) != 2:
        raise ValueError(name + " must contain N and beta")
    return tuple(exact_nonnegative(item, name) for item in value)


def component_totals(components):
    if set(components) != set(ERROR_COMPONENTS):
        missing = sorted(set(ERROR_COMPONENTS) - set(components))
        extra = sorted(set(components) - set(ERROR_COMPONENTS))
        raise ValueError(f"nine-component error schema required; missing={missing}, extra={extra}")
    missing = [name for name in ERROR_COMPONENTS if components[name] is None]
    if missing:
        return None, missing
    pairs = {name: _pair(components[name], name) for name in ERROR_COMPONENTS}
    return tuple(sum((pairs[name][index] for name in ERROR_COMPONENTS), Fraction())
                 for index in range(2)), []


def mutation_gate(results):
    missing = [name for name in REQUIRED_MUTATIONS if results.get(name) is not True]
    return not missing, missing


def existence_decision(residual_upper, components, mutations, *, tolerance=TOLERANCE):
    residual = _pair(residual_upper, "continuous residual upper")
    errors, missing_components = component_totals(components)
    mutations_pass, missing_mutations = mutation_gate(mutations)
    if errors is None or not mutations_pass:
        return {
            "verdict": "OPEN",
            "residual_upper": residual,
            "error_totals": errors,
            "missing_components": missing_components,
            "missing_mutations": missing_mutations,
            "criterion_upper": None,
        }
    criterion = tuple(residual[index] + errors[index] for index in range(2))
    verdict = "EXISTENCE" if all(value <= tolerance for value in criterion) else "OPEN"
    return {
        "verdict": verdict,
        "residual_upper": residual,
        "error_totals": errors,
        "missing_components": [],
        "missing_mutations": [],
        "criterion_upper": criterion,
    }


def nonexistence_decision(obstruction, *, tolerance=TOLERANCE):
    if obstruction is None:
        return {"verdict": "OPEN", "reason": "no class-wide obstruction supplied"}
    required = {"necessary_relation", "entire_declared_class", "separation_lower",
                "error_upper", "optimizer_failure_used"}
    if set(obstruction) != required:
        raise ValueError("complete obstruction contract required")
    separation = exact_nonnegative(obstruction["separation_lower"], "separation lower")
    error = exact_nonnegative(obstruction["error_upper"], "obstruction error")
    valid = (obstruction["necessary_relation"] is True
             and obstruction["entire_declared_class"] is True
             and obstruction["optimizer_failure_used"] is False
             and separation > error)
    return {
        "verdict": "NON_EXISTENCE" if valid else "OPEN",
        "separation_lower": separation,
        "error_upper": error,
        "strict_separation_lower": separation - error,
        "necessary_relation": obstruction["necessary_relation"],
        "entire_declared_class": obstruction["entire_declared_class"],
        "optimizer_failure_used": obstruction["optimizer_failure_used"],
        "tolerance_is_not_the_obstruction_threshold": tolerance,
    }


def serialize_fraction(value):
    if value is None:
        return None
    if isinstance(value, dict):
        return {key: serialize_fraction(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [serialize_fraction(item) for item in value]
    if isinstance(value, Fraction):
        return {"numerator": str(value.numerator), "denominator": str(value.denominator)}
    return value


def build_certificate(*, mode, profile_identity, source_identity, interval, state_law,
                      history_class, numerical_settings, coverage, provenance,
                      residual_upper=None, components=None, mutations=None,
                      obstruction=None):
    if mode not in ("existence", "nonexistence"):
        raise ValueError("certificate mode must be existence or nonexistence")
    if mode == "existence":
        decision = existence_decision(residual_upper, components, mutations)
    else:
        decision = nonexistence_decision(obstruction)
    return serialize_fraction({
        "schema": SCHEMA,
        "verdict": decision["verdict"],
        "mode": mode,
        "profile_identity": profile_identity,
        "source_identity": source_identity,
        "interval": interval,
        "state_law": state_law,
        "history_class": history_class,
        "tolerance": TOLERANCE,
        "decision": decision,
        "error_components": components,
        "mutation_results": mutations,
        "obstruction": obstruction,
        "numerical_settings": numerical_settings,
        "coverage": coverage,
        "provenance": provenance,
        "open_is_not_done": decision["verdict"] == "OPEN",
        "failed_optimization_is_nonexistence": False,
    })
