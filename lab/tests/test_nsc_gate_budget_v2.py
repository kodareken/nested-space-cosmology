"""A numerical field problem must not erase a valid baseline component."""
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import derive_nsc_ks_gate_budget_v2 as B


def baseline():
    return {
        "baseline": {"constraint_order": ["N", "beta"]},
        "partial_error_budget": {
            "constraint_order": ["N", "beta"], "covered_domain_status": "PASS",
            "covered_spectral_regions_action_error_upper": [8.916253820835393e-12, 3.882373567206641e-18],
            "remaining_low_subgap_source_error_bound": None,
        },
    }


def test_reuses_only_covered_baseline_and_does_not_infer_a_complete_budget():
    result = B.assemble(baseline())
    assert result["components"]["covered_regions"]["bound"][0] > 8e-12
    assert result["components"]["covered_regions"]["allocation_excess"][0] > 7e-12
    assert result["components"]["baseline_low_subgap"]["bound"] is None
    assert result["components"]["field_space_time"]["bound"] is None
    assert len(result["missing_components"]) == 8
    assert result["full_error_sum"] is None
    assert result["physical_EXISTENCE_certificate"] is False
    assert result["physical_NONEXISTENCE_certificate"] is False
    assert result["allocation_changed"] is False


@pytest.mark.parametrize("bounds", [[-1e-12, 0], [float("nan"), 0], [0], [0, float("inf")]])
def test_bad_bound_cannot_be_reused(bounds):
    data = baseline()
    data["partial_error_budget"]["covered_spectral_regions_action_error_upper"] = bounds
    with pytest.raises(ValueError, match="bounds required"):
        B.assemble(data)


def test_wrong_scope_or_component_order_is_not_silently_reused():
    data = baseline()
    data["partial_error_budget"]["covered_domain_status"] = "OPEN"
    with pytest.raises(ValueError):
        B.assemble(data)
    data = baseline()
    data["baseline"]["constraint_order"] = ["beta", "N"]
    with pytest.raises(ValueError):
        B.assemble(data)
