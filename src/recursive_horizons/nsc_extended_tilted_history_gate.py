"""Composition gate for the extended tilted/frequency-mixing NSC class."""
from __future__ import annotations


def compose_extended_gate(reference, local, interface, old_homogeneous):
    reference_pass = bool(reference["gate"]["owner_residuals_pass"])
    local_bulk_pass = bool(
        local["gate"]["nodewise_action_forces_pass"]
        and local["gate"]["Euler_boxR_endpoint_ledger_pass"]
        and local["gate"]["weyl_bulk_pass"]
    )
    endpoint_pass = bool(local["gate"]["weyl_free_endpoint_completion_pass"])
    interface_pass = bool(
        interface["gate"]["weighted_interface_variation_pass"]
        and interface["gate"]["frequency_mixing_CAR_witness_pass"]
        and interface["gate"]["channel_block_rank_pass"]
    )
    selector_pass = bool(interface["gate"]["physical_frequency_mixing_kernel_selected"])
    regression_pass = bool(old_homogeneous["gate"]["nonexistence_in_declared_class_proved"])
    composable = reference_pass and local_bulk_pass and interface_pass
    complete = composable and endpoint_pass and selector_pass
    return {
        "components_expose_residuals": composable,
        "reference_pass": reference_pass,
        "local_bulk_pass": local_bulk_pass,
        "weyl_endpoint_completion_pass": endpoint_pass,
        "tilted_interface_construction_pass": interface_pass,
        "physical_kernel_selector_pass": selector_pass,
        "killed_homogeneous_regression_pass": regression_pass,
        "extended_functional_complete": complete,
        "stationarity_solve_authorized": complete,
    }
