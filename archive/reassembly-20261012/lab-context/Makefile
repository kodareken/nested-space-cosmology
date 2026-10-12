PYTHON ?= python3
.DEFAULT_GOAL := help
.PHONY: help verify check test reproduce reproduce-development history clean
help:
	@echo "NSC laboratory: reuse recorded evidence; the approved plan is in the conversation."
	@echo "  make check        check portable links and preserved archive hashes"
	@echo "  make test         run the active test suite explicitly"
	@echo "  make reproduce    explicitly recompute all active scientific records"
	@echo "  make verify       run the integrated scientific verification gate"
	@echo "  make history      explicitly reproduce the frozen scale-closure checkpoint"
	@echo "  make clean        move confirmed caches to recoverable Trash"
check:
	$(PYTHON) -B scripts/check_current.py
test:
	PYTHONPATH=src:. $(PYTHON) -B -m unittest discover -s tests -v
reproduce-development:
	$(PYTHON) -B scripts/check_nsc_compact_interaction.py --check
	$(PYTHON) -B scripts/check_nsc_torsion_uv_map.py --check
	$(PYTHON) -B scripts/check_nsc_flow_compatibility.py --check
	$(PYTHON) -B scripts/check_nsc_charged_sector.py --check
	$(PYTHON) -B scripts/check_nsc_vacuum_charge_matching.py --check
	$(PYTHON) -B scripts/check_nsc_compact_boundary_action.py --check
	$(PYTHON) -B scripts/check_nsc_compact_casimir.py --check
	$(PYTHON) -B scripts/check_nsc_horizon_source.py --check
	$(PYTHON) -B scripts/check_nsc_warped_source.py --check
	$(PYTHON) -B scripts/check_nsc_compact_matching.py --check
	$(PYTHON) -B scripts/check_nsc_gauge_source.py --check
	$(PYTHON) -B scripts/check_nsc_spherical_action.py --check
	$(PYTHON) -B scripts/check_nsc_curvature_eft.py --check
	$(PYTHON) -B scripts/check_nsc_spectral_endpoint.py --check
	$(PYTHON) -B scripts/check_nsc_child_state.py --check
	$(PYTHON) -B scripts/check_nsc_massless_reference.py --check
	$(PYTHON) -B scripts/check_nsc_angular_stress.py --check
	$(PYTHON) -B scripts/check_nsc_unruh_state.py --check
	$(PYTHON) -B scripts/check_nsc_state_regulator.py --check
reproduce: reproduce-development
	$(PYTHON) -B scripts/check_nsc_scale_closure.py --check
	$(PYTHON) -B scripts/check_nsc_regulated_recursion.py --check
	$(PYTHON) -B scripts/check_nsc_radial_spectrum.py --check
	$(PYTHON) -B scripts/check_nsc_geometric_chain.py --check
	$(PYTHON) -B scripts/check_nsc_threshold_response.py --check
	$(PYTHON) -B scripts/check_nsc_boundary_response.py --check
	$(PYTHON) -B scripts/check_nsc_dirac_tetrad.py --check
	$(PYTHON) -B scripts/check_nsc_lorentzian.py --check
	$(PYTHON) -B scripts/check_nsc_smooth_geometry.py --check
	$(PYTHON) -B scripts/check_nsc_covariant_measure.py --check
	$(PYTHON) -B scripts/check_nsc_shape_response.py --check
	$(PYTHON) -B scripts/check_nsc_plateau_conditions.py --check
	$(PYTHON) -B scripts/check_nsc_clock_horizon.py --check
	$(PYTHON) -B scripts/check_nsc_tail_limit.py --check
	$(PYTHON) -B scripts/check_nsc_energy_transfer.py --check
	$(PYTHON) -B scripts/check_nsc_vacuum_work.py --check
	$(PYTHON) -B scripts/check_nsc_observable_bridge.py --check
	$(PYTHON) -B scripts/check_nsc_chiral_boundary.py --check
	$(PYTHON) -B scripts/check_nsc_finite_terms.py --check
	$(PYTHON) -B scripts/check_nsc_covariant_source.py --check
	$(PYTHON) -B scripts/check_nsc_measure_normalization.py --check
	$(PYTHON) -B scripts/check_nsc_influence.py --check
	$(PYTHON) -B scripts/check_nsc_response_matching.py --check
verify: check test reproduce
history:
	$(PYTHON) -B scripts/reproduce_history.py scale-closure
clean:
	$(PYTHON) -B scripts/clean_caches.py
