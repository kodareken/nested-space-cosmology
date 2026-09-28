# Current compact development and repository-level verification routes.

.PHONY: verify-fgc verify-current-development verify-fgc-wave0-algebra verify-fgc-hlt17-components verify-fgc-pro20-ev1-components verify-fgc-pro20-rsrc1-components
.PHONY: fgc-pro20-ev1-frz1 fgc-pro20-ev1-pref1 verify-fgc-pro20-ev1-pref1

# The old verify-fgc recipe stopped before the current TDG8-TDG10 boundary and
# rewrote tracked results. Preserve the familiar name as the complete sealed
# pre-result boundary instead.
verify-fgc: verify-fgc-sf1-foundation

CURRENT_DEVELOPMENT_CERTIFICATES := \
	fgc-tdg9-loc2-pref2 \
	fgc-tdg9-ti2-pref1 \
	fgc-tdg9-ac1-pref1 \
	fgc-tdg9-ur1-pref1 \
	fgc-tdg10-qa1-pref1 \
	fgc-tdg10-qa2-pref1 \
	fgc-tdg11-msel1-frz1 \
	fgc-tdg11-msel1-pref1 \
	fgc-tdg11-imp1 \
	fgc-hlt17-srcq1-rec1-pref1 \
	fgc-def1-stab1-frz1 \
	fgc-def1-stab1-pref1 \
	fgc-sgb1-ctl1-sol1-frz1 \
	fgc-sgb1-ctl1-sol1-pref1 \
	fgc-pro20-ev1-pref1

# Fast current-surface route. It checks compact and bounded synthetic evidence and never
# invokes a historical state-changing runner or opens an ignored campaign.
verify-current-development: $(CURRENT_DEVELOPMENT_CERTIFICATES) verify-fgc-tdg11-msel1-pref1 verify-fgc-tdg11-imp1 verify-fgc-hlt17-srcq1-rec1-pref1 verify-fgc-pro20-ev1-pref1 verify-fgc-wave0-algebra verify-fgc-hlt17-components verify-fgc-pro20-ev1-components verify-fgc-pro20-rsrc1-components
	PYTHONPATH=src:. $(EVOLUTION_PYTHON) -m unittest \
		tests/test_phase_minus1_artifact_catalog.py \
		tests/test_phase_minus1_check_repo_dispatcher.py \
		tests/test_phase_minus1_make_routing.py \
		tests/test_phase_minus1_public_authority.py \
		tests/test_historical_runner_import_boundary.py \
		tests/test_evidence_io.py -v
	$(EVOLUTION_PYTHON) scripts/build_artifact_catalog.py --check
	$(EVOLUTION_PYTHON) scripts/check_repo.py --only-tdg11-imp1

verify-fgc-sf1-foundation: fgc-tdg11-msel1-frz1 fgc-tdg11-msel1-pref1 fgc-tdg11-imp1 fgc-hlt17-srcq1-rec1-pref1 fgc-def1-stab1-frz1 fgc-def1-stab1-pref1 fgc-sgb1-ctl1-sol1-frz1 fgc-sgb1-ctl1-sol1-pref1 fgc-pro20-ev1-pref1

.PHONY: fgc-sgb1-ctl1-sol1-frz1 verify-fgc-sgb1-ctl1-sol1-frz1
.PHONY: fgc-sgb1-ctl1-sol1-pref1 verify-fgc-sgb1-ctl1-sol1-pref1

fgc-sgb1-ctl1-sol1-frz1:
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/reproduce_fgc_sgb1_ctl1_sol1_frz1.py --check

verify-fgc-sgb1-ctl1-sol1-frz1: fgc-sgb1-ctl1-sol1-frz1
	PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. $(PYTHON) -B -m unittest \
		tests/test_fgc_sgb1_ctl1_sol1_frz1_certificate.py \
		tests/test_check_repo_sgb1_ctl1_sol1_frz1.py -v
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/check_repo.py --only-sgb1-ctl1-sol1-frz1

fgc-sgb1-ctl1-sol1-pref1:
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/reproduce_fgc_sgb1_ctl1_sol1_pref1.py --check

verify-fgc-sgb1-ctl1-sol1-pref1: fgc-sgb1-ctl1-sol1-pref1
	PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. $(PYTHON) -B -m unittest \
		tests/test_fgc_sgb1_ctl1_sol1_pref1_binder.py \
		tests/test_check_repo_sgb1_ctl1_sol1_pref1.py -v
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/check_repo.py --only-sgb1-ctl1-sol1-pref1

.PHONY: fgc-def1-stab1-frz1 verify-fgc-def1-stab1-frz1
.PHONY: fgc-def1-stab1-pref1 verify-fgc-def1-stab1-pref1

fgc-def1-stab1-frz1:
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/reproduce_fgc_def1_stab1_frz1.py --check

verify-fgc-def1-stab1-frz1: fgc-def1-stab1-frz1
	PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. $(PYTHON) -B -m unittest \
		tests/test_fgc_def1_stab1_frz1_certificate.py \
		tests/test_check_repo_def1_stab1_frz1.py -v
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/check_repo.py --only-def1-stab1-frz1

fgc-def1-stab1-pref1:
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/reproduce_fgc_def1_stab1_pref1.py --check

verify-fgc-def1-stab1-pref1: fgc-def1-stab1-pref1
	PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. $(PYTHON) -B -m unittest \
		tests/test_fgc_def1_stab1_pref1_binder.py \
		tests/test_check_repo_def1_stab1_pref1.py -v
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/check_repo.py --only-def1-stab1-pref1

# Bounded algebra, interval and synthetic runtime controls, not completed SGB1/DEF1 gates.
verify-fgc-wave0-algebra:
	PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. $(PYTHON) -B -m unittest \
		tests/test_fgc_sgb1_ctl1_constraints.py \
		tests/test_fgc_sgb1_ctl1_time_affinity.py \
		tests/test_fgc_sgb1_ctl1_source.py \
		tests/test_fgc_sgb1_ctl1_family.py \
		tests/test_fgc_sgb1_ctl1_center.py \
		tests/test_fgc_sgb1_ctl1_principal.py \
		tests/test_fgc_sgb1_ctl1_interval_health.py \
		tests/test_fgc_sgb1_ctl1_runtime.py \
		tests/test_fgc_sgb1_ctl1_admission.py \
		tests/test_fgc_sgb1_ctl1_initial_health.py \
		tests/test_fgc_sgb1_ctl1_controls.py \
		tests/test_fgc_sgb1_ctl1_cone.py \
		tests/test_fgc_sgb1_ctl1_continuity.py \
		tests/test_fgc_sgb1_ctl1_trap_refinement.py \
		tests/test_fgc_sgb1_ctl1_trap_refinement2.py \
		tests/test_fgc_sgb1_ctl1_trap_taylor.py \
		tests/test_fgc_sgb1_ctl1_family_principal.py \
		tests/test_fgc_sgb1_ctl1_cell_admission.py \
		tests/test_fgc_sgb1_ctl1_local_symmetrizer.py \
		tests/test_fgc_sgb1_ctl1_trap_barrier.py \
		tests/test_fgc_sgb1_ctl1_sol1.py \
		tests/test_fgc_sgb1_ctl1_sol1_frz1_certificate.py \
		tests/test_check_repo_sgb1_ctl1_sol1_frz1.py \
		tests/test_fgc_sgb1_ctl1_sol1_pref1_binder.py \
		tests/test_check_repo_sgb1_ctl1_sol1_pref1.py \
		tests/test_fgc_def1_stab1.py \
		tests/test_fgc_def1_geometry_error.py \
		tests/test_fgc_def1_stab1_providers.py \
		tests/test_fgc_def1_stab1_qualification.py \
		tests/test_fgc_def1_stab1_freeze_contract.py \
		tests/test_fgc_def1_stab1_frz1_certificate.py \
		tests/test_check_repo_def1_stab1_frz1.py \
		tests/test_fgc_def1_stab1_pref1_binder.py \
		tests/test_check_repo_def1_stab1_pref1.py \
		tests/test_fgc_wave0_instrument_routing.py -v

# C1R1/H17 finite-state and atomic-store components. This is not source,
# origin, resource, binder, or campaign qualification.
verify-fgc-hlt17-components:
	PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. $(EVOLUTION_PYTHON) -B -m unittest \
		tests/test_fgc_tdg11_c1r1_ring.py \
		tests/test_fgc_tdg11_c1r1_enclosure.py \
		tests/test_fgc_tdg11_c1r1_runtime.py \
		tests/test_fgc_hlt17_runtime_member.py \
		tests/test_fgc_hlt17_imp1_cursor.py \
		tests/test_fgc_hlt17_imp1_bridge.py \
		tests/test_fgc_hlt17_member_codec.py \
		tests/test_fgc_hlt17_member_checkpoint.py \
		tests/test_fgc_hlt17_c1r1_integration.py \
		tests/test_fgc_protocol_v19.py \
		tests/test_fgc_hlt17_campaign_store.py -v

# Production-shaped PRO20 persistence, runner policy and prospective authority
# implementation only. No freeze config, live target, source construction,
# namespace, event execution or state advance.
verify-fgc-pro20-ev1-components:
	PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. $(EVOLUTION_PYTHON) -B -m unittest \
		tests/test_fgc_pro20_ev1_protocol.py \
		tests/test_fgc_pro20_ev1_store.py \
		tests/test_fgc_pro20_ev1_runtime.py \
		tests/test_fgc_pro20_ev1_authority.py -v

# Prospective no-store resource isolation only. No seed, namespace, authority,
# parent publisher, status target, or live campaign operation is present.
verify-fgc-pro20-rsrc1-components:
	PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. $(EVOLUTION_PYTHON) -B -m unittest \
		tests/test_fgc_pro20_rsrc1_member.py \
		tests/test_fgc_pro20_rsrc1_attempt.py \
		tests/test_fgc_pro20_rsrc1_isolation.py \
		tests/test_fgc_pro20_rsrc1_seed.py -v

# The one authorized PRO20 run is closed. Compact freeze and independent
# PREF1 verification remain; the live run target is tombstoned.

fgc-pro20-ev1-frz1:
	PYTHONDONTWRITEBYTECODE=1 $(EVOLUTION_PYTHON) -I -B scripts/reproduce_fgc_pro20_ev1_frz1.py --check

fgc-pro20-ev1-pref1:
	PYTHONDONTWRITEBYTECODE=1 $(EVOLUTION_PYTHON) -I -B scripts/reproduce_fgc_pro20_ev1_pref1.py --check

verify-fgc-pro20-ev1-pref1: fgc-pro20-ev1-pref1
	PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. $(EVOLUTION_PYTHON) -B -m unittest \
		tests/test_fgc_pro20_ev1_pref1_certificate.py \
		tests/test_check_repo_pro20_ev1_pref1.py -v

.PHONY: fgc-hlt17-srcq1-rec1-frz1 verify-fgc-hlt17-srcq1-rec1-prelaunch
.PHONY: status-fgc-hlt17-srcq1-rec1
.PHONY: fgc-hlt17-srcq1-rec1-pref1 verify-fgc-hlt17-srcq1-rec1-pref1

FGC_HLT17_SRCQ1_REC1_AUTHORITY_COMMIT ?= $(shell git --no-replace-objects rev-parse --verify HEAD)
FGC_HLT17_SRCQ1_REC1_OUTPUT := $(abspath runs/fgc-2-sf1/hlt17-srcq1-rec1-source-origin-qualification)

fgc-hlt17-srcq1-rec1-frz1:
	PYTHONDONTWRITEBYTECODE=1 $(EVOLUTION_PYTHON) -I -B scripts/reproduce_fgc_hlt17_srcq1_rec1_frz1.py --check

verify-fgc-hlt17-srcq1-rec1-prelaunch: fgc-hlt17-srcq1-rec1-frz1
	PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. $(EVOLUTION_PYTHON) -B -m unittest \
		tests/test_fgc_hlt17_srcq1_rec1_auth1.py \
		tests/test_fgc_hlt17_srcq1_rec1.py -v

status-fgc-hlt17-srcq1-rec1:
	$(EVOLUTION_PYTHON) -I -B scripts/qualify_fgc_hlt17_srcq1_rec1.py \
		--status --authority-commit "$(FGC_HLT17_SRCQ1_REC1_AUTHORITY_COMMIT)"

fgc-hlt17-srcq1-rec1-pref1:
	PYTHONDONTWRITEBYTECODE=1 $(EVOLUTION_PYTHON) -I -B scripts/reproduce_fgc_hlt17_srcq1_rec1_pref1.py --check

verify-fgc-hlt17-srcq1-rec1-pref1: fgc-hlt17-srcq1-rec1-pref1
	PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. $(EVOLUTION_PYTHON) -B -m unittest \
		tests/test_fgc_hlt17_srcq1_rec1_pref1_binder.py \
		tests/test_check_repo_hlt17_srcq1_rec1_pref1.py -v

.PHONY: fgc-tdg11-msel1-frz1 verify-fgc-tdg11-msel1-prelaunch status-fgc-tdg11-msel1
.PHONY: fgc-tdg11-msel1-pref1 verify-fgc-tdg11-msel1-pref1

fgc-tdg11-msel1-frz1:
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/reproduce_fgc_tdg11_msel1_frz1.py --check

verify-fgc-tdg11-msel1-prelaunch: fgc-tdg11-msel1-frz1
	PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. $(EVOLUTION_PYTHON) -B -m unittest \
		tests/test_fgc_tdg11_compensated_rk.py \
		tests/test_fgc_tdg11_rational_complete_c.py \
		tests/test_fgc_tdg11_msel1_reconstruction.py \
		tests/test_fgc_tdg11_msel1_runtime.py \
		tests/test_fgc_tdg11_msel1_contract.py \
		tests/test_fgc_tdg11_msel1_authority.py \
		tests/test_fgc_tdg11_msel1_runner.py \
		tests/test_check_repo_tdg11_msel1_frz1.py -v
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/check_repo.py --only-tdg11-msel1-frz1

# The committed authority verifies this exact HEAD's sole parent and full
# delta. A changed worktree or an existing diagnostic namespace refuses.
FGC_TDG11_MSEL1_AUTHORITY_COMMIT ?= $(shell git --no-replace-objects rev-parse --verify HEAD)

status-fgc-tdg11-msel1:
	$(EVOLUTION_PYTHON) -I -B scripts/run_fgc_tdg11_msel1.py --status --authority-commit "$(FGC_TDG11_MSEL1_AUTHORITY_COMMIT)"

# Compact PREF1 verification is raw/store/shadow/Git-blind. The live
# one-shot construction is never a default Make dependency.
fgc-tdg11-msel1-pref1:
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/reproduce_fgc_tdg11_msel1_pref1.py --check

verify-fgc-tdg11-msel1-pref1: fgc-tdg11-msel1-pref1
	PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. $(EVOLUTION_PYTHON) -B -m unittest \
		tests/test_fgc_tdg11_msel1_pref1_binder.py \
		tests/test_fgc_tdg11_msel1_pref1_protocol.py \
		tests/test_fgc_tdg11_msel1_pref1_reconstruction.py \
		tests/test_fgc_tdg11_msel1_pref1_localization.py \
		tests/test_check_repo_tdg11_msel1_pref1.py -v
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/check_repo.py --only-tdg11-msel1-pref1

.PHONY: fgc-tdg11-imp1 verify-fgc-tdg11-imp1

fgc-tdg11-imp1:
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/reproduce_fgc_tdg11_imp1.py --check

verify-fgc-tdg11-imp1: fgc-tdg11-imp1
	PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. $(EVOLUTION_PYTHON) -B -m unittest \
		tests/test_fgc_tdg11_imp1_enclosure.py \
		tests/test_fgc_tdg11_imp1_ledger.py \
		tests/test_fgc_tdg11_imp1_runtime.py \
		tests/test_fgc_tdg11_imp1_retry_limits.py \
		tests/test_fgc_tdg11_imp1_reproduction.py \
		tests/test_check_repo_tdg11_imp1.py -v
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/check_repo.py --only-tdg11-imp1
