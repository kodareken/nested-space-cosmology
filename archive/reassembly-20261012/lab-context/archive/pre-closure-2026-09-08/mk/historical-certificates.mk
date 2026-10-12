.PHONY: test test-analytic test-evolution-restart test-portable-clean-clone reproduce fgc-action fgc-metric-variation fgc-hyp1-reduction fgc-hyp1-symbol fgc-hyp1-modified-harmonic fgc-hyp1-mhg-reference fgc-hyp1-mhg-implicit fgc-hyp1-mhg-propagation fgc-hyp1-fo1-rc1 fgc-hyp1-dom1-qift1 fgc-hyp1-con1-comp1 fgc-eft0-led1 fgc-hyp1-dom2-mode1 fgc-hyp1-dom3-uhyp1 fgc-hyp1-con2-mprop1 fgc-hyp1-con3-cau1 fgc-hyp1-bnd1-md1 fgc-def0-obs1 fgc-eft1-open1 fgc-src1-nl1 fgc-con4-phy1 fgc-ctr1-reg1 fgc-id1-fam1 fgc-dom4-run1 fgc-hyp2-md1 fgc-bnd2-cp1 fgc-hlt1-mon1 fgc-num1-val1 fgc-run1-sym1 fgc-id0-pref1 fgc-cal0-pref2 fgc-pro4-frz1 fgc-hlt2-mon2 fgc-id2-all1 fgc-cal1-pref3 fgc-pro5-frz1 fgc-hlt3-mon3 fgc-cal2-pref4 fgc-pro6-frz1 fgc-hlt4-mon4 fgc-cal3-pref5 fgc-cal3-pref5-full-replay fgc-pro7-frz1 fgc-hlt5-mon5 fgc-cal4-pref6 fgc-cal4-pref6-full-replay fgc-pro8-frz1 fgc-hlt6-mon6 fgc-cal5-pref7 fgc-pro9-frz1 fgc-hlt7-mon7 fgc-cal6-pref8 fgc-pro10-frz1 fgc-hlt8-mon8 fgc-cal7-pref9 fgc-pro11-frz1 fgc-hlt9-mon9 fgc-cal8-pref10 fgc-src2-cap1 fgc-src2-pref11 fgc-src3 fgc-rsp1-frz1 fgc-rsp1-pref12 fgc-src4-vec1 fgc-pro12-frz1 fgc-hlt10-mon10 fgc-cal9-pref13 fgc-rsp2-frz1 fgc-rsp2-pref14 fgc-pro13-frz1 fgc-hlt11-mon11 run-fgc-pro13-calibration fgc-cal10-pref15 fgc-tdg1-frz1 fgc-tdg1-pref16 fgc-tdg2-frz1 fgc-tdg2-pref17 fgc-tdg3-frz1 fgc-tdg3-pref18 fgc-tdg4-frz1 fgc-tdg4-pref19 fgc-tdg5-frz1 fgc-tdg5-pref20 fgc-tdg5-imp1 fgc-tdg6-frz1 fgc-tdg6-pref21 fgc-tdg6-imp2 fgc-pro14-frz1 fgc-hlt12-mon12 run-fgc-pro14-calibration fgc-cal11-pref22 verify-fgc-rsp2-prelaunch verify-fgc-rsp2 verify-fgc-pro13-prelaunch verify-fgc-pro13-result verify-fgc-tdg1-prelaunch verify-fgc-tdg1-result verify-fgc-tdg2-prelaunch verify-fgc-tdg2-result verify-fgc-tdg3-prelaunch verify-fgc-tdg3-result verify-fgc-tdg4-prelaunch verify-fgc-tdg4-result verify-fgc-tdg5-prelaunch verify-fgc-tdg5-result verify-fgc-tdg5-runtime verify-fgc-tdg6-prelaunch verify-fgc-tdg6-result verify-fgc-tdg6-runtime verify-fgc-pro14-prelaunch verify-fgc-pro14-result verify-fgc verify-fgc-sf1 check data paper verify verify-portable-clean-clone clean

.PHONY: fgc-tdg7-frz1 verify-fgc-tdg7-freeze fgc-tdg7-pref23 verify-fgc-tdg7-pref23 fgc-tdg7-imp3 verify-fgc-tdg7-runtime fgc-pro15-frz1 verify-fgc-pro15-freeze fgc-hlt13-mon13 verify-fgc-hlt13-runtime fgc-pro16-frz1 verify-fgc-pro16-freeze fgc-pro16-pref24 verify-fgc-pro16-pref24 fgc-pro17-frz1 verify-fgc-pro17-freeze fgc-pro17-pref25 verify-fgc-pro17-pref25 fgc-hlt14-mon14 verify-fgc-hlt14-runtime fgc-pro18-frz1 verify-fgc-pro18-freeze fgc-pro18-pref26 verify-fgc-pro18-pref26 fgc-pro18-auth1 verify-fgc-pro18-auth1 fgc-pro18-pref27 verify-fgc-pro18-pref27 fgc-pro19-frz1 verify-fgc-pro19-freeze fgc-hlt16-mon16 verify-fgc-pro19-prelaunch fgc-pro19-sid1-frz1 verify-fgc-pro19-sid1-prelaunch preflight-fgc-pro19-sid1-recovery recover-fgc-pro19-sid1 fgc-pro19-sid2-pref1 verify-fgc-pro19-sid2-postrecovery fgc-pro19-sid3-auth1 verify-fgc-pro19-sid3-prelaunch fgc-pro19-pref28 verify-fgc-pro19-pref28-postattempt fgc-tdg8-run1-auth1 verify-fgc-tdg8-successor-prelaunch status-fgc-tdg8-successor-event1 run-fgc-tdg8-successor-event1 fgc-tdg8-rcv1-auth1 verify-fgc-tdg8-rcv1-prelaunch status-fgc-tdg8-rcv1 recover-fgc-tdg8-rcv1 resume-fgc-tdg8-rcv1 fgc-tdg8-rcv2-auth1 verify-fgc-tdg8-rcv2-prelaunch status-fgc-tdg8-rcv2 run-fgc-tdg8-rcv2 preflight-fgc-pro19-sid3-resume resume-fgc-pro19-sid3-event1 status-fgc-pro19-event1 run-fgc-pro19-event1 verify-fgc-sf1-foundation

.PHONY: fgc-tdg8-rcv3-frz1 verify-fgc-tdg8-rcv3-prelaunch fgc-tdg8-rcv3-pref1 verify-fgc-tdg8-rcv3-pref1 fgc-tdg8-rcv3-auth1 verify-fgc-tdg8-rcv3-event-prelaunch status-fgc-tdg8-rcv3 run-fgc-tdg8-rcv3 fgc-tdg8-rcv3-rec1-auth1 verify-fgc-tdg8-rcv3-rec1-prelaunch status-fgc-tdg8-rcv3-rec1 run-fgc-tdg8-rcv3-rec1 fgc-tdg8-rcv3-pref2 verify-fgc-tdg8-rcv3-pref2

PYTHON ?= /opt/homebrew/bin/python3.14
EVOLUTION_PYTHON ?= /opt/homebrew/Caskroom/miniconda/base/bin/python3
PAPER_PYTHON ?= /Users/admin/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3
override FGC_PRO19_SID3_IMPLEMENTATION_COMMIT := 4d7133c0737524712220f77ee5ae4d9c3550c75a
override FGC_PRO19_SID3_PROCESS_MARKER := FGC-PRO19-SID3-EVENT1
override FGC_PRO19_SID3_SANITIZED_PATH := /opt/homebrew/bin:/usr/bin:/bin:/usr/sbin:/sbin
override FGC_PRO19_SID3_SANITIZED_LOCALE := C.UTF-8
override FGC_TDG8_SUCCESSOR_SANITIZED_PATH := /opt/homebrew/bin:/usr/bin:/bin:/usr/sbin:/sbin
override FGC_TDG8_SUCCESSOR_SANITIZED_LOCALE := C.UTF-8

test: test-analytic test-evolution-restart

test-analytic:
	$(PYTHON) scripts/run_fgc_test_partition.py --partition analytic

test-evolution-restart:
	$(EVOLUTION_PYTHON) scripts/run_fgc_test_partition.py --partition evolution-restart

# Portable CI checks the clean-clone code and compact certificates without
# claiming byte reproduction of either local macOS arithmetic history.
test-portable-clean-clone:
	$(PYTHON) scripts/run_fgc_test_partition.py --partition analytic --portable-clean-clone
	$(PYTHON) scripts/run_fgc_test_partition.py --partition evolution-restart --portable-clean-clone

reproduce:
	$(PYTHON) scripts/reproduce_core.py --output results/core-identities.json
	$(PYTHON) scripts/reproduce_controlled_model.py --output results/controlled-model.json
	$(PYTHON) scripts/reproduce_bao.py --output results/desi-dr2-bao-profile.json
	$(PYTHON) scripts/reproduce_observations.py --output results/gw170817-timing.json
	$(PYTHON) scripts/reproduce_gmf.py --output results/gmf-1-spherical-matching.json
	$(PYTHON) scripts/reproduce_einstein_cartan.py --output results/einstein-cartan-collapse.json
	$(PYTHON) scripts/reproduce_gmf1b_preflight.py --output results/gmf-1b-preflight.json
	$(PYTHON) scripts/reproduce_gmf1b_ecd_symmetry.py --output results/gmf-1b-ecd-symmetry.json
	$(PYTHON) scripts/reproduce_gmf1b_ecd_identity.py --output results/gmf-1b-ecd-identity.json
	$(PYTHON) scripts/reproduce_gmf1b_ecd_torsion.py --output results/gmf-1b-ecd-torsion.json
	$(PYTHON) scripts/reproduce_gmf1b_ecd_kinetic.py --output results/gmf-1b-ecd-kinetic.json
	$(PYTHON) scripts/reproduce_fgc_action.py --output results/fgc-1-action-gate.json
	$(PYTHON) scripts/reproduce_fgc_metric_variation.py --output results/fgc-1-metric-variation.json
	$(PYTHON) scripts/reproduce_fgc_hyp1_reduction.py --output results/fgc-1-hyp1-reduction.json
	$(PYTHON) scripts/reproduce_fgc_hyp1_symbol.py --output results/fgc-1-hyp1-symbol.json
	$(PYTHON) scripts/reproduce_fgc_hyp1_modified_harmonic.py --output results/fgc-1-hyp1-modified-harmonic.json
	$(PYTHON) scripts/reproduce_fgc_hyp1_mhg_reference.py --output results/fgc-1-hyp1-mhg-reference.json
	$(PYTHON) scripts/reproduce_fgc_hyp1_mhg_implicit.py --output results/fgc-1-hyp1-mhg-implicit.json
	$(PYTHON) scripts/reproduce_fgc_hyp1_mhg_propagation.py --output results/fgc-1-hyp1-mhg-propagation.json
	$(PYTHON) scripts/reproduce_fgc_hyp1_fo1_rc1.py --output results/fgc-1-hyp1-fo1-rc1.json
	$(PYTHON) scripts/reproduce_fgc_hyp1_dom1_qift1.py --output results/fgc-1-hyp1-dom1-qift1.json
	$(PYTHON) scripts/reproduce_fgc_hyp1_con1_comp1.py --output results/fgc-1-hyp1-con1-comp1.json
	$(PYTHON) scripts/reproduce_fgc_eft0_led1.py --output results/fgc-1-eft0-led1.json
	$(PYTHON) scripts/reproduce_fgc_hyp1_dom2_mode1.py --output results/fgc-1-hyp1-dom2-mode1.json
	$(PYTHON) scripts/reproduce_fgc_hyp1_dom3_uhyp1.py --output results/fgc-1-hyp1-dom3-uhyp1.json
	$(PYTHON) scripts/reproduce_fgc_hyp1_con2_mprop1.py --output results/fgc-1-hyp1-con2-mprop1.json
	$(PYTHON) scripts/reproduce_fgc_hyp1_con3_cau1.py --output results/fgc-1-hyp1-con3-cau1.json
	$(PYTHON) scripts/reproduce_fgc_hyp1_bnd1_md1.py --output results/fgc-1-hyp1-bnd1-md1.json
	$(PYTHON) scripts/reproduce_fgc_def0_obs1.py --output results/fgc-1-def0-obs1.json
	$(PYTHON) scripts/reproduce_fgc_eft1_open1.py --output results/fgc-1-eft1-open1.json
	$(PYTHON) scripts/reproduce_fgc_id0_pref1.py --output results/fgc-1-id0-pref1.json
	$(PYTHON) scripts/reproduce_fgc_src1_nl1.py --output results/fgc-1-src1-nl1.json
	$(PYTHON) scripts/reproduce_fgc_con4_phy1.py --output results/fgc-1-con4-phy1.json
	$(PYTHON) scripts/reproduce_fgc_ctr1_reg1.py --output results/fgc-1-ctr1-reg1.json
	$(PYTHON) scripts/reproduce_fgc_id1_fam1.py --output results/fgc-1-id1-fam1.json
	$(PYTHON) scripts/reproduce_fgc_dom4_run1.py --output results/fgc-1-dom4-run1.json
	$(PYTHON) scripts/reproduce_fgc_hyp2_md1.py --output results/fgc-1-hyp2-md1.json
	$(PYTHON) scripts/reproduce_fgc_bnd2_cp1.py --output results/fgc-1-bnd2-cp1.json
	$(PYTHON) scripts/reproduce_fgc_hlt1_mon1.py --output results/fgc-1-hlt1-mon1.json
	$(PYTHON) scripts/reproduce_fgc_num1_val1.py --output results/fgc-1-num1-val1.json
	$(PYTHON) scripts/reproduce_fgc_run1_sym1.py --output results/fgc-1-run1-sym1.json
	$(PYTHON) scripts/reproduce_fgc_cal0_pref2.py --output results/fgc-1-cal0-pref2.json
	$(PYTHON) scripts/reproduce_fgc_pro4_frz1.py --output results/fgc-1-pro4-frz1.json
	$(PYTHON) scripts/reproduce_fgc_hlt2_mon2.py --output results/fgc-1-hlt2-mon2.json
	$(PYTHON) scripts/reproduce_fgc_id2_all1.py --output results/fgc-1-id2-all1.json
	$(PYTHON) scripts/reproduce_fgc_cal1_pref3.py --output results/fgc-1-cal1-pref3.json
	$(PYTHON) scripts/reproduce_fgc_pro5_frz1.py --output results/fgc-1-pro5-frz1.json
	$(PYTHON) scripts/reproduce_fgc_hlt3_mon3.py --check --output results/fgc-1-hlt3-mon3.json
	$(PYTHON) scripts/reproduce_fgc_cal2_pref4.py --check --output results/fgc-1-cal2-pref4.json
	$(PYTHON) scripts/reproduce_fgc_pro6_frz1.py --output results/fgc-1-pro6-frz1.json
	$(PYTHON) scripts/reproduce_fgc_hlt4_mon4.py --check --output results/fgc-1-hlt4-mon4.json
	$(PYTHON) scripts/reproduce_fgc_cal3_pref5.py --check --output results/fgc-1-cal3-pref5.json
	$(PYTHON) scripts/reproduce_fgc_pro7_frz1.py --output results/fgc-1-pro7-frz1.json
	$(PYTHON) scripts/reproduce_fgc_hlt5_mon5.py --check --output results/fgc-1-hlt5-mon5.json
	$(PYTHON) scripts/reproduce_fgc_cal4_pref6.py --check --output results/fgc-1-cal4-pref6.json
	$(PYTHON) scripts/reproduce_fgc_pro8_frz1.py --output results/fgc-1-pro8-frz1.json
	$(PYTHON) scripts/reproduce_fgc_hlt6_mon6.py --check --output results/fgc-1-hlt6-mon6.json
	$(PYTHON) scripts/reproduce_nested_spectrum.py --output results/nested-gradient-spectrum.json

fgc-action:
	$(PYTHON) scripts/reproduce_fgc_action.py --output results/fgc-1-action-gate.json

fgc-metric-variation: fgc-action
	$(PYTHON) scripts/reproduce_fgc_metric_variation.py --output results/fgc-1-metric-variation.json

fgc-hyp1-reduction: fgc-metric-variation
	$(PYTHON) scripts/reproduce_fgc_hyp1_reduction.py --output results/fgc-1-hyp1-reduction.json

fgc-hyp1-symbol: fgc-hyp1-reduction
	$(PYTHON) scripts/reproduce_fgc_hyp1_symbol.py --output results/fgc-1-hyp1-symbol.json

fgc-hyp1-modified-harmonic: fgc-hyp1-symbol
	$(PYTHON) scripts/reproduce_fgc_hyp1_modified_harmonic.py --output results/fgc-1-hyp1-modified-harmonic.json

fgc-hyp1-mhg-reference: fgc-hyp1-modified-harmonic
	$(PYTHON) scripts/reproduce_fgc_hyp1_mhg_reference.py --output results/fgc-1-hyp1-mhg-reference.json

fgc-hyp1-mhg-implicit: fgc-hyp1-mhg-reference
	$(PYTHON) scripts/reproduce_fgc_hyp1_mhg_implicit.py --output results/fgc-1-hyp1-mhg-implicit.json

fgc-hyp1-mhg-propagation: fgc-hyp1-mhg-implicit
	$(PYTHON) scripts/reproduce_fgc_hyp1_mhg_propagation.py --output results/fgc-1-hyp1-mhg-propagation.json

fgc-hyp1-fo1-rc1: fgc-hyp1-mhg-propagation
	$(PYTHON) scripts/reproduce_fgc_hyp1_fo1_rc1.py --output results/fgc-1-hyp1-fo1-rc1.json

fgc-hyp1-dom1-qift1: fgc-hyp1-fo1-rc1
	$(PYTHON) scripts/reproduce_fgc_hyp1_dom1_qift1.py --output results/fgc-1-hyp1-dom1-qift1.json

fgc-hyp1-con1-comp1: fgc-hyp1-dom1-qift1
	$(PYTHON) scripts/reproduce_fgc_hyp1_con1_comp1.py --output results/fgc-1-hyp1-con1-comp1.json

fgc-eft0-led1: fgc-hyp1-con1-comp1
	$(PYTHON) scripts/reproduce_fgc_eft0_led1.py --output results/fgc-1-eft0-led1.json

fgc-hyp1-dom2-mode1: fgc-hyp1-con1-comp1
	$(PYTHON) scripts/reproduce_fgc_hyp1_dom2_mode1.py --output results/fgc-1-hyp1-dom2-mode1.json

# This exact interval certificate is intentionally the expensive FGC target.
fgc-hyp1-dom3-uhyp1: fgc-hyp1-dom2-mode1
	$(PYTHON) scripts/reproduce_fgc_hyp1_dom3_uhyp1.py --output results/fgc-1-hyp1-dom3-uhyp1.json

fgc-hyp1-con2-mprop1: fgc-hyp1-con1-comp1
	$(PYTHON) scripts/reproduce_fgc_hyp1_con2_mprop1.py --output results/fgc-1-hyp1-con2-mprop1.json

fgc-hyp1-con3-cau1: fgc-hyp1-con2-mprop1
	$(PYTHON) scripts/reproduce_fgc_hyp1_con3_cau1.py --output results/fgc-1-hyp1-con3-cau1.json

fgc-hyp1-bnd1-md1: fgc-hyp1-dom3-uhyp1 fgc-hyp1-con2-mprop1
	$(PYTHON) scripts/reproduce_fgc_hyp1_bnd1_md1.py --output results/fgc-1-hyp1-bnd1-md1.json

fgc-def0-obs1: fgc-hyp1-reduction
	$(PYTHON) scripts/reproduce_fgc_def0_obs1.py --output results/fgc-1-def0-obs1.json

fgc-eft1-open1: fgc-eft0-led1 fgc-hyp1-bnd1-md1 fgc-hyp1-con3-cau1 fgc-def0-obs1
	$(PYTHON) scripts/reproduce_fgc_eft1_open1.py --output results/fgc-1-eft1-open1.json

# This focused floating certificate consumes the existing exact predecessor
# records without rerunning DOM3; the outer repository audit checks their
# implementation-hash freshness.
fgc-src1-nl1:
	$(PYTHON) scripts/reproduce_fgc_src1_nl1.py --output results/fgc-1-src1-nl1.json

# Like SRC1, CON4 consumes the sealed exact predecessor records without
# rerunning their expensive proofs during focused constraint development.
fgc-con4-phy1:
	$(PYTHON) scripts/reproduce_fgc_con4_phy1.py --output results/fgc-1-con4-phy1.json

# CTR1 rederives the affected reference descendants in exact centre series;
# it consumes but does not rerun the sealed annular predecessor certificates.
fgc-ctr1-reg1: fgc-src1-nl1 fgc-con4-phy1
	$(PYTHON) scripts/reproduce_fgc_ctr1_reg1.py --output results/fgc-1-ctr1-reg1.json

# ID0 audits immutable PROTO1 history and therefore precedes, rather than
# depends on, active PROTO2 RUN1 composition.
fgc-id0-pref1:
	$(PYTHON) scripts/reproduce_fgc_id0_pref1.py --output results/fgc-1-id0-pref1.json

# ID1 consumes the repaired protocol and the exact physical-constraint and
# regular-centre certificates. It constructs initial hypersurfaces only.
fgc-id1-fam1: fgc-id0-pref1 fgc-con4-phy1 fgc-ctr1-reg1
	$(PYTHON) scripts/reproduce_fgc_id1_fam1.py --output results/fgc-1-id1-fam1.json

# DOM4 bridges the constraint-solved ID1 slice to complete REF1 second jets
# and certifies only a nonzero classical state domain, not a trajectory tube.
fgc-dom4-run1: fgc-src1-nl1 fgc-con4-phy1 fgc-ctr1-reg1 fgc-id1-fam1
	$(PYTHON) scripts/reproduce_fgc_dom4_run1.py --output results/fgc-1-dom4-run1.json

# HYP2 applies the full covariant 24-characteristic all-covector early-kill
# test to the exact DOM4 domain definition. It remains classical, not EFT.
fgc-hyp2-md1: fgc-dom4-run1
	$(PYTHON) scripts/reproduce_fgc_hyp2_md1.py --output results/fgc-1-hyp2-md1.json

# BND2 protects only the retained measurement region by an evolving all-cone
# causal ledger. It does not claim a nonlinear constraint-preserving IBVP.
fgc-bnd2-cp1: fgc-con4-phy1 fgc-id1-fam1 fgc-dom4-run1 fgc-hyp2-md1
	$(PYTHON) scripts/reproduce_fgc_bnd2_cp1.py --output results/fgc-1-bnd2-cp1.json

# HLT1 makes every classical domain, constraint, scale, alias, and boundary
# stop transactional before stage acceptance. It is not an evolution result.
fgc-hlt1-mon1: fgc-src1-nl1 fgc-con4-phy1 fgc-dom4-run1 fgc-hyp2-md1 fgc-bnd2-cp1
	$(PYTHON) scripts/reproduce_fgc_hlt1_mon1.py --output results/fgc-1-hlt1-mon1.json

# NUM1 validates both numerical methods, the unchanged-equation batched source,
# atomic HLT1/BND2 composition, restart, centre limits, and affine measurement
# on outcome-blind controls. It opens no FGC-QR holdout result.
fgc-num1-val1: fgc-src1-nl1 fgc-con4-phy1 fgc-ctr1-reg1 fgc-id1-fam1 fgc-dom4-run1 fgc-hyp2-md1 fgc-bnd2-cp1 fgc-hlt1-mon1 fgc-def0-obs1
	$(PYTHON) scripts/reproduce_fgc_num1_val1.py --output results/fgc-1-num1-val1.json

# This target composes existing hash-bound evidence and intentionally does not
# rerun the expensive DOM3 interval proof during focused RUN1 development.
fgc-run1-sym1: fgc-id0-pref1 fgc-src1-nl1 fgc-con4-phy1 fgc-ctr1-reg1 fgc-id1-fam1 fgc-dom4-run1 fgc-hyp2-md1 fgc-bnd2-cp1 fgc-hlt1-mon1 fgc-num1-val1
	$(PYTHON) scripts/reproduce_fgc_run1_sym1.py --output results/fgc-1-run1-sym1.json

# CAL0 evaluates PROTO3's declared GR-0 input with HLT1's own estimator and
# audits branch-specific stop ownership. It certifies a pre-holdout contract
# obstruction, not a candidate trajectory or mechanism result.
fgc-cal0-pref2: fgc-run1-sym1
	$(PYTHON) scripts/reproduce_fgc_cal0_pref2.py --output results/fgc-1-cal0-pref2.json

# PRO4-FRZ1 binds the premise-only PROTO4 overlay to immutable PROTO3/CAL0
# Git blobs. It freezes no candidate outcome and opens no holdout execution.
fgc-pro4-frz1: fgc-cal0-pref2
	$(PYTHON) scripts/reproduce_fgc_pro4_frz1.py --output results/fgc-1-pro4-frz1.json

# HLT2 implements PROTO4's branch-owned and convergence-aware successor
# admission on static/synthetic controls. It opens no calibration trajectory.
fgc-hlt2-mon2: fgc-pro4-frz1
	$(PYTHON) scripts/reproduce_fgc_hlt2_mon2.py --output results/fgc-1-hlt2-mon2.json

# ID2 evaluates every possible calibration amplitude against all five future
# FGC-QR input modifiers before reading a trajectory. It freezes only the
# statically eligible GR-0 candidate set, not a collapse outcome.
fgc-id2-all1: fgc-hlt2-mon2
	$(PYTHON) scripts/reproduce_fgc_id2_all1.py --output results/fgc-1-id2-all1.json

# CAL1-PREF3 composes the ID2 continuum data with each native semidiscrete
# method at t=0. It diagnoses PROTO4's contract and reads no trajectory.
fgc-cal1-pref3: fgc-id2-all1
	$(PYTHON) scripts/reproduce_fgc_cal1_pref3.py --output results/fgc-1-cal1-pref3.json

# PRO5-FRZ1 freezes CAL1's premise-only repair against immutable commit
# 0781d7d. It does not implement the runtime compositor or authorize a run.
fgc-pro5-frz1: fgc-cal1-pref3
	$(PYTHON) scripts/reproduce_fgc_pro5_frz1.py --output results/fgc-1-pro5-frz1.json

# HLT3 composes the real branch-specific GR-0 transaction, freezes all twelve
# fresh inputs, and authorizes only CAL2. The canonical pre-launch observation
# remains verifiable after launch through the hash-bound campaign manifest.
fgc-hlt3-mon3: fgc-pro5-frz1
	$(PYTHON) scripts/reproduce_fgc_hlt3_mon3.py --check --output results/fgc-1-hlt3-mon3.json

# CAL2-PREF4 consumes the immutable completed PROTO5 GR-0 campaign and proves
# that its threshold-adjacent stop belongs to an unaccepted trial step. It is
# a protocol diagnosis, not a calibration or mechanism result.
fgc-cal2-pref4: fgc-hlt3-mon3
	$(PYTHON) scripts/reproduce_fgc_cal2_pref4.py --check --output results/fgc-1-cal2-pref4.json

# PRO6-FRZ1 freezes CAL2's source-stop ownership repair against immutable
# checkpoint 47f4c83. It does not implement HLT4 or authorize a trajectory.
fgc-pro6-frz1: fgc-cal2-pref4
	$(PYTHON) scripts/reproduce_fgc_pro6_frz1.py --output results/fgc-1-pro6-frz1.json

# HLT4 implements only PROTO6's accepted-versus-internal-trial source-stop
# ownership, freezes the twelve unchanged inputs, and authorizes one fresh
# GR-0 recalibration. It does not open either candidate branch.
fgc-hlt4-mon4: fgc-pro6-frz1
	$(PYTHON) scripts/reproduce_fgc_hlt4_mon4.py --check --output results/fgc-1-hlt4-mon4.json

# CAL3-PREF5 binds the completed PROTO6 campaign and deterministically replays
# the finest terminal member once when generated. The normal target performs
# fast immutable/hash/semantic validation; the explicit full target repeats
# the roughly four-minute trajectory reconstruction.
fgc-cal3-pref5: fgc-hlt4-mon4
	$(PYTHON) scripts/reproduce_fgc_cal3_pref5.py --check --output results/fgc-1-cal3-pref5.json

fgc-cal3-pref5-full-replay: fgc-hlt4-mon4
	$(PYTHON) scripts/reproduce_fgc_cal3_pref5.py --check --replay --output results/fgc-1-cal3-pref5.json

# PRO7-FRZ1 freezes the general accepted-state versus wholly-unaccepted
# proposal transaction rule against immutable CAL3 evidence. It authorizes no
# runtime or trajectory.
fgc-pro7-frz1: fgc-cal3-pref5
	$(PYTHON) scripts/reproduce_fgc_pro7_frz1.py --output results/fgc-1-pro7-frz1.json

# HLT5 implements only PROTO7's general unaccepted-proposal ownership,
# freezes the twelve unchanged inputs, and authorizes one fresh GR-0 campaign.
# It does not open either candidate branch.
fgc-hlt5-mon5: fgc-pro7-frz1
	$(PYTHON) scripts/reproduce_fgc_hlt5_mon5.py --check --output results/fgc-1-hlt5-mon5.json

# CAL4-PREF6 binds the completed PROTO7 campaign and its deterministic
# amplitude-5/2 replay. The normal target validates immutable raw/cache hashes;
# the explicit full target repeats the roughly six-minute reconstruction.
fgc-cal4-pref6: fgc-hlt5-mon5
	$(PYTHON) scripts/reproduce_fgc_cal4_pref6.py --check --output results/fgc-1-cal4-pref6.json

fgc-cal4-pref6-full-replay: fgc-hlt5-mon5
	$(PYTHON) scripts/reproduce_fgc_cal4_pref6.py --check --replay --output results/fgc-1-cal4-pref6.json

# PRO8-FRZ1 freezes CAL4's common-event diagnostic ownership correction
# against immutable PROTO7/CAL4 evidence. It authorizes no runtime or
# trajectory and changes no numerical threshold.
fgc-pro8-frz1: fgc-cal4-pref6
	$(PYTHON) scripts/reproduce_fgc_pro8_frz1.py --output results/fgc-1-pro8-frz1.json

# HLT6 implements only PROTO8's common-event diagnostic ownership,
# reconstructs the unchanged twelve inputs, and authorizes one fresh GR-0
# campaign. It reads no PROTO8 trajectory and opens no candidate branch.
fgc-hlt6-mon6: fgc-pro8-frz1
	$(PYTHON) scripts/reproduce_fgc_hlt6_mon6.py --check --output results/fgc-1-hlt6-mon6.json

# CAL5-PREF7 binds the completed PROTO8 campaign and reduces its two early
# GR-0 stops to a resolution-ladder adequacy question. It changes no threshold
# and authorizes no successor runtime or candidate branch.
fgc-cal5-pref7: fgc-hlt6-mon6
	$(PYTHON) scripts/reproduce_fgc_cal5_pref7.py --check --output results/fgc-1-cal5-pref7.json

# PRO9-FRZ1 freezes CAL5's one-level resolution-ladder successor against the
# immutable PROTO8 campaign. It implements no runtime and creates no namespace.
fgc-pro9-frz1: fgc-cal5-pref7
	$(PYTHON) scripts/reproduce_fgc_pro9_frz1.py --output results/fgc-1-pro9-frz1.json

# HLT7 constructs the new ladder and binds the immutable runner, but the
# unchanged nested-tail rule vetoes all four t=0 cases. No campaign is opened.
fgc-hlt7-mon7: fgc-pro9-frz1
	$(PYTHON) scripts/reproduce_fgc_hlt7_mon7.py --output results/fgc-1-hlt7-mon7.json

# CAL6-PREF8 preserves HLT7's raw veto while testing whether that ratio is
# conditioned below the proper-grid map's measured non-idempotence. It opens
# no trajectory and freezes no successor protocol.
fgc-cal6-pref8: fgc-hlt7-mon7
	$(PYTHON) scripts/reproduce_fgc_cal6_pref8.py --output results/fgc-1-cal6-pref8.json

# PRO10-FRZ1 froze only CAL6's conditioning-guarded interpretation of a failed
# finest-pair ratio.  Its generator is intentionally pre-launch-only because
# the empty-namespace observation cannot be recreated after the campaign.  At
# this post-campaign boundary, prove the canonical freeze through HLT8's
# immutable-commit lineage without mutating the hash-bound generator.
fgc-pro10-frz1: fgc-cal6-pref8
	$(PYTHON) -m unittest tests/test_fgc_pro10_frz1_reproduction.py -v

# HLT8 implements the frozen PROTO10 compositor, reconstructs all twelve
# unchanged physical inputs, attacks every guarded-saturation premise, and
# records its authorization of only one fresh GR-0 calibration. Its canonical
# check consumes
# the stored pre-launch namespace evidence and remains valid post-campaign.
fgc-hlt8-mon8:
	$(PYTHON) scripts/reproduce_fgc_hlt8_mon8.py --check --output results/fgc-1-hlt8-mon8.json

# CAL7 binds the completed PROTO10 campaign and localizes its repeated
# RK4-8193 stop to a centre-adjacent binary64 epsilon/r^2 source floor.  Its
# well-balanced reference-state calculation is a t=0 control, not a runtime.
fgc-cal7-pref9: fgc-hlt8-mon8
	$(PYTHON) scripts/reproduce_fgc_cal7_pref9.py --check --output results/fgc-1-cal7-pref9.json

# PRO11-FRZ1's pre-launch absence observation is temporal. HLT9 now proves its
# exact canonical blob from immutable commit b53c0f1, so the same target remains
# valid after an authorized calibration has populated the proto11 namespace.
fgc-pro11-frz1: fgc-cal7-pref9
	$(PYTHON) -m unittest tests/test_fgc_pro11_frz1_reproduction.py -v

# HLT9 implements and attacks only PROTO11's reference-state numerical map,
# rebuilds all twelve inputs, and authorizes one fresh GR-0 calibration. It
# advances no trajectory and opens no candidate or physical claim.
fgc-hlt9-mon9: fgc-pro11-frz1
	$(PYTHON) scripts/reproduce_fgc_hlt9_mon9.py --check --output results/fgc-1-hlt9-mon9.json

# CAL8 binds the completed PROTO11 campaign, proves both amplitudes reached one
# evolved common event, and separates an affine source-arithmetic wall from an
# independent direct coarse derivative-spectrum veto. It opens no candidate.
fgc-cal8-pref10: fgc-hlt9-mon9
	$(PYTHON) scripts/reproduce_fgc_cal8_pref10.py --check --output results/fgc-1-cal8-pref10.json

# SRC2 freezes and attacks only the replay/capture instrument. The cheap target
# validates lineage, arithmetic controls, and nonclaims without launching the
# expensive evolved replay or creating its ignored raw namespace.
fgc-src2-cap1: fgc-cal8-pref10
	$(PYTHON) scripts/run_fgc_src2_affine_capture.py --check
	$(PYTHON) -m unittest tests/test_fgc_src2_affine_arithmetic.py tests/test_fgc_src2_capture_runner.py -v

# PREF11 reduces the completed ignored CAP1 bundle to one compact, tracked
# point fixture and exact-rational cancellation diagnosis. It never reruns the
# expensive replay and does not freeze a successor runtime.
fgc-src2-pref11: fgc-src2-cap1
	$(PYTHON) scripts/reproduce_fgc_src2_pref11.py --check --output results/fgc-1-src2-pref11.json
	$(PYTHON) -m unittest tests/test_fgc_src2_exact_oracle.py tests/test_fgc_src2_pref11_reproduction.py -v

# SRC3 changes only the binary64 evaluation graph of the unchanged GR-0 REF1
# rows. The optional CAP1 raw fixture is recomputed when locally present; a
# clean clone still reproduces the compact and independent exact controls.
fgc-src3: fgc-src2-pref11
	$(PYTHON) scripts/reproduce_fgc_src3.py --check --output results/fgc-1-src3.json
	$(PYTHON) -m unittest tests/test_fgc_src3_reference_balanced_source.py tests/test_fgc_src3_reproduction.py -v

# RSP1 prospectively freezes one amplitude-three, two-method,
# three-resolution study. Its authorization publishes the failed t=0 target
# ratios but reads no evolved endpoint and opens no candidate branch.
fgc-rsp1-frz1: fgc-src3
	$(PYTHON) scripts/reproduce_fgc_rsp1_frz1.py --check --output results/fgc-1-rsp1-frz1.json
	$(PYTHON) -m unittest tests/test_fgc_rsp1_resolution_runtime.py tests/test_fgc_rsp1_frz1_reproduction.py tests/test_fgc_rsp1_runner.py -v

# PREF12 binds the completed ignored RSP1 namespace, restores all six terminal
# states from the checkpoint, and recomputes the endpoint before clearing only
# the frozen amplitude-three target veto for successor design.
fgc-rsp1-pref12: fgc-rsp1-frz1
	$(PYTHON) scripts/reproduce_fgc_rsp1_pref12.py --check --output results/fgc-1-rsp1-pref12.json
	$(PYTHON) -m unittest tests/test_fgc_rsp1_campaign_diagnosis.py tests/test_fgc_rsp1_pref12_reproduction.py -v

# SRC4 preserves SRC3's unchanged equations and strict source gate while
# contracting only the four-dimensional tensor indices through NumPy. Its
# equivalence certificate opens no trajectory and does not freeze PROTO12.
fgc-src4-vec1: fgc-rsp1-pref12
	$(PYTHON) scripts/reproduce_fgc_src4_vec1.py --check --output results/fgc-1-src4-vec1.json
	$(PYTHON) -m unittest tests/test_fgc_src4_vectorized_reference_source.py tests/test_fgc_src4_vec1_reproduction.py -v

# PRO12's absent-namespace observation is temporal pre-launch evidence. HLT10
# now proves the exact canonical freeze from immutable commit a427801 and
# hash-matches the current artifact after the authorized campaign populated
# the proto12 namespace, as well as testing every non-temporal contract.
fgc-pro12-frz1: fgc-src4-vec1
	$(PYTHON) -m unittest tests/test_fgc_protocol_v12.py tests/test_fgc_spectral_sensitivity_v12.py tests/test_fgc_pro12_frz1_reproduction.py -v

# HLT10 installs PROTO12's batched SRC4 evaluator and pairwise
# resolved-or-diagnostically-saturated admission in the live GR-0 runtime. It
# authorizes only a fresh PROTO12 calibration; every candidate and physical
# claim remains closed.
fgc-hlt10-mon10: fgc-pro12-frz1
	$(PYTHON) scripts/reproduce_fgc_hlt10_mon10.py --check --output results/fgc-1-hlt10-mon10.json
	$(PYTHON) -m unittest tests/test_fgc_proto12_runtime.py tests/test_fgc_hlt10_mon10_reproduction.py tests/test_fgc_gr0_campaign_runner_v12.py -v

# CAL9/PREF13 binds the completed ignored PROTO12 campaign, restores its
# terminal checkpoint, and localizes the sole late comparator veto. It changes
# no threshold and opens no candidate branch.
fgc-cal9-pref13: fgc-hlt10-mon10
	$(PYTHON) scripts/reproduce_fgc_cal9_pref13.py --check --output results/fgc-1-cal9-pref13.json
	$(PYTHON) -m unittest tests/test_fgc_cal9_proto12_campaign_diagnosis.py tests/test_fgc_cal9_pref13_reproduction.py -v

# RSP2 freezes one and only one genuinely finer SSPRK3 member.  This target
# reconstructs its two CAL9 predecessors and authorizes the run, but never
# creates the raw namespace or advances the trajectory.
fgc-rsp2-frz1: fgc-cal9-pref13
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_rsp2_frz1.py --check --output results/fgc-1-rsp2-frz1.json
	$(EVOLUTION_PYTHON) -m unittest tests/test_fgc_rsp2_constraint_runtime.py tests/test_fgc_rsp2_frz1_reproduction.py tests/test_fgc_rsp2_runner.py -v

verify-fgc-rsp2-prelaunch: fgc-rsp2-frz1
	$(PYTHON) scripts/check_repo.py

# PREF14 consumes the immutable prelaunch commit and post-run hashes. It must
# not rerun FRZ1's historical namespace-absence observation or the mutable
# floating replay chain after execution; its reproducer restores and checks
# the required CAL9/RSP2 checkpoints directly.
fgc-rsp2-pref14:
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_rsp2_pref14.py --check --output results/fgc-1-rsp2-pref14.json
	$(EVOLUTION_PYTHON) -m unittest tests/test_fgc_rsp2_pref14_diagnosis.py tests/test_fgc_rsp2_pref14_reproduction.py -v

verify-fgc-rsp2: fgc-rsp2-pref14
	$(PYTHON) scripts/check_repo.py

# PRO13 prospectively binds the six exact t=23/16 restart payloads and the
# method-owned resolution ladders. It creates neither successor namespace and
# advances no state. Once HLT11 is launched, later verification must consume
# the immutable PRO13 commit rather than replay this temporal absence check.
fgc-pro13-frz1:
	$(PYTHON) scripts/reproduce_fgc_cal10_pref15.py --check-pro13-freeze
	$(PYTHON) -m unittest tests/test_fgc_protocol_v13.py tests/test_fgc_pro13_frz1_reproduction.py -v

# HLT11 restores every complete restart payload and binds the method-owned
# runtime without advancing beyond t=23/16. Its namespace evidence is temporal;
# the canonical record retains the historical observation after launch.
fgc-hlt11-mon11: fgc-pro13-frz1
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_hlt11_mon11.py --check --output results/fgc-1-hlt11-mon11.json
	$(EVOLUTION_PYTHON) -m unittest tests/test_fgc_proto13_runtime.py tests/test_fgc_hlt11_mon11_reproduction.py tests/test_fgc_gr0_campaign_runner_v13.py -v

verify-fgc-pro13-prelaunch: fgc-hlt11-mon11
	$(PYTHON) scripts/check_repo.py

# PREF15 binds the completed PROTO13 terminal checkpoint and independently
# recomputes its final spatial, constraint, trapped, and temporal admissions.
# It never resumes the terminal checkpoint or promotes the temporal stop to a
# candidate-action or physical result. A clean clone may omit the ignored raw
# bundle; a partial bundle fails closed.
fgc-cal10-pref15:
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_cal10_pref15.py --check --output results/fgc-1-cal10-pref15.json
	$(EVOLUTION_PYTHON) -m unittest tests/test_fgc_cal10_proto13_campaign_diagnosis.py tests/test_fgc_cal10_pref15_reproduction.py -v

verify-fgc-pro13-result: fgc-cal10-pref15
	$(PYTHON) scripts/check_repo.py

# TDG1 freezes a no-trajectory diagnosis after a synthetic theorem-control
# proves that normalized raw temporal tail fractions cannot generally measure
# spatial convergence of a nonzero temporal continuum signal. It hashes but
# does not load the terminal checkpoint arrays.
fgc-tdg1-frz1: fgc-cal10-pref15
	$(PYTHON) scripts/reproduce_fgc_tdg1_frz1.py --verify --output results/fgc-1-tdg1-frz1.json
	$(PYTHON) -m unittest tests/test_fgc_tdg1_temporal_diagnosis.py tests/test_fgc_tdg1_frz1_reproduction.py -v

verify-fgc-tdg1-prelaunch: fgc-tdg1-frz1
	$(PYTHON) scripts/check_repo.py

# PREF16 consumes the exact Accelerate-owned terminal histories read-only and
# binds TDG1's outcome. It advances no state and defines no replacement gate.
fgc-tdg1-pref16: fgc-tdg1-frz1
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_tdg1_pref16.py --check --output results/fgc-1-tdg1-pref16.json
	$(EVOLUTION_PYTHON) -m unittest tests/test_fgc_tdg1_pref16_reproduction.py -v

verify-fgc-tdg1-result: fgc-tdg1-pref16
	$(PYTHON) scripts/check_repo.py

# TDG2 freezes an absolute, unnormalized tail-power discriminator after TDG1
# has bound the normalized-ratio defect. It verifies the source checkpoint hash
# but may not load a real history array until this freeze is committed.
fgc-tdg2-frz1: fgc-tdg1-pref16
	$(PYTHON) scripts/reproduce_fgc_tdg2_frz1.py --verify --output results/fgc-1-tdg2-frz1.json
	$(PYTHON) -m unittest tests/test_fgc_tdg2_absolute_tail_diagnosis.py tests/test_fgc_tdg2_frz1_reproduction.py -v

verify-fgc-tdg2-prelaunch: fgc-tdg2-frz1
	$(PYTHON) scripts/check_repo.py

# PREF17 consumes the immutable Accelerate-owned terminal histories read-only,
# binds all 576 signal ladders, and preserves the mixed outcome. It defines no
# replacement temporal admission and authorizes no trajectory.
fgc-tdg2-pref17: fgc-tdg2-frz1
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_tdg2_pref17.py --check --output results/fgc-1-tdg2-pref17.json
	$(EVOLUTION_PYTHON) -m unittest tests/test_fgc_tdg2_pref17_reproduction.py -v

verify-fgc-tdg2-result: fgc-tdg2-pref17
	$(PYTHON) scripts/check_repo.py

# TDG3 prospectively freezes two native-grid, no-resampling estimators after
# PREF17 localized the unresolved majority to the interpolation enclosure. It
# may inspect the compact PREF17 result but loads no terminal history array.
fgc-tdg3-frz1: fgc-tdg2-pref17
	$(PYTHON) scripts/reproduce_fgc_tdg3_frz1.py --verify --output results/fgc-1-tdg3-frz1.json
	$(PYTHON) -m unittest tests/test_fgc_tdg3_native_tail_diagnosis.py tests/test_fgc_tdg3_frz1_reproduction.py -v

verify-fgc-tdg3-prelaunch: fgc-tdg3-frz1
	$(PYTHON) scripts/check_repo.py

# PREF18 consumes the immutable Accelerate-owned histories read-only, applies
# both frozen native estimators, and binds the mixed result. It advances no
# state and defines no replacement temporal admission.
fgc-tdg3-pref18: fgc-tdg3-frz1
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_tdg3_pref18.py --check --output results/fgc-1-tdg3-pref18.json
	$(EVOLUTION_PYTHON) -m unittest tests/test_fgc_tdg3_pref18_reproduction.py -v

verify-fgc-tdg3-result: fgc-tdg3-pref18
	$(PYTHON) scripts/check_repo.py

# TDG4 freezes only the finite-sampling identifiability theorem and exact
# no-history controls. It reads PREF18's compact record, never the raw history,
# and authorizes no replacement gate or trajectory.
fgc-tdg4-frz1: fgc-tdg3-pref18
	$(PYTHON) scripts/reproduce_fgc_tdg4_frz1.py --verify --output results/fgc-1-tdg4-frz1.json
	$(PYTHON) -m unittest tests/test_fgc_tdg4_sampling_identifiability.py tests/test_fgc_tdg4_frz1_reproduction.py -v

verify-fgc-tdg4-prelaunch: fgc-tdg4-frz1
	$(PYTHON) scripts/check_repo.py

# PREF19 independently executes the frozen theorem without loading a terminal
# history. It retires the present sampled rule only as continuum admission and
# authorizes replacement-gate design, not acceptance or a trajectory.
fgc-tdg4-pref19: fgc-tdg4-frz1
	$(PYTHON) scripts/reproduce_fgc_tdg4_pref19.py --verify --output results/fgc-1-tdg4-pref19.json
	$(PYTHON) -m unittest tests/test_fgc_tdg4_sampling_identifiability_theorem.py tests/test_fgc_tdg4_pref19_reproduction.py -v

verify-fgc-tdg4-result: fgc-tdg4-pref19
	$(PYTHON) scripts/check_repo.py

# TDG5 prospectively freezes one sufficient replacement architecture selected
# from TDG4. It reads no raw history, defines no acceptance threshold, and
# authorizes neither runtime implementation nor a new trajectory.
fgc-tdg5-frz1: fgc-tdg4-pref19
	$(PYTHON) scripts/reproduce_fgc_tdg5_frz1.py --verify --output results/fgc-1-tdg5-frz1.json
	$(PYTHON) -m unittest tests/test_fgc_tdg5_stage_complete_refinement.py tests/test_fgc_tdg5_frz1_reproduction.py -v

verify-fgc-tdg5-prelaunch: fgc-tdg5-frz1
	$(PYTHON) scripts/check_repo.py

# PREF20 independently executes TDG5's frozen finite-dimensional theorem,
# checks the immutable engine's endpoint-RHS records, and authorizes only
# implementation and synthetic validation of the paired runtime.
fgc-tdg5-pref20: fgc-tdg5-frz1
	$(PYTHON) scripts/reproduce_fgc_tdg5_pref20.py --verify --output results/fgc-1-tdg5-pref20.json
	$(PYTHON) -m unittest tests/test_fgc_tdg5_stage_complete_refinement_theorem.py tests/test_fgc_tdg5_pref20_reproduction.py -v

verify-fgc-tdg5-result: fgc-tdg5-pref20
	$(PYTHON) scripts/check_repo.py

# IMP1 implements the authorized two-phase runtime pair and continuous
# binary64 envelope on synthetic controls only. It freezes no thresholds and
# cannot open PROTO14, GR-0 recalibration, or a candidate branch.
fgc-tdg5-imp1: fgc-tdg5-pref20
	$(PYTHON) scripts/reproduce_fgc_tdg5_imp1.py --verify --output results/fgc-1-tdg5-imp1.json
	$(PYTHON) -m unittest tests/test_fgc_tdg5_stage_complete_refinement_runtime.py tests/test_fgc_tdg5_imp1_reproduction.py -v

verify-fgc-tdg5-runtime: fgc-tdg5-imp1
	$(PYTHON) scripts/check_repo.py

# TDG6 freezes the replacement threshold and production-admission semantics
# without loading CAL10 or implementing the production compositor. Only the
# separately committed TDG6 binder may authorize that implementation.
fgc-tdg6-frz1: fgc-tdg5-imp1
	$(PYTHON) scripts/reproduce_fgc_tdg6_frz1.py --verify --output results/fgc-1-tdg6-frz1.json
	$(PYTHON) -m unittest tests/test_fgc_tdg6_temporal_admission_design.py tests/test_fgc_tdg6_frz1_reproduction.py -v

verify-fgc-tdg6-prelaunch: fgc-tdg6-frz1
	$(PYTHON) scripts/check_repo.py

fgc-tdg6-pref21: fgc-tdg6-frz1
	$(PYTHON) scripts/reproduce_fgc_tdg6_pref21.py --verify --output results/fgc-1-tdg6-pref21.json
	$(PYTHON) -m unittest tests/test_fgc_tdg6_temporal_admission_theorem.py tests/test_fgc_tdg6_pref21_reproduction.py -v

verify-fgc-tdg6-result: fgc-tdg6-pref21
	$(PYTHON) scripts/check_repo.py

# IMP2 implements PREF21's authorized production-shaped compositor and proves
# it only on synthetic controls. It authorizes a prospective PROTO14 freeze,
# but it cannot advance a production state or open any candidate branch.
fgc-tdg6-imp2: fgc-tdg6-pref21
	$(PYTHON) scripts/reproduce_fgc_tdg6_imp2.py --verify --output results/fgc-1-tdg6-imp2.json
	$(PYTHON) -m unittest tests/test_fgc_tdg6_temporal_admission_runtime.py tests/test_fgc_tdg6_imp2_reproduction.py -v

verify-fgc-tdg6-runtime: fgc-tdg6-imp2
	$(PYTHON) scripts/check_repo.py

# PRO14 and HLT12 are prospective, prelaunch records.  Their live namespace
# absence assertions are not rerun after launch; CAL11 owns postlaunch proof.
# The historical targets remain available only for the prelaunch boundary.
fgc-pro14-frz1: fgc-tdg6-imp2
	$(PYTHON) -m unittest tests/test_fgc_protocol_v14.py -v
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_pro14_frz1.py --check
	$(EVOLUTION_PYTHON) -m unittest tests/test_fgc_pro14_frz1_reproduction.py -v

# HLT12 wires the already-qualified TDG6 compositor into the production
# PROTO14 restart path, attacks the integration, and advances no production
# state.  Its namespace-absence observation is sealed before the command below
# may launch the one authorized fresh GR-0 continuation.
fgc-hlt12-mon12: fgc-pro14-frz1
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_hlt12_mon12.py --check --output results/fgc-1-hlt12-mon12.json
	$(EVOLUTION_PYTHON) -m unittest tests/test_fgc_proto14_runtime.py tests/test_fgc_hlt12_mon12_reproduction.py tests/test_fgc_gr0_campaign_runner_v14.py -v

verify-fgc-pro14-prelaunch: fgc-hlt12-mon12
	$(PYTHON) scripts/check_repo.py

# CAL11/PREF22 consumes the already-terminal PROTO14 bundle read-only.  It
# binds the checkpoint-owned invalid-runtime result and must never run or
# resume the historical campaign.
fgc-cal11-pref22:
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_cal11_pref22.py --check --output results/fgc-1-cal11-pref22.json
	$(EVOLUTION_PYTHON) -m unittest tests/test_fgc_cal11_pref22_reproduction.py -v

# TDG7 freezes a stage-safe exact-binary64 shared-lattice diagnosis of CAL11's
# pre-shadow binary64 failure.  Its prospective G=8Q selection makes actual
# RK4/SSPRK3 c={0,1/2,1} stage times and endpoints exact; it does not repair,
# mutate, resume, or launch PROTO14.
fgc-tdg7-frz1: fgc-cal11-pref22
	$(PYTHON) scripts/reproduce_fgc_tdg7_frz1.py --verify --output results/fgc-1-tdg7-frz1.json
	$(PYTHON) -m unittest tests/test_fgc_tdg7_binary64_subdivision_lattice.py tests/test_fgc_tdg7_frz1_reproduction.py -v

verify-fgc-tdg7-freeze: fgc-tdg7-frz1
	$(PYTHON) scripts/check_repo.py

# PREF23 independently binds the stage-safe lattice and authorizes only a
# later runtime-repair implementation plus synthetic qualification.  It never
# resumes PROTO14 or creates a calibration/candidate outcome.
fgc-tdg7-pref23: fgc-tdg7-frz1
	$(PYTHON) scripts/reproduce_fgc_tdg7_pref23.py --verify --output results/fgc-1-tdg7-pref23.json
	$(PYTHON) -m unittest tests/test_fgc_tdg7_stage_safe_lattice_theorem.py tests/test_fgc_tdg7_pref23_reproduction.py -v

verify-fgc-tdg7-pref23: fgc-tdg7-pref23
	$(PYTHON) scripts/check_repo.py

# IMP3 implements PREF23's authorized stage-safe runtime repair and qualifies
# it solely with deterministic synthetic controls. It does not load CAL11 or
# PROTO14 history, create a namespace, or advance a calibration trajectory.
fgc-tdg7-imp3: fgc-tdg7-pref23
	$(PYTHON) scripts/reproduce_fgc_tdg7_imp3.py --verify --output results/fgc-1-tdg7-imp3.json
	$(PYTHON) -m unittest tests/test_fgc_tdg7_stage_safe_runtime.py tests/test_fgc_tdg7_imp3_reproduction.py -v

verify-fgc-tdg7-runtime: fgc-tdg7-imp3
	$(PYTHON) scripts/check_repo.py

# PRO15 prospectively freezes the campaign-owned FRESH_READY/RETRY_PENDING
# cursor and crash-recovery contract required after IMP3.  It implements no
# cursor runtime, creates no namespace, and authorizes no calibration.
fgc-pro15-frz1: fgc-tdg7-imp3
	$(PYTHON) scripts/reproduce_fgc_pro15_frz1.py --verify --output results/fgc-1-pro15-frz1.json
	$(PYTHON) -m unittest tests/test_fgc_protocol_v15.py tests/test_fgc_pro15_frz1_reproduction.py -v

verify-fgc-pro15-freeze: fgc-pro15-frz1
	$(PYTHON) scripts/check_repo.py

# HLT13 implements and synthetically qualifies PROTO15's durable campaign
# cursor. It reads no historical campaign data and creates neither output root.
fgc-hlt13-mon13:
	$(PYTHON) scripts/reproduce_fgc_hlt13_mon13.py --verify --output results/fgc-1-hlt13-mon13.json
	$(PYTHON) -m unittest tests/test_fgc_proto15_runtime.py tests/test_fgc_hlt13_mon13_reproduction.py -v

verify-fgc-hlt13-runtime: fgc-hlt13-mon13
	$(PYTHON) scripts/check_repo.py --only-hlt13-mon13

# PRO16 freezes the missing successful common-event transition and the external
# trusted-genesis contract. It implements no runner and creates no namespace.
fgc-pro16-frz1: fgc-hlt13-mon13
	$(PYTHON) scripts/reproduce_fgc_pro16_frz1.py --verify --output results/fgc-1-pro16-frz1.json
	$(PYTHON) -m unittest tests/test_fgc_protocol_v16.py tests/test_fgc_pro16_frz1_reproduction.py -v

verify-fgc-pro16-freeze: fgc-pro16-frz1
	$(PYTHON) scripts/check_repo.py --only-pro16-frz1

# PREF24 independently proves the receipt graph acyclic, but diagnoses the
# successor/genesis as under-specified. It authorizes only a PROTO17 freeze.
fgc-pro16-pref24: fgc-pro16-frz1
	$(PYTHON) scripts/reproduce_fgc_pro16_pref24.py --verify --output results/fgc-1-pro16-pref24.json
	$(PYTHON) -m unittest tests/test_fgc_proto16_transition_theorem.py tests/test_fgc_pro16_pref24_reproduction.py -v

verify-fgc-pro16-pref24: fgc-pro16-pref24
	$(PYTHON) scripts/check_repo.py --only-pro16-pref24

# PROTO17 freezes the exact pure transition/genesis construction identified as
# absent by PREF24.  It remains side-effect-free and creates no namespace.
fgc-pro17-frz1: fgc-pro16-pref24
	$(PYTHON) scripts/reproduce_fgc_pro17_frz1.py --verify --output results/fgc-1-pro17-frz1.json
	$(PYTHON) -m unittest tests/test_fgc_protocol_v17.py tests/test_fgc_proto17_pure_construction.py tests/test_fgc_pro17_frz1_reproduction.py -v

verify-fgc-pro17-freeze: fgc-pro17-frz1
	$(PYTHON) scripts/check_repo.py --only-pro17-frz1

# PREF25 independently reconstructs the sealed pure mapping without importing
# either PROTO17 implementation module. It authorizes only HLT14 implementation
# plus synthetic qualification; no namespace or trajectory is opened.
fgc-pro17-pref25: fgc-pro17-frz1
	$(PYTHON) scripts/reproduce_fgc_pro17_pref25.py --verify --output results/fgc-1-pro17-pref25.json
	$(PYTHON) -m unittest tests/test_fgc_proto17_construction_theorem.py tests/test_fgc_pro17_pref25_reproduction.py -v

verify-fgc-pro17-pref25: fgc-pro17-pref25
	$(PYTHON) scripts/check_repo.py --only-pro17-pref25

# HLT14 implements and synthetically qualifies PROTO17's authenticated
# six-bundle genesis and one-common-event durable runtime. It opens no real
# campaign input, production namespace, trajectory, calibration, or candidate.
fgc-hlt14-mon14: fgc-pro17-pref25
	$(PYTHON) scripts/reproduce_fgc_hlt14_mon14.py --verify --output results/fgc-1-hlt14-mon14.json
	$(PYTHON) -m unittest tests/test_fgc_proto17_hlt14_inputs.py tests/test_fgc_proto17_hlt14_runtime.py tests/test_fgc_hlt14_mon14_reproduction.py -v

verify-fgc-hlt14-runtime: fgc-hlt14-mon14
	$(PYTHON) scripts/check_repo.py --only-hlt14-mon14

# PRO18 freezes only the compact production-prelaunch overlay.  Ordinary
# verification never opens campaign archives, re-observes temporal namespace
# absence, creates a namespace, or runs a trajectory.
fgc-pro18-frz1: fgc-hlt14-mon14
	$(PYTHON) scripts/reproduce_fgc_pro18_frz1.py --verify
	$(PYTHON) -m unittest tests/test_fgc_pro18_frz1_reproduction.py -v

verify-fgc-pro18-freeze: fgc-pro18-frz1
	$(PYTHON) scripts/check_repo.py --only-pro18-frz1

# PREF26 is intentionally raw-dependent and read-only: it opens only the two
# declared historical source containers and journals, never a future PROTO17
# output root. It is excluded from portable/default aggregate verification.
fgc-pro18-pref26: fgc-pro18-frz1
	$(PYTHON) scripts/reproduce_fgc_pro18_pref26.py --verify
	$(PYTHON) -m unittest tests/test_fgc_proto18_production_inputs.py tests/test_fgc_proto18_historical_replay.py tests/test_fgc_proto18_pref26_binder.py tests/test_fgc_pro18_pref26_reproduction.py -v

verify-fgc-pro18-pref26: fgc-pro18-pref26
	$(PYTHON) scripts/check_repo.py --only-pro18-pref26

# AUTH1 consumes only PREF26's already-committed compact result. It constructs
# one calibration-only GenesisSpec and pins the prospective launch source and
# evolution environment, but performs no Git launch authentication, raw
# reimport, future-root observation, namespace creation, or state advance.
fgc-pro18-auth1: fgc-pro18-frz1
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_pro18_auth1.py --verify
	$(EVOLUTION_PYTHON) -m unittest tests/test_fgc_proto18_auth1_inputs.py tests/test_fgc_proto18_authority.py tests/test_fgc_pro18_auth1_reproduction.py tests/test_fgc_proto17_hlt15_runtime.py tests/test_fgc_proto17_hlt15_runner.py -v

verify-fgc-pro18-auth1: fgc-pro18-auth1
	$(PYTHON) scripts/check_repo.py --only-pro18-auth1

# PREF27's explicit historical target binds the real HLT15 generation-zero
# store and runs only bounded temporary no-adoption/staging-residue controls.
# It never advances the calibration state. Ordinary aggregate verification
# consumes only the compact tracked certificate and does not replay these
# one-time temporal or refusal observations after authorized progression.
fgc-pro18-pref27: fgc-pro18-auth1
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_pro18_pref27.py --verify
	$(EVOLUTION_PYTHON) -m unittest tests/test_fgc_proto18_pref27_binder.py tests/test_check_repo_pro18_pref27.py tests/test_fgc_pro18_pref27_reproduction.py -v

# Explicit one-time historical live audit. It is intentionally not the
# postlaunch-safe compact PREF27 path used by ordinary aggregate verification.
verify-fgc-pro18-pref27:
	$(PYTHON) scripts/check_repo.py --only-pro18-pref27

# PRO19 freezes only the first GR-0 edge from sealed generation zero.  It has
# no HLT16 runtime, raw-array import, output namespace, or state-advance API.
fgc-pro19-frz1:
	$(PYTHON) scripts/reproduce_fgc_pro19_frz1.py --verify
	$(PYTHON) -m unittest tests/test_fgc_proto19_progression_freeze.py -v

verify-fgc-pro19-freeze: fgc-pro19-frz1
	$(PYTHON) scripts/check_repo.py --only-pro19-frz1

# MON16 qualifies only the durable state/payload/journal/checkpoint machinery,
# recovery semantics, status surface, and bounded first-event entrypoint.  It
# does not open or advance the authenticated calibration store.
fgc-hlt16-mon16: fgc-pro19-frz1
	$(EVOLUTION_PYTHON) scripts/check_repo.py --only-hlt16-mon16
	PYTHONPATH=src $(EVOLUTION_PYTHON) -m unittest tests/test_hlt16_lifecycle.py tests/test_hlt16_progression_attempt.py tests/test_fgc_hlt16_member_codec.py tests/test_fgc_hlt16_state_store.py tests/test_fgc_hlt16_campaign_schema.py tests/test_fgc_hlt16_campaign_store.py tests/test_fgc_hlt16_campaign_recovery.py tests/test_fgc_hlt16_campaign_runtime.py tests/test_fgc_proto19_progression_contract.py tests/test_fgc_proto19_progression_inputs.py tests/test_fgc_proto19_launch_authority.py tests/test_fgc_proto19_launch_manifest.py tests/test_fgc_pro19_event1_status.py tests/test_fgc_pro19_event1_runner.py tests/test_fgc_hlt16_mon16_artifact.py tests/test_check_repo_hlt16_prelaunch.py -v

# Compact prelaunch verification only.  No trajectory command is reachable
# from this target or from any ordinary verifier.
verify-fgc-pro19-prelaunch: fgc-hlt16-mon16
	$(PYTHON) scripts/check_repo.py --only-pro19-prelaunch

# CFL1 preserves the failed first attempt as an implementation-contract
# diagnosis and authorizes only the exact generation-six continuation.  This
# target is read-only: it does not take over the stale lease or advance state.
fgc-pro19-cfl1-frz1:
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_pro19_cfl1_frz1.py --verify
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_hlt16_mon16.py --verify
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_pro19_launch_manifest.py --verify
	PYTHONPATH=src $(EVOLUTION_PYTHON) -m unittest \
		tests.test_fgc_proto19_cfl_continuation_authority \
		tests.test_fgc_hlt16_campaign_runtime.HLT16CampaignRuntimeDispatchTests \
		tests.test_fgc_pro19_event1_runner.Proto19Event1RunnerTests.test_continuation_preflight_authenticates_exact_checkpoint_without_writer \
		tests.test_fgc_pro19_event1_runner.Proto19Event1RunnerTests.test_corrected_resume_uses_persisted_shells_and_original_plan_only \
		tests.test_fgc_pro19_event1_runner.Proto19Event1RunnerTests.test_continuation_refuses_store_movement_after_anchor_capture \
		tests.test_fgc_pro19_event1_runner.Proto19Event1RunnerTests.test_corrected_resume_never_reopens_generation_zero -v

verify-fgc-pro19-cfl1-prelaunch: fgc-pro19-cfl1-frz1
	$(PYTHON) scripts/check_repo.py --only-pro19-cfl1-prelaunch

# SID1 authenticates exactly the pre-recovery generation-seven suffix and the
# finite recovery it permits.  Its live-anchor observation is intentionally a
# one-time pre-recovery action: ordinary verification only checks the compact
# tracked certificate so it remains valid after the recovery checkpoint exists.
fgc-pro19-sid1-frz1: fgc-hlt16-mon16
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_pro19_sid1_frz1.py --verify-compact
	PYTHONPATH=src $(EVOLUTION_PYTHON) -m unittest tests/test_check_repo_pro19_sid1_prelaunch.py -v

verify-fgc-pro19-sid1-prelaunch: fgc-pro19-sid1-frz1
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_pro19_sid1_frz1.py --verify
	PYTHONPATH=src $(EVOLUTION_PYTHON) -m unittest tests/test_fgc_proto19_sid1_authority.py tests/test_fgc_pro19_event1_runner.py -v
	$(PYTHON) scripts/check_repo.py --only-pro19-sid1-prelaunch

# Explicit provenance-only preflight. The historical recovery mutation is
# tombstoned in closed-live-targets.mk.
preflight-fgc-pro19-sid1-recovery:
	$(EVOLUTION_PYTHON) scripts/run_fgc_pro19_event1.py --sid1-preflight

# SID2 independently binds the already-published metadata-only generation-8
# recovery. The ordinary target is compact and remains valid after lawful
# progression; the explicit post-recovery target owns the one-time live audit.
fgc-pro19-sid2-pref1: fgc-pro19-sid1-frz1
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_pro19_sid2_pref1.py --verify-compact
	PYTHONPATH=src $(EVOLUTION_PYTHON) -m unittest tests/test_check_repo_pro19_sid2_postrecovery.py -v

verify-fgc-pro19-sid2-postrecovery: fgc-pro19-sid2-pref1
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_pro19_sid2_pref1.py --verify
	PYTHONPATH=src $(EVOLUTION_PYTHON) -m unittest tests/test_fgc_proto19_sid2_binder.py tests/test_check_repo_pro19_sid2_postrecovery.py -v
	$(EVOLUTION_PYTHON) scripts/check_repo.py --only-pro19-sid2-postrecovery

# SID3 compactly binds the exact A-owned source/import image and immutable
# generation-eight inputs. This target and its focused repository verifier are
# store-blind: the live campaign is reachable only through the explicit
# committed-stdin targets below.
fgc-pro19-sid3-auth1: fgc-pro19-sid2-pref1
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_pro19_sid3_auth1.py --verify
	PYTHONPATH=src $(EVOLUTION_PYTHON) -m unittest \
		tests/test_fgc_proto19_execution_closure.py \
		tests/test_fgc_proto19_gr0_static_factory.py \
		tests/test_fgc_proto19_resume_authority.py \
		tests/test_fgc_pro19_sid3_bootstrap.py \
		tests/test_fgc_pro19_sid3_runner.py \
		tests/test_check_repo_pro19_sid3_prelaunch.py -v

verify-fgc-pro19-sid3-prelaunch: fgc-pro19-sid3-auth1
	$(EVOLUTION_PYTHON) scripts/check_repo.py --only-pro19-sid3-prelaunch

# PREF28 compactly binds the generation-nine invalid terminal and unchanged
# accepted member map. Ordinary verification is store-blind; only the explicit
# post-attempt target below re-authenticates the live 56-leaf terminal store.
fgc-pro19-pref28:
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_pro19_pref28.py --verify-compact
	PYTHONPATH=src $(EVOLUTION_PYTHON) -m unittest tests/test_fgc_proto19_pref28_binder.py tests/test_check_repo_pro19_pref28.py -v

verify-fgc-pro19-pref28-postattempt: fgc-pro19-pref28
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_pro19_pref28.py --verify
	$(EVOLUTION_PYTHON) scripts/check_repo.py --only-pro19-pref28-postattempt

# TDG8 freezes exactly one read-only replay-evidence diagnosis and one
# conditional representation-only repair. It is store-blind and authorizes no
# successor execution.
fgc-tdg8-frz1:
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_tdg8_frz1.py --verify
	PYTHONPATH=src $(EVOLUTION_PYTHON) -m unittest tests/test_fgc_tdg8_replay_freeze.py -v

fgc-tdg8-diag1: fgc-tdg8-frz1
	$(EVOLUTION_PYTHON) scripts/diagnose_fgc_tdg8_replay.py --verify-compact
	PYTHONPATH=src $(EVOLUTION_PYTHON) -m unittest \
		tests/test_hlt16_lifecycle.py \
		tests/test_hlt16_progression_attempt.py -v

# Store-blind compact authority for exactly one fresh GR-0 successor event.
# It never re-observes destination absence and never opens either campaign.
fgc-tdg8-run1-auth1:
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_tdg8_run1_auth1.py --verify-compact

verify-fgc-tdg8-successor-prelaunch: fgc-tdg8-run1-auth1
	PYTHONPATH=src $(EVOLUTION_PYTHON) -m unittest \
		tests/test_fgc_tdg8_successor_authority.py \
		tests/test_fgc_tdg8_successor_runtime.py \
		tests/test_fgc_tdg8_successor_runner.py -v
	$(EVOLUTION_PYTHON) scripts/check_repo.py --only-tdg8-successor-prelaunch

# Read-only status for the fresh successor namespace. HEAD is resolved at
# invocation time; the runner independently authenticates that exact clean
# committed image and the compact typed authority.
status-fgc-tdg8-successor-event1:
	@set -eu; \
	authority_commit="$$(env -i PATH="$(FGC_TDG8_SUCCESSOR_SANITIZED_PATH)" LC_ALL="$(FGC_TDG8_SUCCESSOR_SANITIZED_LOCALE)" \
		git --no-replace-objects --no-optional-locks \
		-c core.fsmonitor=false -c core.untrackedCache=false \
		-C "$(CURDIR)" rev-parse --verify 'HEAD^{commit}')"; \
	env -i PATH="$(FGC_TDG8_SUCCESSOR_SANITIZED_PATH)" LC_ALL="$(FGC_TDG8_SUCCESSOR_SANITIZED_LOCALE)" \
		"$(EVOLUTION_PYTHON)" -I -B scripts/run_fgc_tdg8_successor_event.py \
		--authority-commit "$$authority_commit" --status

# Compact, store-blind authority for the exact failed TDG8 retry edge.  The
# one-time --write action is performed before its immutable prelaunch commit;
# every ordinary verifier below consumes only the frozen compact result.
fgc-tdg8-rcv1-auth1:
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_tdg8_rcv1_auth1.py --verify-compact

verify-fgc-tdg8-rcv1-prelaunch: fgc-tdg8-rcv1-auth1
	PYTHONPATH=src $(EVOLUTION_PYTHON) -m unittest \
		tests/test_fgc_hlt16_campaign_runtime.py \
		tests/test_fgc_tdg8_retry_recovery_authority.py \
		tests/test_fgc_tdg8_retry_recovery_runner.py -v
	$(EVOLUTION_PYTHON) scripts/check_repo.py --only-tdg8-rcv1-prelaunch

# All three commands authenticate the same committed image. Status is
# read-only; recover publishes only the unique metadata successor; resume
# continues the original GR-0 event from that exact recovered ancestor.
status-fgc-tdg8-rcv1:
	@set -eu; \
	authority_commit="$$(env -i PATH="$(FGC_TDG8_SUCCESSOR_SANITIZED_PATH)" LC_ALL="$(FGC_TDG8_SUCCESSOR_SANITIZED_LOCALE)" \
		git --no-replace-objects --no-optional-locks \
		-c core.fsmonitor=false -c core.untrackedCache=false \
		-C "$(CURDIR)" rev-parse --verify 'HEAD^{commit}')"; \
	env -i PATH="$(FGC_TDG8_SUCCESSOR_SANITIZED_PATH)" LC_ALL="$(FGC_TDG8_SUCCESSOR_SANITIZED_LOCALE)" \
		"$(EVOLUTION_PYTHON)" -I -B scripts/run_fgc_tdg8_retry_recovery.py \
		--authority-commit "$$authority_commit" --status

# One prospective authority for the bounded TDG8 temporal-retry relation.
# Its compact verifier is store-blind; only the one-time result writer and the
# explicit status/run commands inspect the authenticated campaign.
fgc-tdg8-rcv2-auth1:
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_tdg8_rcv2_auth1.py --verify-compact

verify-fgc-tdg8-rcv2-prelaunch: fgc-tdg8-rcv2-auth1
	PYTHONPATH=src $(EVOLUTION_PYTHON) -m unittest \
		tests/test_fgc_hlt16_campaign_recovery.py \
		tests/test_fgc_tdg8_bounded_retry_authority.py \
		tests/test_fgc_tdg8_bounded_retry_runner.py -v
	$(EVOLUTION_PYTHON) scripts/check_repo.py --only-tdg8-rcv2-prelaunch

status-fgc-tdg8-rcv2:
	@set -eu; \
	authority_commit="$$(env -i PATH="$(FGC_TDG8_SUCCESSOR_SANITIZED_PATH)" LC_ALL="$(FGC_TDG8_SUCCESSOR_SANITIZED_LOCALE)" \
		git --no-replace-objects --no-optional-locks \
		-c core.fsmonitor=false -c core.untrackedCache=false \
		-C "$(CURDIR)" rev-parse --verify 'HEAD^{commit}')"; \
	env -i PATH="$(FGC_TDG8_SUCCESSOR_SANITIZED_PATH)" LC_ALL="$(FGC_TDG8_SUCCESSOR_SANITIZED_LOCALE)" \
		"$(EVOLUTION_PYTHON)" -I -B scripts/run_fgc_tdg8_bounded_retry.py \
		--authority-commit "$$authority_commit" --status

# Store-blind RCV3 projection freeze.  This consumes only the tracked compact
# certificate; it never repeats the one-time destination-absence observation
# and never invokes the separately authorized bootstrap installer.
fgc-tdg8-rcv3-frz1:
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_tdg8_rcv3_frz1.py --verify-compact

verify-fgc-tdg8-rcv3-prelaunch: fgc-tdg8-rcv3-frz1
	PYTHONPATH=src $(EVOLUTION_PYTHON) -m unittest \
		tests/test_fgc_tdg8_rcv3_freeze.py \
		tests/test_fgc_tdg8_rcv3_fork_runtime.py \
		tests/test_check_repo_tdg8_rcv3_prelaunch.py -v
	$(EVOLUTION_PYTHON) scripts/check_repo.py --only-tdg8-rcv3-prelaunch

# Post-install compact evidence. Ordinary verification never opens either run
# store. The explicit verification target below performs the independent,
# read-only live binding and its focused adversarial suite.
fgc-tdg8-rcv3-pref1:
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_tdg8_rcv3_pref1.py --verify-compact

verify-fgc-tdg8-rcv3-pref1: fgc-tdg8-rcv3-pref1
	PYTHONPATH=src $(EVOLUTION_PYTHON) -m unittest \
		tests/test_fgc_tdg8_rcv3_pref1_binder.py \
		tests/test_check_repo_tdg8_rcv3_pref1.py -v
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_tdg8_rcv3_pref1.py --verify
	$(EVOLUTION_PYTHON) scripts/check_repo.py --only-tdg8-rcv3-pref1

# Prospective same-event authority. Ordinary verification consumes compact
# bytes only and never opens or advances the installed projection.
fgc-tdg8-rcv3-auth1:
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_tdg8_rcv3_auth1.py --verify-compact

verify-fgc-tdg8-rcv3-event-prelaunch: fgc-tdg8-rcv3-auth1
	PYTHONPATH=src $(EVOLUTION_PYTHON) -m unittest \
		tests/test_fgc_tdg8_rcv3_execution_authority.py \
		tests/test_fgc_tdg8_rcv3_event_runner.py \
		tests/test_check_repo_tdg8_rcv3_auth1.py -v
	$(EVOLUTION_PYTHON) scripts/check_repo.py --only-tdg8-rcv3-auth1

status-fgc-tdg8-rcv3:
	@set -eu; \
	authority_commit="$$(env -i PATH="$(FGC_TDG8_SUCCESSOR_SANITIZED_PATH)" LC_ALL="$(FGC_TDG8_SUCCESSOR_SANITIZED_LOCALE)" \
		git --no-replace-objects --no-optional-locks \
		-c core.fsmonitor=false -c core.untrackedCache=false \
		-C "$(CURDIR)" rev-parse --verify 'HEAD^{commit}')"; \
	env -i PATH="$(FGC_TDG8_SUCCESSOR_SANITIZED_PATH)" LC_ALL="$(FGC_TDG8_SUCCESSOR_SANITIZED_LOCALE)" \
		"$(EVOLUTION_PYTHON)" -I -B scripts/run_fgc_tdg8_rcv3_event.py \
		--authority-commit "$$authority_commit" --status

# Prospective REC1 authority for the already-durable generation-nine plus
# sequence-11/12 recovery suffix. Ordinary verification is store-blind: the
# explicit runner alone may take over the stale lease, publish the exact frozen
# generation-ten checkpoint, verify it before a proposal, and remain in event 23.
fgc-tdg8-rcv3-rec1-auth1:
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_tdg8_rcv3_rec1_auth1.py --verify-compact

verify-fgc-tdg8-rcv3-rec1-prelaunch: fgc-tdg8-rcv3-rec1-auth1
	PYTHONPATH=src $(EVOLUTION_PYTHON) -m unittest \
		tests/test_fgc_tdg8_rcv3_rec1_authority.py \
		tests/test_fgc_tdg8_rcv3_rec1_event_runner.py \
		tests/test_check_repo_tdg8_rcv3_rec1_auth1.py -v
	$(EVOLUTION_PYTHON) scripts/check_repo.py --only-tdg8-rcv3-rec1-auth1

status-fgc-tdg8-rcv3-rec1:
	@set -eu; \
	authority_commit="$$(env -i PATH="$(FGC_TDG8_SUCCESSOR_SANITIZED_PATH)" LC_ALL="$(FGC_TDG8_SUCCESSOR_SANITIZED_LOCALE)" \
		git --no-replace-objects --no-optional-locks \
		-c core.fsmonitor=false -c core.untrackedCache=false \
		-C "$(CURDIR)" rev-parse --verify 'HEAD^{commit}')"; \
	env -i PATH="$(FGC_TDG8_SUCCESSOR_SANITIZED_PATH)" LC_ALL="$(FGC_TDG8_SUCCESSOR_SANITIZED_LOCALE)" \
		"$(EVOLUTION_PYTHON)" -I -B scripts/run_fgc_tdg8_rcv3_rec1_event.py \
		--authority-commit "$$authority_commit" --status

# Outcome-neutral terminal evidence. Ordinary and foundation verification read
# only the tracked compact record; the explicit verification target performs
# the independent read-only store reconstruction and payload decode.
fgc-tdg8-rcv3-pref2:
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_tdg8_rcv3_pref2.py --verify-compact

verify-fgc-tdg8-rcv3-pref2: fgc-tdg8-rcv3-pref2
	PYTHONPATH=src $(EVOLUTION_PYTHON) -m unittest \
		tests/test_fgc_tdg8_rcv3_pref2_binder.py \
		tests/test_check_repo_tdg8_rcv3_pref2.py -v
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_tdg8_rcv3_pref2.py --verify
	$(EVOLUTION_PYTHON) scripts/check_repo.py --only-tdg8-rcv3-pref2

.PHONY: fgc-tdg9-ar1-auth1 verify-fgc-tdg9-ar1-prelaunch status-fgc-tdg9-ar1 run-fgc-tdg9-ar1 fgc-tdg9-ar1-pref1 verify-fgc-tdg9-ar1-pref1

# Prospective, compact authority for one read-only exact-arithmetic diagnostic
# over RCV3 retry widths 3, 4, and 5. It binds the exact environment, selection,
# and 17-path single-successor committed image. Ordinary verification never
# opens the campaign or diagnostic namespace and never replays a proposal.
fgc-tdg9-ar1-auth1:
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_tdg9_ar1_auth1.py --verify-compact

verify-fgc-tdg9-ar1-prelaunch: fgc-tdg9-ar1-auth1
	PYTHONPATH=src $(EVOLUTION_PYTHON) -m unittest \
		tests/test_fgc_tdg9_exact_temporal_arithmetic.py \
		tests/test_fgc_tdg9_ar1_authority.py \
		tests/test_fgc_tdg9_ar1_runner.py \
		tests/test_check_repo_tdg9_ar1_auth1.py -v
	$(EVOLUTION_PYTHON) scripts/check_repo.py --only-tdg9-ar1-auth1

status-fgc-tdg9-ar1:
	@set -eu; \
	authority_commit="$$(env -i PATH="$(FGC_TDG8_SUCCESSOR_SANITIZED_PATH)" LC_ALL="$(FGC_TDG8_SUCCESSOR_SANITIZED_LOCALE)" \
		git --no-replace-objects --no-optional-locks \
		-c core.fsmonitor=false -c core.untrackedCache=false \
		-C "$(CURDIR)" rev-parse --verify 'HEAD^{commit}')"; \
	env -i PATH="$(FGC_TDG8_SUCCESSOR_SANITIZED_PATH)" LC_ALL="$(FGC_TDG8_SUCCESSOR_SANITIZED_LOCALE)" \
		"$(EVOLUTION_PYTHON)" -I -B scripts/run_fgc_tdg9_ar1.py \
		--authority-commit "$$authority_commit" --status

# Outcome-neutral post-run binder. Ordinary verification reads only the
# compact tracked record; --live is an explicit, bounded 108-call replay.
fgc-tdg9-ar1-pref1:
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_tdg9_ar1_pref1.py

verify-fgc-tdg9-ar1-pref1: fgc-tdg9-ar1-pref1
	PYTHONPATH=src $(EVOLUTION_PYTHON) -m unittest \
		tests/test_fgc_tdg9_ar1_pref1_binder.py \
		tests/test_check_repo_tdg9_ar1_pref1.py -v
	$(EVOLUTION_PYTHON) scripts/check_repo.py --only-tdg9-ar1-pref1

.PHONY: fgc-tdg9-loc1-frz1 verify-fgc-tdg9-loc1-prelaunch run-fgc-tdg9-loc1

# Prospective compact authority for exact localization of PREF1's ten failed
# occurrences. Ordinary verification is store-blind and never executes LOC1.
fgc-tdg9-loc1-frz1:
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_tdg9_loc1_frz1.py --verify-compact

verify-fgc-tdg9-loc1-prelaunch: fgc-tdg9-loc1-frz1
	PYTHONPATH=src $(EVOLUTION_PYTHON) -m unittest \
		tests/test_fgc_tdg9_local_extrema.py \
		tests/test_fgc_tdg9_loc1_authority.py \
		tests/test_fgc_tdg9_loc1_runner.py \
		tests/test_check_repo_tdg9_loc1_frz1.py -v
	$(EVOLUTION_PYTHON) scripts/check_repo.py --only-tdg9-loc1-frz1

.PHONY: fgc-tdg9-loc1-err1 fgc-tdg9-loc2-frz1 verify-fgc-tdg9-loc2-prelaunch run-fgc-tdg9-loc2 fgc-tdg9-loc2-pref2 verify-fgc-tdg9-loc2-pref2 fgc-tdg9-ti1-frz1 verify-fgc-tdg9-ti1-compact verify-fgc-tdg9-ti1-prelaunch status-fgc-tdg9-ti1 run-fgc-tdg9-ti1 fgc-tdg9-ti2-frz1 verify-fgc-tdg9-ti2-prelaunch status-fgc-tdg9-ti2 run-fgc-tdg9-ti2 fgc-tdg9-ti2-pref1 verify-fgc-tdg9-ti2-pref1 fgc-tdg9-ac1-frz1 verify-fgc-tdg9-ac1-prelaunch status-fgc-tdg9-ac1 run-fgc-tdg9-ac1 fgc-tdg9-ac1-pref1 verify-fgc-tdg9-ac1-pref1 fgc-tdg9-ur1-frz1 verify-fgc-tdg9-ur1-prelaunch status-fgc-tdg9-ur1 run-fgc-tdg9-ur1 fgc-tdg9-ur1-pref1 verify-fgc-tdg9-ur1-pref1 fgc-tdg10-qa1-frz1 verify-fgc-tdg10-qa1-prelaunch status-fgc-tdg10-qa1 run-fgc-tdg10-qa1 fgc-tdg10-qa1-pref1 verify-fgc-tdg10-qa1-pref1 fgc-tdg10-qa2-frz1 verify-fgc-tdg10-qa2-prelaunch status-fgc-tdg10-qa2 run-fgc-tdg10-qa2 fgc-tdg10-qa2-pref1 verify-fgc-tdg10-qa2-pref1

fgc-tdg9-loc1-err1:
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_tdg9_loc1_err1.py --verify-compact

# Prospective LOC2 authority. Ordinary verification is store-blind and never
# executes the localization diagnostic.
fgc-tdg9-loc2-frz1: fgc-tdg9-loc1-err1
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_tdg9_loc2_frz1.py --verify-compact

verify-fgc-tdg9-loc2-prelaunch: fgc-tdg9-loc2-frz1
	PYTHONPATH=src:. $(EVOLUTION_PYTHON) -m unittest \
		tests/test_fgc_tdg9_loc1_err1.py \
		tests/test_fgc_tdg9_local_extrema_v2.py \
		tests/test_fgc_tdg9_loc2_authority.py \
		tests/test_fgc_tdg9_loc2_runner.py \
		tests/test_check_repo_tdg9_loc2_frz1.py -v
	$(EVOLUTION_PYTHON) scripts/check_repo.py --only-tdg9-loc2-frz1

# Outcome-neutral post-run binder. Ordinary verification reads only the
# tracked compact certificate; --live is the explicit one-time 60-route replay.
fgc-tdg9-loc2-pref2:
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_tdg9_loc2_pref2.py

verify-fgc-tdg9-loc2-pref2: fgc-tdg9-loc2-pref2
	PYTHONPATH=src $(EVOLUTION_PYTHON) -m unittest \
		tests/test_fgc_tdg9_loc2_pref2_binder.py \
		tests/test_check_repo_tdg9_loc2_pref2.py -v
	$(EVOLUTION_PYTHON) scripts/check_repo.py --only-tdg9-loc2-pref2

# Historical TI1 compact record. The authority attempt stopped invalid before
# the first shadow proposal; status/run now fail closed through its retirement
# shim while compact verification preserves the immutable evidence.
fgc-tdg9-ti1-frz1: fgc-tdg9-loc2-pref2
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_tdg9_ti1_frz1.py --verify-compact

verify-fgc-tdg9-ti1-compact: fgc-tdg9-ti1-frz1
	PYTHONPATH=src:. $(EVOLUTION_PYTHON) -m unittest \
		tests/test_fgc_tdg9_ti1_authority.py \
		tests/test_fgc_tdg9_ti1_runner.py \
		tests/test_check_repo_tdg9_ti1_frz1.py -v
	$(EVOLUTION_PYTHON) scripts/check_repo.py --only-tdg9-ti1-frz1

verify-fgc-tdg9-ti1-prelaunch: verify-fgc-tdg9-ti1-compact

# Historical compatibility target. The retirement shim always exits 2.
status-fgc-tdg9-ti1:
	@set -eu; \
	authority_commit="$$(env -i PATH="$$PATH" HOME="$$HOME" LC_ALL=C LANG=C \
		git --no-replace-objects --no-optional-locks -c core.fsmonitor=false \
		-c core.untrackedCache=false -C "$(CURDIR)" rev-parse --verify 'HEAD^{commit}')"; \
	"$(EVOLUTION_PYTHON)" -I -B scripts/run_fgc_tdg9_ti1.py \
		--authority-commit "$$authority_commit" --status

# Prospective fingerprint-only recovery. Ordinary verification consumes the
# tracked compact certificate and never constructs an SSPRK3 shadow proposal.
fgc-tdg9-ti2-frz1: verify-fgc-tdg9-ti1-compact
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_tdg9_ti2_frz1.py --verify-compact

verify-fgc-tdg9-ti2-prelaunch: fgc-tdg9-ti2-frz1
	PYTHONPATH=src:. $(EVOLUTION_PYTHON) -m unittest \
		tests/test_fgc_tdg9_ti2_authority.py \
		tests/test_fgc_tdg9_ti2_runner.py \
		tests/test_check_repo_tdg9_ti2_frz1.py -v
	$(EVOLUTION_PYTHON) scripts/check_repo.py --only-tdg9-ti2-frz1

status-fgc-tdg9-ti2:
	@set -eu; \
	authority_commit="$$(env -i PATH="$$PATH" HOME="$$HOME" LC_ALL=C LANG=C \
		git --no-replace-objects --no-optional-locks -c core.fsmonitor=false \
		-c core.untrackedCache=false -C "$(CURDIR)" rev-parse --verify 'HEAD^{commit}')"; \
	"$(EVOLUTION_PYTHON)" -I -B scripts/run_fgc_tdg9_ti2.py \
		--authority-commit "$$authority_commit" --status

# Outcome-neutral post-run binder. Ordinary verification reads only the
# compact certificate; --live is the explicit one-time 60-route shadow replay.
fgc-tdg9-ti2-pref1:
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_tdg9_ti2_pref1.py

verify-fgc-tdg9-ti2-pref1: fgc-tdg9-ti2-pref1
	PYTHONPATH=src:. $(EVOLUTION_PYTHON) -m unittest \
		tests/test_fgc_tdg9_ti2_pref1_binder.py \
		tests/test_check_repo_tdg9_ti2_pref1.py -v
	$(EVOLUTION_PYTHON) scripts/check_repo.py --only-tdg9-ti2-pref1

# Prospective retry-3 all-channel production-classifier audit. Ordinary
# verification consumes the tracked compact certificate and never constructs
# an SSPRK3 shadow proposal.
fgc-tdg9-ac1-frz1: fgc-tdg9-ti2-pref1
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_tdg9_ac1_frz1.py --verify-compact

verify-fgc-tdg9-ac1-prelaunch: fgc-tdg9-ac1-frz1
	PYTHONPATH=src:. $(EVOLUTION_PYTHON) -m unittest \
		tests/test_fgc_tdg9_ac1_authority.py \
		tests/test_fgc_tdg9_ac1_runner.py \
		tests/test_check_repo_tdg9_ac1_frz1.py -v
	$(EVOLUTION_PYTHON) scripts/check_repo.py --only-tdg9-ac1-frz1

status-fgc-tdg9-ac1:
	@set -eu; \
	authority_commit="$$(env -i PATH="$$PATH" HOME="$$HOME" LC_ALL=C LANG=C \
		git --no-replace-objects --no-optional-locks -c core.fsmonitor=false \
		-c core.untrackedCache=false -C "$(CURDIR)" rev-parse --verify 'HEAD^{commit}')"; \
	"$(EVOLUTION_PYTHON)" -I -B scripts/run_fgc_tdg9_ac1.py \
		--authority-commit "$$authority_commit" --status

# Outcome-neutral post-run binder. Ordinary verification reads only the
# compact certificate; --live is the explicit one-time raw/store bind.
fgc-tdg9-ac1-pref1:
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_tdg9_ac1_pref1.py

verify-fgc-tdg9-ac1-pref1: fgc-tdg9-ac1-pref1
	PYTHONPATH=src:. $(EVOLUTION_PYTHON) -m unittest \
		tests/test_fgc_tdg9_ac1_pref1_binder.py \
		tests/test_check_repo_tdg9_ac1_pref1.py -v
	$(EVOLUTION_PYTHON) scripts/check_repo.py --only-tdg9-ac1-pref1

# Prospective retry-3 u:R envelope-owner diagnostic. Ordinary verification
# consumes the tracked compact certificate and never constructs an SSPRK3
# shadow proposal.
fgc-tdg9-ur1-frz1: fgc-tdg9-ac1-pref1
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_tdg9_ur1_frz1.py --verify-compact

verify-fgc-tdg9-ur1-prelaunch: fgc-tdg9-ur1-frz1
	PYTHONPATH=src:. $(EVOLUTION_PYTHON) -m unittest \
		tests/test_fgc_tdg9_ur1_authority.py \
		tests/test_fgc_tdg9_ur1_runner.py \
		tests/test_check_repo_tdg9_ur1_frz1.py -v
	$(EVOLUTION_PYTHON) scripts/check_repo.py --only-tdg9-ur1-frz1

status-fgc-tdg9-ur1:
	@set -eu; \
	authority_commit="$$(env -i PATH="$$PATH" HOME="$$HOME" LC_ALL=C LANG=C \
		git --no-replace-objects --no-optional-locks -c core.fsmonitor=false \
		-c core.untrackedCache=false -C "$(CURDIR)" rev-parse --verify 'HEAD^{commit}')"; \
	"$(EVOLUTION_PYTHON)" -I -B scripts/run_fgc_tdg9_ur1.py \
		--authority-commit "$$authority_commit" --status

# Outcome-neutral post-run binder. Ordinary verification reads only the
# compact certificate; --live is the explicit one-time raw/store bind.
fgc-tdg9-ur1-pref1:
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_tdg9_ur1_pref1.py

verify-fgc-tdg9-ur1-pref1: fgc-tdg9-ur1-pref1
	PYTHONPATH=src:. $(EVOLUTION_PYTHON) -m unittest \
		tests/test_fgc_tdg9_ur1_pref1_binder.py \
		tests/test_check_repo_tdg9_ur1_pref1.py -v
	$(EVOLUTION_PYTHON) scripts/check_repo.py --only-tdg9-ur1-pref1

# Prospective retry-3 exact complete-C all-channel qualification. Compact
# verification is store-blind and never constructs an SSPRK3 shadow proposal.
fgc-tdg10-qa1-frz1: fgc-tdg9-ur1-pref1
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_tdg10_qa1_frz1.py --verify-compact

verify-fgc-tdg10-qa1-prelaunch: fgc-tdg10-qa1-frz1
	PYTHONPATH=src:. $(EVOLUTION_PYTHON) -m unittest \
		tests/test_fgc_tdg10_exact_complete_c_admission.py \
		tests/test_fgc_tdg10_exact_complete_c_runtime.py \
		tests/test_fgc_tdg10_qa1_authority.py \
		tests/test_fgc_tdg10_qa1_runner.py \
		tests/test_check_repo_tdg10_qa1_frz1.py -v
	$(EVOLUTION_PYTHON) scripts/check_repo.py --only-tdg10-qa1-frz1

status-fgc-tdg10-qa1:
	@set -eu; \
	authority_commit="$$(env -i PATH="$$PATH" HOME="$$HOME" LC_ALL=C LANG=C \
		git --no-replace-objects --no-optional-locks -c core.fsmonitor=false \
		-c core.untrackedCache=false -C "$(CURDIR)" rev-parse --verify 'HEAD^{commit}')"; \
	"$(EVOLUTION_PYTHON)" -I -B scripts/run_fgc_tdg10_qa1.py \
		--authority-commit "$$authority_commit" --status

# Independently bound retry-3 exact complete-C terminal. Ordinary verification
# consumes only the compact result and never invokes the seven-shadow live bind.
fgc-tdg10-qa1-pref1:
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_tdg10_qa1_pref1.py

verify-fgc-tdg10-qa1-pref1: fgc-tdg10-qa1-pref1
	PYTHONPATH=src:. $(EVOLUTION_PYTHON) -m unittest \
		tests/test_fgc_tdg10_qa1_pref1_binder.py \
		tests/test_check_repo_tdg10_qa1_pref1.py -v
	$(EVOLUTION_PYTHON) scripts/check_repo.py --only-tdg10-qa1-pref1

# Prospective retry-4/5 SSPRK3-on-inherited-SBP4 exact complete-C two-width
# robustness qualification. Compact verification is store-blind and never
# constructs an SSPRK3 shadow proposal. QA1 licenses this method-design freeze,
# never old-member adoption.
fgc-tdg10-qa2-frz1: fgc-tdg10-qa1-pref1
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_tdg10_qa2_frz1.py --verify-compact

verify-fgc-tdg10-qa2-prelaunch: fgc-tdg10-qa2-frz1
	PYTHONPATH=src:. $(EVOLUTION_PYTHON) -m unittest \
		tests/test_fgc_tdg10_qa2_authority.py \
		tests/test_fgc_tdg10_qa2_runner.py \
		tests/test_check_repo_tdg10_qa2_frz1.py -v
	$(EVOLUTION_PYTHON) scripts/check_repo.py --only-tdg10-qa2-frz1

status-fgc-tdg10-qa2:
	@set -eu; \
	authority_commit="$$(env -i PATH="$$PATH" HOME="$$HOME" LC_ALL=C LANG=C \
		git --no-replace-objects --no-optional-locks -c core.fsmonitor=false \
		-c core.untrackedCache=false -C "$(CURDIR)" rev-parse --verify 'HEAD^{commit}')"; \
	"$(EVOLUTION_PYTHON)" -I -B scripts/run_fgc_tdg10_qa2.py \
		--authority-commit "$$authority_commit" --status

# Independently bound QA2 retry-4/5 nonpass terminal. Ordinary verification
# consumes only the compact result. It never invokes --live, the QA2 runner,
# the raw namespace, or the campaign store.
fgc-tdg10-qa2-pref1:
	$(EVOLUTION_PYTHON) scripts/reproduce_fgc_tdg10_qa2_pref1.py

verify-fgc-tdg10-qa2-pref1: fgc-tdg10-qa2-pref1
	PYTHONPATH=src:. $(EVOLUTION_PYTHON) -m unittest \
		tests/test_fgc_tdg10_qa2_pref1_binder.py \
		tests/test_check_repo_tdg10_qa2_pref1.py -v
	$(EVOLUTION_PYTHON) scripts/check_repo.py --only-tdg10-qa2-pref1

# Explicit read-only live preflight. The bootstrap source comes from immutable
# commit A, never from the worktree. Both outer Git and Python receive a clean
# environment; the bootstrap then independently requires clean HEAD == C,
# A -> C, the exact metadata-only delta, and the complete guarded closure.
preflight-fgc-pro19-sid3-resume:
	@set -eu; \
	authority_commit="$$(env -i PATH="$(FGC_PRO19_SID3_SANITIZED_PATH)" LC_ALL="$(FGC_PRO19_SID3_SANITIZED_LOCALE)" \
		git --no-replace-objects --no-optional-locks \
		-c core.fsmonitor=false -c core.untrackedCache=false \
		-C "$(CURDIR)" rev-parse --verify 'HEAD^{commit}')"; \
	env -i PATH="$(FGC_PRO19_SID3_SANITIZED_PATH)" LC_ALL="$(FGC_PRO19_SID3_SANITIZED_LOCALE)" \
		git --no-replace-objects --no-optional-locks \
		-c core.fsmonitor=false -c core.untrackedCache=false \
		-C "$(CURDIR)" show \
		"$(FGC_PRO19_SID3_IMPLEMENTATION_COMMIT):scripts/bootstrap_fgc_pro19_sid3.py" | \
	env -i PATH="$(FGC_PRO19_SID3_SANITIZED_PATH)" LC_ALL="$(FGC_PRO19_SID3_SANITIZED_LOCALE)" \
		"$(EVOLUTION_PYTHON)" -I -B - \
		--root "$(CURDIR)" \
		--implementation-commit "$(FGC_PRO19_SID3_IMPLEMENTATION_COMMIT)" \
		--authority-commit "$$authority_commit" \
		--process-marker "$(FGC_PRO19_SID3_PROCESS_MARKER)" \
		--preflight

# Explicit heavy resume. It is deliberately absent from every ordinary test,
# check, compact verifier, and foundation dependency. The guarded runner
# repeats the complete live preflight immediately before writer acquisition.
preflight-fgc-pro19-cfl1-resume:
	$(EVOLUTION_PYTHON) scripts/run_fgc_pro19_event1.py --resume-preflight

status-fgc-pro19-event1:
	$(EVOLUTION_PYTHON) scripts/status_fgc_pro19_event1.py

verify-fgc-pro14-result: fgc-cal11-pref22
	$(PYTHON) scripts/check_repo.py

# Focused SF1 integration boundary. HLT13/HLT14 qualify bounded synthetic
# runtime contracts; PRO18/PREF26/AUTH1 now seal compact production inputs and
# one commit-independent calibration GenesisSpec. HLT15/PREF27 bind the
# authenticated origin; HLT16 qualifies the durable runtime. The first real
# attempt advanced five finite RK4-2049 steps; SID1/SID2 bind the no-advance
# generation-eight recovery. SID3 adds the committed GR-0 resume-image
# authority; the terminal source is preserved through generation ten. RCV3
# freezes and PREF1 independently authenticates the exact generation-nine
# projection. AUTH1 authorized the unchanged same-event GR-0 continuation;
# REC1 froze the exact metadata-only recovery forced by its first durable
# rejection suffix. PREF2 now binds its numerical temporal-retry-exhaustion
# terminal while preserving the accepted physical state. Common-event
# AR1/PREF1 binds a finite sampled exact-discrete obstruction. LOC1/ERR1
# localizes the v1 boundary-root defect; LOC2/PREF2 independently authenticates
# all 60 repaired primary-v2 localizations and the ten component-ownership
# reductions without selecting a temporal remedy. Completion, calibration,
# candidate execution, and physics remain outside this foundation target.
# FGC-1-TDG10-QA1-FRZ1/PREF1 and FGC-1-TDG10-QA2-FRZ1/PREF1 enter the
# foundation only through compact evidence. The seven-shadow and two-width
# runs remain explicit commands.
verify-fgc-sf1-foundation: fgc-cal10-pref15 fgc-tdg1-frz1 fgc-tdg1-pref16 fgc-tdg2-frz1 fgc-tdg2-pref17 fgc-tdg3-frz1 fgc-tdg3-pref18 fgc-tdg4-frz1 fgc-tdg4-pref19 fgc-tdg5-frz1 fgc-tdg5-pref20 fgc-tdg5-imp1 fgc-tdg6-frz1 fgc-tdg6-pref21 fgc-tdg6-imp2 fgc-cal11-pref22 fgc-tdg7-frz1 fgc-tdg7-pref23 fgc-tdg7-imp3 fgc-pro15-frz1 fgc-hlt13-mon13 fgc-pro16-frz1 fgc-pro16-pref24 fgc-pro17-frz1 fgc-pro17-pref25 fgc-hlt14-mon14 fgc-pro18-frz1 fgc-pro18-auth1 fgc-pro19-frz1 fgc-hlt16-mon16 fgc-pro19-sid1-frz1 fgc-pro19-sid2-pref1 fgc-pro19-sid3-auth1 fgc-pro19-pref28 fgc-tdg8-run1-auth1 fgc-tdg8-rcv1-auth1 fgc-tdg8-rcv2-auth1 fgc-tdg8-rcv3-frz1 fgc-tdg8-rcv3-pref1 fgc-tdg8-rcv3-auth1 fgc-tdg8-rcv3-rec1-auth1 fgc-tdg8-rcv3-pref2 fgc-tdg9-ar1-auth1 fgc-tdg9-ar1-pref1 fgc-tdg9-loc1-frz1 fgc-tdg9-loc2-frz1 fgc-tdg9-loc2-pref2 fgc-tdg9-ti2-frz1 fgc-tdg9-ti2-pref1 fgc-tdg9-ac1-frz1 fgc-tdg9-ac1-pref1 fgc-tdg9-ur1-frz1 fgc-tdg9-ur1-pref1 fgc-tdg10-qa1-frz1 fgc-tdg10-qa1-pref1 fgc-tdg10-qa2-frz1 fgc-tdg10-qa2-pref1
	$(PYTHON) -m unittest tests/test_fgc_evolution_nonlinear_source.py tests/test_fgc_src1_nl1_reproduction.py tests/test_fgc_constraint_system.py tests/test_fgc_con4_phy1_reproduction.py tests/test_fgc_regular_center.py tests/test_fgc_ctr1_reg1_reproduction.py tests/test_fgc_initial_data_preflight.py tests/test_fgc_id0_pref1_reproduction.py tests/test_fgc_initial_data_family.py tests/test_fgc_id1_fam1_reproduction.py tests/test_fgc_initial_state_bridge.py tests/test_fgc_covariant_principal_health.py tests/test_fgc_multidirectional_health.py tests/test_fgc_dom4_run1_reproduction.py tests/test_fgc_hyp2_md1_reproduction.py tests/test_fgc_boundary_domain.py tests/test_fgc_bnd2_cp1_reproduction.py tests/test_fgc_health_monitor.py tests/test_fgc_hlt1_mon1_reproduction.py tests/test_fgc_numerical_engine.py tests/test_fgc_runtime_transaction.py tests/test_fgc_vectorized_source.py tests/test_fgc_numerical_controls.py tests/test_fgc_affine_measurement.py tests/test_fgc_num1_val1_reproduction.py tests/test_fgc_run1_authorization.py tests/test_fgc_run1_sym1_reproduction.py tests/test_fgc_gr0_direct_source.py tests/test_fgc_gr0_calibration.py tests/test_fgc_cal0_pref2_reproduction.py tests/test_fgc_protocol_v4.py tests/test_fgc_proto4_admission.py tests/test_fgc_hlt2_mon2_reproduction.py tests/test_fgc_static_initial_admission.py tests/test_fgc_id2_all1_reproduction.py tests/test_fgc_calibration_runtime.py tests/test_fgc_cal1_pref3_reproduction.py tests/test_fgc_protocol_v5.py tests/test_fgc_proto5_runtime.py tests/test_fgc_gr0_campaign_runner.py tests/test_fgc_cal2_source_diagnosis.py tests/test_fgc_cal2_pref4_reproduction.py tests/test_fgc_protocol_v6.py tests/test_fgc_proto6_runtime.py tests/test_fgc_gr0_campaign_runner_v6.py tests/test_fgc_hlt4_mon4_reproduction.py tests/test_fgc_cal3_pref5_reproduction.py tests/test_fgc_protocol_v7.py tests/test_fgc_proto7_runtime.py tests/test_fgc_gr0_campaign_runner_v7.py tests/test_fgc_hlt5_mon5_reproduction.py tests/test_fgc_cal4_common_event_diagnosis.py tests/test_fgc_cal4_pref6_reproduction.py tests/test_fgc_protocol_v8.py tests/test_fgc_proto8_runtime.py tests/test_fgc_gr0_campaign_runner_v8.py tests/test_fgc_hlt6_mon6_reproduction.py tests/test_fgc_cal5_resolution_diagnosis.py tests/test_fgc_cal5_pref7_reproduction.py tests/test_fgc_protocol_v9.py tests/test_fgc_pro9_frz1_reproduction.py tests/test_fgc_proto9_runtime.py tests/test_fgc_gr0_campaign_runner_v9.py tests/test_fgc_hlt7_mon7_reproduction.py tests/test_fgc_spectral_sensitivity.py tests/test_fgc_cal6_pref8_reproduction.py tests/test_fgc_protocol_v10.py tests/test_fgc_pro10_frz1_reproduction.py tests/test_fgc_proto10_runtime.py tests/test_fgc_hlt8_mon8_reproduction.py tests/test_fgc_gr0_campaign_runner_v10.py tests/test_fgc_cal7_center_roundoff_diagnosis.py tests/test_fgc_cal7_pref9_reproduction.py tests/test_fgc_protocol_v11.py tests/test_fgc_pro11_frz1_reproduction.py tests/test_fgc_proto11_runtime.py tests/test_fgc_hlt9_mon9_reproduction.py tests/test_fgc_gr0_campaign_runner_v11.py tests/test_fgc_cal8_proto11_campaign_diagnosis.py tests/test_fgc_cal8_pref10_reproduction.py tests/test_fgc_cal9_proto12_campaign_diagnosis.py tests/test_fgc_cal9_pref13_reproduction.py -v
	$(EVOLUTION_PYTHON) scripts/check_repo.py

# Reserved for the eventual result-bearing outer proof.  A green foundation
# must never be mistaken for DEF1/ROB1 completion.
verify-fgc-sf1:
	@echo "FGC-2-SF1 is incomplete: use verify-fgc-sf1-foundation for the sealed pre-result boundary" >&2
	@false

check: reproduce paper
	$(EVOLUTION_PYTHON) scripts/check_repo.py

data:
	$(PYTHON) scripts/verify_data.py

paper:
	$(PAPER_PYTHON) scripts/build_paper.py

# The test suite already recomputes every canonical record in memory and the
# repository audit verifies the tracked bytes and implementation hashes.  Do
# not call `reproduce` again here: that would duplicate the exact DOM3 proof.
# The outer audit inherits the pinned evolution runtime because SID3 verifies
# the exact NumPy package, native extension, and distribution metadata image.
verify: test paper data
	$(EVOLUTION_PYTHON) scripts/check_repo.py

verify-portable-clean-clone: test-portable-clean-clone paper data
	$(PYTHON) scripts/check_repo.py

clean:
	$(PAPER_PYTHON) scripts/build_paper.py --clean
