PYTHON ?= python3
HISTORY_PYTHON ?= $(PYTHON)
JOBS ?= auto

.PHONY: help install verify reproduce reproduce-development reproduce-exact paper paper-check test check

help:
	@echo "Nested-Space Cosmology public commands"
	@echo "  make install          install the pinned project and paper dependencies"
	@echo "  make check            validate manifests, claims, paths, and links"
	@echo "  make test             run affected atlas, preservation and publication tests"
	@echo "  make demonstrate      display authenticated results without recomputing"
	@echo "  make demonstrate-recompute  explicitly rerun the selected demonstrations"
	@echo "  make reproduce*       retired campaign routes; select an explicit historical context"
	@echo "  make paper            rebuild the tracked paper PDF"
	@echo "  make paper-check      rebuild twice and verify deterministic bytes"
	@echo "  make draft-check      authenticate the focused OPEN draft and its evidence"
	@echo "  make finite-check     authenticate the measured finite-realization article and complete evidence graph"
	@echo "  make history-check    replay original editions with the selected HISTORY_PYTHON profile"
	@echo "  make verify           validate the current paper, atlas and preserved editions"

install:
	$(PYTHON) -m pip install -e '.[paper,dev]'

check:
	$(PYTHON) scripts/check_publication.py
	$(PYTHON) scripts/check_lab_snapshot.py
	$(PYTHON) scripts/check_claim_bindings.py

test:
	$(PYTHON) scripts/test_reassembly.py

reproduce reproduce-development reproduce-exact:
	@echo "This campaign route is retired: merged current sources are not the frozen release context."
	@echo "Use make history-check to authenticate preserved results without restarting the campaigns."
	@echo "A campaign rerun needs a fresh admission and budget, then an explicitly selected context from scripts/replay_history.py prepare."
	@echo "See docs/reproducing.md for the historical source boundary."
	@$(PYTHON) -c "raise SystemExit(2)"

paper:
	$(PYTHON) scripts/build_paper.py

paper-check:
	$(PYTHON) scripts/build_paper.py --check

verify: check test paper-check history-check

.PHONY: history-check
history-check:
	$(HISTORY_PYTHON) scripts/replay_history.py all

.PHONY: finite-paper finite-check
finite-paper:
	@echo "This edition is frozen. Use scripts/replay_history.py prepare and its explicit historical context for rebuild experiments."

finite-check:
	$(PYTHON) scripts/verify_finite_regeneration.py

.PHONY: draft draft-check
draft:
	@echo "This edition is frozen. Use scripts/replay_history.py prepare and its explicit historical context for rebuild experiments."

draft-check:
	$(PYTHON) scripts/verify_local_gate_draft.py

.PHONY: demonstrate demonstrate-recompute
demonstrate:
	$(PYTHON) -B scripts/demonstrate.py

demonstrate-recompute:
	$(PYTHON) -B scripts/demonstrate.py --recompute
