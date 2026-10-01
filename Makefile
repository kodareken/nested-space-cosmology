PYTHON ?= python3
JOBS ?= auto

.PHONY: help install verify reproduce reproduce-development reproduce-exact paper paper-check test check

help:
	@echo "Nested-Space Cosmology public commands"
	@echo "  make install          install the pinned project and paper dependencies"
	@echo "  make check            validate manifests, claims, paths, and links"
	@echo "  make test             run focused publication tests"
	@echo "  make demonstrate      display authenticated results without recomputing"
	@echo "  make demonstrate-recompute  explicitly rerun the selected demonstrations"
	@echo "  make reproduce        recompute the released chain and development evidence"
	@echo "  make reproduce-development  check subsequent development records"
	@echo "  make reproduce-exact  require byte-identical recomputation"
	@echo "  make paper            rebuild the tracked paper PDF"
	@echo "  make paper-check      rebuild twice and verify deterministic bytes"
	@echo "  make draft-check      authenticate the focused OPEN draft and its evidence"
	@echo "  make finite-check     authenticate the measured finite-realization article and complete evidence graph"
	@echo "  make finite-paper     rebuild the finite article twice from its pinned science commit"
	@echo "  make draft            rebuild the focused draft twice (requires pdflatex/bibtex)"
	@echo "  make verify           validate locked evidence, tests, and both papers"

install:
	$(PYTHON) -m pip install -e '.[paper,dev]'

check:
	$(PYTHON) scripts/check_publication.py
	$(PYTHON) scripts/check_lab_snapshot.py

test:
	$(PYTHON) -m pytest -q

reproduce-development:
	$(PYTHON) scripts/reproduce_public_results.py --mode portable --jobs $(JOBS) --only results/development/compact-interaction.json,results/development/torsion-uv-map.json,results/development/flow-compatibility.json,results/development/charged-sector.json,results/development/vacuum-charge-matching.json,results/development/compact-boundary-action.json,results/development/compact-casimir.json,results/development/horizon-source.json,results/development/warped-source.json,results/development/compact-matching.json,results/development/gauge-source.json,results/development/spherical-action.json,results/development/curvature-eft.json,results/development/spectral-endpoint.json,results/development/child-state.json,results/development/massless-reference.json,results/development/angular-stress.json,results/development/unruh-state.json,results/development/state-regulator.json

reproduce:
	$(PYTHON) scripts/reproduce_public_results.py --mode portable --jobs $(JOBS)

reproduce-exact:
	$(PYTHON) scripts/reproduce_public_results.py --mode exact --jobs $(JOBS)

paper:
	$(PYTHON) scripts/build_paper.py

paper-check:
	$(PYTHON) scripts/build_paper.py --check

verify: check test paper-check draft-check finite-check

.PHONY: finite-paper finite-check
finite-paper:
	$(PYTHON) scripts/build_finite_regeneration.py --science-commit 8a52257fc2829c73d485c7b4c6e310baadfe3c38

finite-check:
	$(PYTHON) scripts/verify_finite_regeneration.py

.PHONY: draft draft-check
draft:
	$(PYTHON) scripts/build_local_gate_draft.py --check

draft-check:
	$(PYTHON) scripts/verify_local_gate_draft.py

.PHONY: demonstrate demonstrate-recompute
demonstrate:
	$(PYTHON) -B scripts/demonstrate.py

demonstrate-recompute:
	$(PYTHON) -B scripts/demonstrate.py --recompute
