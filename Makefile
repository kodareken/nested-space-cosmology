PYTHON ?= python3
JOBS ?= auto

.PHONY: help install verify reproduce reproduce-development reproduce-exact paper paper-check test check

help:
	@echo "Nested-Space Cosmology public commands"
	@echo "  make install          install the pinned project and paper dependencies"
	@echo "  make check            validate manifests, claims, paths, and links"
	@echo "  make test             run focused publication tests"
	@echo "  make reproduce        recompute the released chain and development evidence"
	@echo "  make reproduce-development  check subsequent development records"
	@echo "  make reproduce-exact  require byte-identical recomputation"
	@echo "  make paper            rebuild the tracked paper PDF"
	@echo "  make paper-check      rebuild twice and verify deterministic bytes"
	@echo "  make verify           run the integrated public-repository gate"

install:
	$(PYTHON) -m pip install -e '.[paper]'

check:
	$(PYTHON) scripts/check_publication.py

test:
	$(PYTHON) -m unittest discover -s tests -v

reproduce-development:
	$(PYTHON) scripts/check_nsc_compact_interaction.py --check
	$(PYTHON) scripts/check_nsc_torsion_uv_map.py --check

reproduce: reproduce-development
	$(PYTHON) scripts/reproduce_public_results.py --mode portable --jobs $(JOBS)

reproduce-exact:
	$(PYTHON) scripts/reproduce_public_results.py --mode exact --jobs $(JOBS)

paper:
	$(PYTHON) scripts/build_paper.py

paper-check:
	$(PYTHON) scripts/build_paper.py --check

verify: check test reproduce paper-check

.PHONY: demonstrate
demonstrate:
	$(PYTHON) -B scripts/demonstrate.py
