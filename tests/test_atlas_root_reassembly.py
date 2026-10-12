"""Current root producers cannot impersonate the frozen pre-migration source."""
import importlib
from pathlib import Path

import pytest


@pytest.mark.parametrize("version", range(1, 7))
def test_current_producer_has_one_root_and_rejects_old_source_binding(version):
    suffix = "" if version == 1 else f"_v{version}"
    producer = importlib.import_module("derive_nsc_cosmological_atlas" + suffix)
    root = Path(__file__).resolve().parents[1]
    assert producer.LAB == producer.REPO == root
    sources = producer.source_bindings()
    assert sources
    assert not any(path.startswith("lab/") for path in sources)
    assert all((root / path).is_file() for path in sources)
    with pytest.raises(ValueError, match="source binding mismatch"):
        producer.check(root / f"results/development/nsc-cosmological-atlas-v{version}")
