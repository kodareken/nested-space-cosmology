"""One root package, with explicit Windows limitations for historical tests."""
from pathlib import Path
import sys
import warnings
import pytest

ROOT = Path(__file__).resolve().parents[1]
for _path in (str(ROOT), str(ROOT / "src"), str(ROOT / "scripts")):
    if _path not in sys.path:
        sys.path.insert(0, _path)

from recursive_horizons import nsc_posix  # noqa: E402

_installed = nsc_posix.install()
if _installed:
    warnings.warn(
        "POSIX stand-ins registered for: " + ", ".join(_installed)
        + ". CPU budget guards do not fire and advisory run locks are not"
        " enforced on this platform; see recursive_horizons.nsc_posix.",
        RuntimeWarning,
        stacklevel=1,
    )


@pytest.fixture(scope="session")
def original_middle_correction():
    """Authenticate the original witness independently of fresh adaptive paths."""
    import hashlib
    import json
    import derive_nsc_vacuum_source_correction as producer
    record_bytes = (ROOT / producer.OUTPUT).read_bytes()
    assert hashlib.sha256(record_bytes).hexdigest() == (
        "3e5e31cd670288eac26e8cb139d3b24c886230962560dd0813e76595c8772fb5")
    prior = json.loads(record_bytes)
    replayed, _ = producer.calculate(
        (prior, (ROOT / producer.PAYLOAD).read_bytes()),
        source_context="329ef809e466a86ed4a290ee32c19485f85add9d")
    return replayed
