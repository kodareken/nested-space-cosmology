"""Make the recursive_horizons package importable for every lab test.

Several modules under lab/src are only reachable when that directory is on
the path, and tests that import them otherwise fail at collection for a
reason that has nothing to do with the test.

Nineteen laboratory files also import `resource` or `fcntl`, which are
POSIX-only and absent on Windows. The stand-ins in
`recursive_horizons.nsc_posix` are registered here so those files import and
their tests can run. The stand-ins are documented there, including what they
do not do.

`lab/` itself is on the path as well, because two test modules import from
`tests.test_nsc_ks_difference_envelope`. Without it those two fail at
collection for a reason that has nothing to do with what they check.
"""
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
for _path in (str(LAB), str(LAB / "src"), str(LAB / "scripts")):
    if _path not in sys.path:
        sys.path.insert(0, _path)

from recursive_horizons import nsc_posix  # noqa: E402

_installed = nsc_posix.install()
if _installed:
    import warnings

    warnings.warn(
        "POSIX stand-ins registered for: " + ", ".join(_installed)
        + ". CPU budget guards do not fire and advisory run locks are not"
        " enforced on this platform; see recursive_horizons.nsc_posix.",
        RuntimeWarning,
        stacklevel=1,
    )
