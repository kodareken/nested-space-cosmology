"""Stand-ins for two POSIX-only modules, for running this laboratory on Windows.

The checkout was developed on macOS. Nineteen files import `resource` for CPU
budgets and resident-memory readings, and `fcntl` for advisory run locks.
Neither module exists on Windows, so those files fail at import for a reason
that has nothing to do with what they compute.

What the stand-ins actually do, stated plainly so nothing is hidden:

`resource.getrusage(RUSAGE_SELF)` reports CPU time through
`time.process_time` and reports zero resident memory. A zero memory reading is
not a measurement; callers that compare it against a threshold must treat the
absence as "not measured" rather than as "small".

`resource.getrlimit` and `setrlimit` do nothing. A CPU budget guard therefore
does not fire. Callers relying on it must enforce their own wall-clock limit.

`fcntl.flock` does nothing. Advisory run locks are not enforced, so two
processes could collide over one run directory. On a single-machine Windows
run this has not bitten, but it is a real difference and not a silent
improvement.

Nothing here is imported by the package itself. `install()` registers the
stand-ins in `sys.modules` only when the real modules are absent, which is how
the laboratory test suite reaches those nineteen files on Windows.
"""

from __future__ import annotations

from pathlib import Path

import sys
import time
import types

#: Recorded so a caller, or a reader of a record, can tell which came from here.
STANDIN = True


class _Usage:
    """The fields callers read from `resource.getrusage`."""

    __slots__ = ("ru_utime", "ru_stime", "ru_maxrss")

    def __init__(self, cpu_seconds: float) -> None:
        self.ru_utime = cpu_seconds
        self.ru_stime = 0.0
        # Not measured on this platform. Zero means unknown, not small.
        self.ru_maxrss = 0


def build_resource() -> types.ModuleType:
    module = types.ModuleType("resource")
    module.RUSAGE_SELF = 0
    module.RLIMIT_CPU = 0
    module.RLIM_INFINITY = -1
    module.STANDIN = True

    def getrusage(_who: int = 0) -> _Usage:
        return _Usage(time.process_time())

    def getrlimit(_which: int = 0):
        return (-1, -1)

    def setrlimit(_which: int = 0, _limits=None) -> None:
        return None

    module.getrusage = getrusage
    module.getrlimit = getrlimit
    module.setrlimit = setrlimit
    return module


def build_fcntl() -> types.ModuleType:
    module = types.ModuleType("fcntl")
    module.LOCK_EX = 2
    module.LOCK_SH = 1
    module.LOCK_NB = 4
    module.LOCK_UN = 8
    module.STANDIN = True

    def flock(_fd: int, _operation: int) -> None:
        """No-op. See the module docstring: locks are not enforced here."""
        return None

    module.flock = flock
    return module


def install() -> list[str]:
    """Register a stand-in for each absent module. Returns what was installed."""
    installed = []
    try:
        import resource  # noqa: F401
    except ImportError:
        sys.modules["resource"] = build_resource()
        installed.append("resource")
    try:
        import fcntl  # noqa: F401
    except ImportError:
        sys.modules["fcntl"] = build_fcntl()
        installed.append("fcntl")
    return installed


if __name__ == "__main__":
    print("installed:", install() or "nothing (both modules are present)")

ROOT = Path(__file__).resolve().parents[2]


def repo_path(*parts: str) -> Path:
    """Resolve a current input inside the single canonical checkout."""
    path = ROOT.joinpath(*parts).resolve()
    if not path.is_relative_to(ROOT):
        raise ValueError("path escapes the repository")
    return path
