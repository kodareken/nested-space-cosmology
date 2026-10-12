#!/usr/bin/env python3
"""Stable CLI and compatibility facade for modular repository checks."""

from __future__ import annotations

import sys
from pathlib import Path
from types import ModuleType

if __package__:
    from .repo_checks import _shared as _shared_module
    from .repo_checks import calibration as _calibration_module
    from .repo_checks import catalog as _catalog_module
    from .repo_checks import core as _core_module
    from .repo_checks import progression as _progression_module
    from .repo_checks import public as _public_module
    from .repo_checks import temporal_diagnostics as _temporal_module
    from .repo_checks import temporal_selection as _selection_module
    from .repo_checks.core import *  # noqa: F401,F403
    from .repo_checks.core import __all__ as _CORE_ALL
else:
    _SCRIPT_DIRECTORY = str(Path(__file__).resolve().parent)
    if _SCRIPT_DIRECTORY not in sys.path:
        sys.path.insert(0, _SCRIPT_DIRECTORY)
    from repo_checks import _shared as _shared_module  # type: ignore[no-redef]
    from repo_checks import calibration as _calibration_module  # type: ignore[no-redef]
    from repo_checks import catalog as _catalog_module  # type: ignore[no-redef]
    from repo_checks import core as _core_module  # type: ignore[no-redef]
    from repo_checks import progression as _progression_module  # type: ignore[no-redef]
    from repo_checks import public as _public_module  # type: ignore[no-redef]
    from repo_checks import temporal_diagnostics as _temporal_module  # type: ignore[no-redef]
    from repo_checks import temporal_selection as _selection_module  # type: ignore[no-redef]
    from repo_checks.core import *  # type: ignore[no-redef]  # noqa: F401,F403
    from repo_checks.core import __all__ as _CORE_ALL  # type: ignore[no-redef]

__all__ = _CORE_ALL

_DOMAIN_MODULES = (
    _shared_module,
    _public_module,
    _calibration_module,
    _progression_module,
    _temporal_module,
    _selection_module,
    _catalog_module,
    _core_module,
)


class _CompatibilityFacade(ModuleType):
    """Keep legacy monkey-patches effective across the split modules.

    Historical focused tests import ``scripts.check_repo`` and patch private
    helpers or ``REPOSITORY`` on that module.  The old monolith made those
    patches affect every check.  Propagating assignments to each domain that
    exports the same name preserves that compatibility without moving the
    implementation back into this dispatcher.
    """

    def __setattr__(self, name: str, value: object) -> None:
        super().__setattr__(name, value)
        if name.startswith("_CompatibilityFacade"):
            return
        for module in _DOMAIN_MODULES:
            if name in module.__dict__:
                setattr(module, name, value)


sys.modules[__name__].__class__ = _CompatibilityFacade


if __name__ == "__main__":
    raise SystemExit(main())
