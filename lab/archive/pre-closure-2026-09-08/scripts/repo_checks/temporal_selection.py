"""Current TDG11 compact routing, independent of live diagnostic state."""

from __future__ import annotations

import sys

from ._shared import REPOSITORY


def _check_tdg11_msel1_frz1_result() -> list[str]:
    previous = list(sys.path)
    try:
        sys.path.insert(0, str(REPOSITORY / "src"))
        from recursive_horizons.fgc.evolution.tdg11_msel1_authority import validate_compact_bundle

        validate_compact_bundle(REPOSITORY)
    except (OSError, ValueError, RuntimeError) as error:
        return [f"TDG11-MSEL1 compact freeze differs: {error}"]
    finally:
        sys.path[:] = previous
    return []


def _check_tdg11_msel1_pref1_result() -> list[str]:
    previous = list(sys.path)
    try:
        sys.path.insert(0, str(REPOSITORY / "src"))
        from recursive_horizons.fgc.evolution.tdg11_msel1_pref1_binder import (
            verify_compact,
        )

        verify_compact(REPOSITORY)
    except (OSError, ValueError, RuntimeError) as error:
        return [f"TDG11-MSEL1-PREF1 compact result differs: {error}"]
    finally:
        sys.path[:] = previous
    return []


def _check_tdg11_imp1_result() -> list[str]:
    previous = list(sys.path)
    try:
        sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]
        from scripts.reproduce_fgc_tdg11_imp1 import verify_compact

        verify_compact(REPOSITORY)
    except (OSError, ValueError, RuntimeError) as error:
        return [f"TDG11-IMP1 compact implementation result differs: {error}"]
    finally:
        sys.path[:] = previous
    return []


def _check_hlt17_srcq1_rec1_pref1_result() -> list[str]:
    previous = list(sys.path)
    try:
        sys.path.insert(0, str(REPOSITORY / "src"))
        from recursive_horizons.fgc.evolution.hlt17_srcq1_rec1_pref1_binder import (
            verify_compact,
        )

        verify_compact(REPOSITORY)
    except (OSError, ValueError, RuntimeError) as error:
        return [f"HLT17-SRCQ1-REC1-PREF1 compact result differs: {error}"]
    finally:
        sys.path[:] = previous
    return []


def _check_pro20_ev1_pref1_result() -> list[str]:
    previous = list(sys.path)
    try:
        sys.path.insert(0, str(REPOSITORY / "src"))
        from recursive_horizons.fgc.evolution.pro20_ev1_pref1_certificate import (
            verify_compact,
        )

        verify_compact(REPOSITORY)
    except (OSError, ValueError, RuntimeError) as error:
        return [f"PRO20-EV1-PREF1 compact result differs: {error}"]
    finally:
        sys.path[:] = previous
    return []


__all__ = [
    "_check_tdg11_msel1_frz1_result",
    "_check_tdg11_msel1_pref1_result",
    "_check_tdg11_imp1_result",
    "_check_hlt17_srcq1_rec1_pref1_result",
    "_check_pro20_ev1_pref1_result",
]
