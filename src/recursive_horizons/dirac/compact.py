"""Use the retained compact-matching owner without duplicating its numerics."""
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys


def _matching_owner():
    path = Path(__file__).resolve().parents[3] / "scripts/check_nsc_compact_matching.py"
    name = "check_nsc_compact_matching"
    loaded = sys.modules.get(name)
    if loaded is not None and Path(loaded.__file__).resolve() == path:
        return loaded
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"compact-matching owner unavailable: {path}")
    owner = module_from_spec(spec)
    sys.modules[name] = owner
    spec.loader.exec_module(owner)
    return owner


_owner = _matching_owner()
proper_ratio = _owner.proper_ratio
zero_limit = _owner.zero_limit
matched_weight = _owner.matched_weight
spectral_moments = _owner.spectral_moments
local_coefficients = _owner.local_coefficients
flat_moments = _owner.flat_moments
static_basis_coefficients = _owner.static_basis_coefficients
