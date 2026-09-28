"""Metadata recovery cannot authorize a changed calculation or invalid binding."""
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"scripts"))
from derive_nsc_ks_value_metadata_successor import (
    corrected_metadata_for_validation, unchanged_numerical_code)


def test_only_the_writer_may_change():
    old = 'import numpy\n\ndef save_family():\n    pass\n\ndef evolve():\n    return 1\n'
    new = old.replace('    pass', '    return {"evaluation_identity": "bound"}')
    unchanged_numerical_code(old, new)
    with pytest.raises(ValueError, match="numerical code changed"):
        unchanged_numerical_code(old, new.replace('return 1', 'return 2'))
    with pytest.raises(ValueError, match="numerical code changed"):
        unchanged_numerical_code(old, new.replace('import numpy', 'import scipy'))


def test_only_the_redundant_missing_identity_is_inferred_without_mutating_inputs():
    record = {"evaluation_identity": "x", "binding_digest": "x"}
    arrays = {"metadata": {"binding_digest": "x"}}
    result = corrected_metadata_for_validation(record, arrays)
    assert result["metadata"]["evaluation_identity"] == "x"
    assert "evaluation_identity" not in arrays["metadata"]
    with pytest.raises(ValueError, match="binding identity differ"):
        corrected_metadata_for_validation(record, {"metadata": {"binding_digest": "y"}})
    with pytest.raises(ValueError, match="omitted identity"):
        corrected_metadata_for_validation(record, result)
