"""Current publication version must not relabel the frozen foundation."""
import importlib.util
import json
from pathlib import Path
import tomllib

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("current_publication_provenance", ROOT / "scripts/build_publication_provenance.py")
provenance = importlib.util.module_from_spec(spec)
spec.loader.exec_module(provenance)


def test_current_package_version_and_frozen_foundation_have_distinct_owners(tmp_path):
    (tmp_path / "paper").mkdir()
    (tmp_path / "pyproject.toml").write_text('[project]\nversion = "0.27.0"\n')
    (tmp_path / "paper/metadata.json").write_text(json.dumps({"version": "0.26.0"}))
    assert provenance.publication_versions(tmp_path) == {
        "public_version": "0.27.0", "foundation_version": "0.26.0"}


def test_public_import_version_matches_the_current_project_version():
    package_spec = importlib.util.spec_from_file_location(
        "current_public_core_version", ROOT / "src/recursive_horizons/__init__.py")
    package = importlib.util.module_from_spec(package_spec)
    package_spec.loader.exec_module(package)
    project = tomllib.loads((ROOT / "pyproject.toml").read_text())
    assert package.__version__ == project["project"]["version"]
