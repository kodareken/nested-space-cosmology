#!/usr/bin/env python3
"""Stage only the pinned, hash-checked public allowlist from a research checkout."""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
from pathlib import Path, PurePosixPath
import subprocess

ROOT = Path(__file__).resolve().parents[1]
SPEC = ROOT / 'results/release-spec.json'


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def safe_relative(relative: str) -> str:
    path = PurePosixPath(relative)
    fixture = (len(path.parts) == 5 and path.parts[:3] == ('tests', 'fixtures', 'source-checkpoints')
               and len(path.parts[3]) == 40 and all(c in '0123456789abcdef' for c in path.parts[3])
               and path.name == 'pyproject.toml')
    native = relative == 'src/recursive_horizons/_angular_transport.cpp'
    if (path.is_absolute() or '..' in path.parts or len(path.parts) < 2
            or path.parts[0] not in {'results', 'scripts', 'src', 'tests', 'docs'}
            or (path.suffix not in {'.json', '.py', '.md'} and not fixture and not native)
            or any(part.startswith('.') for part in path.parts)):
        raise ValueError(f'path is outside the curated scientific allowlist: {relative}')
    return relative


def function_source(data: bytes, name: str) -> bytes:
    source = data.decode('utf-8')
    node = next(node for node in ast.parse(source).body
                if isinstance(node, ast.FunctionDef) and node.name == name)
    return ''.join(source.splitlines(keepends=True)[node.lineno-1:node.end_lineno]).encode()


def stage(source: Path, destination: Path, spec: dict, *, check: bool = False) -> dict:
    commit = spec['source_commit']
    resolved = subprocess.check_output(['git', 'rev-parse', '--verify', commit + '^{commit}'],
                                       cwd=source, text=True).strip()
    if resolved != commit:
        raise ValueError('release source must be pinned by its full commit hash')
    def blob(relative, target=None):
        if relative == 'pyproject.toml':
            if target is None or not safe_relative(target).endswith('/pyproject.toml'):
                raise ValueError('laboratory configuration requires an immutable source fixture')
        else:
            safe_relative(relative)
            if target is not None and target != relative:
                raise ValueError('only the declared laboratory configuration can be relocated')
        return subprocess.check_output(['git', 'show', f'{commit}:{relative}'], cwd=source)
    prepared = []
    seen = set()
    for entry in spec['import_files'] + spec['retained_byte_identical_files']:
        relative = safe_relative(entry['path'])
        if relative in seen:
            raise ValueError(f'duplicate import path: {relative}')
        seen.add(relative)
        data = blob(entry.get('source_path', relative), relative)
        if digest(data) != entry['sha256']:
            raise ValueError(f'pinned source hash mismatch: {relative}')
        target = destination / relative
        if target.exists() and target.read_bytes() != data:
            raise ValueError(f'refusing to overwrite different public bytes: {relative}')
        if entry in spec['retained_byte_identical_files'] and not target.is_file():
            raise ValueError(f'retained public dependency is absent: {relative}')
        prepared.append((target, data))
    comparator = spec['retained_comparator']
    comparator_path = destination / safe_relative(comparator['path'])
    if digest(comparator_path.read_bytes()) != comparator['sha256']:
        raise ValueError('public comparison helper changed')
    for data in (comparator_path.read_bytes(), blob(comparator['path'])):
        if digest(function_source(data, comparator['function'])) != comparator['function_sha256']:
            raise ValueError('public and laboratory compare functions are not byte-identical')
    preserved = [
        entry
        for key, value in spec.items()
        if key.startswith('preserved_') and key.endswith('_scientific_files')
        for entry in value
    ]
    for entry in preserved:
        if digest((destination / safe_relative(entry['path'])).read_bytes()) != entry['sha256']:
            raise ValueError(f'previous scientific bytes changed: {entry["path"]}')
    # Validate the entire transaction before performing its reversible writes.
    created = 0
    if not check:
        for target, data in prepared:
            if not target.exists():
                target.parent.mkdir(parents=True, exist_ok=True)
                with target.open('xb') as stream:
                    stream.write(data)
                created += 1
    return {'source_commit': commit, 'validated_files': len(prepared),
            'created_files': created, 'check_only': check}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--destination', type=Path, default=ROOT)
    parser.add_argument('--spec', type=Path, default=SPEC)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    print(json.dumps(stage(args.source, args.destination,
                           json.loads(args.spec.read_text()), check=args.check), indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
