"""Read-only authentication of an explicitly pinned implementation AST graph.

No historical source is imported or executed, and output artifacts are excluded
by nsc_frozen_sources. Current implementation graphs retain their strict owner.
"""
from functools import lru_cache
from hashlib import sha256
from pathlib import Path

@lru_cache(maxsize=4)
def original_source_graph(root, owners, context):
    """Authenticate the original AST dependency graph without executing it."""
    import ast
    import subprocess
    from pathlib import PurePosixPath
    from .nsc_frozen_sources import authenticate_source, source_sha256
    root = Path(root)
    pending, visited, hashes, source_bytes = list(owners), set(), {}, {}
    # Validate the explicit pin before inventorying its immutable source tree.
    source_sha256(root, owners[0], context=context)
    available = set(subprocess.check_output(
        ['git', '-C', str(root), 'ls-tree', '-r', '--name-only', context],
        stderr=subprocess.PIPE, timeout=30,
    ).decode('utf-8').splitlines())

    def read(name):
        if name not in source_bytes:
            expected = source_sha256(root, name, context=context)
            source_bytes[name] = authenticate_source(root, name, expected, context=context)
        return source_bytes[name]

    def add_module(parts):
        if not parts:
            return
        if parts[0] == 'recursive_horizons':
            base = PurePosixPath('src')
        elif len(parts) == 1 and parts[0].startswith(('derive_', 'check_', 'validate_')):
            base = PurePosixPath('scripts')
        else:
            return
        for candidate in (base.joinpath(*parts).with_suffix('.py'),
                          base.joinpath(*parts, '__init__.py')):
            if candidate.as_posix() in available:
                pending.append(candidate.as_posix())
                return

    while pending:
        name = pending.pop()
        if name in visited:
            continue
        visited.add(name)
        raw = read(name)
        path = PurePosixPath(name)
        if path.suffix != '.py':
            hashes[name] = sha256(raw).hexdigest()
            continue
        tree = ast.parse(raw.decode('utf-8'), filename=name)
        if path.name == '__init__.py' and all(
                (isinstance(item, ast.Expr) and isinstance(item.value, ast.Constant)
                 and isinstance(item.value.value, str))
                or (isinstance(item, ast.Assign)
                    and all(isinstance(target, ast.Name) and target.id == '__version__'
                            for target in item.targets)
                    and isinstance(item.value, ast.Constant) and isinstance(item.value.value, str))
                for item in tree.body):
            continue
        hashes[name] = sha256(raw).hexdigest()
        package = path.parent.relative_to('src').parts if name.startswith('src/') else ()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    add_module(tuple(alias.name.split('.')))
            elif isinstance(node, ast.ImportFrom):
                if node.level:
                    if node.level > len(package):
                        continue
                    prefix = package[:len(package)-node.level+1]
                else:
                    prefix = ()
                parts = prefix + tuple(node.module.split('.') if node.module else ())
                add_module(parts)
                for alias in node.names:
                    if alias.name != '*':
                        add_module(parts + tuple(alias.name.split('.')))
    return hashes

