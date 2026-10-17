"""Source provenance for preregistered evidence without freezing the tree.

Two separate guarantees replace the old "hash a few files in the working
tree" locks:

* **Historical locks** are verified against the git commit they were made
  at (``docs/evidence/lock_provenance.json``). A completed campaign stays
  provable forever while the working tree keeps evolving.
* **Running campaigns** lock the *full import closure* of their entry point
  at launch (:func:`closure_digest`), so no imported module can change
  unnoticed mid-campaign, and nothing is frozen when no campaign runs.
"""

from __future__ import annotations

import ast
import hashlib
import json
import re
import subprocess
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
LOCK_PROVENANCE = REPO_ROOT / 'docs' / 'evidence' / 'lock_provenance.json'


def _module_path(module: str, root: Path) -> Path | None:
    """Resolve a repo-internal dotted module name to its source file."""
    base = root.joinpath(*module.split('.'))
    for candidate in (base.with_suffix('.py'), base / '__init__.py'):
        if candidate.is_file():
            return candidate
    return None


def _imported_modules(path: Path, root: Path) -> set[str]:
    """Return dotted names imported anywhere in a file (lazy imports too)."""
    tree = ast.parse(path.read_text(), filename=str(path))
    package = '.'.join(path.relative_to(root).with_suffix('').parts[:-1])
    names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                parts = package.split('.') if package else []
                parts = parts[: len(parts) - (node.level - 1)]
                base = '.'.join(parts + ([node.module] if node.module else []))
            else:
                base = node.module or ''
            names.add(base)
            names.update(f'{base}.{alias.name}' for alias in node.names)
    return names


def source_closure(entry_points: list[str | Path], root: str | Path = REPO_ROOT) -> list[Path]:
    """Return every repo-internal source file reachable by import.

    Imports are followed statically (including function-local imports), so
    the closure is a superset of what any code path can execute.

    Args:
        entry_points: Entry-point files, relative to ``root`` or absolute.
        root: Repository root that contains the ``pipeline`` package.

    Returns:
        Sorted list of absolute source paths in the closure.
    """
    root = Path(root).resolve()
    pending = [(root / Path(entry)).resolve() for entry in entry_points]
    seen: set[Path] = set()
    while pending:
        path = pending.pop()
        if path in seen or not path.is_file():
            continue
        seen.add(path)
        for module in _imported_modules(path, root):
            if module.split('.')[0] != 'pipeline':
                continue
            parts = module.split('.')
            # Every package __init__ on the way runs on import, too.
            for depth in range(1, len(parts) + 1):
                resolved = _module_path('.'.join(parts[:depth]), root)
                if resolved is not None and resolved not in seen:
                    pending.append(resolved)
    return sorted(seen)


def closure_digest(entry_points: list[str | Path], root: str | Path = REPO_ROOT) -> str:
    """Return a checkout-independent SHA-256 of an entry point's import closure.

    Args:
        entry_points: Entry-point files, relative to ``root`` or absolute.
        root: Repository root that contains the ``pipeline`` package.

    Returns:
        Hex digest over relative paths and CRLF-normalized contents.
    """
    root = Path(root).resolve()
    digest = hashlib.sha256()
    for path in source_closure(entry_points, root):
        digest.update(path.relative_to(root).as_posix().encode())
        digest.update(b'\0')
        digest.update(path.read_bytes().replace(b'\r\n', b'\n'))
        digest.update(b'\0')
    return digest.hexdigest()


def git_commit(root: str | Path = REPO_ROOT) -> str | None:
    """Return the current commit, or None outside a git checkout.

    Args:
        root: Repository root.

    Returns:
        Full commit hash or None.
    """
    result = subprocess.run(
        ['git', 'rev-parse', 'HEAD'], cwd=root, capture_output=True, text=True
    )
    return result.stdout.strip() if result.returncode == 0 else None


def commit_present(commit: str, root: str | Path = REPO_ROOT) -> bool:
    """Return True when ``commit`` is present in the local git object database.

    Args:
        commit: Commit hash or tag to check.
        root: Repository root.

    Returns:
        True when ``git cat-file -e <commit>^{commit}`` succeeds.
    """
    result = subprocess.run(
        ['git', 'cat-file', '-e', f'{commit}^{{commit}}'],
        cwd=root,
        capture_output=True,
    )
    return result.returncode == 0


def _show(root: Path, commit: str, path: str) -> bytes | None:
    result = subprocess.run(
        ['git', 'show', f'{commit}:{path}'], cwd=root, capture_output=True
    )
    return result.stdout if result.returncode == 0 else None


def _legacy_policy_digest(root: Path, commit: str) -> str | None:
    """Recompute the retired six-file policy digest exactly as of ``commit``."""
    runner = _show(root, commit, 'run_divide_conquer_pilot.py')
    if runner is None:
        return None
    match = re.search(rb'POLICY_SOURCE_PATHS\s*=\s*\((.*?)\n\)', runner, re.S)
    if match is None:
        return None
    digest = hashlib.sha256()
    for raw in re.findall(rb"Path\('([^']+)'\)", match.group(1)):
        relative = raw.decode()
        content = _show(root, commit, relative)
        if content is None:
            return None
        digest.update(relative.encode())
        digest.update(b'\0')
        digest.update(content.replace(b'\r\n', b'\n'))
        digest.update(b'\0')
    return digest.hexdigest()


def verify_lock_at_commit(
    lock_path: str | Path, commit: str, root: str | Path = REPO_ROOT
) -> list[str]:
    """Check a historical lock's source hashes against one git commit.

    Supports the three historical lock shapes: a ``sha256`` file map, a
    single ``implementation_sha256`` for ``molten_neb.py``, and the retired
    six-file ``policy_source_sha256`` digest.

    Args:
        lock_path: Lock JSON document.
        commit: Commit-ish the lock claims to describe.
        root: Repository root (a git checkout).

    Returns:
        Mismatch descriptions; empty when the lock verifies.
    """
    root = Path(root)
    lock = json.loads(Path(lock_path).read_text())
    errors: list[str] = []
    if lock.get('sha256'):
        for relative, expected in lock['sha256'].items():
            content = _show(root, commit, relative)
            if content is None:
                errors.append(f'{relative}: missing at {commit}')
            elif hashlib.sha256(content).hexdigest() != expected:
                errors.append(f'{relative}: changed at {commit}')
    elif lock.get('implementation_sha256'):
        content = _show(root, commit, 'pipeline/evidence/molten_neb.py')
        if content is None or (
            hashlib.sha256(content).hexdigest() != lock['implementation_sha256']
        ):
            errors.append(f'molten_neb.py does not match at {commit}')
    elif lock.get('policy_source_sha256'):
        if _legacy_policy_digest(root, commit) != lock['policy_source_sha256']:
            errors.append(f'policy sources do not match at {commit}')
    else:
        errors.append('lock records no source hashes')
    return errors
