#!/usr/bin/env python3
"""Contracts that keep maintained source APIs documented and reviewable."""

import ast
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

PLACEHOLDER_PHRASES = (
    "Input value for ", "Value supplying ", "The value produced by ",
    "The resulting ", "Computed result in the representation",
)


def source_paths():
    """Return every maintained implementation file covered by this contract."""
    paths = list((REPO_ROOT / "pipeline").rglob("*.py"))
    paths.extend(REPO_ROOT.glob("*.py"))
    return sorted(paths)


def public_apis(tree):
    """Yield public top-level APIs and public methods, excluding local helpers."""
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)) \
                and not node.name.startswith("_"):
            yield node
            if isinstance(node, ast.ClassDef):
                yield from (
                    method for method in node.body
                    if isinstance(method, (ast.FunctionDef, ast.AsyncFunctionDef))
                    and not method.name.startswith("_"))


def test_every_source_module_explains_its_purpose():
    missing = []
    for path in source_paths():
        if not ast.get_docstring(ast.parse(path.read_text())):
            missing.append(str(path.relative_to(REPO_ROOT)))
    assert not missing, f"source modules missing purpose docstrings: {missing}"


def test_every_public_api_explains_its_contract():
    missing = []
    placeholders = []
    for path in source_paths():
        tree = ast.parse(path.read_text())
        for node in public_apis(tree):
            doc = ast.get_docstring(node)
            identity = f"{path.relative_to(REPO_ROOT)}:{node.lineno}:{node.name}"
            if not doc:
                missing.append(identity)
            elif any(phrase in doc for phrase in PLACEHOLDER_PHRASES):
                placeholders.append(identity)
    assert not missing, f"public APIs missing docstrings: {missing}"
    assert not placeholders, f"public APIs contain placeholder docs: {placeholders}"


def test_public_callable_documents_parameters_and_returns():
    """Require named parameter roles and explicit non-None return semantics."""
    missing_parameters = []
    missing_returns = []
    for path in source_paths():
        tree = ast.parse(path.read_text())
        for node in public_apis(tree):
            if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            doc = ast.get_docstring(node) or ""
            arguments = [
                argument.arg for argument in (
                    node.args.posonlyargs + node.args.args + node.args.kwonlyargs)
                if argument.arg not in {"self", "cls"}
            ]
            identity = f"{path.relative_to(REPO_ROOT)}:{node.lineno}:{node.name}"
            absent = [name for name in arguments if name not in doc]
            if absent:
                missing_parameters.append((identity, absent))
            annotation = ast.unparse(node.returns) if node.returns else ""
            if annotation not in {"", "None", "NoneType"} and \
                    "Returns:" not in doc and "Return\n" not in doc:
                missing_returns.append(identity)
    assert not missing_parameters, missing_parameters
    assert not missing_returns, missing_returns


TESTS = (
    test_every_source_module_explains_its_purpose,
    test_every_public_api_explains_its_contract,
    test_public_callable_documents_parameters_and_returns,
)


def main():
    for test in TESTS:
        test()
        print("PASS", test.__name__)
    print(f"{len(TESTS)}/{len(TESTS)} source-documentation contracts passed")


if __name__ == "__main__":
    main()
