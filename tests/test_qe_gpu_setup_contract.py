#!/usr/bin/env python3
"""Ensure documented production QE setup cannot silently install a CPU build."""

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def test_qe_gpu_installer_is_blackwell_pinned_and_fail_closed():
    script = (ROOT / "scripts" / "install_qe_gpu.sh").read_text()
    for required in (
        "NVHPC_VERSION=25.5",
        "CUDA_VERSION=12.9",
        "CUDA_ARCH=120",
        "QE_VERSION=7.5",
        "--with-cuda-cc=\"$CUDA_ARCH\"",
        "expected Blackwell compute capability 12.0",
        "sha256sum --check --status",
        "-D__CUDA",
        "GPU acceleration is ACTIVE",
        "OPAL_PREFIX",
    ):
        assert required in script


def test_python_requirements_do_not_recommend_cpu_qe_package():
    requirements = (ROOT / "requirements.txt").read_text().lower()
    assert "conda install -c conda-forge qe" not in requirements
    assert "scripts/install_qe_gpu.sh" in requirements


def test_production_qe_execution_requires_gpu_by_default():
    from pipeline.validation.qe_workflows import QEExecutionConfig

    assert QEExecutionConfig().require_gpu is True
    try:
        QEExecutionConfig(require_gpu=False)
    except TypeError:
        pass
    else:
        raise AssertionError('require_gpu must not be caller-configurable')


def test_qe_resolution_rejects_path_and_requires_activation():
    """A generic same-named executable must never become production QE."""
    import os
    import tempfile
    from unittest.mock import patch

    from pipeline.simulation.executables import resolve_qe_executable

    with tempfile.TemporaryDirectory() as tmp:
        fake = Path(tmp) / 'pw.x'
        fake.write_text('#!/bin/sh\nexit 0\n')
        fake.chmod(0o755)
        with patch.dict(os.environ, {'PATH': tmp}, clear=True):
            try:
                resolve_qe_executable('pw.x')
            except RuntimeError as exc:
                assert 'PW_X is not set' in str(exc)
            else:
                raise AssertionError('generic PATH QE must be rejected')
        with patch.dict(os.environ, {'PW_X': str(fake)}, clear=True):
            assert resolve_qe_executable('pw.x') == str(fake.resolve())


def test_executable_resolution_preserves_multicall_symlink_name():
    import tempfile

    from pipeline.simulation.executables import resolve_executable

    with tempfile.TemporaryDirectory() as tmp:
        implementation = Path(tmp) / 'env.sh'
        implementation.write_text('#!/bin/sh\nexit 0\n')
        implementation.chmod(0o755)
        launcher = Path(tmp) / 'mpirun'
        launcher.symlink_to(implementation.name)
        assert resolve_executable(str(launcher)) == str(launcher.absolute())


def main():
    tests = (
        test_qe_gpu_installer_is_blackwell_pinned_and_fail_closed,
        test_python_requirements_do_not_recommend_cpu_qe_package,
        test_production_qe_execution_requires_gpu_by_default,
        test_qe_resolution_rejects_path_and_requires_activation,
        test_executable_resolution_preserves_multicall_symlink_name,
    )
    for test in tests:
        test()
        print("PASS", test.__name__)
    print(f"{len(tests)}/{len(tests)} GPU-QE setup contracts passed")


if __name__ == "__main__":
    main()
