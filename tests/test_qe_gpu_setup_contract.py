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


def test_neb_gpu_evidence_requires_every_child_image_log():
    import tempfile

    from pipeline.validation.qe_workflows import neb_gpu_accelerated

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        input_path = root / 'candidate.neb.in'
        output_path = root / 'candidate.neb.out'
        input_path.write_text(
            "&PATH num_of_images=3 /\n&CONTROL prefix='candidate' /\n")
        output_path.write_text('NEB parent output without a GPU banner\n')
        for index in range(1, 4):
            child = root / 'tmp' / f'candidate_{index}' / 'PW.out'
            child.parent.mkdir(parents=True)
            child.write_text('GPU acceleration is ACTIVE.\n')
        assert neb_gpu_accelerated(input_path, output_path)
        (root / 'tmp/candidate_2/PW.out').write_text('CPU-only output\n')
        assert not neb_gpu_accelerated(input_path, output_path)


def test_single_rank_qe_still_uses_matching_mpi_launcher():
    import os
    import tempfile
    from unittest.mock import patch

    from pipeline.validation.qe_workflows import (
        QEExecutionConfig, build_qe_command)

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        pw, mpi = root / 'pw.x', root / 'mpirun'
        for executable in (pw, mpi):
            executable.write_text('#!/bin/sh\nexit 0\n')
            executable.chmod(0o755)
        with patch.dict(os.environ, {'PW_X': str(pw), 'MPIEXEC': str(mpi)}):
            command = build_qe_command(
                'pw.x', str(root / 'input.in'), QEExecutionConfig())
    assert command[:3] == [str(mpi), '-np', '1']


def main():
    tests = (
        test_qe_gpu_installer_is_blackwell_pinned_and_fail_closed,
        test_python_requirements_do_not_recommend_cpu_qe_package,
        test_production_qe_execution_requires_gpu_by_default,
        test_qe_resolution_rejects_path_and_requires_activation,
        test_executable_resolution_preserves_multicall_symlink_name,
        test_neb_gpu_evidence_requires_every_child_image_log,
        test_single_rank_qe_still_uses_matching_mpi_launcher,
    )
    for test in tests:
        test()
        print("PASS", test.__name__)
    print(f"{len(tests)}/{len(tests)} GPU-QE setup contracts passed")


if __name__ == "__main__":
    main()
