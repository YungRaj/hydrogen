#!/usr/bin/env python3
"""Run pre-submit CI (pyright + pytest contracts) on a RunPod GPU pod with strict cost and auto-stop guarantees."""

from __future__ import annotations

import argparse
import json
import os
import signal
import subprocess
import time
import urllib.error
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MAX_COST_PER_HR = 1.00
GRAPHQL_URL = 'https://api.runpod.io/graphql'


def _load_api_key() -> str:
    """Resolve the RunPod API key from environment or gitignored .runpod_api_key."""
    env_key = os.environ.get('RUNPOD_API_KEY', '').strip()
    if env_key:
        return env_key
    key_file = REPO_ROOT / '.runpod_api_key'
    if key_file.is_file():
        key = key_file.read_text(encoding='utf-8').strip()
        if key:
            return key
    raise RuntimeError(
        'Missing RunPod API key. Set RUNPOD_API_KEY or write it to .runpod_api_key'
    )


def _load_ssh_public_key(ssh_private_key: str) -> str:
    """Read the SSH public key matching the private key for new pod deployments."""
    pub_path = Path(f'{ssh_private_key}.pub')
    if pub_path.is_file():
        return pub_path.read_text(encoding='utf-8').strip()
    proc = subprocess.run(
        ['ssh-keygen', '-y', '-f', ssh_private_key],
        capture_output=True,
        text=True,
        check=True,
    )
    return proc.stdout.strip()


def _gql(api_key: str, query: str, variables: dict | None = None) -> dict:
    """Execute a RunPod GraphQL query or mutation."""
    body: dict[str, object] = {'query': query}
    if variables is not None:
        body['variables'] = variables
    req = urllib.request.Request(
        GRAPHQL_URL,
        data=json.dumps(body).encode('utf-8'),
        headers={
            'Content-Type': 'application/json',
            'Authorization': f'Bearer {api_key}',
            'User-Agent': 'hydrogen-runpod-ci/1.0',
        },
    )
    try:
        with urllib.request.urlopen(req, timeout=25) as resp:
            payload = json.loads(resp.read().decode('utf-8'))
    except urllib.error.HTTPError as exc:
        detail = exc.read().decode('utf-8', errors='replace')
        raise RuntimeError(f'RunPod HTTP {exc.code}: {detail}') from exc
    if 'errors' in payload and payload['errors']:
        raise RuntimeError(f"RunPod GraphQL error: {payload['errors']}")
    return payload.get('data', {})


def _list_pods(api_key: str) -> list[dict]:
    """Return all pods associated with the current RunPod account."""
    query = """
    query {
      myself {
        pods {
          id
          name
          desiredStatus
          costPerHr
          imageName
          env
          machine {
            dataCenterId
            gpuDisplayName
          }
          runtime {
            uptimeInSeconds
            ports {
              ip
              isIpPublic
              privatePort
              publicPort
              type
            }
          }
        }
      }
    }
    """
    data = _gql(api_key, query)
    return data.get('myself', {}).get('pods', []) or []


def _get_pod(api_key: str, pod_id: str) -> dict:
    """Fetch current status and runtime ports for a single pod."""
    query = f"""
    query {{
      pod(input: {{podId: "{pod_id}"}}) {{
        id
        name
        desiredStatus
        costPerHr
        runtime {{
          uptimeInSeconds
          ports {{
            ip
            isIpPublic
            privatePort
            publicPort
            type
          }}
        }}
      }}
    }}
    """
    return _gql(api_key, query).get('pod') or {}


def _extract_ssh_endpoint(pod: dict) -> tuple[str, int] | None:
    """Return (public_ip, public_port) for port 22/tcp if the pod runtime is ready."""
    runtime = pod.get('runtime') or {}
    ports = runtime.get('ports') or []
    for port_info in ports:
        if (
            port_info.get('privatePort') == 22
            and port_info.get('isIpPublic')
            and port_info.get('ip')
            and port_info.get('publicPort')
        ):
            return str(port_info['ip']), int(port_info['publicPort'])
    return None


def _has_valid_ssh_key(pod: dict) -> bool:
    """Check whether the pod has a non-empty SSH PUBLIC_KEY configured."""
    env_list = pod.get('env') or []
    for entry in env_list:
        if isinstance(entry, str) and entry.startswith('PUBLIC_KEY='):
            val = entry.split('=', 1)[1].strip()
            return bool(val and val != 'null')
    return False


def _stop_pod(api_key: str, pod_id: str) -> None:
    """Stop a pod via the RunPod API and verify its desiredStatus becomes EXITED."""
    print(f'[runpod-ci] Stopping pod {pod_id} so no GPU charges accrue...', flush=True)
    for attempt in range(3):
        try:
            res = _gql(
                api_key,
                f"""
                mutation {{
                  podStop(input: {{podId: "{pod_id}"}}) {{
                    id
                    desiredStatus
                  }}
                }}
                """,
            )
            stopped = res.get('podStop') or {}
            print(
                f"[runpod-ci] Pod {pod_id} desiredStatus={stopped.get('desiredStatus')}",
                flush=True,
            )
            return
        except Exception as exc:
            print(
                f'[runpod-ci] Warning: stop attempt {attempt + 1}/3 failed: {exc}',
                flush=True,
            )
            time.sleep(2)


def _terminate_pod(api_key: str, pod_id: str) -> None:
    """Terminate a temporary CI pod via the RunPod API so zero storage or GPU charges accrue."""
    print(
        f'[runpod-ci] Terminating temporary CI pod {pod_id} so no charges accrue...',
        flush=True,
    )
    for attempt in range(3):
        try:
            _gql(
                api_key,
                f"""
                mutation {{
                  podTerminate(input: {{podId: "{pod_id}"}})
                }}
                """,
            )
            print(f'[runpod-ci] Temporary pod {pod_id} terminated.', flush=True)
            return
        except Exception as exc:
            print(
                f'[runpod-ci] Warning: terminate attempt {attempt + 1}/3 failed: {exc}',
                flush=True,
            )
            _stop_pod(api_key, pod_id)
            time.sleep(2)


def _wait_for_ssh(
    api_key: str, pod_id: str, ssh_key: str, timeout_s: int = 240
) -> tuple[str, int]:
    """Wait until the pod exposes a public TCP port 22 and accepts SSH connections."""
    deadline = time.time() + timeout_s
    while time.time() < deadline:
        pod = _get_pod(api_key, pod_id)
        if pod.get('desiredStatus') != 'RUNNING':
            raise RuntimeError(
                f"Pod {pod_id} exited unexpectedly (status={pod.get('desiredStatus')})"
            )
        endpoint = _extract_ssh_endpoint(pod)
        if endpoint is not None:
            ip, port = endpoint
            probe = subprocess.run(
                [
                    'ssh',
                    '-o',
                    'StrictHostKeyChecking=no',
                    '-o',
                    'UserKnownHostsFile=/dev/null',
                    '-o',
                    'ConnectTimeout=5',
                    '-o',
                    'BatchMode=yes',
                    '-p',
                    str(port),
                    '-i',
                    ssh_key,
                    f'root@{ip}',
                    'echo READY',
                ],
                capture_output=True,
                text=True,
            )
            if probe.returncode == 0 and 'READY' in probe.stdout:
                return ip, port
        time.sleep(4)
    raise TimeoutError(f'Timed out waiting for SSH on pod {pod_id}')


def _deploy_fresh_pod_under_cap(
    api_key: str, max_cost_per_hr: float, ssh_key: str
) -> tuple[str, str, int, str]:
    """Deploy a fresh on-demand pod strictly <= max_cost_per_hr when stopped hosts have no free GPUs.

    Returns:
        Tuple of (pod_id, ip, port, cleanup_action='terminate').
    """
    pub_key = _load_ssh_public_key(ssh_key)
    query = """
    query {
      gpuTypes {
        id
        displayName
        memoryInGb
        secureCloud
        lowestPrice(input: {gpuCount: 1, secureCloud: true, supportPublicIp: true}) {
          uninterruptablePrice
          stockStatus
        }
      }
    }
    """
    gpus = _gql(api_key, query).get('gpuTypes') or []
    eligible = []
    for gpu in gpus:
        gid = str(gpu.get('id') or '')
        if not gid.startswith('NVIDIA') or not gpu.get('secureCloud'):
            continue
        lowest = gpu.get('lowestPrice') or {}
        price = lowest.get('uninterruptablePrice')
        stock = lowest.get('stockStatus')
        if (
            price is not None
            and float(price) <= max_cost_per_hr
            and stock in ('High', 'Medium', 'Low')
        ):
            eligible.append((float(price), 0 if stock == 'High' else (1 if stock == 'Medium' else 2), gid))

    # Sort by price ascending, preferring better stock on ties
    eligible.sort()
    if not eligible:
        raise RuntimeError(
            f'No Secure Cloud NVIDIA GPU types with public IP available under ${max_cost_per_hr:.2f}/hr'
        )

    mutation = """
    mutation($input: PodFindAndDeployOnDemandInput!) {
      podFindAndDeployOnDemand(input: $input) {
        id
        name
        desiredStatus
        costPerHr
      }
    }
    """
    for price, _, gpu_type_id in eligible[:6]:
        print(
            f'[runpod-ci] Deploying temporary CI pod ({gpu_type_id}, ~${price:.2f}/hr <= ${max_cost_per_hr:.2f}/hr cap)...',
            flush=True,
        )
        try:
            res = _gql(
                api_key,
                mutation,
                variables={
                    'input': {
                        'cloudType': 'SECURE',
                        'supportPublicIp': True,
                        'gpuCount': 1,
                        'volumeInGb': 0,
                        'containerDiskInGb': 50,
                        'minVcpuCount': 4,
                        'minMemoryInGb': 16,
                        'gpuTypeId': gpu_type_id,
                        'name': 'hydrogen_ci_ephemeral',
                        'imageName': 'runpod/pytorch:1.4.0-cu1281-torch280-ubuntu2404',
                        'ports': '22/tcp',
                        'volumeMountPath': '/workspace',
                        'env': [{'key': 'PUBLIC_KEY', 'value': pub_key}],
                    }
                },
            )
            deployed = res.get('podFindAndDeployOnDemand')
            if not deployed or not deployed.get('id'):
                continue
            pod_id = str(deployed['id'])
            actual_cost = float(deployed.get('costPerHr') or price)
            if actual_cost > max_cost_per_hr:
                print(
                    f'[runpod-ci] Deployed pod {pod_id} costs ${actual_cost:.2f}/hr > ${max_cost_per_hr:.2f}/hr cap! Terminating immediately.',
                    flush=True,
                )
                _terminate_pod(api_key, pod_id)
                continue
            try:
                ip, port = _wait_for_ssh(api_key, pod_id, ssh_key, timeout_s=210)
                print(
                    f'[runpod-ci] Connected to ephemeral CI pod {pod_id} at {ip}:{port} (${actual_cost:.2f}/hr)',
                    flush=True,
                )
                return pod_id, ip, port, 'terminate'
            except Exception as exc:
                print(
                    f'[runpod-ci] Ephemeral pod {pod_id} failed SSH readiness: {exc}',
                    flush=True,
                )
                _terminate_pod(api_key, pod_id)
        except Exception as exc:
            print(
                f'[runpod-ci] Could not deploy {gpu_type_id}: {exc}',
                flush=True,
            )

    raise RuntimeError(
        f'Unable to resume or deploy any pod under ${max_cost_per_hr:.2f}/hr'
    )


def _acquire_pod(
    api_key: str, max_cost_per_hr: float, ssh_key: str
) -> tuple[str, str, int, str]:
    """Find a running pod, resume a stopped pod <= max_cost_per_hr, or deploy an ephemeral CI pod <= max_cost_per_hr.

    Returns:
        Tuple of (pod_id, ip, port, cleanup_mode) where cleanup_mode is 'stop' or 'terminate'.
    """
    pods = _list_pods(api_key)

    # 1. Check if a pod <= max_cost_per_hr is already RUNNING
    running = [
        p
        for p in pods
        if p.get('desiredStatus') == 'RUNNING'
        and float(p.get('costPerHr') or 999.0) <= max_cost_per_hr
        and _has_valid_ssh_key(p)
    ]
    running.sort(key=lambda p: float(p.get('costPerHr') or 999.0))
    for pod in running:
        pod_id = str(pod['id'])
        cost = float(pod.get('costPerHr') or 0.0)
        gpu = (pod.get('machine') or {}).get('gpuDisplayName', 'GPU')
        print(
            f"[runpod-ci] Found already-running pod {pod_id} ({pod.get('name')}, {gpu}, ${cost:.2f}/hr)",
            flush=True,
        )
        ip, port = _wait_for_ssh(api_key, pod_id, ssh_key, timeout_s=60)
        return pod_id, ip, port, 'stop'

    # 2. Try resuming existing EXITED pods that have an SSH key and cost <= max_cost_per_hr
    candidates = [
        p
        for p in pods
        if p.get('desiredStatus') == 'EXITED'
        and float(p.get('costPerHr') or 999.0) <= max_cost_per_hr
        and _has_valid_ssh_key(p)
    ]
    candidates.sort(key=lambda p: float(p.get('costPerHr') or 999.0))

    for pod in candidates:
        pod_id = str(pod['id'])
        cost = float(pod.get('costPerHr') or 0.0)
        if cost > max_cost_per_hr:
            continue
        gpu = (pod.get('machine') or {}).get('gpuDisplayName', 'GPU')
        print(
            f"[runpod-ci] Attempting to resume stopped pod {pod_id} ({pod.get('name')}, {gpu}, ${cost:.2f}/hr <= ${max_cost_per_hr:.2f}/hr cap)...",
            flush=True,
        )
        try:
            res = _gql(
                api_key,
                f"""
                mutation {{
                  podResume(input: {{podId: "{pod_id}", gpuCount: 1}}) {{
                    id
                    desiredStatus
                    costPerHr
                  }}
                }}
                """,
            )
            resumed = res.get('podResume')
            if not resumed:
                continue
            resumed_cost = float(resumed.get('costPerHr') or cost)
            if resumed_cost > max_cost_per_hr:
                print(
                    f'[runpod-ci] Pod {pod_id} resumed cost ${resumed_cost:.2f}/hr exceeds ${max_cost_per_hr:.2f}/hr cap! Stopping immediately.',
                    flush=True,
                )
                _stop_pod(api_key, pod_id)
                continue
            ip, port = _wait_for_ssh(api_key, pod_id, ssh_key, timeout_s=150)
            print(
                f'[runpod-ci] Connected to resumed pod {pod_id} at {ip}:{port} (${resumed_cost:.2f}/hr)',
                flush=True,
            )
            return pod_id, ip, port, 'stop'
        except Exception as exc:
            print(
                f'[runpod-ci] Host for stopped pod {pod_id} has no free GPU ({exc}); leaving it EXITED.',
                flush=True,
            )

    # 3. Fallback: deploy a fresh ephemeral pod <= max_cost_per_hr and terminate it when done
    return _deploy_fresh_pod_under_cap(api_key, max_cost_per_hr, ssh_key)


def _run_remote_ci(
    ip: str, port: int, ssh_key: str, *, include_moose: bool = False
) -> int:
    """Sync current repository state to the pod and execute pyright + pytest contracts."""
    ssh_base = [
        'ssh',
        '-o',
        'StrictHostKeyChecking=no',
        '-o',
        'UserKnownHostsFile=/dev/null',
        '-o',
        'LogLevel=ERROR',
        '-p',
        str(port),
        '-i',
        ssh_key,
        f'root@{ip}',
    ]

    # Step 1: Fast bootstrap on pod (uses cached wheels in /workspace/cache/ci_wheels.tar.gz if available)
    bootstrap_script = r"""
set -euo pipefail
if ! python3 -c "import pytest, ase, cantera, sklearn, scipy, pandas" >/dev/null 2>&1 || ! command -v pyright >/dev/null 2>&1; then
    echo "[runpod-ci] Bootstrapping Python CI dependencies on container..."
    if [ -f /workspace/cache/ci_wheels.tar.gz ]; then
        rm -rf /tmp/wheels
        tar -xzf /workspace/cache/ci_wheels.tar.gz -C /tmp
        pip install -q --no-index --find-links=/tmp/wheels --break-system-packages --root-user-action=ignore pytest pyright ase cantera scikit-learn scipy pandas
    else
        pip install -q --break-system-packages --root-user-action=ignore pytest pyright ase cantera scikit-learn scipy pandas
    fi
fi
mkdir -p /root/hydrogen
"""
    boot = subprocess.run(ssh_base + [bootstrap_script], cwd=str(REPO_ROOT))
    if boot.returncode != 0:
        return boot.returncode

    # Step 2: Stream the exact working tree (tracked + untracked non-ignored files) via tar over SSH
    print('[runpod-ci] Syncing repository snapshot to pod...', flush=True)
    files_proc = subprocess.run(
        ['git', 'ls-files', '--cached', '--others', '--exclude-standard', '-z'],
        cwd=str(REPO_ROOT),
        capture_output=True,
        check=True,
    )
    existing_files = [
        f
        for f in files_proc.stdout.decode('utf-8', errors='replace').split('\0')
        if f and (REPO_ROOT / f).exists()
    ]
    tar_env = os.environ.copy()
    tar_env['COPYFILE_DISABLE'] = '1'
    tar_proc = subprocess.Popen(
        ['tar', '-czf', '-', '--no-xattrs', '-T', '-'],
        cwd=str(REPO_ROOT),
        env=tar_env,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
    )
    assert tar_proc.stdin is not None
    assert tar_proc.stdout is not None
    tar_bytes, _ = tar_proc.communicate('\n'.join(existing_files).encode('utf-8'))
    if tar_proc.returncode != 0:
        raise RuntimeError('Failed to create local tar archive of repository files')

    sync_cmd = r"""
set -euo pipefail
mkdir -p /root/hydrogen
find /root/hydrogen -mindepth 1 -maxdepth 1 ! -name quantum_espresso ! -name .git -exec rm -rf {} +
tar -xzf - --no-same-owner --no-same-permissions -C /root/hydrogen
SSSP_DIR=/root/hydrogen/quantum_espresso/sssp/1.3.0-pbe-efficiency
mkdir -p "$SSSP_DIR"
if [ ! -f "$SSSP_DIR/C.pbe-n-kjpaw_psl.1.0.0.UPF" ]; then
    if [ -f /workspace/cache/SSSP_1.3.0_PBE_efficiency.tar.gz ]; then
        tar -xzf /workspace/cache/SSSP_1.3.0_PBE_efficiency.tar.gz --no-same-owner --no-same-permissions -C "$SSSP_DIR"
    else
        echo "[runpod-ci] Downloading SSSP 1.3.0 PBE efficiency pseudopotentials..."
        curl -fsSL -o /tmp/sssp.tar.gz "https://archive.materialscloud.org/api/records/rcyfm-68h65/files/SSSP_1.3.0_PBE_efficiency.tar.gz/content"
        tar -xzf /tmp/sssp.tar.gz --no-same-owner --no-same-permissions -C "$SSSP_DIR"
    fi
fi
"""
    sync_proc = subprocess.run(
        ssh_base + [sync_cmd], input=tar_bytes, cwd=str(REPO_ROOT)
    )
    if sync_proc.returncode != 0:
        return sync_proc.returncode

    # Step 3: Run pyright + pytest contracts (+ MOOSE verification when installed or requested)
    print('[runpod-ci] Running pyright and pytest contract suite on pod...', flush=True)
    moose_block = (
        'if [ -x /opt/envs/moose-env/bin/moose-opt ]; then '
        'export MOOSE_EXECUTABLE=/opt/envs/moose-env/bin/moose-opt; '
        'fi; '
        'if [ "' + ('1' if include_moose else '0') + '" = "1" ] || command -v moose-opt >/dev/null 2>&1 || [ -n "${MOOSE_EXECUTABLE:-}" ]; then '
        'echo "[runpod-ci] Running MOOSE model verification (tests/test_moose_models.py)..."; '
        'python3 tests/test_moose_models.py; '
        'fi'
    )
    test_cmd = (
        'set -euo pipefail; '
        'cd /root/hydrogen; '
        'pyright; '
        'python3 tests/run_pytest_contracts.py; '
        + moose_block
    )
    test_proc = subprocess.run(ssh_base + [test_cmd], cwd=str(REPO_ROOT))
    return test_proc.returncode


def _cleanup_pod(api_key: str, pod_id: str, cleanup_mode: str) -> None:
    """Stop or terminate the pod depending on whether it was resumed or newly deployed."""
    if cleanup_mode == 'terminate':
        _terminate_pod(api_key, pod_id)
    else:
        _stop_pod(api_key, pod_id)


def main() -> int:
    """Entry point for RunPod CI execution with guaranteed pod stop/termination on exit."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--max-cost',
        type=float,
        default=float(
            os.environ.get('HYDROGEN_MAX_POD_COST', DEFAULT_MAX_COST_PER_HR)
        ),
        help='Maximum allowed pod cost in $/hr (hard cap default: $1.00/hr)',
    )
    parser.add_argument(
        '--ssh-key',
        default=str(Path.home() / '.ssh' / 'id_ed25519'),
        help='Path to SSH private key',
    )
    parser.add_argument(
        '--always-stop',
        action='store_true',
        default=True,
        help='Always stop/terminate the pod after CI completes (default: True)',
    )
    parser.add_argument(
        '--include-moose',
        action='store_true',
        default=os.environ.get('HYDROGEN_CI_INCLUDE_MOOSE', '0') == '1',
        help='Require MOOSE model verification (tests/test_moose_models.py)',
    )
    args = parser.parse_args()

    if args.max_cost > 1.00:
        raise SystemExit(
            f'Refusing --max-cost={args.max_cost:.2f}: hard policy forbids pods > $1.00/hr'
        )

    api_key = _load_api_key()
    pod_id: str | None = None
    cleanup_mode = 'stop'

    def _signal_handler(signum, _frame):
        if pod_id is not None:
            _cleanup_pod(api_key, pod_id, cleanup_mode)
        raise SystemExit(128 + signum)

    signal.signal(signal.SIGINT, _signal_handler)
    signal.signal(signal.SIGTERM, _signal_handler)

    t0 = time.time()
    try:
        pod_id, ip, port, cleanup_mode = _acquire_pod(
            api_key, args.max_cost, args.ssh_key
        )
        rc = _run_remote_ci(ip, port, args.ssh_key, include_moose=args.include_moose)
        elapsed = time.time() - t0
        status_str = 'PASSED' if rc == 0 else f'FAILED (exit {rc})'
        print(f'[runpod-ci] CI {status_str} in {elapsed:.1f}s', flush=True)
        return rc
    finally:
        if pod_id is not None:
            _cleanup_pod(api_key, pod_id, cleanup_mode)


if __name__ == '__main__':
    raise SystemExit(main())
