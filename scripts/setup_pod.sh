#!/usr/bin/env bash
# Idempotent GPU-pod bootstrap: conda envs, repository, SSSP, GPU Quantum ESPRESSO.
#
# Network volumes such as RunPod's /workspace (geesefs) cannot hold
# executables, so conda, environments, the repository and QE live on the
# container disk; only caches, model weights and results persist on the volume.
#
# Usage: bash scripts/setup_pod.sh            (everything)
#        SKIP_QE=1 SKIP_MOOSE=1 bash scripts/setup_pod.sh
set -euo pipefail

PERSIST=${HYDROGEN_PERSIST_ROOT:-/workspace}
CONDA_ROOT=${HYDROGEN_CONDA_ROOT:-/opt/miniforge3}
ENV_ROOT=${HYDROGEN_ENV_ROOT:-/opt/envs}
REPO=${HYDROGEN_REPO:-/root/hydrogen}
QE_ROOT=${HYDROGEN_GPU_QE_ROOT:-/opt/qe-gpu}
SSSP_RECORD=https://archive.materialscloud.org/api/records/rcyfm-68h65/files
SSSP_ARCHIVE=SSSP_1.3.0_PBE_efficiency.tar.gz

export HF_HOME=${HF_HOME:-$PERSIST/hf}
export XDG_CACHE_HOME=${XDG_CACHE_HOME:-$PERSIST/cache}
mkdir -p "$HF_HOME" "$XDG_CACHE_HOME" "$ENV_ROOT"

log() { printf '\n==> %s\n' "$*"; }

log "System tools (nvtop for live GPU monitoring)"
if ! command -v nvtop >/dev/null 2>&1 && command -v apt-get >/dev/null 2>&1; then
    apt-get update -qq && DEBIAN_FRONTEND=noninteractive apt-get install -y -qq nvtop tmux ||
        echo "nvtop install skipped (no apt access)"
fi

log "Miniforge"
if [ ! -x "$CONDA_ROOT/bin/conda" ]; then
    curl -fsSL -o /tmp/miniforge.sh \
        https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
    bash /tmp/miniforge.sh -b -p "$CONDA_ROOT"
fi
# shellcheck disable=SC1091
source "$CONDA_ROOT/bin/activate"

log "Repository"
if [ ! -d "$REPO/.git" ]; then
    git clone https://github.com/YungRaj/hydrogen "$REPO"
fi
cd "$REPO"
git pull --ff-only || echo "git pull skipped (local changes or detached head)"
git log -1 --oneline || true

create_env() {  # name, environment file
    if [ ! -d "$ENV_ROOT/$1" ]; then
        conda env create -q -p "$ENV_ROOT/$1" -f "$2"
    else
        conda env update -q -p "$ENV_ROOT/$1" -f "$2"
    fi
}

log "fairchem-env (screening, reactors, tests)"
create_env fairchem-env environments/environment-fairchem.yml

if [ -z "${SKIP_MOOSE:-}" ]; then
    log "moose-env (NTEC / electrochemical models)"
    create_env moose-env environments/environment-moose.yml
fi

log "SSSP 1.3.0 PBE efficiency pseudopotentials"
SSSP_DIR=$REPO/quantum_espresso/sssp/1.3.0-pbe-efficiency
if ! "$ENV_ROOT/fairchem-env/bin/python" - <<'EOF'
import json, sys
from pipeline.validation.qe_workflows import SSSP_MANIFEST, verify_sssp
sys.exit(0 if verify_sssp(list(json.load(open(SSSP_MANIFEST))))['valid'] else 1)
EOF
then
    curl -fsSL -o "$XDG_CACHE_HOME/$SSSP_ARCHIVE" "$SSSP_RECORD/$SSSP_ARCHIVE/content"
    rm -rf /tmp/sssp && mkdir -p /tmp/sssp
    tar -xzf "$XDG_CACHE_HOME/$SSSP_ARCHIVE" -C /tmp/sssp
    find /tmp/sssp \( -iname '*.upf' \) -exec cp {} "$SSSP_DIR/" \;
fi
"$ENV_ROOT/fairchem-env/bin/python" - <<'EOF'
import json
from pipeline.validation.qe_workflows import SSSP_MANIFEST, verify_sssp
result = verify_sssp(list(json.load(open(SSSP_MANIFEST))))
assert result['valid'], result['errors']
print(f"SSSP verified: {len(result['records'])} elements (md5-pinned)")
EOF

if [ -z "${SKIP_QE:-}" ]; then
    log "GPU Quantum ESPRESSO (NVHPC build for the visible GPU)"
    # The NVHPC archive is extracted with permissions/ownership, which network
    # volumes (geesefs) reject, so the QE build cache stays on the local disk.
    XDG_CACHE_HOME=${HYDROGEN_LOCAL_CACHE:-$HOME/.cache} \
        HYDROGEN_GPU_QE_ROOT=$QE_ROOT bash scripts/install_qe_gpu.sh
fi

log "Smoke checks"
"$ENV_ROOT/fairchem-env/bin/python" -c "import torch, fairchem.core as f; print('torch', torch.__version__, 'fairchem', f.__version__, torch.cuda.get_device_name(0))"
[ -n "${SKIP_MOOSE:-}" ] || "$ENV_ROOT/moose-env/bin/moose-opt" --version | head -1 || true

cat <<EOF

SETUP_DONE. Per-shell environment:
  source $CONDA_ROOT/bin/activate
  export HF_HOME=$HF_HOME
  export MOOSE_EXECUTABLE=$ENV_ROOT/moose-env/bin/moose-opt
  source $QE_ROOT/qe-*/activate.sh     # pw.x, neb.x, MPIEXEC
EOF
