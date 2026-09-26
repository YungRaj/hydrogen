#!/usr/bin/env bash
# Install a pinned Blackwell-capable NVIDIA HPC SDK and GPU Quantum ESPRESSO.

set -euo pipefail

NVHPC_VERSION=25.5
NVHPC_RELEASE=255
CUDA_VERSION=12.9
CUDA_ARCH=120
QE_VERSION=7.5
NVHPC_SHA256=32ed33840bf38448bca346eff5efeb3f80796faa863e9ae0934d5745f2a561be
QE_SHA256=7e1f7a9a21b63192f5135218bee20a5321b66582e4756536681b76e9c59b3cc8

INSTALL_ROOT=${HYDROGEN_GPU_QE_ROOT:-${XDG_DATA_HOME:-$HOME/.local/share}/hydrogen/qe-gpu}
CACHE_ROOT=${XDG_CACHE_HOME:-$HOME/.cache}/hydrogen/qe-gpu
NVHPC_ROOT=$INSTALL_ROOT/nvhpc
QE_PREFIX=$INSTALL_ROOT/qe-$QE_VERSION-nvhpc$NVHPC_VERSION-cc$CUDA_ARCH
JOBS=${QE_BUILD_JOBS:-$(getconf _NPROCESSORS_ONLN)}
NVHPC_ARCHIVE=nvhpc_2025_${NVHPC_RELEASE}_Linux_x86_64_cuda_${CUDA_VERSION}.tar.gz
NVHPC_DIR=nvhpc_2025_${NVHPC_RELEASE}_Linux_x86_64_cuda_${CUDA_VERSION}
NVHPC_URL=https://developer.download.nvidia.com/hpc-sdk/${NVHPC_VERSION}/${NVHPC_ARCHIVE}
QE_ARCHIVE=qe-${QE_VERSION}.tar.gz
QE_URL=https://gitlab.com/QEF/q-e/-/archive/qe-${QE_VERSION}/q-e-qe-${QE_VERSION}.tar.gz

die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
command -v curl >/dev/null || die 'curl is required'
command -v tar >/dev/null || die 'tar is required'
command -v nvidia-smi >/dev/null || die 'an NVIDIA driver is required'
test "$(uname -s)" = Linux || die 'GPU Quantum ESPRESSO setup supports Linux only'
test "$(uname -m)" = x86_64 || die 'this pinned SDK artifact requires x86_64'

mapfile -t capabilities < <(nvidia-smi --query-gpu=compute_cap --format=csv,noheader)
test "${#capabilities[@]}" -gt 0 || die 'no NVIDIA GPU is visible'
for capability in "${capabilities[@]}"; do
    test "$capability" = 12.0 || die "expected Blackwell compute capability 12.0, found $capability"
done

mkdir -p "$CACHE_ROOT" "$INSTALL_ROOT"
download() {
    local url=$1 target=$2 checksum=$3
    if ! printf '%s  %s\n' "$checksum" "$target" | sha256sum --check --status 2>/dev/null; then
        curl --fail --location --continue-at - --output "$target" "$url"
    fi
    printf '%s  %s\n' "$checksum" "$target" | sha256sum --check --status ||
        die "checksum verification failed for $target"
}

download "$NVHPC_URL" "$CACHE_ROOT/$NVHPC_ARCHIVE" "$NVHPC_SHA256"
if test ! -x "$NVHPC_ROOT/Linux_x86_64/$NVHPC_VERSION/compilers/bin/nvfortran"; then
    rm -rf "$CACHE_ROOT/$NVHPC_DIR"
    tar -xzf "$CACHE_ROOT/$NVHPC_ARCHIVE" -C "$CACHE_ROOT"
    NVHPC_SILENT=true NVHPC_INSTALL_TYPE=single NVHPC_INSTALL_DIR="$NVHPC_ROOT" \
        "$CACHE_ROOT/$NVHPC_DIR/install"
fi

SDK=$NVHPC_ROOT/Linux_x86_64/$NVHPC_VERSION
MPIROOT=$SDK/comm_libs/$CUDA_VERSION/hpcx/hpcx-2.22.1/ompi
export PATH="$SDK/comm_libs/mpi/bin:$SDK/compilers/bin:$PATH"
export LD_LIBRARY_PATH="$SDK/comm_libs/mpi/lib:$SDK/compilers/lib:$SDK/cuda/$CUDA_VERSION/lib64:${LD_LIBRARY_PATH:-}"
export OPAL_PREFIX="$MPIROOT"
"$SDK/compilers/bin/nvaccelinfo" | grep -q 'Default Target:.*cc120' ||
    die 'NVIDIA HPC SDK does not detect a cc120 target'

download "$QE_URL" "$CACHE_ROOT/$QE_ARCHIVE" "$QE_SHA256"
SOURCE=$CACHE_ROOT/q-e-qe-$QE_VERSION
rm -rf "$SOURCE"
tar -xzf "$CACHE_ROOT/$QE_ARCHIVE" -C "$CACHE_ROOT"
(
    cd "$SOURCE"
    export FC=nvfortran F90=nvfortran F77=nvfortran CC=nvc MPIF90=mpif90
    ./configure --prefix="$QE_PREFIX" \
        --with-cuda="$SDK/cuda/$CUDA_VERSION" \
        --with-cuda-runtime="$CUDA_VERSION" --with-cuda-cc="$CUDA_ARCH" \
        --with-cuda-mpi=no --enable-openmp --with-scalapack=no
    grep -q -- '-D__CUDA' make.inc || die 'QE configure did not enable CUDA'
    make -j"$JOBS" pw neb
    mkdir -p "$QE_PREFIX/bin"
    install -m 0755 bin/pw.x bin/neb.x "$QE_PREFIX/bin/"
)

cat > "$QE_PREFIX/activate.sh" <<EOF
export PATH="$QE_PREFIX/bin:$SDK/comm_libs/mpi/bin:\$PATH"
export LD_LIBRARY_PATH="$SDK/comm_libs/mpi/lib:$SDK/compilers/lib:$SDK/cuda/$CUDA_VERSION/lib64:\${LD_LIBRARY_PATH:-}"
export OPAL_PREFIX="$MPIROOT"
export PW_X="$QE_PREFIX/bin/pw.x"
export NEB_X="$QE_PREFIX/bin/neb.x"
export MPIEXEC="$MPIROOT/bin/mpirun"
EOF

printf 'GPU QE installed. Before running the pipeline:\n  source %s\n' "$QE_PREFIX/activate.sh"
printf 'Runtime output must contain: GPU acceleration is ACTIVE.\n'
