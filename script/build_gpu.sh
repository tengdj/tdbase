#!/usr/bin/env bash
set -euo pipefail
tdbase_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$tdbase_root/script/env_cuda.sh"
cmake -S "$tdbase_root/src" -B "$tdbase_root/build/gpu" -G Ninja \
    -DUSE_GPU=ON -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CUDA_ARCHITECTURES="${CUDA_ARCHITECTURES:-native}"
cmake --build "$tdbase_root/build/gpu" --parallel "${BUILD_JOBS:-2}"
