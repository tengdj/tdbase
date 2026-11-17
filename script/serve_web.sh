#!/usr/bin/env bash
# 在服务器上启动局域网查询页面；同一项目只运行一个服务实例。
set -euo pipefail

tdbase_web_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$tdbase_web_root"
mkdir -p build/deployment

# Keep the LAN query service on the second physical GPU by default.  CUDA
# renumbers the selected device as logical device 0 inside this process; an
# explicit CUDA_VISIBLE_DEVICES value can still override it for maintenance.
export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-1}"

exec 9>build/deployment/frontend.lock
if ! flock -n 9; then
    echo "TDBase 网页服务已在运行。"
    exit 0
fi

echo "$$" > build/deployment/frontend.pid
exec python3 -u frontend/server.py \
    --gpu-build build/gpu --cpu-build build/cpu \
    --host "${TDBASE_WEB_HOST:-0.0.0.0}" \
    --port "${TDBASE_WEB_PORT:-8765}"
