#!/usr/bin/env python3
"""Local web UI/API for the unified TDBase + FARM executable.

Run from WSL:
    python3 frontend/server.py --build build/gpu --port 8765
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import time
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlparse


ROOT = Path(__file__).resolve().parents[1]
FRONTEND = ROOT / "frontend"
DEFAULT_GPU = ROOT / "build" / "gpu" / "tdbase"
DEFAULT_CPU = ROOT / "build" / "cpu" / "tdbase"
MAX_OUTPUT = 2_000_000


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="TDBase/FARM local web UI")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--build", type=Path, help="directory containing tdbase")
    parser.add_argument("--gpu-build", type=Path, help="GPU build directory (overrides --build for GPU)")
    parser.add_argument("--cpu-build", type=Path, help="CPU build directory (overrides --build for CPU)")
    parser.add_argument("--timeout", type=int, default=600, help="query timeout in seconds")
    return parser.parse_args()


ARGS = parse_args()
BUILD_OVERRIDE = ARGS.build.resolve() if ARGS.build else None
GPU_BUILD_OVERRIDE = ARGS.gpu_build.resolve() if ARGS.gpu_build else None
CPU_BUILD_OVERRIDE = ARGS.cpu_build.resolve() if ARGS.cpu_build else None


def executable(backend: str) -> Path:
    override = GPU_BUILD_OVERRIDE if backend == "gpu" else CPU_BUILD_OVERRIDE
    if override:
        candidate = override / "tdbase"
    elif BUILD_OVERRIDE and backend == "gpu":
        # A directory explicitly named cpu is a CPU-only selection.
        candidate = (BUILD_OVERRIDE / "tdbase.gpu-missing"
                     if "cpu" in BUILD_OVERRIDE.name.lower()
                     else BUILD_OVERRIDE / "tdbase")
    elif BUILD_OVERRIDE and backend == "cpu":
        # A GPU build contains the same CPU-capable unified executable, but
        # prefer the separately built CPU artifact when it exists.
        candidate = (ROOT / "build" / "cpu" / "tdbase") if (ROOT / "build" / "cpu" / "tdbase").is_file() else BUILD_OVERRIDE / "tdbase"
    else:
        candidate = DEFAULT_GPU if backend == "gpu" else DEFAULT_CPU
    return candidate.resolve()


def clean_path(value: object, name: str, required: bool = True) -> str:
    text = str(value or "").strip()
    if required and not text:
        raise ValueError(f"{name} 不能为空")
    if "\x00" in text:
        raise ValueError(f"{name} 包含非法字符")
    return text


def integer(value: object, name: str, minimum: int, maximum: int | None = None, default: int | None = None) -> int:
    if value in (None, ""):
        if default is not None:
            return default
        raise ValueError(f"{name} 不能为空")
    try:
        number = int(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} 必须是整数") from error
    if number < minimum or (maximum is not None and number > maximum):
        limit = f"{minimum}..{maximum}" if maximum is not None else f">={minimum}"
        raise ValueError(f"{name} 必须在 {limit} 范围内")
    return number


def build_command(payload: dict) -> tuple[list[str], Path, dict]:
    dimension = payload.get("dimension", "2d")
    backend = payload.get("backend", "cpu")
    query = payload.get("query", "intersect")
    if dimension not in {"2d", "3d"}:
        raise ValueError("dimension 必须是 2d 或 3d")
    if backend not in {"cpu", "gpu"}:
        raise ValueError("backend 必须是 cpu 或 gpu")
    allowed = {"intersect", "intersection", "within"} if dimension == "2d" else {"intersect", "within", "nn"}
    if query not in allowed:
        raise ValueError(f"{dimension} 不支持 {query} 查询")

    binary = executable(backend)
    if not binary.is_file():
        raise FileNotFoundError(f"找不到 {backend.upper()} 版本: {binary}")

    source = clean_path(payload.get("source"), "源数据")
    # Both engines support self joins. FARM's within query and all TDBase
    # queries accept an omitted target and reuse the source dataset.
    target_required = dimension == "2d" and query != "within"
    target = clean_path(payload.get("target"), "目标数据", target_required)
    command = [str(binary)]

    if dimension == "2d":
        command += ["farm", query, "-s", source]
        if target:
            command += ["-t", target]
        if backend == "gpu":
            command.append("-g")
        if payload.get("hierarchy"):
            command.append("-h")
        if payload.get("approximation"):
            command.append("-a")
        command += ["-n", str(integer(payload.get("threads"), "线程数", 1, 256, 2))]
        command += ["-v", str(integer(payload.get("vpr"), "二维栅格粒度", 1, 10000, 10))]
        command += ["-w", str(integer(payload.get("bitwidth"), "位宽", 2, 8, 4))]
        if query == "within":
            command += ["-d", str(integer(payload.get("distance"), "距离", 0, 10_000_000, 10))]
    else:
        command += ["join", "--dim", "3", "-q", query, "--tile1", source]
        if target:
            command += ["--tile2", target]
        target_object = payload.get("target_object")
        if target_object not in (None, ""):
            command += ["--target", str(integer(target_object, "目标物体 ID", 0, 2_147_483_647))]
        if backend == "gpu":
            command.append("-g")
        command += ["--cn", str(integer(payload.get("compute_threads"), "几何计算线程数", 1, 256, 2))]
        lods = str(payload.get("lods", "100")).strip()
        if not re.fullmatch(r"(?:20|40|60|80|100)(?:\s*,\s*(?:20|40|60|80|100))*", lods):
            raise ValueError("LOD 只能由 20、40、60、80、100 组成，用逗号分隔")
        for lod in dict.fromkeys(part.strip() for part in lods.split(",")):
            command += ["-l", lod]
        if query == "within":
            command += ["-w", str(integer(payload.get("distance"), "距离", 0, 10_000_000, 100))]
        elif query == "nn":
            command += ["-k", str(integer(payload.get("k"), "K", 1, 100_000, 1))]

    if payload.get("print_result", True):
        command.append("-p")
    return command, binary, {"dimension": dimension, "backend": backend, "query": query}


def run_query(payload: dict) -> dict:
    command, binary, meta = build_command(payload)
    started = time.monotonic()
    try:
        process = subprocess.run(
            command,
            cwd=ROOT,
            capture_output=True,
            text=True,
            errors="replace",
            timeout=ARGS.timeout,
            env=os.environ.copy(),
        )
    except subprocess.TimeoutExpired as error:
        raise TimeoutError(f"查询超过 {ARGS.timeout} 秒，已终止") from error
    elapsed = round(time.monotonic() - started, 3)
    stdout = process.stdout[-MAX_OUTPUT:]
    stderr = process.stderr[-MAX_OUTPUT:]
    count_match = re.findall(r"count-found:\s*(\d+)", stdout + stderr)
    return {
        **meta,
        "ok": process.returncode == 0,
        "exit_code": process.returncode,
        "duration_seconds": elapsed,
        "command": command,
        "binary": str(binary),
        "count": int(count_match[-1]) if count_match else None,
        "stdout": stdout,
        "stderr": stderr,
    }


class Handler(BaseHTTPRequestHandler):
    server_version = "TDBaseWeb/1.0"

    def log_message(self, format: str, *args: object) -> None:
        print(f"[{self.log_date_time_string()}] {format % args}")

    def send_json(self, value: dict, status: int = HTTPStatus.OK) -> None:
        data = json.dumps(value, ensure_ascii=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(data)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self) -> None:
        parsed = urlparse(self.path)
        if parsed.path == "/api/health":
            self.send_json({
                "ok": True,
                "gpu_binary": str(executable("gpu")),
                "gpu_available": executable("gpu").is_file(),
                "cpu_binary": str(executable("cpu")),
                "cpu_available": executable("cpu").is_file(),
            })
            return
        if parsed.path in {"/", "/index.html"}:
            self.serve_file(FRONTEND / "index.html", "text/html; charset=utf-8")
            return
        if parsed.path in {"/app.js", "/styles.css"}:
            mime = "text/javascript; charset=utf-8" if parsed.path.endswith("js") else "text/css; charset=utf-8"
            self.serve_file(FRONTEND / parsed.path.lstrip("/"), mime)
            return
        self.send_json({"error": "Not found"}, HTTPStatus.NOT_FOUND)

    def serve_file(self, path: Path, mime: str) -> None:
        try:
            data = path.read_bytes()
        except FileNotFoundError:
            self.send_json({"error": "Frontend file not found"}, HTTPStatus.NOT_FOUND)
            return
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", mime)
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def do_POST(self) -> None:
        if urlparse(self.path).path != "/api/query":
            self.send_json({"error": "Not found"}, HTTPStatus.NOT_FOUND)
            return
        try:
            size = int(self.headers.get("Content-Length", "0"))
            if size > 100_000:
                raise ValueError("请求过大")
            payload = json.loads(self.rfile.read(size).decode("utf-8"))
            result = run_query(payload)
            self.send_json(result, HTTPStatus.OK if result["ok"] else HTTPStatus.BAD_REQUEST)
        except (ValueError, FileNotFoundError, TimeoutError, json.JSONDecodeError) as error:
            self.send_json({"ok": False, "error": str(error)}, HTTPStatus.BAD_REQUEST)
        except Exception as error:  # Keep the local UI responsive on process errors.
            self.send_json({"ok": False, "error": f"服务器错误: {error}"}, HTTPStatus.INTERNAL_SERVER_ERROR)


def main() -> None:
    print(f"TDBase Web UI: http://{ARGS.host}:{ARGS.port}")
    print(f"Project root: {ROOT}")
    if BUILD_OVERRIDE:
        print(f"Build override: {BUILD_OVERRIDE}")
    if GPU_BUILD_OVERRIDE:
        print(f"GPU build override: {GPU_BUILD_OVERRIDE}")
    if CPU_BUILD_OVERRIDE:
        print(f"CPU build override: {CPU_BUILD_OVERRIDE}")
    server = ThreadingHTTPServer((ARGS.host, ARGS.port), Handler)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopping TDBase Web UI")
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
