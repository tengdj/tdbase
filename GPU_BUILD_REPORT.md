# GPU 环境与编译验证

本文件记录第一轮 GPU 环境与验证；后续统一入口、面积数值导出和完整集成
结果见 [FARM_INTEGRATION_REPORT.md](FARM_INTEGRATION_REPORT.md)。

验证日期：2026-11-05。使用 WSL `Ubuntu-22.04`，项目目录
`/home/moli/tdfarm/tdbase`；本轮修改已同步到 `D:\tdfarm\tdbase`。

## 环境

| 项目 | 实际环境 |
| --- | --- |
| GPU | NVIDIA GeForce RTX 4080，16376 MiB，计算能力 8.9 |
| Windows NVIDIA 驱动 | 560.94，WSL 可正常识别 GPU |
| CUDA Toolkit | 已安装 12.6，nvcc V12.6.85 |
| CUDA 路径 | `/usr/local/cuda` → CUDA 12.6 |
| 主机编译器 | GCC/G++ 11.4.0 |
| 构建 | CMake + Ninja，Release，CUDA C++17，sm_89 |

已有 Toolkit 和驱动满足编译与运行要求，无需重复安装。新增
`script/env_cuda.sh` 和 `script/build_gpu.sh`；WSL 的 `~/.profile` 已加载
环境脚本，原配置备份为 `~/.profile.before-tdbase-cuda`。原有 `.bashrc`
中的 CUDA 配置保留。没有安装 Linux NVIDIA 显卡驱动。

## 构建结果

GPU 全量构建成功，输出在 `/home/moli/tdfarm/tdbase/build/gpu/`，可执行文件
也复制到 `D:\tdfarm\tdbase\build\gpu\`（均为 Linux ELF，需通过 WSL 运行）。

- 三维：`tdbase`，运行查询时用 `-g` 选择 GPU。
- 二维 GPU：`intersect`、`intersection`、`within_polygon`。
- 二维 CPU：`intersect_cpu`、`intersection_cpu`、`within_polygon_cpu`。
- 数据工具：`dump_polygons`、`filter_polygons`；另有原有 `test` 目标。

修正了 CMake 写死的 `/usr/bin/nvcc` 和 GPU 架构，改用 CUDA Toolkit
发现机制与导入库；`USE_GPU` 改为可显式覆盖的选项，并限定到相应目标，
使 CPU 和 GPU 程序能够同时编译。为 FARM CUDA 内核启用了扩展 lambda
支持，并修正 pinned-memory 释放函数中的错误宏名。

实际运行还发现 GPU `within_polygon` 的非层级模式会访问空层级元数据。
现将平面栅格序列化为单层元数据，修复了 2/4/8 位编码下的非法显存访问。

## 验证结果

- CPU 回归：**53/53 通过**，使用本次启用 GPU 的 Release 构建产物。
- GPU 功能测试：**55/55 通过**，启用 `CUDA_LAUNCH_BLOCKING=1` 同步报告内核错误。
- 二维测试覆盖相交、包含、不相交、交集面积执行、距离阈值、自连接、
  近似模式、2/4/8 位编码，以及平面和层级栅格。
- 三维测试覆盖渐进 LOD 和 LOD100 距离查询；日志确认进入 GPU
  距离计算，GPU 和 CPU 输出的对象对一致。
- `cuobjdump` 确认可执行文件包含 sm_89 内核；动态链接解析到 WSL
  `libcuda.so.1` 和 CUDA 12 的 `libcudart.so.12`。

验证范围：交集面积命令没有输出数值面积，因此只验证 GPU 执行成功，
未完成数值面积精度断言。`filter_polygons` 和原有 `test` 仅完成编译。
未进行大规模数据性能或完整几何精度验证。

额外尝试了 Compute Sanitizer 显存检查，但 Windows WDDM 调试接口未启用，
工具报告 `Failed to initialize WDDM debugger interface` 和 `Device not supported`。
因此本轮**没有完成 Sanitizer 检查**；这与上述正常 GPU 查询成功运行分别记录。
没有为可选检查修改 Windows 调试或注册表设置。

原始构建日志、逐项测试日志和 `results.json` 位于两份项目的
`build/gpu-validation/`。FARM 当前两个固定 scratch buffer 合计占用 8 GiB
显存，测试按顺序执行。

## 重复编译与验证

在 WSL 中运行：

```bash
cd /home/moli/tdfarm/tdbase
bash script/build_gpu.sh
python3 tests/cpu_smoke.py build/gpu . build/gpu-validation/cpu
python3 tests/gpu_smoke.py build/gpu . build/gpu-validation/gpu
```

从 Windows PowerShell 重新编译：

```powershell
wsl -d Ubuntu-22.04 -- bash /home/moli/tdfarm/tdbase/script/build_gpu.sh
```
