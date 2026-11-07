# CPU 编译与冒烟测试（2026-11-05）

本文件记录第一轮 CPU 验证；后续统一入口与面积数值验证见
[FARM_INTEGRATION_REPORT.md](FARM_INTEGRATION_REPORT.md)。

本次在 WSL 中清除 `USE_GPU` 环境变量，以 GCC 11.4、CMake/Ninja、
CGAL 5.6.1、Debug 配置编译，保留断言。未编译 CUDA 源文件。

## 结果

- 全工程编译通过：`tdbase`、`test`、`intersect_cpu`、`intersection_cpu`、
  `within_polygon_cpu`、`dump_polygons`、`filter_polygons`。
- `tests/cpu_smoke.py`：53/53 项通过。
- 二维：WKT 转二进制、相交/分离/包含、1 km 与 2 km 距离阈值、
  自连接、2/4/8 bit 位宽、开启/关闭层级栅格、多线程及近似模式。
- 交集面积命令验证正常执行；原命令不输出面积数值，本次未做面积数值精度验收。
- 三维：使用仓库自带 nuclei/vessel 数据，各取首个对象，验证 KNN，
  以及渐进式和 LOD 100 距离查询。渐进式查询实际经过 LOD 20/40 解码、
  数据打包及 CPU 几何计算，在 LOD 40 提前确定结果；两种方式均返回 `(0, 0)`。
- `filter_polygons` 和原 `test` 为编译检查，未运行其原有依赖特定输入的逻辑。

## 本次修复

1. `src/include/himesh.h` 按 CGAL 版本选择 5.x/6.x 的 AABB 头文件和类型名称。
2. `src/farm2d/geometry/farm.cpp` 的距离查询在未开启 `-h` 时使用基础栅格，
   避免访问空的层级数组及层级数量减一时的下溢。相关平面/层级模式都已回归。

## 文件位置

- D 盘与 WSL 的上述源代码修改已同步。
- 已生成的 Linux 可执行文件：`/home/moli/tdfarm/tdbase/build/cpu/`。
- 本次 CMake 构建目录：`/tmp/tdfarm-tdbase-cpu-20261105/`。
- 测试结果和逐项日志：`D:\tdfarm\tdbase\build\cpu-validation\`，
  汇总为 `results.json`。

## 复现

在 WSL 中执行；以下命令会在独立目录构建 CPU 版本：

```bash
cd /home/moli/tdfarm/tdbase
env -u USE_GPU cmake -S src -B build/cpu-cmake -G Ninja -DCMAKE_BUILD_TYPE=Debug
cmake --build build/cpu-cmake -j2
python3 tests/cpu_smoke.py build/cpu-cmake .
```

测试仅依赖 Python 标准库。运行已生成的程序可直接使用：

```bash
./build/cpu/intersect_cpu --help
./build/cpu/tdbase
```
