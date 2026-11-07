# 在 TDBase 中使用 FARM

在 WSL 中执行以下命令。主程序已经直接链接二维 FARM 模块，复制一个
`tdbase` 可执行文件即可使用二维、三维命令，不需要旁边放置 FARM 辅助程序。
系统仍需满足该可执行文件的动态库依赖。

```bash
cd /home/moli/tdfarm/tdbase
source script/env_cuda.sh
./build/gpu/tdbase farm --help
```

## 数据转换与过滤

输入是每行一个闭合、无洞的二维 `POLYGON` WKT。例如：

```text
POLYGON ((0 0,1 0,1 1,0 1,0 0))
```

```bash
./build/gpu/tdbase farm convert source.wkt source.dt
./build/gpu/tdbase farm convert target.wkt target.dt
./build/gpu/tdbase farm filter input.wkt filtered.dt
```

`filter` 延续 FARM 原来的经纬度过滤规则：保留所有顶点都满足
`-180 < x < 180` 且 `-90 < y < 90` 的多边形，并输出二进制数据。
空输入或全部被过滤时也会生成有效的空数据集。

FARM 的 `.dt` 是二维专用格式，与三维 Tile `.dt` 不通用；统一程序会校验
输入格式。`MULTIPOLYGON` 和带内环的多边形尚未实现正确的 FARM 内核处理，
转换入口会明确拒绝，避免原始解析器静默丢失内环。

## 二维查询

统一命令默认使用 CPU，加 `-g` 使用 GPU：

```bash
# 是否相交，输出匹配对象对
./build/gpu/tdbase farm intersect -s source.dt -t target.dt -h -p

# GPU 相交
./build/gpu/tdbase farm intersect -s source.dt -t target.dt -h -g -p

# 交集面积；每行 source_id、target_id、area
./build/gpu/tdbase farm intersection -s source.dt -t target.dt -h -g -p

# 距离不超过 2 公里
./build/gpu/tdbase farm within -s source.dt -t target.dt -d 2 -h -g -p

# 自连接，排除对象与自身，保留两个方向的匹配
./build/gpu/tdbase farm within -s source.dt -d 2 -h -g -p

# 近似查询
./build/gpu/tdbase farm intersect -s source.dt -t target.dt -h -a -g -p

# 输出到 TSV 或 JSON 文件
./build/gpu/tdbase farm intersection -s source.dt -t target.dt -o areas.tsv
./build/gpu/tdbase farm intersect -s source.dt -t target.dt --format json -o pairs.json
```

相交和距离查询的 TSV 是两列：`source_id<TAB>target_id`。
面积查询的 TSV 是三列，包含所有经过包围盒筛选的候选对象对，包括零面积。
编号从 0 开始，分别对应各自数据文件的对象顺序；线程加载顺序不会改变编号。
行按源编号、目标编号排序。面积以输入坐标单位的平方表示，经纬度数据下
不是平方公里。公共接口已把 FARM 内部的两倍面积累加值换算为实际面积。

运行日志和 `candidate-pairs`、`count-found` 写入标准错误；标准输出用于
`-p` 请求的数据。对面积查询，`count-found` 为面积大于零的候选对数量。
JSON 格式包含 `count`、`candidates`、`matches`，每条记录包含 `source`、
`target`，面积查询另含 `area`。

## 统一的 join 命令

```bash
# 二维：也可以使用 --tile1/--tile2，与三维路径参数风格一致
./build/gpu/tdbase join --dim 2 -q intersect --tile1 source.dt --tile2 target.dt -g -p
./build/gpu/tdbase join --dim 2 -q intersection --tile1 source.dt --tile2 target.dt -p
./build/gpu/tdbase join --dim 2 -q within --tile1 source.dt --tile2 target.dt -d 2 -g -p

# 三维：原有参数语义保留，--dim 3 可省略
./build/gpu/tdbase join --dim 3 -q within --tile1 source3d.dt --tile2 target3d.dt -w 100 -l 100 -g -p
```

二维 `-w` 表示填充率编码位宽，距离使用 `-d` 或 `--within_dist`；
三维 `-w` 保留原来的距离阈值含义。二维不提供 KNN；三维的查询类型仍为
`intersect`、`within`、`nn`。

## FARM 参数

| 参数 | 含义 | 默认值 |
| --- | --- | --- |
| `-s`, `--source` | 源二维数据；也可用 `--tile1` | 必填 |
| `-t`, `--target` | 目标二维数据；也可用 `--tile2` | within 可省略 |
| `-g`, `--gpu` | 使用 GPU | CPU |
| `--cpu` | 显式使用 CPU | — |
| `-h`, `--hierarchy` | 层级栅格 | 关闭 |
| `-a`, `--approximation` | 近似模式 | 关闭 |
| `-n`, `--threads` | CPU 工作线程数，正整数 | 系统 CPU 数 |
| `-v`, `--vpr` | 每栅格顶点数参数 | 10 |
| `-w`, `--bitwidth` | 填充率编码位宽，2～8 | 4 |
| `-d`, `--within_distance` | 距离阈值，非负整数，公里 | 10 |
| `-b`, `--batch_size` | GPU 每批候选对数量 | 自动 |
| `-m`, `--merge_threshold` | GPU 合并阈值，0～1 | 0.9 |
| `-c`, `--approx_confidence` | 近似置信阈值，0～1 | 0.9 |
| `-l`, `--NLow` | 最低层参数，正整数 | 1 |
| `-u`, `--unroll_size` | GPU 精化展开长度 | 16 |
| `-p`, `--print_result` | 标准输出匹配行 | 关闭 |
| `-o`, `--output` | 结果文件 | 无 |
| `--format` | `tsv` 或 `json` | tsv |
| `--help` | 帮助；`-h` 是层级栅格 | — |

兼容名称 `tdbase 2d ...`、`tdbase farm within_polygon ...`、
`tdbase farm dump_polygons ...`、`tdbase farm filter_polygons ...` 也可使用。
原来的 `intersect[_cpu]`、`intersection[_cpu]`、`within_polygon[_cpu]`、
`dump_polygons`、`filter_polygons` 程序仍可运行，它们调用同一套库。
独立 GPU 程序默认使用 GPU；`_cpu` 程序默认使用 CPU。

## C++ 接口与构建

公共接口为 `src/farm2d/api.h`，命名空间 `tdbase::farm2d`，库目标为
`tdbase_farm2d`。接口提供 `execute`、`convert_wkt`、`filter_wkt`，
不暴露二维/三维同名的内部几何类型。

```cpp
#include "farm2d/api.h"

tdbase::farm2d::Options options;
options.source = "source.dt";
options.target = "target.dt";
options.query = tdbase::farm2d::Query::Intersection;
options.backend = tdbase::farm2d::Backend::CPU;
options.hierarchy = true;
auto result = tdbase::farm2d::execute(options);
// result.count、result.candidates、result.matches[i].area
```

在本项目的 CMake 中用 `target_link_libraries(your_target PRIVATE tdbase_farm2d)`。
`Options::collect_results=false` 可以只获取计数，避免保存结果行。
`gpu_available()` 表示构建是否包含 CUDA 支持；GPU 实际可用性仍由 CUDA
运行时检查。接口的输入错误通过 C++ 异常报告，命令行失败返回非零退出码。

```bash
# GPU 构建
bash script/build_gpu.sh

# 纯 CPU 构建（即使 shell 已 export USE_GPU=1，也会明确禁用）
cmake -S src -B build/unified-cpu -G Ninja -DUSE_GPU=OFF -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTING=ON
cmake --build build/unified-cpu --parallel 2

# 进程内接口和统一命令测试
ctest --test-dir build/gpu --output-on-failure
python3 tests/unified_cli.py build/gpu . build/unified-validation/gpu --gpu
```

`make -C src USE_GPU=1` 和 `make -C src/farm2d USE_GPU=1` 也会转到同一份
CMake 构建，默认目录为 `build/cmake`。通过 `-DFARM_PROFILE_QUERY=ON`
开启原有详细查询分析统计。

FARM SQL 示例、辅助脚本和许可证保留在 `src/farm2d/sql`、
`src/farm2d/scripts` 和 `src/farm2d/LICENSE`。SQL 仍在对应数据库中运行。
GPU 内核保留两个 4 GiB scratch buffer 的设计，在 RTX 4080 上应顺序执行。
