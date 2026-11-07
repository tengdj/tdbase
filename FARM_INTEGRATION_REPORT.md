# FARM 与 TDBase 统一集成验收

日期：2026-11-05。Windows 源码目录 `D:\tdfarm\tdbase`，WSL 源码目录
`/home/moli/tdfarm/tdbase`。原始 `D:\tdfarm\FARM` 仓库保留。

## 最终实现

- FARM 核心、栅格、索引、存储和 CUDA 实现位于 `src/farm2d`，已迁入
  `tdbase::farm2d` 命名空间，头文件保护宏也已隔离。
- 新增 `tdbase_farm2d` 共享使用的静态库目标；`tdbase` 直接链接二维实现。
  不通过 shell、子进程或同目录辅助程序转发查询。
- `tdbase farm convert/filter/intersect/intersection/within` 提供二维命令。
  `tdbase join --dim 2` 与 `--dim 3` 提供统一维度选择；省略维度保持原三维行为。
- 所有原 FARM 查询参数仍可使用；统一入口支持 CPU/GPU 显式选择。
  原独立二维程序成为同一库的兼容入口。
- 新增 C++ 公共接口 `src/farm2d/api.h`，提供配置、查询、转换、过滤和结构化结果。
- 支持匹配对象对、实际交集面积、TSV/JSON，以及标准输出和文件导出。
  输出编号对应源/目标文件各自的位置；GPU 的全局内部编号会换回本地编号。
- 修正多线程分块加载造成的对象顺序不确定；释放查询持有的对象、GPU
  和 pinned-memory 缓冲，支持同一进程内重复调用。
- 输入路径、参数、二维二进制格式和 WKT 进行入口校验；过滤器不再依赖硬编码路径。
  显式拒绝内环，避免原始 FARM 读取逻辑静默丢失孔洞。
- Make 和 CMake 指向同一套实现。SQL、辅助脚本、许可证均保留。

## 构建与验证

| 验证项 | 结果 |
| --- | --- |
| CUDA 12.6、RTX 4080 sm_89、Release 构建 | 通过 |
| 显式 `USE_GPU=OFF`、Debug 构建 | 通过 |
| GPU 构建的统一命令测试 | 91/91 |
| 纯 CPU 构建的统一命令测试 | 64/64 |
| GPU 构建的进程内 C++ API CTest | 1/1 |
| 纯 CPU 构建的进程内 C++ API CTest | 1/1 |
| 原 CPU 命令和三维回归 | 53/53 |
| 原 GPU 命令和三维 CPU/GPU 对照 | 55/55 |

统一命令测试把 `tdbase` 复制到只有该程序的独立目录，再从不同的当前目录
运行，PATH 不包含构建目录，以验证不依赖 FARM 辅助可执行文件。

二维测试覆盖 CPU/GPU、平面/层级栅格、2/4/8 位编码、相交、包含、分离、
距离阈值、自连接、近似模式、数据转换、过滤、空数据、TSV/JSON、文件导出、
错误参数、错误维度与错误数据格式。数值面积断言包括部分重叠和完全包含的
矩形，期望面积均为 0.25。14 万零 1 个对象的数据用于跨越 10 MiB 加载任务
边界，验证多线程加载后仍输出正确的文件位置编号。

C++ 测试在同一个进程内重复转换、调用查询、切换 CPU/GPU，并检查
面积、距离结果与对象编号。原三维渐进 LOD 和 LOD100 距离回归仍通过；
GPU 日志确认执行三维 GPU 距离计算，并与 CPU 输出对象对一致。

结果与日志保存于两份项目的 `build/unified-validation/`，包括各项
`results.json`、构建日志、CTest 日志和二进制 SHA256 清单。

## 可执行文件

- GPU 版本：`/home/moli/tdfarm/tdbase/build/gpu/tdbase`。
- 纯 CPU Debug 构建：`/home/moli/tdfarm/tdbase/build/unified-cpu/tdbase`。
- 上述产物另复制到两份项目的 `build/cpu/` 和 Windows 的 `build/gpu/`，
  可执行文件均为 Linux ELF，应通过 WSL 运行。
- 原来的开发测试程序改名为 `tdbase_test`，避免与 CTest 的 `test` 保留目标冲突。

完整使用方法和接口说明见 [FARM_USAGE.md](FARM_USAGE.md)。

## 验证边界

本轮完成 FARM 的模块、命令、结果和构建集成。输入几何沿用 FARM 算子的
能力范围，目前公开入口支持无洞的 `POLYGON`；没有新增 MULTIPOLYGON 或
内环算法。二维数据不会自动转成三维 Tile。原三维实验性相交问题也未在
本轮改写。

数值测试覆盖了上述代表性几何，不等于任意复杂几何的完整精度证明或性能
基准。Compute Sanitizer 仍受 Windows WDDM 调试接口限制，未作为通过项。
GPU 算子保留原有两个固定 4 GiB scratch buffer 的设计。
