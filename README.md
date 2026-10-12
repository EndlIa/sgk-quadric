# sgk-quadric

基于海星 SGGK 1.7.11 的二次曲面识别与转换组件:判定 B-rep 面能否精确等价为平面 / 圆柱 / 圆锥 / 球 / 圆环,
能等价则拟合解析参数、核验几何精度并把该面的底层曲面替换为解析面,不能等价则拒绝并保留原面。

## 1. 简介

NURBS 是 CAD 内核的通用曲面表示:解析面、扫描面与自由曲面统一写成样条形式,内核因而只需维护一套求值与存储路径。
代价是面一旦被样条化,其求交、投影与布尔运算都必须走样条算法,开销显著高于解析算法。本组件从样条面自身的几何微分量
判定其形状,在可能时恢复解析表示,使下游运算重新走解析路径。


## 2. 目录

```text
sgk-quadric/
├─ CMakeLists.txt          # 2 个静态库 + 14 个可执行目标
├─ cmake/FindSGGK.cmake    # SDK 查找模块
├─ build.bat               # 配置与构建入口
├─ include/sgkquadric/     # 契约头
├─ src/                    # 实现
├─ apps/
│  ├─ checks/              # 7 个自检程序
│  ├─ harness/             # 语料与基准:run_corpus / bench_perf / batch_stability
│  └─ probes/              # 冒烟与探针:hello / probe_writeback / probe_speed
├─ tools/                  # gen_corpus
└─ data/                   # 语料清单与原始数据记录
```

## 3. 如何构建

### 3.1 前置条件

| 项 | 要求 |
|---|---|
| 平台 | Windows x64 |
| 编译器 | MSVC |
| 构建工具 | CMake ≥ 3.20 与 Ninja(随 VS 的 CMake 组件安装) |
| SGGK SDK |`sggk.lic` 必须位于 SDK 根目录 |

SDK 只提供 MSVC ABI 的导入库与 DLL,故必须使用 MSVC 工具链。

### 3.2 构建

```bat
build.bat
```

### 3.3 运行

14 个可执行目标输出到 `build/`。

| 程序 | 用途 | 通过判据 |
|---|---|---|
| `hello` | SDK 连通性冒烟 | 退出码 0,打印版本 |
| `test_math` | 数学核心自检 | `failed=0`(30 项) |
| `test_io` | `.sgt` / STEP 往返、损坏输入拒绝、版本字段 | `failed=0`(18 项) |
| `test_detect` | 判型:五类正例、裁剪块、退化点、负例拒绝 | `failed=0`(14 项) |
| `test_fit` | 拟合:五类参数恢复、精修收敛、退化拒绝 | `failed=0`(98 项) |
| `test_pipeline` | 端到端:判定→拟合→核验→写回 | `failed=0`(33 项) |
| `test_exceptions` | 异常与越界输入 | `failed=0`(16 项) |
| `test_g1_tangent` | 相切邻面 G1 证据 | `failed=0`(14 项) |
| `gen_corpus` | 生成带真值的合成语料 | 输出 `[CORPUS]` 统计,无 `[SKIP]` / `[WARN]` |
| `run_corpus` | 合成集批量判定与混淆表 | `accuracy=100.00%`、`loadFailed=0` |
| `bench_perf` | 判定/流水线耗时与转换前后配对操作 | 两行 `[ASSERT]` 均为 `true` |
| `batch_stability` | ≥1000 个不同曲面的连续处理与跨遍一致性 | `inconsistent=0`、`caseExceptions=0`、`loadFailed=0` |
| `probe_writeback` | 写回自洽探针(NURBS↔解析往返、pcurve) | 见 `[SUMMARY] (1)…(4)` 行的 `true` / `false` 判定 |
| `probe_speed` | 投影核探针(为 `bench_perf` 的投影替代方案提供对拍) | 观察三项加速比 |

任一自检程序以 `[SUMMARY] total=<断言数> failed=<n>` 结尾,`failed=0` 且退出码 0 为通过。
全量复跑顺序:

```bat
build.bat
cd build
gen_corpus.exe --seed 20261008 --cases-per-type 250
test_math.exe & test_io.exe & test_detect.exe & test_fit.exe
test_pipeline.exe & test_exceptions.exe & test_g1_tangent.exe
run_corpus.exe
bench_perf.exe --repeat 20
batch_stability.exe --min-unique-cases 1000 --passes 2
```
