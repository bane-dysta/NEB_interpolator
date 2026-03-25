```{=latex}
\begin{titlepage}
\centering
\vspace*{2.8cm}
{\Huge\bfseries NEB interpo 软件说明书 \par}
\vspace{0.9cm}
{\Large 面向分子几何分析、RMSD 对齐与反应路径插值/NEB 优化的命令行工具集 \par}
\vspace{2.2cm}
\begin{tabular}{rl}
产品名称： & NEB interpo \\
版本： & 1.0 \\
文档类型： & 软件说明书 \\
语言： & 中文 \\
\end{tabular}
\vfill
{\large NEB interpo\par}
\vspace*{1.2cm}
\end{titlepage}

\tableofcontents
\newpage
```

# 概述

## 产品定位

NEB_interpo 是一组围绕分子三维结构处理构建的命令行程序，覆盖以下三类工作：

- 单个结构的几何分析与结构编辑；
- 两个结构之间的 RMSD 刚体对齐；
- 反应路径的端点插值与 NEB（Nudged Elastic Band）迭代更新。

该工具集面向以 XYZ 几何文件为中心组织的计算化学与分子建模工作流，适用于路径初猜构造、端点对齐、原子索引分析、结构局部导出以及外部量化程序驱动的 NEB 迭代。

## 适用场景

NEB_interpo 适用于以下场景：

- 在两个端点结构之间快速生成一组中间图像，用于反应路径初猜；
- 在进入 NEB、GSM 或其他路径算法之前，对端点进行 RMSD 刚体对齐；
- 对 XYZ 结构执行距离、角度、二面角、平面夹角和几何中心分析；
- 在 NEB 迭代中调用外部程序提供梯度或力；
- 在没有外部量化引擎时，使用内置的距离矩阵惩罚模型对路径进行几何平滑示例。

## 程序组成

当前版本以以下程序为正式命令入口；其中部分后端程序依赖构建条件：

| 程序 | 作用 |
| --- | --- |
| `xyzgeom` | XYZ 几何分析器，支持交互式菜单与少量命令行操作。 |
| `calc_rmsd_xyz` | 两个结构之间的 RMSD 计算与刚体对齐程序。该程序在 Fortran 与 BLAS/LAPACK 可用时构建。 |
| `neb_interpolator` | 路径插值与 NEB 迭代程序。 |
| `allign_lapack` | RMSD 对齐后端，可作为自动对齐回退后端之一。该程序在 BLAS/LAPACK 可用时构建。 |
| `allign_eigen` | RMSD 对齐后端，可作为自动对齐回退后端之一。 |

其中，`allign_lapack` 与 `allign_eigen` 主要作为对齐后端被自动调用；面向用户的主入口为 `xyzgeom`、`calc_rmsd_xyz` 与 `neb_interpolator`。

## 核心能力

- 读取标准单帧 XYZ 结构，并执行常见几何分析；
- 使用刚体旋转和平移对移动结构进行 RMSD 最小化对齐；
- 生成 LIC、LIIC 与 DM 三类中间路径；
- 以 LIC 或 LIIC 路径为初始猜测，执行 NEB 迭代；
- 在 NEB 迭代中调用外部引擎读取梯度或力；
- 自动查找可用的 RMSD 对齐后端；
- 输出独立 XYZ 文件序列或单个多帧轨迹文件。

## 核心对象

| 对象 | 说明 |
| --- | --- |
| endpoint | 路径端点结构，即初始结构与终态结构。 |
| image | 路径中的单个图像；`nimages` 仅统计中间图像，不包含两个端点。 |
| path | 从初始端点、若干中间图像到终态端点组成的完整路径。 |
| alignment backend | 用于刚体对齐的后端可执行程序，按查找顺序自动选择。 |
| external engine | 在 NEB 周期中为每个中间图像提供梯度或力的外部程序。 |

# 系统组成与构建

## 目录结构

源码包的典型结构如下：

```text
NEB_interpo/
  CMakeLists.txt
  README.md
  src/
    xyzgeom.cpp
    calc_rmsd_xyz.f90
    neb_interpolator.cpp
    neb_driver.{h,cpp}
    rmsd_align.{h,cpp}
    neb_core.{h,cpp}
    neb_engine.{h,cpp}
    internal_ic.{h,cpp}
    geometry.h
    util.h
    zmat/
  tests/
  examples/
  third_party/
```

## 构建要求

### 编译器与工具

| 组件 | 要求 |
| --- | --- |
| C++ 编译器 | 支持 C++17 |
| Fortran 编译器 | 用于构建 `calc_rmsd_xyz`；若缺失，则该程序不生成 |
| CMake | 3.10 或更高 |
| Python 3 | 仅用于部分测试脚本 |

### 线性代数依赖

| 依赖 | 用途 |
| --- | --- |
| BLAS/LAPACK | `calc_rmsd_xyz` 与 `allign_lapack` |
| Intel MKL | 可选；启用时优先作为 BLAS/LAPACK 后端 |
| Eigen | 源码包内附带，用于 `allign_eigen` |

### CMake 选项

当前版本公开使用的构建选项如下：

| 选项 | 默认值 | 说明 |
| --- | --- | --- |
| `CMAKE_BUILD_TYPE` | `Release` | 构建类型；未显式给出时默认为 `Release`。 |
| `USE_MKL` | `ON` | 若检测到可用 MKL，则使用 MKL；否则自动回退到系统 BLAS/LAPACK。 |

## 构建步骤

```bash
mkdir build
cd build
cmake ..
make -j4
```

构建完成后，主要可执行文件位于：

```text
build/bin/
```

## 安装

```bash
sudo make install
```

默认安装内容包括已成功构建的程序与相关文件，通常包括：

- `bin/` 下的可执行程序；
- `include/` 下的头文件；
- `share/NEB_interpo/examples/` 下的示例 XYZ 文件。

## 常用构建目标

| 目标 | 说明 |
| --- | --- |
| `all` | 构建全部可执行程序 |
| `xyzgeom` | 仅构建 `xyzgeom` |
| `calc_rmsd_xyz` | 仅构建 `calc_rmsd_xyz`；若缺失 Fortran 或 BLAS/LAPACK，则该目标不可用 |
| `neb_interpolator` | 仅构建 `neb_interpolator` |
| `test_all` | 运行全部自定义测试目标 |
| `usage` | 打印可用目标与示例命令 |
| `install` | 执行安装 |

# 输入数据与基本规则

## XYZ 输入格式

`xyzgeom` 与 `neb_interpolator` 使用相同的 XYZ 读取逻辑。当前支持的输入格式为标准单帧 XYZ 文本：

```text
<N>
<Comment line>
<Symbol1> <x1> <y1> <z1>
<Symbol2> <x2> <y2> <z2>
...
```

使用规则如下：

- 第 1 行必须为原子数；
- 第 2 行为注释行，可为空；
- 后续必须为逐原子的元素符号与三列坐标；
- 当前版本以标准四列 XYZ 为正式输入形式；额外附加列不属于正式输入范围；
- 读取逻辑按单帧 XYZ 处理；若文件中包含多帧内容，仅首帧参与读取；
- `xyzgeom` 与 `neb_interpolator` 要求输入内容与 XYZ 语义一致，文档统一使用 `.xyz` 作为扩展名。

## Gaussian 输入格式

`calc_rmsd_xyz` 当前支持以下两类文件作为输入：

- `.xyz`
- `.gjf`

两种格式可以混合使用，例如：

```bash
calc_rmsd_xyz ref.gjf mobile.xyz
```

当前版本不将 `.com` 作为 `calc_rmsd_xyz` 的正式输入扩展名。

## 结构兼容性规则

### `neb_interpolator` 与 `xyzgeom` 路径插值

在路径插值与 NEB 处理中，两个端点必须满足以下条件：

- 原子数完全一致；
- 原子顺序完全一致；
- 元素符号逐项完全一致。

当前版本的兼容性检查以符号文本逐项比较为准。因此，元素符号的大小写与书写形式应在两个端点中保持一致。

### `calc_rmsd_xyz`

`calc_rmsd_xyz` 不执行化学对应关系推断。两个输入文件之间的原子对应关系由用户负责保证。

## Gview 风格原子选择语法

`xyzgeom` 的交互式几何分析使用 Gview 风格的原子选择字符串。当前支持的语法如下：

```text
1
1,2,5
1-4
1,3,6-10
```

规则如下：

- 索引为 1 起始；
- 区间 `a-b` 为闭区间；
- 支持半角逗号与全角逗号；
- 距离、角度、二面角、向量等有顺序含义的操作保留输入顺序；
- “导出一组原子”“镜像平面定义”等集合型操作会自动去重并按索引排序。

# 对齐后端查找规则

`xyzgeom` 与 `neb_interpolator` 在需要执行 RMSD 对齐时，会自动选择可用后端。当前查找规则如下：

1. 若显式给出路径且该路径可执行，则直接使用该路径；
2. 在 `PATH` 中依次查找：
   - 指定默认名或显式路径对应的文件名；
   - `calc_rmsd_xyz`；
   - `allign_lapack`；
   - `allign_eigen`；
3. 在当前可执行程序所在目录中按相同顺序查找；
4. 若以上均未命中，则保留默认路径字符串，并在实际调用时由系统报告执行失败。

`neb_interpolator` 的 `-r`/`--rmsd-exec` 选项用于指定优先路径。

# 命令行工具

# `xyzgeom`

## 概述

`xyzgeom` 是面向单个 XYZ 结构的几何分析器。该程序支持两种运行方式：

- 交互式菜单模式；
- 命令行单次处理模式。

## 运行方式

### 交互式模式

当仅提供一个 XYZ 文件时，程序进入交互式模式：

```bash
xyzgeom molecule.xyz
```

### 命令行模式

命令行模式的正式语法如下：

```bash
xyzgeom <xyz_file> --print
xyzgeom <xyz_file> --swap <selection>
xyzgeom <xyz_file> --mirror <selection>
```

程序要求输入文件始终位于第一个位置参数；命令行选项不采用 GNU 风格的“纯选项入口”。

## 命令行选项

| 选项 | 参数 | 说明 |
| --- | --- | --- |
| `--print` | 无 | 将当前结构以 XYZ 形式打印到标准输出。 |
| `--swap` | `<selection>` | 交换两个原子。选择字符串必须解析为恰好两个原子。 |
| `--mirror` | `<selection>` | 将整个分子关于由三个原子定义的平面镜像。选择字符串必须解析为恰好三个原子。 |

### 命令行输出行为

- `--print` 直接输出当前结构；
- `--swap` 与 `--mirror` 在内存中完成变换后，将变换后的 XYZ 输出到标准输出；
- 命令行模式不自动覆盖原文件；
- 命令行模式不自动生成单独输出文件，除非将标准输出重定向到文件。

## 交互式菜单

当前版本的交互式菜单包含以下条目：

| 菜单项 | 功能 |
| --- | --- |
| `-1` | 将当前结构打印到屏幕 |
| `1` | 计算两个原子之间的距离与方向向量 |
| `2` | 交换两个原子 |
| `3` | 计算键角与对应平面法向量 |
| `4` | 计算两个平面之间的夹角 |
| `5` | 关于选定平面镜像整个分子 |
| `6` | 计算二面角 |
| `7` | 计算一组原子的几何中心 |
| `8` | 导出选定原子为新的 XYZ 文件 |
| `9` | 查找给定半径内的原子并导出 |
| `10` | 导出当前结构为新的 XYZ 文件 |
| `11` | 读取新的 XYZ 文件 |
| `12` | 切换距离单位显示（Å / Bohr） |
| `13` | 与第二个 XYZ 文件进行 RMSD 对齐 |
| `14` | 与第二个 XYZ 文件进行路径插值（LIC 或 LIIC） |
| `0` | 退出 |

## 功能说明

### 距离与向量

- 输入为两个原子索引，顺序为 `i,j`；
- 输出距离以及方向向量 `i -> j`；
- 距离单位与向量单位受当前显示单位控制。

### 键角

- 输入为三个原子 `i,j,k`，其中 `j` 为顶点；
- 输出键角（度）以及由两条键向量张成平面的单位法向量；
- 当任一边长接近零时，程序报错并终止该操作。

### 平面夹角

- 分别输入两个平面，各由 3 个原子定义；
- 输出为两个平面之间的锐角；
- 当前实现对法向量取绝对点积，因此结果范围为 `0°` 到 `90°`。

### 二面角

- 输入为四个按顺序排列的原子 `i,j,k,l`；
- 输出有符号二面角，范围为 `-180°` 到 `180°`；
- 当任一中间平面退化时，程序报错。

### 几何中心

- 计算所选原子坐标的算术平均值；
- 当前实现为几何中心，不进行质量加权；
- 输出单位受当前显示单位控制。

### 镜像

- 由三个原子定义镜像平面；
- 镜像操作作用于当前结构中的全部原子；
- 定义平面的三个原子若共线或重复，操作失败。

### 导出选定原子

- 输入一组原子选择；
- 生成新的 XYZ 文件；
- 输出文件仅包含选定原子，索引按升序写出。

### 半径搜索

- 选择一个中心原子与搜索半径；
- 程序输出命中的原子索引列表（Gview 风格）；
- 随后写出一个新的 XYZ 文件；
- 导出文件中，中心原子始终位于首行，随后为命中的邻近原子。

### 与第二个 XYZ 对齐

- 参考文件固定不动；
- 当前已加载结构作为移动结构执行刚体对齐；
- 对齐成功后，结果只加载到内存，不自动覆盖磁盘上的原始输入文件；
- 若对齐失败，当前结构保持不变。

### 与第二个 XYZ 插值

- 当前已加载结构作为初始端点；
- 第二个 XYZ 文件作为终态端点；
- 可选择 `LIC` 或 `LIIC`；
- 默认中间图像数为 `5`；
- 默认输出前缀为 `intrp_`；
- 输出模式固定为独立 XYZ 文件序列；
- 该功能内部调用与 `neb_interpolator` 相同的路径驱动器。

当用户选择 `LIIC` 时，若 LIIC 初始化失败而触发了 DM 或 LIC 回退，程序会给出警告并继续输出实际生成的路径。

## 单位行为

`xyzgeom` 的单位切换仅影响以下内容：

- 屏幕上输出的距离、向量与几何中心；
- 半径搜索时输入的半径值解释；
- 当前菜单标题中的单位显示。

当前结构在内存中始终保持原始坐标数值。导出到 XYZ 文件时，程序不执行坐标单位换算。

## 当前行为与限制

- `xyzgeom` 的正式命令行入口仅包括 `--print`、`--swap` 与 `--mirror`；
- `--help` 不是正式选项。若将 `--help` 作为第一个位置参数之后的文件名使用，程序不会进入帮助模式；
- 对未知命令行选项，当前版本不输出专门的错误提示；
- 路径插值菜单只提供 `LIC` 与 `LIIC`，不提供 `DM`、`NEB` 或 `multiframe` 输出模式；
- 路径插值菜单默认启用端点对齐，且不提供关闭对齐的交互式选项；
- 输入必须为单帧 XYZ；
- 菜单中的结构操作在内存中进行，导出前不会自动写回原文件。

# `calc_rmsd_xyz`

## 概述

`calc_rmsd_xyz` 用于计算两个结构之间的 RMSD，并将第二个输入文件刚体对齐到第一个输入文件。该程序接受 `.xyz` 与 `.gjf` 输入，并支持对选定原子范围进行拟合。

## 用法

```bash
calc_rmsd_xyz <file1> <file2>
calc_rmsd_xyz <file1> <file2> <range1> <range2>
```

示例：

```bash
calc_rmsd_xyz ref.xyz mobile.xyz
calc_rmsd_xyz ref.gjf mobile.xyz
calc_rmsd_xyz ref.xyz mobile.xyz 1-10 1-10
```

## 参数说明

| 参数 | 说明 |
| --- | --- |
| `file1` | 参考结构，保持不动。支持 `.xyz` 与 `.gjf`。 |
| `file2` | 移动结构，将被旋转和平移后对齐到 `file1`。支持 `.xyz` 与 `.gjf`。 |
| `range1` | 参考结构的原子范围，格式为 `i-j`。 |
| `range2` | 移动结构的原子范围，格式为 `i-j`。 |

### 范围规则

- 范围为 1 起始闭区间；
- `range1` 与 `range2` 必须包含相同数量的原子；
- 当指定范围时，仅使用该范围内的原子计算拟合变换与 RMSD；
- 范围之外的原子会应用同一刚体变换一起移动。

## 输入输出行为

### 标准输出

程序在标准输出打印 RMSD 值，例如：

```text
RMSD =     0.000000
```

### 输出文件

程序不会修改第一个输入文件。第二个输入文件会生成一个新的对齐结果文件：

| `file2` 类型 | 输出文件 |
| --- | --- |
| `mobile.xyz` | `mobile_new.xyz` |
| `mobile.gjf` | `mobile_new.gjf` |

输出文件位于第二个输入文件所在目录。

## 当前行为与限制

- 当前正式支持的扩展名仅为 `.xyz` 与 `.gjf`；
- 参考结构固定，移动结构为第二个输入文件；
- 程序不执行元素映射、重排或化学等价判断；
- 原子对应关系由用户保证；
- 程序要求参数个数为 2 或 4；
- 范围语法只接受单个 `i-j` 形式，不接受列表式选择；
- 输出文件通过为第二个输入文件添加 `_new` 后缀生成，不会默认覆盖原文件。

# `neb_interpolator`

## 概述

`neb_interpolator` 用于在两个 XYZ 端点之间生成路径，并可进一步执行 NEB 迭代。该程序既可以只做几何插值，也可以在每个 NEB 周期中调用外部引擎提供梯度或力。

## 用法

```bash
neb_interpolator [options] <initial.xyz> <final.xyz>
```

## 方法总览

| 方法 | 说明 |
| --- | --- |
| `lic` | 笛卡尔坐标线性插值。 |
| `liic` | 基于内坐标的线性插值。 |
| `dm` | 基于距离矩阵误差最小化的插值。 |
| `neb` | 使用 LIC 作为初始路径执行 NEB。 |
| `neb-liic` | 使用 LIIC 作为初始路径执行 NEB；必要时按 LIIC -> DM -> LIC 回退。 |

## 基本选项

| 选项 | 默认值 | 说明 |
| --- | --- | --- |
| `-n`, `--nimages NUM` | `5` | 中间图像数量，不包含两个端点。 |
| `-m`, `--method METHOD` | `neb` | 路径方法；可取 `lic`、`liic`、`dm`、`neb`、`neb-liic`。 |
| `-p`, `--prefix PREFIX` | 空字符串 | 输出文件名前缀。 |
| `-o`, `--output MODE` | `separate` | 输出模式；`separate` 或 `multiframe`。 |
| `-s`, `--step STEP` | `0.0001` | NEB 坐标更新步长。仅对 `neb` 与 `neb-liic` 生效。 |
| `-c`, `--conv THRESHOLD` | `0.01` | NEB 收敛阈值。仅对 `neb` 与 `neb-liic` 生效。 |
| `-i`, `--maxiter ITER` | `10000` | NEB 最大循环次数。仅对 `neb` 与 `neb-liic` 生效。 |
| `-a`, `--align` | 开启 | 使用 RMSD 后端对齐端点。 |
| `--no-align` | 关闭 | 禁用端点对齐。 |
| `-r`, `--rmsd-exec PATH` | 自动查找 | 指定 RMSD 对齐后端优先路径。 |
| `-h`, `--help` | 无 | 打印帮助并退出。 |

## LIIC 选项

以下选项在 `liic` 与 `neb-liic` 中生效：

| 选项 | 默认值 | 说明 |
| --- | --- | --- |
| `--bond-factor F` | `1.25` | 建立键连接图时，键长截断系数。 |
| `--fd-step H` | `1e-4` | 构造原始 B 矩阵的有限差分步长。 |
| `--ev-thresh T` | `1e-3` | DLC 特征值筛选阈值。 |
| `--liic-maxiter N` | `50` | 内坐标回代最大迭代次数。 |
| `--liic-tol T` | `1e-4` | 原始内坐标 RMS 残差收敛阈值。 |
| `--liic-damp D` | `1e-8` | 回代时加入的阻尼项。 |
| `--liic-max-step S` | `0.20` | 单次回代允许的最大笛卡尔位移。 |
| `--liic-verbose V` | `0` | LIIC 详细级别：`0` 静默，`1` 每图像，`2` 每迭代。 |

## DM 选项

以下选项在 `dm` 中生效，也作为 LIIC 失败后的回退参数：

| 选项 | 默认值 | 说明 |
| --- | --- | --- |
| `--dm-maxiter N` | `800` | DM 最大迭代次数。 |
| `--dm-step S` | `5e-3` | DM 梯度下降步长。 |
| `--dm-tol T` | `1e-3` | 距离矩阵 RMS 误差阈值。 |
| `--dm-max-step S` | `0.20` | 单次迭代允许的最大笛卡尔位移。 |
| `--no-dm-fallback` | 关闭 | 禁用 LIIC -> DM 回退，改为 LIIC 失败后直接回退到 LIC。 |

## 输出模式

### `separate`

当 `-o separate` 时，程序输出独立 XYZ 文件序列：

```text
<prefix>00.xyz
<prefix>01.xyz
...
<prefix>NN.xyz
```

规则如下：

- `00.xyz` 为初始端点；
- `01.xyz` 到 `0N.xyz`/`NN.xyz` 为中间图像；
- 最后一个文件为终态端点；
- 总文件数为 `nimages + 2`；
- 文件输出到当前工作目录，或输出前缀中指定的现有目录。

### `multiframe`

当 `-o multiframe` 时，程序写出：

```text
<prefix>trajectory.xyz
```

该文件按多帧 XYZ 形式依次写出：

1. 初始端点；
2. 各中间图像；
3. 终态端点。

## 处理流程

`neb_interpolator` 的标准执行流程如下：

1. 读取初始端点与终态端点；
2. 检查原子数与原子符号顺序是否一致；
3. 若启用对齐，则尝试将终态端点刚体对齐到初始端点；
4. 根据所选方法生成初始路径：
   - `lic` 直接使用笛卡尔线性插值；
   - `liic` 采用内坐标插值，失败时可回退到 DM 或 LIC；
   - `dm` 采用距离矩阵插值，失败时回退到 LIC；
   - `neb` 使用 LIC 作为起始路径；
   - `neb-liic` 使用 LIIC 作为起始路径，并继承 LIIC 的回退链；
5. 若为 `neb` 或 `neb-liic`，执行 NEB 迭代；
6. 将结果写出为独立 XYZ 文件或多帧轨迹文件。

## 端点对齐规则

当启用对齐时，程序行为如下：

- 初始端点作为参考结构；
- 终态端点作为移动结构；
- 程序先计算对齐前 RMSD；
- 若后端成功生成对齐结果，则比较对齐前后 RMSD；
- 仅当对齐后 RMSD 更小，才采用对齐后的终态结构；
- 若对齐失败、无法读取对齐结果，或对齐后 RMSD 未改善，则保留原始终态结构继续计算。

对齐过程使用临时文件，不会修改用户提供的端点文件。

## 方法说明

### LIC

LIC 直接对每个原子坐标在笛卡尔空间进行线性插值。该方法速度最快，不依赖键连接或内部坐标构造，但在大幅构象变化或刚体旋转未完全消除的情况下，中间图像可能不够平滑。

### LIIC

LIIC 基于内部坐标插值。当前实现包含以下步骤：

- 根据端点结构的并集键图建立原始键、角与二面角集合；
- 使用有限差分构造原始 B 矩阵；
- 对 `G = B B^T` 做特征分解并选择 DLC 子空间；
- 在 DLC 空间中线性插值，并回代到笛卡尔坐标；
- 必要时执行 Kabsch 对齐和步长截断。

LIIC 更适合需要保持内部几何连续性的路径初猜。若 LIIC 构造失败，程序会按配置回退到 DM 或 LIC。

### DM

DM 通过最小化当前图像的距离矩阵与目标距离矩阵之间的误差来构造路径。该方法适合在 LIIC 失效时提供较平滑的备用路径。

### NEB 与 NEB-LIIC

当前 NEB 实现遵循以下规则：

- 两个端点固定不动，仅更新中间图像；
- 切向量采用全局形式 `R(i+1) - R(i-1)` 并归一化；
- 实力项取真实力在切向量垂直方向上的分量；
- 弹簧力取相邻图像间距差在切向量方向上的分量；
- 每个周期用固定步长 `step` 更新坐标；
- 收敛判据为所有中间图像全部原子上的最大 NEB 力模长。

当前版本不包含以下功能：

- climbing image；
- 能量加权切向量；
- 自适应步长；
- 端点能量或每周期能量输出。

## 外部引擎模式

外部引擎仅对 `neb` 与 `neb-liic` 生效。其作用是在每个 NEB 周期为各中间图像提供梯度或力。

### 选项

| 选项 | 默认值 | 说明 |
| --- | --- | --- |
| `--engine-cmd CMD` | 未启用 | 启用外部引擎，并以 `CMD` 作为调用命令。 |
| `--engine-in FILE` | `neb_engine_in.dat` | 引擎输入文件名。 |
| `--engine-out FILE` | `neb_engine_out.dat` | 引擎输出文件名。 |
| `--engine-units U` | `Angstrom` | 写入引擎输入/输出头部的单位标签。当前仅作标签使用，不执行单位换算。 |
| `--engine-output-is-force` | 关闭 | 将引擎输出解释为力；默认解释为梯度。 |
| `--engine-spring K` | `1.0` | 外部引擎模式下的 NEB 弹簧常数。 |
| `--engine-run-every N` | `1` | 每 `N` 个周期运行一次引擎；跳过的周期复用上一次的结果。 |
| `--engine-keep-cycle-files` | 关闭 | 为每个周期保留独立的引擎输入/输出文件。 |

### 调用规则

若命令字符串中包含占位符，程序会执行替换：

| 占位符 | 含义 |
| --- | --- |
| `{in}` | 引擎输入文件路径 |
| `{out}` | 引擎输出文件路径 |
| `{cycle}` | 当前 NEB 周期编号 |

若命令中缺少 `{in}` 或 `{out}`，程序会自动在命令末尾追加缺失的文件路径；若两个占位符都不存在，则自动追加 `<infile> <outfile>`。

同时，程序会设置以下环境变量：

| 变量 | 含义 |
| --- | --- |
| `NEB_ENGINE_IN` | 当前周期的引擎输入文件路径 |
| `NEB_ENGINE_OUT` | 当前周期的引擎输出文件路径 |
| `NEB_CYCLE` | 当前周期编号 |
| `NEB_NATOMS` | 原子数 |
| `NEB_NIMAGES` | 中间图像数 |

### 输入文件格式

引擎输入文件由多个 XYZ 块组成，每个块对应一个中间图像：

```text
<natoms>
image=<i> units=<U> type=xyz cycle=<cycle>
<sym> <x> <y> <z>
...
```

规则如下：

- `image` 从 `1` 开始编号；
- `units` 为标签；
- `type` 固定为 `xyz`；
- `cycle` 为当前 NEB 周期编号。

### 输出文件格式

外部引擎必须为每个中间图像输出一个向量块：

```text
<natoms>
image=<i> units=<U> type=gradient
<sym> <vx> <vy> <vz>
...
```

或：

```text
<natoms>
image=<i> units=<U> type=force
<sym> <fx> <fy> <fz>
...
```

当前解析规则如下：

- 输出块数量必须与中间图像数量一致；
- `image=` 可省略；省略时按文件中块的顺序依次分配给第 1、2、3... 个中间图像；
- 若指定了 `image=`，则不允许重复；
- 每个块的原子数必须与当前图像一致；
- 每个块的元素符号顺序必须与当前图像一致；
- 默认情况下，输出向量解释为梯度，程序内部按 `force = -gradient` 转换；
- 当指定 `--engine-output-is-force` 时，输出向量按力直接使用。

### 默认引擎

当未提供 `--engine-cmd` 时，`neb` 与 `neb-liic` 使用内置 `DistancePenaltyEngine`。该引擎根据当前图像距离矩阵与线性目标距离矩阵的偏差生成虚拟力。

该模式的作用是：

- 在无外部量化引擎时验证程序流程；
- 对路径做几何平滑示例；
- 为开发与测试提供无外部依赖的 NEB 驱动。

该模式不代表真实势能面上的 NEB 计算结果。

## 结果写出与退出行为

### 插值方法

`lic`、`liic` 与 `dm` 在成功生成路径后直接写出结果文件并退出。

### NEB 方法

`neb` 与 `neb-liic` 在完成最多 `maxiter` 个周期后写出结果文件。当前实现中：

- 当 `max_force < conv` 时，程序会报告收敛并提前结束；
- 当达到 `maxiter` 但未满足收敛阈值时，程序仍会正常写出最后一轮图像并返回成功状态；
- 程序不会因为“未收敛”而自动返回失败码。

# 示例

## 1. 直接生成 LIC 路径

```bash
neb_interpolator -n 10 -m lic initial.xyz final.xyz
```

输出：

```text
00.xyz
01.xyz
...
11.xyz
```

## 2. 生成 LIIC 初猜并写出多帧轨迹

```bash
neb_interpolator --no-align -o multiframe -n 8 -m liic -p liic_ initial.xyz final.xyz
```

输出：

```text
liic_trajectory.xyz
```

## 3. 使用外部引擎执行 NEB

```bash
neb_interpolator -m neb -n 8 \
  --engine-cmd "python3 engine.py {in} {out}" \
  --engine-units AU \
  initial.xyz final.xyz
```

## 4. 对齐第二个结构并生成对齐结果

```bash
calc_rmsd_xyz ref.xyz mobile.xyz
```

输出：

```text
RMSD = ...
mobile_new.xyz
```

## 5. 使用选择范围进行局部拟合

```bash
calc_rmsd_xyz ref.xyz mobile.xyz 1-20 1-20
```

## 6. 以命令行方式打印或修改结构

```bash
xyzgeom molecule.xyz --print
xyzgeom molecule.xyz --swap 1,3
xyzgeom molecule.xyz --mirror 1-3
```

# 当前行为与限制

## 输入与格式

- `xyzgeom` 与 `neb_interpolator` 按单帧 XYZ 读取输入；
- `neb_interpolator` 不读取 Gaussian 输入文件；
- `calc_rmsd_xyz` 当前只正式支持 `.xyz` 与 `.gjf`；
- `neb_interpolator` 的端点兼容性检查要求原子数、原子顺序和元素符号完全一致；
- 多帧 XYZ 可作为输出格式，但再次作为 `xyzgeom` 或 `neb_interpolator` 输入时，仅首帧会被读取。

## 对齐

- 自动对齐是“尽量执行”的步骤，而不是强制条件；
- 当对齐失败时，`neb_interpolator` 与 `xyzgeom` 会继续使用原始终态结构；
- 对齐仅执行刚体平移与旋转，不做原子重排。

## LIIC / DM

- LIIC 与 DM 依赖元素符号能够映射到共价半径数据；
- 当元素符号无法识别时，LIIC/DM 初始化会返回错误；
- `liic` 与 `neb-liic` 可能在成功返回的同时实际采用 DM 或 LIC 回退路径。

## NEB

- 默认 `neb`/`neb-liic` 在未给出外部引擎时使用虚拟距离矩阵惩罚力；
- 当前 NEB 实现不输出能量曲线；
- 当前 NEB 实现不包含 climbing image 与自适应更新策略；
- 达到 `maxiter` 而未收敛时，程序仍返回成功并写出最后结果。

## 外部引擎

- `--engine-cmd` 仅对 `neb` 与 `neb-liic` 生效；对纯插值方法无效；
- `--engine-run-every N` 在 `N>1` 时会复用上一个实际调用周期的向量结果；
- `--engine-units` 仅写入标签，不进行任何单位转换；
- 外部引擎必须对每个中间图像都返回一个完整向量块。

## `xyzgeom`

- 交互模式中的单位切换只影响显示与半径输入解释，不改变实际存储坐标；
- `xyzgeom` 命令行模式不自动写文件；
- `xyzgeom` 当前没有正式的 `--help` 入口；
- 路径插值菜单不提供 `DM` 与 `NEB`。

# 附：推荐工作流

## 端点预处理

1. 使用 `xyzgeom` 或 `calc_rmsd_xyz` 检查端点结构的一致性；
2. 必要时先执行刚体对齐；
3. 确认两个端点原子顺序完全对应。

## 路径初猜生成

- 几何变化较小：优先使用 `lic`；
- 需要更平滑的内部几何连续性：优先使用 `liic`；
- `liic` 失败时，可保留默认 DM 回退。

## 真实 NEB 计算

1. 使用 `liic` 或 `neb-liic` 生成更合理的初始路径；
2. 提供外部引擎命令，确保每个周期返回完整梯度或力；
3. 根据外部程序数值特性调节 `--step`、`--conv`、`--maxiter` 与 `--engine-spring`；
4. 仅在明确接受数值近似时使用 `--engine-run-every N`。
