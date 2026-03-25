*README由GPT 5.2，可查看帮助；彩虹屁勿信。*

# 分子几何工具

用于分子结构分析、RMSD 对齐，以及反应路径插值/NEB 优化的一组小工具。

## 概述

本项目包含三个主要程序：

- **xyzgeom**：交互式/命令行分子几何分析器（距离/角度/二面角/导出/镜像/对齐/插值等）。
- **calc_rmsd_xyz**：Fortran/LAPACK 实现的 RMSD 计算与结构对齐（生成 `*_new.xyz`）。
- **neb_interpolator**：反应路径插值与 NEB（Nudged Elastic Band）优化的命令行工具。

## 功能特性

### xyzgeom

- 计算距离、角度和二面角
- 原子操作（交换、镜像、导出）
- RMSD 对齐（调用 `calc_rmsd_xyz`；与 `neb_interpolator` 使用同一套“查找 + 调用”逻辑）
- 与第二个 XYZ 的 LIC/NEB 路径生成
- 支持 Angstrom/Bohr 单位
- 兼容 Gview 风格的原子选择语法

### calc_rmsd_xyz

基于 jxzou 的 RMSD 程序微调，增加 `.xyz` 支持。

- Kabsch/SVD 对齐
- 支持 `.xyz` / `.gjf`
- 支持原子范围选择的部分对齐
- 输出对齐后的 `mobile_new.xyz`

### neb_interpolator

- 插值/优化方法：LIC / LIIC / DM / NEB / NEB-LIIC
  - **LIC**：Cartesian linear interpolation（笛卡尔坐标线性插值）
  - **LIIC**：internal-coordinate interpolation（基于 DLC/Delocalized Internal Coordinates 的内坐标插值）
  - **DM**：distance-matrix interpolation（距离矩阵插值）
  - **NEB**：NEB 优化（默认用 LIC 初始化路径）
  - **NEB-LIIC**：NEB 优化（用 LIIC 初始化路径）
- LIIC 自动回退：LIIC → DM → LIC
- 插值/NEB 前可选端点对齐（默认开启；用 `calc_rmsd_xyz`）
- 支持外部引擎模式（NEB/NEB-LIIC 每个 cycle 调用外部程序返回梯度/力）
- **重要**：未指定外部引擎时，NEB/NEB-LIIC 会使用 `DistancePenaltyEngine`（距离矩阵惩罚的“虚拟力”）进行路径平滑示例，并不代表真实势能面上的 NEB。
- 自动查找 `calc_rmsd_xyz`：优先 PATH，其次同目录，最后回退到 `./calc_rmsd_xyz`

## 系统要求

### 编译器/工具链

- **C++17** 编译器（GCC/Clang/ICC 均可）
- Fortran 编译器（gfortran / ifort 等）
- CMake ≥ 3.10

### 库依赖

- BLAS/LAPACK（系统库或 Intel MKL）

### 可选：Intel MKL

CMake 提供 `USE_MKL` 选项：

- `-DUSE_MKL=ON`（默认）：若 `MKLROOT` 可用且库文件存在，则启用 MKL
- `-DUSE_MKL=OFF`：使用系统 BLAS/LAPACK

## 构建与安装

```bash
mkdir build
cd build
cmake ..
make -j4

# 可选：系统级安装
sudo make install
```

常用目标：

```bash
make test_all
make usage
```

## 目录结构

```
NEB_interpo/
├── CMakeLists.txt
├── README.md
├── src/
│   ├── xyzgeom.cpp
│   ├── calc_rmsd_xyz.f90
│   ├── neb_interpolator.cpp      # CLI
│   ├── neb_driver.{h,cpp}        # 共享库层 API（xyzgeom/neb_interpolator 统一调用）
│   ├── rmsd_align.{h,cpp}        # RMSD 对齐封装（查找 + 调用 calc_rmsd_xyz）
│   ├── internal_ic.{h,cpp}       # LIIC/DM 后端
│   ├── neb_engine.{h,cpp}
│   ├── neb_core.{h,cpp}
│   └── zmat/
│       ├── covalent_radii.{h,cpp}
│       ├── atom.h
│       ├── zmat.f90
│       └── zmatrix_generator.cpp
└── tests/
```

> `src/zmat/` 下除 `covalent_radii` 外还包含一些实验性/未接入主程序的代码文件（如 `zmat.f90`、`zmatrix_generator.cpp`）。

## 使用方法

### xyzgeom

交互模式：

```bash
./bin/xyzgeom molecule.xyz
```

核心菜单（节选）：

- 13. 与第二个 XYZ 文件对齐（RMSD，对齐工具为 `calc_rmsd_xyz`）
- 14. 与第二个 XYZ 进行路径生成

命令行模式示例：

```bash
# 打印结构
./bin/xyzgeom molecule.xyz --print

# 交换原子 1 和 3
./bin/xyzgeom molecule.xyz --swap 1,3

# 通过原子 1,2,3 定义的平面镜像
./bin/xyzgeom molecule.xyz --mirror 1-3
```

### calc_rmsd_xyz

```bash
# 基本对齐
./bin/calc_rmsd_xyz reference.xyz mobile.xyz

# 带原子选择的对齐
./bin/calc_rmsd_xyz reference.xyz mobile.xyz 1-10 1-10

# 混合格式
./bin/calc_rmsd_xyz reference.gjf mobile.xyz
```

### neb_interpolator

基本用法：

```bash
# 默认 NEB（LIC 初始化）
./bin/neb_interpolator initial.xyz final.xyz

# LIC 插值
./bin/neb_interpolator -n 10 -m lic initial.xyz final.xyz

# LIIC 插值（推荐作为高质量初始路径）
./bin/neb_interpolator -n 10 -m liic initial.xyz final.xyz

# DM 插值
./bin/neb_interpolator -m dm -n 8 initial.xyz final.xyz

# NEB-LIIC（用 LIIC 初始化再做 NEB）
./bin/neb_interpolator -n 10 -m neb-liic initial.xyz final.xyz

# 关闭端点对齐
./bin/neb_interpolator --no-align initial.xyz final.xyz

# 多帧输出（trajectory.xyz）
./bin/neb_interpolator -o multiframe initial.xyz final.xyz
```

#### 常用选项

- `-n, --nimages NUM`：中间图像数量（默认 5）
- `-m, --method METHOD`：`lic | liic | dm | neb | neb-liic`（默认 `neb`）
- `-p, --prefix PREFIX`：输出前缀
- `-o, --output MODE`：`separate | multiframe`
- `--no-align`：禁用端点对齐
- `-r, --rmsd-exec PATH`：手动指定 `calc_rmsd_xyz`

LIIC 参数（`-m liic` / `-m neb-liic`）：

- `--bond-factor F`
- `--fd-step H`
- `--ev-thresh T`
- `--liic-maxiter N`
- `--liic-tol T`
- `--liic-damp D`
- `--liic-max-step S`

DM 参数：

- `--dm-maxiter N`
- `--dm-step S`
- `--dm-tol T`
- `--dm-max-step S`
- `--no-dm-fallback`：禁用 LIIC → DM 回退（LIIC 失败直接回退到 LIC）

## 外部引擎模式（External engine）

外部引擎仅对 **NEB / NEB-LIIC** 生效：每个 NEB cycle 调用外部程序返回梯度/力。

如果你**没有**提供 `--engine-cmd`，则会自动回退到 `DistancePenaltyEngine`（虚拟势能/力），用于路径平滑/测试。
要在真实势能面上做 NEB，请提供外部引擎并输出梯度/力。

示例：

```bash
./bin/neb_interpolator -m neb -n 8 \
  --engine-cmd "python engine.py {in} {out}" \
  --engine-units AU \
  initial.xyz final.xyz
```

（详细格式说明请直接运行 `./bin/neb_interpolator --help` 查看。）

### 重要提示：`--engine-run-every`

`--engine-run-every N`（N>1）会在中间 cycle **复用上一次外部引擎返回的 forces/gradients**（即“陈旧梯度”），这会改变数值行为，可能显著影响收敛性。
建议仅用于调试/性能试验，不建议用于生产计算。

### 兼容性别名（可能弃用）

代码里保留了一些历史别名（供旧脚本使用），未来可能移除：

- `--engine-vector gradient|force`（选择外部引擎输出向量的语义：梯度或力；等价于设置/不设置 `--engine-output-is-force`）
- `--engine-every`（≈ `--engine-run-every`）
- `--engine-keep-files`（≈ `--engine-keep-cycle-files`）

这些别名默认不在 `--help` 中列出；以 `--help` 中列出的正式参数为准。

## 故障排除

### 1) MKL 相关

```bash
# 强制不用 MKL
cmake .. -DUSE_MKL=OFF -DBLA_VENDOR=Generic
```

### 2) 找不到 calc_rmsd_xyz

`neb_interpolator` 和 `xyzgeom` 会按以下顺序查找：

1. PATH 中的 `calc_rmsd_xyz`
2. 与当前可执行文件同目录的 `calc_rmsd_xyz`
3. 回退到 `./calc_rmsd_xyz`

也可以显式指定：

```bash
./bin/neb_interpolator -r /path/to/calc_rmsd_xyz initial.xyz final.xyz
```

### 3) LIIC 插值失败

- 调大 `--bond-factor`（如 1.30）以包含更多键
- 调小 `--ev-thresh`（如 5e-4）以保留更多 DLC 模式
- 增加 `--liic-maxiter` 迭代次数
- 若 LIIC 失败，默认会回退到 DM，再回退到 LIC

---

*最后更新：2026-02-24*
