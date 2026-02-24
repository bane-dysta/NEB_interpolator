*README由Claude 4.0生成，可以查看帮助，彩红屁勿信*

# 分子几何工具

用于分子结构分析、对齐和反应路径插值的综合工具套件。

## 概述

本项目提供三个集成的分子几何操作工具：
- **xyzgeom**: 增强版分子几何分析器，具有14+种分析功能
- **calc_rmsd_xyz**: 使用Fortran/LAPACK的高精度RMSD计算器和结构对齐工具
- **neb_interpolator**: 用于反应路径的NEB（弹性带法）和LIIC插值

## 功能特性

### XYZ几何分析器 (xyzgeom)
- 计算距离、角度和二面角
- 原子操作（交换、镜像、导出）
- 使用RMSD进行结构对齐
- NEB插值集成
- 支持埃（Angstrom）和玻尔（Bohr）单位
- 兼容Gview的原子选择语法

### RMSD计算器 (calc_rmsd_xyz)

RMSD计算器是基于jxzou的[RMSD程序](https://gitlab.com/jxzou/rmsd)的微调版本，增加了对xyz文件的支持。

- 使用SVD分解进行精确结构对齐
- 支持.xyz和.gjf格式
- 支持原子范围选择的部分结构对齐
- 最优旋转矩阵计算

### NEB插值器
- **多种插值方法**：
  - **LIIC** (Linear Internal Coordinate Interpolation): 线性内坐标插值，快速简单
  - **IIC** (Internal Coordinate Interpolation): 基于DLC（Delocalized Internal Coordinates）的高质量内坐标插值
  - **DM** (Distance Matrix): 距离矩阵插值，适用于复杂体系
  - **NEB**: 传统弹性带法优化
  - **NEB-IIC**: 使用IIC初始化的NEB，结合两者优势
- **智能回退机制**：IIC失败时自动回退到DM，DM失败时回退到LIIC
- 插值前自动结构对齐（使用Fortran RMSD）
- 可配置收敛参数和算法选项
- 自动查找calc_rmsd_xyz可执行文件（PATH或同目录）

## 系统要求

### 编译器要求
- 兼容C++11的编译器（GCC 4.8+, Intel C++ 2013+）
- Fortran编译器（gfortran 4.8+, ifort 2013+）
- CMake 3.10+

### 库依赖
- BLAS/LAPACK（系统库或Intel MKL）
- POSIX系统用于文件操作

### 可选依赖
- Intel MKL用于优化的线性代数运算

## 安装

### 1. 克隆或下载项目
```bash
git clone <repository-url>
cd MolecularGeometryTools
```

### 2. 目录结构
```
MolecularGeometryTools/
├── CMakeLists.txt
├── README.md
├── src/
│   ├── xyzgeom.cpp           # 主几何分析器
│   ├── calc_rmsd_xyz.f90     # Fortran RMSD计算器
│   ├── neb_interpolator.cpp  # 独立NEB工具
│   └── neb_interpolator.h    # NEB库头文件
└── build/                    # 构建目录（需要创建）
```

### 3. 构建项目
```bash
mkdir build
cd build
cmake ..
make -j4
```

### 4. （可选）系统级安装
```bash
sudo make install
```

### 5. 验证安装
```bash
make test_all  # 运行所有测试
make usage     # 显示使用信息
```

## 使用方法

### XYZ几何分析器

#### 交互模式
```bash
./bin/xyzgeom molecule.xyz
```

类Multiwfn的交互菜单选项：
- 1.  选择2个原子，计算距离与向量
- 2.  选择2个原子，交换其索引位置，更新内存(命令行模式直接输出新xyz)
- 3.  选择3个原子，计算键角并拟合平面法向量
- 4.  选择两组3个原子，计算两个平面夹角
- 5.  选择3个原子定义平面，将分子关于平面进行对称，更新内存(命令行模式直接输出新xyz)
- 6.  选择4个原子，计算二面角
- 7.  选择n个原子，计算几何中心
- 8.  选择n个原子，导出新xyz文件
- 9.  选择1个原子，找到半径n埃内所有原子，输出新xyz文件与Gview风格索引
- 10. 将内存中的数据导出xyz文件
- 11. 更换内存中的xyz文件
- 12. 主界面切换计算数学量时输出单位Bohr/Angstrom
- 13. 与第二个XYZ文件对齐（RMSD）
- 14. 与第二个XYZ进行NEB插值

#### 命令行模式
```bash
# 打印结构
./bin/xyzgeom molecule.xyz --print

# 交换原子1和3
./bin/xyzgeom molecule.xyz --swap 1,3

# 通过原子1,2,3定义的平面镜像
./bin/xyzgeom molecule.xyz --mirror 1-3
```

### RMSD计算器

本程序是基于jxzou的[RMSD程序](https://gitlab.com/jxzou/rmsd)的微调版本，增加了对xyz文件的支持。

```bash
# 基本对齐
./bin/calc_rmsd_xyz reference.xyz mobile.xyz

# 带原子选择的对齐
./bin/calc_rmsd_xyz reference.xyz mobile.xyz 1-10 1-10

# 混合格式支持
./bin/calc_rmsd_xyz reference.gjf mobile.xyz
```

输出：创建对齐结构的`mobile_new.xyz`文件

### NEB插值器

#### 基本用法

```bash
# 基本NEB插值，生成5个图像
./bin/neb_interpolator initial.xyz final.xyz

# LIIC插值，生成10个图像（快速简单）
./bin/neb_interpolator -n 10 -m liic initial.xyz final.xyz

# IIC插值（高质量内坐标插值）
./bin/neb_interpolator -n 10 -m iic initial.xyz final.xyz

# 距离矩阵插值（适用于复杂体系）
./bin/neb_interpolator -m dm -n 8 initial.xyz final.xyz

# NEB-IIC（使用IIC初始化的NEB，推荐用于复杂反应）
./bin/neb_interpolator -n 10 -m neb-iic initial.xyz final.xyz

# 自定义NEB参数
./bin/neb_interpolator -n 7 -m neb -s 0.0002 -c 0.005 -i 5000 initial.xyz final.xyz

# 不进行对齐
./bin/neb_interpolator --no-align initial.xyz final.xyz

# 多帧输出（用于VMD/PyMOL可视化）
./bin/neb_interpolator -o multiframe initial.xyz final.xyz
```

#### 选项说明

**基本选项**：
- `-n, --nimages NUM`: 中间图像数量（默认：5）
- `-m, --method METHOD`: 插值方法：`liic` | `iic` | `dm` | `neb` | `neb-iic`（默认：`neb`）
- `-p, --prefix PREFIX`: 输出文件名前缀
- `-o, --output MODE`: 输出模式：`separate`（分别文件）或 `multiframe`（单文件轨迹）
- `-s, --step STEP`: NEB步长（默认：0.0001）
- `-c, --conv THRESHOLD`: 收敛阈值（默认：0.01）
- `-i, --maxiter ITER`: 最大迭代次数（默认：10000）
- `-a, --align`: 启用结构对齐（默认）
- `--no-align`: 禁用结构对齐
- `-r, --rmsd-exec PATH`: calc_rmsd_xyz可执行文件路径（默认自动查找）

**IIC选项**（用于 `-m iic` 或 `-m neb-iic`）：
- `--bond-factor F`: 键截断因子，用于确定键连接（默认：1.25）
- `--fd-step H`: B矩阵有限差分步长（默认：1e-4 Å）
- `--ev-thresh T`: DLC坐标选择的本征值阈值（默认：1e-3）
- `--iic-maxiter N`: 最大回变换迭代次数（默认：50）
- `--iic-tol T`: 原始残差RMS容差（默认：1e-4）
- `--iic-damp D`: 本征值阻尼（默认：1e-8）
- `--iic-max-step S`: 每次迭代最大笛卡尔位移（默认：0.20 Å）

**DM选项**（用于 `-m dm` 或作为IIC回退）：
- `--dm-maxiter N`: 最大DM迭代次数（默认：800）
- `--dm-step S`: DM梯度下降步长（默认：5e-3）
- `--dm-tol T`: RMS距离误差容差（默认：1e-3 Å）
- `--dm-max-step S`: 每次迭代最大笛卡尔位移（默认：0.20 Å）
- `--no-dm-fallback`: 禁用DM回退（IIC失败时直接回退到LIIC）

**外部引擎模式（External engine）**（用于 `-m neb` 或 `-m neb-iic`）：

当你希望每个NEB cycle 调用外部程序（量化化学/机器学习势能等）来计算梯度/力时，使用外部引擎模式。

- `--engine-cmd CMD`: 启用外部引擎并在每个cycle运行CMD。
  - 若 CMD 包含占位符 `{in}`/`{out}`/`{cycle}`，会自动替换。
  - 若不包含占位符，程序会在命令末尾追加：`<infile> <outfile>`。
- `--engine-in FILE`: 引擎输入文件名（默认：`neb_engine_in.dat`）
- `--engine-out FILE`: 引擎输出文件名（默认：`neb_engine_out.dat`）
- `--engine-units U`: 写入文件头的单位标签（默认：`Angstrom`；也可设为 `AU`）
- `--engine-vector TYPE`: 引擎输出向量类型：`gradient` 或 `force`（默认：`gradient`）
- `--engine-spring K`: 外部引擎NEB的弹簧常数（默认：1.0）
- `--engine-every N`: 每N个cycle调用一次引擎（默认：1，即每cycle）
- `--engine-keep-files`: 保留每个cycle的I/O文件（文件名会加 `_cycle####` 后缀）

**外部引擎I/O文件格式**（与扩展XYZ类似，按image分块）：

1) 程序写给引擎的输入（`type=xyz`）：

```text
3
image=1 units=AU type=xyz cycle=0
O   0.00000000  0.00000000  0.00000000
H   0.75860200  0.00000000  0.50428400
H   0.75860200  0.00000000 -0.50428400
...
```

2) 引擎返回给程序的输出（`type=gradient` 或 `type=force`）：

```text
3
image=1 units=AU type=gradient
O   0.01  -0.02  0.00
H  -0.01   0.01  0.00
H   0.00   0.01  0.00
...
```

> 注意：若 `--engine-vector gradient`，程序会将梯度取负号作为力；若 `--engine-vector force` 则直接当作力使用。

**引擎可用环境变量**（方便脚本读取，无需解析命令行参数）：
- `NEB_ENGINE_IN`：输入文件路径
- `NEB_ENGINE_OUT`：输出文件路径
- `NEB_CYCLE`：当前cycle编号
- `NEB_NIMAGES`：中间图像数
- `NEB_NATOMS`：原子数

**示例**：

```bash
./bin/neb_interpolator -m neb -n 8 \
  --engine-cmd "python engine.py {in} {out}" \
  --engine-units AU \
  initial.xyz final.xyz
```

#### 方法选择建议

- **LIIC** (`-m liic`): 快速简单，适合简单体系或快速预览
- **IIC** (`-m iic`): 高质量内坐标插值，适合大多数反应路径，自动回退到DM/LIIC
- **DM** (`-m dm`): 距离矩阵方法，适合复杂体系或IIC失败的情况
- **NEB** (`-m neb`): 传统弹性带法，需要好的初始猜测
- **NEB-IIC** (`-m neb-iic`): **推荐**，使用IIC生成高质量初始路径，然后用NEB优化

## 示例工作流程

### 1. 反应路径分析
```bash
# 加载初始结构
./bin/xyzgeom reactant.xyz

# 在交互模式中：
# - 选择选项13与product.xyz对齐
# - 选择选项14进行NEB插值
# - 输入图像数量（如10）
# - 选择NEB方法
# 输出：neb_00.xyz到neb_11.xyz
```

### 2. 结构比较
```bash
# 计算两个构象间的RMSD
./bin/calc_rmsd_xyz conformer1.xyz conformer2.xyz

# 对齐和分析
./bin/xyzgeom conformer1.xyz
# 选择选项13，输入conformer2.xyz
# 选择选项1测量特定距离
```

### 3. 高质量反应路径插值
```bash
# 使用IIC方法生成高质量初始路径
./bin/neb_interpolator -n 10 -m iic reactant.xyz product.xyz

# 使用NEB-IIC结合方法
./bin/neb_interpolator -n 10 -m neb-iic \
    --bond-factor 1.2 --ev-thresh 1e-3 \
    reactant.xyz product.xyz

# 复杂体系使用DM方法
./bin/neb_interpolator -m dm -n 8 \
    --dm-maxiter 1000 --dm-tol 5e-4 \
    complex_initial.xyz complex_final.xyz
```

### 4. 批处理
```bash
#!/bin/bash
# 处理多个结构对
for i in {1..10}; do
    ./bin/neb_interpolator structure${i}_initial.xyz structure${i}_final.xyz \
        -n 5 -p pathway${i}_ -m neb
done
```

## 原子选择语法

工具支持兼容Gview的选择语法：
- 单个原子：`1,3,5`
- 范围：`1-10`
- 混合：`1-3,5,7-10`
- 支持全角逗号：`1，2，3`

## 输出文件

- **XYZ文件**: 标准XYZ格式，包含原子坐标
- **_new后缀**: calc_rmsd_xyz产生的对齐/修改结构
- **编号序列**: 00.xyz（初始）、01.xyz...N.xyz（中间）、N+1.xyz（最终）
- **trajectory.xyz**: 用于VMD/PyMOL可视化的多帧XYZ

## 性能提示

1. **使用Intel MKL**: 构建前设置`MKLROOT`环境变量
2. **方法选择**:
   - 简单体系：使用`-m liic`快速预览
   - 一般反应：使用`-m iic`获得高质量路径
   - 复杂反应：使用`-m neb-iic`结合IIC和NEB的优势
   - 大体系或IIC失败：使用`-m dm`
3. **IIC参数调优**:
   - 增大`--bond-factor`（如1.3）以包含更多键连接
   - 减小`--ev-thresh`（如5e-4）以包含更多DLC坐标
   - 增加`--iic-maxiter`（如100）以提高收敛性
4. **收敛性**: 调整步长和阈值以获得更快/更准确的结果
5. **对齐**: 预对齐结构以获得更好的插值质量（默认启用）
6. **初始路径**: IIC生成的初始路径通常比LIIC质量更高，适合后续NEB优化

## 故障排除

### 常见问题

1. **MKL编译错误**
   ```bash
   # 改用系统BLAS/LAPACK
   cmake .. -DUSE_MKL=OFF
   ```

2. **找不到calc_rmsd_xyz**
   ```bash
   # neb_interpolator会自动查找：
   # 1. PATH中的calc_rmsd_xyz
   # 2. neb_interpolator同目录下的calc_rmsd_xyz
   # 3. 如果都找不到，使用默认路径./calc_rmsd_xyz
   
   # 也可以手动指定路径
   ./bin/neb_interpolator -r /path/to/calc_rmsd_xyz initial.xyz final.xyz
   ```

3. **NEB收敛问题**
   - 减小步长：`-s 0.00001`
   - 增加最大迭代次数：`-i 20000`
   - 使用IIC生成更好的初始路径：`-m neb-iic`
   - 或先尝试IIC：`-m iic`，然后手动优化

4. **IIC插值失败或质量不佳**
   - 调整键因子：`--bond-factor 1.3`（增大以包含更多键）
   - 减小有限差分步长：`--fd-step 5e-5`
   - 调整本征值阈值：`--ev-thresh 5e-4`（减小以包含更多DLC坐标）
   - 增加回变换迭代：`--iic-maxiter 100`
   - 如果IIC失败，程序会自动回退到DM或LIIC

5. **重排路径逆天**
   - 多因原子顺序、朝向不同导致
   - 使用gview为QST2设计的的connection功能找到合理的原子顺序
   - 使用rmsd模块进行旋转对齐
   - 使用主功能2交换位置异常的原子
   - 使用主功能5进行镜像操作

## 算法详情

### RMSD对齐
使用带SVD分解的Kabsch算法找到最小化两个结构间RMSD的最优旋转矩阵。

### NEB（弹性带法）算法细节

#### 算法原理
弹性带法（Nudged Elastic Band, NEB）是一种寻找反应路径的经典方法，通过在初始和最终结构之间插入一系列中间图像，并用弹簧连接这些图像形成"弹性带"，通过优化使其收敛到最小能量路径。

#### 实现细节

**1. 初始化**
- 使用线性内坐标插值（LIIC）生成初始中间图像
- 中间图像数量可配置（默认5个）
- 每个中间图像的位置通过线性插值确定：
  ```
  r_img = r_initial + factor × (r_final - r_initial)
  其中 factor = (img + 1) / (num_images + 1)
  ```

**2. 力的计算**
每个中间图像受到两种力的作用：

**弹簧力（Spring Force）**
- 沿着反应路径切线方向的弹簧力，保持图像间的均匀分布
- 计算公式：`F_spring = k × [(r_next - r_current) - (r_current - r_prev)]`
- 弹簧常数 k = 1.0
- 只保留沿切线方向的分量，避免影响垂直于路径的优化

**垂直力（Perpendicular Force）**
- 基于距离矩阵保持的内坐标约束力
- 目标：维持合理的分子几何，防止产生不现实的中间结构
- 每个图像的目标距离矩阵通过线性插值确定：
  ```
  d_target(i,j) = d_initial(i,j) + factor × [d_final(i,j) - d_initial(i,j)]
  ```
- 垂直力移除沿切线方向的分量，确保不干扰反应路径优化

**3. 切线向量计算**
- 对每个中间图像计算局部切线向量
- 切线方向定义为相邻图像间的位移差：
  ```
  tangent = (r_next - r_current) - (r_current - r_prev)
  ```
- 归一化处理确保数值稳定性

**4. 位置更新**
- 使用简单的欧拉积分更新原子坐标：
  ```
  r_new = r_old + step_size × F_total
  ```
- 默认步长：0.0001
- 自适应步长可根据收敛情况调整

**5. 收敛判据**
- 监控所有图像上所有原子的最大力
- 当最大力小于收敛阈值时停止优化（默认：0.01）
- 最大迭代次数限制（默认：10000）

#### 算法特点

**优势：**
- 能够找到连接反应物和产物的最小能量路径
- 通过距离矩阵约束保持合理的分子几何
- 自动对齐功能提高插值质量
- 支持大分子体系

**局限性：**
- 需要预先知道反应的初始和最终状态
- 对于复杂反应可能需要调整参数
- 不包含真实的势能面信息（仅基于几何约束）

**参数调优建议：**
- 收敛困难时：减小步长（如 0.00001）
- 需要更高精度时：降低收敛阈值（如 0.005）
- 复杂体系：增加中间图像数量
- 快速预览：先用LIIC方法

### 内坐标插值方法

#### LIIC (Linear Internal Coordinate Interpolation)
线性内坐标插值，直接在内坐标空间进行线性插值，然后转换回笛卡尔坐标。快速但可能产生不合理的中间结构。

#### IIC (Internal Coordinate Interpolation)
基于DLC（Delocalized Internal Coordinates）的高质量内坐标插值方法，灵感来自GSM（Growing String Method）：

1. **键连接识别**：基于共价半径和键因子自动识别分子键连接
2. **原始内坐标构建**：生成键长、键角、二面角等原始内坐标
3. **DLC基构建**：
   - 通过有限差分计算原始B矩阵（Bp）
   - 构建Gram矩阵 G = Bp × Bp^T
   - 通过Jacobi对角化获得DLC基（本征值大于阈值的本征向量）
4. **内坐标插值**：在DLC空间中线性插值
5. **回变换**：迭代优化将内坐标变化转换回笛卡尔坐标
   - 使用阻尼最小二乘：dx = Bp^T × Σ(v_i × (v_i·dq) / (λ_i + damp))
   - 每次迭代后进行Kabsch对齐以消除平移/旋转
   - 支持线搜索以稳定收敛

**特点**：
- 保持键长、键角、二面角的合理性
- 避免原子重叠和不合理的中间结构
- 自动处理二面角周期性（包装到[-π, π]）
- 失败时自动回退到DM或LIIC

#### DM (Distance Matrix Interpolation)
距离矩阵插值方法：
- 在距离矩阵空间进行插值
- 使用梯度下降优化满足目标距离矩阵
- 适合复杂体系或IIC失败的情况

#### 回退机制
IIC方法具有智能回退机制：
1. 首先尝试IIC插值
2. 如果IIC失败，自动回退到DM插值
3. 如果DM也失败，回退到简单的LIIC
4. 可通过`--no-dm-fallback`禁用DM回退，直接使用LIIC

### 距离矩阵保持
在插值过程中维持内坐标以防止不现实的中间结构。通过约束原子间距离的变化，确保中间图像保持合理的化学键长和几何构型。

## 引用

如果您在研究中使用这些工具，请引用：
```
分子几何工具
[姓名/机构]
[年份]
```

## 贡献

欢迎贡献！请在项目仓库中提交拉取请求或问题。

---
*最后更新：26.2.11*
