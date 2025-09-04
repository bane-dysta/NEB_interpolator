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
- 内坐标线性插值（LIIC）
- 弹性带法（NEB）优化
- 插值前自动结构对齐
- 可配置收敛参数

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

```bash
# 基本NEB插值，生成5个图像
./bin/neb_interpolator initial.xyz final.xyz

# LIIC插值，生成10个图像
./bin/neb_interpolator -n 10 -m liic initial.xyz final.xyz

# 自定义参数
./bin/neb_interpolator -n 7 -m neb -s 0.0002 -c 0.005 -i 5000 initial.xyz final.xyz

# 不进行对齐
./bin/neb_interpolator --no-align initial.xyz final.xyz

# 多帧输出
./bin/neb_interpolator -o multiframe initial.xyz final.xyz
```

选项说明：
- `-n, --nimages NUM`: 中间图像数量（默认：5）
- `-m, --method METHOD`: 插值方法：liic或neb（默认：neb）
- `-p, --prefix PREFIX`: 输出文件名前缀
- `-o, --output MODE`: 输出模式：separate或multiframe
- `-s, --step STEP`: NEB步长（默认：0.0001）
- `-c, --conv THRESHOLD`: 收敛阈值（默认：0.01）
- `-i, --maxiter ITER`: 最大迭代次数（默认：10000）
- `-a, --align`: 启用结构对齐（默认）
- `--no-align`: 禁用结构对齐
- `-r, --rmsd-exec PATH`: calc_rmsd_xyz可执行文件路径

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

### 3. 批处理
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
2. **收敛性**: 调整步长和阈值以获得更快/更准确的结果
3. **大体系**: 使用LIIC作为初始猜测，然后用NEB精化
4. **对齐**: 预对齐结构以获得更好的插值质量

## 故障排除

### 常见问题

1. **MKL编译错误**
   ```bash
   # 改用系统BLAS/LAPACK
   cmake .. -DUSE_MKL=OFF
   ```

2. **找不到calc_rmsd_xyz**
   ```bash
   # 在xyzgeom中指定完整路径
   ./bin/xyzgeom molecule.xyz
   # 提示对齐时，确保calc_rmsd_xyz在PATH中或同一目录下
   ```

3. **NEB收敛问题**
   - 减小步长：`-s 0.00001`
   - 增加最大迭代次数：`-i 20000`
   - 先尝试LIIC：`-m liic`

4. **重排路径逆天**
   - 多因原子顺序、朝向不同导致
   - 使用gview为QST2设计的的connection功能找到合理的原子顺序
   - 使用rmsd模块进行旋转对齐
   - 使用主功能2交换位置异常的原子
   - 使用主功能5进行镜像操作

## 算法详情

### RMSD对齐
使用带SVD分解的Kabsch算法找到最小化两个结构间RMSD的最优旋转矩阵。

### NEB方法
实现爬坡图像NEB，包含：
- 沿反应坐标的弹簧力
- 距离保持的垂直力
- 自适应步长优化

### 距离矩阵保持
在插值过程中维持内坐标以防止不现实的中间结构。

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
*最后更新：25.9.4*
