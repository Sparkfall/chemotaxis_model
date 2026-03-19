# README.md - 趋化性模拟建模项目

## 项目简介

本项目是一个微重力下细菌趋化行为模拟模型，用于模拟大肠杆菌在正常重力（1g）和微重力（~0g）条件下的空间运动和底物降解行为。

## 项目结构

```
chemotaxis_model/
├── config.py           # 全局参数配置（所有物理参数集中定义）
├── bacteria.py         # 细菌数据结构定义
├── environment.py      # 模块A：环境场模块
├── signaling.py        # 模块B：细菌信号通路模块
├── movement.py         # 模块C：细菌运动模块
├── metabolism.py       # 模块D：代谢与消耗模块
├── main.py             # 模块E：主控模块
├── visualization.py    # 模块E：可视化模块
├── utils.py            # 公用工具函数
├── test_modules.py     # 单元测试脚本
├── __init__.py         # 包初始化文件
└── README.md           # 项目说明文档
```

## 模块说明

### 模块A：环境场模块（environment.py）

- 底物颗粒位置与浓度场计算
- 扩散方程求解（有限差分法）
- 重力对流项（可开关）

### 模块B：细菌信号通路模块（signaling.py）

- 简化Barkai-Leibler模型
- 输入：当前位置浓度 → 输出：翻转概率
- 适应动力学
- 支持三种模式：野生型/无趋化/趋化增强

### 模块C：细菌运动模块（movement.py）

- Run-and-tumble随机游走
- 重力沉降（可开关）
- 边界条件处理（镜面反射）

### 模块D：代谢与消耗模块（metabolism.py）

- 单细胞营养消耗（Michaelis-Menten动力学）
- 底物浓度场更新

### 模块E：主控与可视化模块（main.py + visualization.py）

- 初始化所有模块
- 时间步进循环
- 数据记录与输出
- 图表生成

## 使用方法

### 1. 运行完整仿真（六组对比实验）

```bash
python main.py --mode all --output results
```

这将运行以下六组条件：

1. 有趋化 + 1g
2. 无趋化 + 1g
3. 增强趋化 + 1g
4. 有趋化 + 微重力
5. 无趋化 + 微重力
6. 增强趋化 + 微重力

### 2. 运行单组条件仿真

```bash
# 野生型 + 重力
python main.py --mode wild_type --gravity --output results

# 无趋化 + 微重力
python main.py --mode no_chemotaxis --no-gravity --output results

# 增强趋化 + 重力
python main.py --mode enhanced --gravity --output results
```

### 3. 自定义参数

```bash
# 自定义仿真时间和细菌数量
python main.py --mode wild_type --gravity --time 1800 --bacteria 1000 --output results
```

### 4. 运行单元测试

```bash
python test_modules.py
```

## 参数配置

所有物理参数集中在 `config.py` 中定义，包括：

- **仿真空间参数**：域大小、网格分辨率、时间步长
- **细菌物理参数**：游速、run时间、翻转角度等
- **信号通路参数**：解离常数、适应时间、Hill系数等
- **环境参数**：重力开关、粘度、温度
- **底物颗粒参数**：半径、位置、表面浓度
- **代谢参数**：消耗速率、生长产率

修改 `config.py` 中的参数即可调整仿真行为。

## 输出结果

运行仿真后，将在输出目录生成以下图表：

1. `density_comparison.png` - 底物附近菌密度随时间变化对比
2. `consumption_comparison.png` - 累计营养消耗随时间变化对比
3. `concentration_field.png` - 最终时刻浓度场热力图
4. `final_dist_*.png` - 各条件的最终菌群空间分布图

## 物理模型参考

- Berg & Brown (1972) Nature 239:500-504 - 细菌运动参数
- Barkai & Leibler (1997) Nature 387:913-917 - 信号通路模型
- Cluzel et al. (2000) Science 287:1652-1655 - 翻转概率Hill函数
- Sourjik & Berg (2002) PNAS 99:123-127 - 受体参数
- Segall et al. (1986) PNAS 83:8987 - 适应时间
- Lendenmann et al. (1996) Microbiology 142:1131 - 葡萄糖消耗
- Crank J (1975) "The Mathematics of Diffusion" - 扩散方程

## 依赖项

- Python 3.7+
- NumPy
- Matplotlib
- SciPy

安装依赖：

```bash
pip install numpy matplotlib scipy
```

## 验证标准

模型完成后，需要通过以下验证：

1. **无趋化基线验证**：菌群空间分布保持均匀，均方位移符合扩散关系
2. **趋化聚集验证**：底物颗粒附近菌密度显著高于背景
3. **精确适应验证**：恒定浓度下受体活性在~20秒内回到基线（误差<5%）

## 作者

Chemotaxis Modeling Team

## 版本

1.0.0
