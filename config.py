# config.py —— 所有参数均需标注文献来源
"""
趋化性模拟全局参数配置文件
微重力下细菌趋化行为模拟模型

所有物理参数集中在此文件中定义，任何模块不得硬编码物理常数。
"""

import numpy as np

# === 仿真空间 ===
DOMAIN_SIZE = 2000e-6        # 仿真域边长，单位 m（2000 μm）
GRID_RESOLUTION = 100        # 网格分辨率（100×100×100）
DT = 0.1                     # 时间步长，单位 s
TOTAL_TIME = 30            # 总仿真时间，单位 s（1小时）

# === 细菌物理参数 ===
SWIM_SPEED = 25e-6           # 游速，单位 m/s（25 μm/s）
                             # 来源：Berg & Brown (1972) Nature 239:500-504
MEAN_RUN_TIME = 1.0          # 平均run时间，单位 s
                             # 来源：Berg & Brown (1972)
TUMBLE_ANGLE_MEAN = 62       # 翻转角度均值，单位 度
                             # 来源：Berg & Brown (1972)
TUMBLE_ANGLE_STD = 26        # 翻转角度标准差，单位 度
                             # 来源：Berg & Brown (1972)
TUMBLE_DURATION = 0.1        # 翻转持续时间，单位 s
CELL_RADIUS = 0.5e-6         # 细胞半径，单位 m
SEDIMENTATION_SPEED = 0.5e-6 # 重力沉降速度，单位 m/s
                             # 来源：Cluzel et al. (2000) Science 287:1652

# === 趋化信号通路参数 ===
K_D = 2.5e-6                 # Tar受体对天冬氨酸的解离常数，单位 M
                             # 来源：Sourjik & Berg (2002) PNAS 99:123-127
K_A = K_D / 35               # 半最大响应浓度（适应放大后）
                             # 来源：Sourjik & Berg (2002) 灵敏度放大35倍
ADAPTATION_TIME = 4.0        # 适应时间常数，单位 s
                             # 来源：Segall et al. (1986) PNAS 83:8987
CW_BIAS_BASELINE = 0.35      # 基线CW偏好（翻转概率基线）
                             # 来源：Cluzel et al. (2000) Science 287:1652
HILL_COEFF = 10.3            # CheY-P对马达的Hill系数
                             # 来源：Cluzel et al. (2000)
CHEYP_HALF = 3.1e-6          # CheY-P半最大浓度，单位 M
                             # 来源：Cluzel et al. (2000)

# === 信号通路ODE参数 ===
TAU_A = 0.1                  # 受体活性响应时间常数，单位 s
A_HALF = 0.5                 # Hill函数半最大活性

# === 葡萄糖相关参数 ===
# 注意：项目中趋化引诱物为葡萄糖，通过PTS系统介导
# 以下用天冬氨酸/Tar系统的参数作为初始近似
# 后续可根据文献调整为葡萄糖特异性参数
GLUCOSE_DIFFUSION_COEFF = 6.7e-10   # 葡萄糖在水中扩散系数，单位 m²/s
                                     # 来源：CRC Handbook of Chemistry and Physics
ATTRACTANT_TYPE = "glucose"
GLUCOSE_KM = 2e-6            # 葡萄糖转运Michaelis常数，单位 M
                             # 来源：Lendenmann et al. (1996)

# === 环境参数 ===
GRAVITY_ON = True            # 重力开关：True=1g，False=微重力
VISCOSITY = 1e-3             # 培养基粘度，单位 Pa·s（水的粘度）
TEMPERATURE = 310            # 温度，单位 K（37°C）

# === 底物颗粒参数 ===
SUBSTRATE_RADIUS = 50e-6     # 底物颗粒半径，单位 m（50 μm）
SUBSTRATE_POSITIONS = [       # 底物颗粒中心坐标（可配置多个）
    (1000e-6, 1000e-6, 1000e-6),  # 域中心放一个颗粒
]
SUBSTRATE_SURFACE_CONC = 10e-3  # 颗粒表面引诱物浓度，单位 M（10 mM）

# === 细菌群体参数 ===
N_BACTERIA = 500             # 模拟细菌数量
INITIAL_DISTRIBUTION = "uniform"  # 初始分布："uniform" 或 "cluster"

# === 代谢参数 ===
CONSUMPTION_RATE = 1e-17     # 单细胞葡萄糖消耗速率，单位 mol/s
                             # 来源：Lendenmann et al. (1996) Microbiology 142:1131
GROWTH_YIELD = 0.5           # 生长产率（无量纲）
DEATH_RATE = 0.0             # 死亡速率，单位 1/s（简化模型中可设为0）

# === 数值方法参数 ===
CONSUMPTION_THRESHOLD = 1e-7  # 营养耗尽阈值，单位 M
ENHANCED_KD_FACTOR = 0.2      # 增强模式下K_D的缩放因子（1/5）
CONVECTION_VELOCITY = 1e-5    # 自然对流速度量级，单位 m/s

# === 记录参数 ===
RECORD_INTERVAL = 1.0        # 数据记录间隔，单位 s
SUBSTRATE_PROXIMITY_RADIUS = 200e-6  # 底物附近区域半径，单位 m（200 μm）
