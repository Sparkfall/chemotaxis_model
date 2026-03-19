# signaling.py —— 模块B：细菌信号通路模块
"""
细菌信号通路模块（模块B）

实现简化的Barkai-Leibler趋化信号通路模型。
输入：当前位置浓度 → 输出：翻转概率

物理模型（简化版Barkai-Leibler）：
使用两变量ODE系统描述受体活性 a(t) 和甲基化水平 m(t)：

    da/dt = (1/τ_a) * [f(m, c) - a]
    dm/dt = (1/τ_m) * [a₀ - a]

其中：
- f(m, c) = m / (m + (1 + c/K_D)·(1-m)) 是受体活性对甲基化和浓度的响应函数
- τ_a ≈ 0.1 s 是受体活性的快速响应时间常数
- τ_m ≈ ADAPTATION_TIME (4 s) 是甲基化适应的慢时间常数
- a₀ = CW_BIAS_BASELINE 是适应目标活性
- c 是当前位置的引诱物浓度
- K_D 是受体解离常数

来源：Barkai & Leibler (1997) Nature 387:913-917; Tu et al. (2008) PNAS 105:14855-14860.

受体活性通过Hill函数映射为翻转概率：
    P_tumble = a^H / (a^H + a_{1/2}^H)

来源：Cluzel et al. (2000) Science 287:1652-1655.
"""

import numpy as np
import config


class SignalingPathway:
    """
    细菌趋化信号通路类。
    
    实现简化的Barkai-Leibler模型，支持三种趋化模式：
    - "wild_type": 野生型（正常趋化）
    - "no_chemotaxis": 无趋化（翻转概率恒定）
    - "enhanced": 趋化增强（灵敏度提高）
    """

    def __init__(self, cfg=config, chemotaxis_mode="wild_type"):
        """
        初始化信号通路。
        
        参数：
            cfg: 配置模块（默认为config）
            chemotaxis_mode: 趋化模式
                - "wild_type": 野生型
                - "no_chemotaxis": 无趋化
                - "enhanced": 趋化增强
        """
        self.cfg = cfg
        self.mode = chemotaxis_mode
        
        # 信号通路参数
        self.k_d = cfg.K_D
        self.tau_a = cfg.TAU_A  # 受体活性响应时间常数
        self.tau_m = cfg.ADAPTATION_TIME  # 适应时间常数
        self.a0 = cfg.CW_BIAS_BASELINE  # 适应目标活性
        self.hill_coeff = cfg.HILL_COEFF
        self.a_half = cfg.A_HALF  # Hill函数半最大活性
        
        # 增强模式：降低K_D（提高灵敏度）
        if self.mode == "enhanced":
            self.k_d *= cfg.ENHANCED_KD_FACTOR

    def update(self, bacterium, current_concentration, dt):
        """
        更新一个细菌的信号通路状态和翻转概率。
        
        直接修改bacterium对象的tumble_probability等属性。
        
        参数：
            bacterium: Bacterium对象
            current_concentration: 当前位置的引诱物浓度（单位M）
            dt: 时间步长（单位s）
        """
        # 无趋化模式：翻转概率恒定
        if self.mode == "no_chemotaxis":
            bacterium.set_tumble_probability(self.cfg.CW_BIAS_BASELINE)
            return
        
        # 获取当前状态
        a = bacterium.receptor_activity
        m = bacterium.methylation
        c = current_concentration
        
        # 计算受体响应对甲基化和浓度的响应函数 f(m, c)
        # f(m, c) = m / (m + (1 + c/K_D) * (1 - m))
        # 来源：Barkai-Leibler模型
        denominator = m + (1 + c / self.k_d) * (1 - m)
        if denominator > 0:
            f_m_c = m / denominator
        else:
            f_m_c = 0.5
        
        # 用欧拉法更新ODE
        # da/dt = (f(m,c) - a) / tau_a
        # dm/dt = (a0 - a) / tau_m
        da_dt = (f_m_c - a) / self.tau_a
        dm_dt = (self.a0 - a) / self.tau_m
        
        a_new = a + da_dt * dt
        m_new = m + dm_dt * dt
        
        # 裁剪到有效范围（防止数值溢出）
        a_new = np.clip(a_new, 0.01, 0.99)
        m_new = np.clip(m_new, 0.01, 0.99)
        
        # 更新细菌状态
        bacterium.receptor_activity = a_new
        bacterium.methylation = m_new
        
        # 将受体活性通过Hill函数映射为翻转概率
        # P = a^H / (a^H + a_half^H)
        # 来源：Cluzel et al. (2000) Science 287:1652-1655
        a_pow = a_new ** self.hill_coeff
        a_half_pow = self.a_half ** self.hill_coeff
        tumble_prob = a_pow / (a_pow + a_half_pow)
        
        bacterium.set_tumble_probability(tumble_prob)
