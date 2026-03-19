# metabolism.py —— 模块D：代谢与消耗模块
"""
代谢与消耗模块（模块D）

处理单细胞的营养消耗和状态更新。

物理模型：

单细胞的葡萄糖消耗速率遵循Michaelis-Menten动力学：
    q = q_max * c / (c + K_m)

其中：
- q_max = CONSUMPTION_RATE（最大消耗速率）
- c 是局部葡萄糖浓度
- K_m ≈ 2 μM 是葡萄糖转运的Michaelis常数

来源：Lendenmann et al. (1996) Microbiology 142:1131-1140;
      Natarajan & Bhatt (2002) Biochem Eng J 11:193-199.

当局部浓度低于阈值（0.1 μM）时，消耗速率降为0（营养耗尽）。
"""

import numpy as np
import config


class MetabolismEngine:
    """
    细菌代谢引擎类。
    
    处理单细胞的营养消耗，遵循Michaelis-Menten动力学。
    本版本不模拟细胞生长和分裂，细菌数量保持恒定。
    """

    def __init__(self, cfg=config):
        """
        初始化代谢引擎。
        
        参数：
            cfg: 配置模块（默认为config）
        """
        self.cfg = cfg
        self.max_consumption_rate = cfg.CONSUMPTION_RATE
        self.k_m = cfg.GLUCOSE_KM
        self.consumption_threshold = cfg.CONSUMPTION_THRESHOLD

    def step(self, bacterium, local_concentration, dt):
        """
        处理一个细菌的营养消耗和状态更新。
        
        参数：
            bacterium: Bacterium对象
            local_concentration: 局部葡萄糖浓度（单位M）
            dt: 时间步长（单位s）
        
        返回：
            该细菌在此时间步消耗的营养量（单位mol）
        """
        # 如果细菌已死亡，不消耗营养
        if not bacterium.alive:
            return 0.0
        
        # 检查营养是否耗尽
        if local_concentration < self.consumption_threshold:
            return 0.0
        
        # Michaelis-Menten动力学计算消耗速率
        # q = q_max * c / (c + K_m)
        # 来源：Lendenmann et al. (1996)
        consumption_rate = (self.max_consumption_rate * local_concentration / 
                           (local_concentration + self.k_m))
        
        # 计算消耗量
        consumption = consumption_rate * dt
        
        # 更新细菌累计消耗量
        bacterium.nutrient_consumed += consumption
        
        return consumption
