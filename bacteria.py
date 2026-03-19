# bacteria.py —— 单个细菌的数据结构
"""
细菌数据结构定义

单个细菌的状态。所有模块通过读写这个对象来交换数据。
"""

import numpy as np
from config import CW_BIAS_BASELINE


def random_unit_vector():
    """生成三维空间中均匀分布的随机单位方向向量。"""
    vec = np.random.randn(3)
    return vec / np.linalg.norm(vec)


class Bacterium:
    """单个细菌的状态。所有模块通过读写这个对象来交换数据。"""

    def __init__(self, position, bacterium_id):
        # === 空间状态 ===
        self.id = bacterium_id
        self.position = np.array(position, dtype=float)  # [x, y, z]，单位 m
        self.direction = random_unit_vector()              # 单位方向向量

        # === 运动状态 ===
        self.state = "run"           # "run" 或 "tumble"
        self.state_timer = 0.0       # 当前状态已持续时间，单位 s

        # === 信号通路状态 ===
        self.receptor_activity = 0.5  # 受体活性 a(t)，范围 [0, 1]
        self.methylation = 0.5        # 甲基化水平 m(t)，范围 [0, 1]
        self.conc_history = []        # 最近若干步的浓度记录
        self.tumble_probability = CW_BIAS_BASELINE  # 当前翻转概率

        # === 代谢状态 ===
        self.alive = True
        self.nutrient_consumed = 0.0  # 累计消耗营养量

    def get_position(self):
        """返回位置坐标，供模块A查询浓度用。"""
        return self.position.copy()

    def set_tumble_probability(self, p):
        """由模块B调用，设置翻转概率。"""
        self.tumble_probability = np.clip(p, 0.0, 1.0)
