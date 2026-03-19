"""
趋化性模拟建模包

微重力下细菌趋化行为模拟模型

模块结构：
- config.py: 全局参数配置
- bacteria.py: 细菌数据结构定义
- environment.py: 环境场模块（模块A）
- signaling.py: 细菌信号通路模块（模块B）
- movement.py: 细菌运动模块（模块C）
- metabolism.py: 代谢与消耗模块（模块D）
- visualization.py: 可视化模块（模块E）
- main.py: 主控模块（模块E）
- utils.py: 公用工具函数
"""

__version__ = "1.0.0"
__author__ = "Chemotaxis Modeling Team"

from .config import *
from .bacteria import Bacterium
from .environment import Environment
from .signaling import SignalingPathway
from .movement import MovementEngine
from .metabolism import MetabolismEngine
from .visualization import Visualizer
