# movement.py —— 模块C：细菌运动模块
"""
细菌运动模块（模块C）

实现细菌的Run-and-Tumble随机游走运动。
包含重力沉降（可开关）和边界条件处理。

物理模型：

run阶段：
    细菌沿当前方向以恒定速度 v = SWIM_SPEED 做直线运动。
    x(t+dt) = x(t) + v * d * dt
    其中 d 是单位方向向量。

run持续时间：
    服从指数分布，每个时间步以概率 P_tumble * dt / MEAN_RUN_TIME 决定是否翻转。
    来源：Berg & Brown (1972) Nature 239:500-504.

tumble阶段：
    持续TUMBLE_DURATION秒，菌在原地旋转。
    新方向相对旧方向的偏转角 θ 服从正态分布 N(TUMBLE_ANGLE_MEAN, TUMBLE_ANGLE_STD²)。
    偏转的方位角 φ 在 [0, 2π) 上均匀分布。
    来源：Berg & Brown (1972); Saragosti et al. (2012) PLoS ONE 7:e35412.

重力沉降（仅GRAVITY_ON=True时）：
    每个时间步，位置的z分量额外减少 SEDIMENTATION_SPEED * dt。
    来源：估算值，基于Stokes沉降公式对1μm球体在水中的计算。
"""

import numpy as np
import config
from utils import rodrigues_rotation


class MovementEngine:
    """
    细菌运动引擎类。
    
    实现Run-and-Tumble随机游走运动，包括：
    - run阶段：直线运动
    - tumble阶段：原地旋转并改变方向
    - 重力沉降（可选）
    - 边界反射
    """

    def __init__(self, cfg=config):
        """
        初始化运动引擎。
        
        参数：
            cfg: 配置模块（默认为config）
        """
        self.cfg = cfg
        self.swim_speed = cfg.SWIM_SPEED
        self.mean_run_time = cfg.MEAN_RUN_TIME
        self.tumble_angle_mean = cfg.TUMBLE_ANGLE_MEAN
        self.tumble_angle_std = cfg.TUMBLE_ANGLE_STD
        self.tumble_duration = cfg.TUMBLE_DURATION
        self.gravity_on = cfg.GRAVITY_ON
        self.sedimentation_speed = cfg.SEDIMENTATION_SPEED
        self.domain_size = cfg.DOMAIN_SIZE

    def step(self, bacterium, dt):
        """
        执行一个细菌的运动更新。
        
        直接修改bacterium对象的position、direction、state等属性。
        包含重力沉降（如果GRAVITY_ON=True）和边界反射。
        
        参数：
            bacterium: Bacterium对象
            dt: 时间步长（单位s）
        """
        if bacterium.state == "run":
            self._run_step(bacterium, dt)
        elif bacterium.state == "tumble":
            self._tumble_step(bacterium, dt)
        
        # 边界反射处理
        self._apply_boundary_reflection(bacterium)

    def _run_step(self, bacterium, dt):
        """
        执行run阶段的运动。
        
        - 位置更新：沿当前方向直线运动
        - 重力沉降（如果开启）
        - 判断是否翻转
        """
        # 位置更新：position += SWIM_SPEED * direction * dt
        bacterium.position += self.swim_speed * bacterium.direction * dt
        
        # 重力沉降（仅z方向）
        if self.gravity_on:
            bacterium.position[2] -= self.sedimentation_speed * dt
        
        # 更新状态计时器
        bacterium.state_timer += dt
        
        # 判断是否翻转
        # 翻转概率：P_tumble * dt / MEAN_RUN_TIME
        # 来源：Berg & Brown (1972) 指数分布run时间
        flip_prob = bacterium.tumble_probability * dt / self.mean_run_time
        
        if np.random.random() < flip_prob:
            # 切换到tumble状态
            bacterium.state = "tumble"
            bacterium.state_timer = 0.0

    def _tumble_step(self, bacterium, dt):
        """
        执行tumble阶段的运动。
        
        - 位置不变（原地旋转）
        - 计时器增加
        - 完成后计算新方向
        """
        # 位置不变（原地旋转）
        bacterium.state_timer += dt
        
        # 检查是否完成tumble
        if bacterium.state_timer >= self.tumble_duration:
            # 计算新方向
            new_direction = self._calculate_new_direction(bacterium.direction)
            bacterium.direction = new_direction
            
            # 切换回run状态
            bacterium.state = "run"
            bacterium.state_timer = 0.0

    def _calculate_new_direction(self, old_direction):
        """
        计算tumble后的新方向。
        
        新方向相对旧方向的偏转角 θ 服从正态分布。
        偏转的方位角 φ 在 [0, 2π) 上均匀分布。
        
        使用Rodrigues旋转公式实现。
        
        参数：
            old_direction: 旧方向向量（单位向量）
        
        返回：
            新方向向量（单位向量）
        """
        # 采样偏转角（度转弧度）
        theta_deg = np.random.normal(self.tumble_angle_mean, self.tumble_angle_std)
        theta = np.radians(theta_deg)
        
        # 采样方位角（均匀分布）
        phi = np.random.uniform(0, 2 * np.pi)
        
        # 构建旋转轴（垂直于旧方向的随机方向）
        # 首先找到一个垂直于old_direction的向量
        if abs(old_direction[2]) < 0.9:
            # 如果z分量不大，用(0,0,1)叉乘得到垂直向量
            perp = np.cross(old_direction, [0, 0, 1])
        else:
            # 如果z分量大，用(1,0,0)叉乘
            perp = np.cross(old_direction, [1, 0, 0])
        
        perp = perp / np.linalg.norm(perp)
        
        # 根据方位角phi旋转perp，得到最终的旋转轴
        rotation_axis = (perp * np.cos(phi) + 
                        np.cross(old_direction, perp) * np.sin(phi))
        
        # 使用Rodrigues公式旋转
        new_direction = rodrigues_rotation(old_direction, rotation_axis, theta)
        
        # 确保是单位向量
        new_direction = new_direction / np.linalg.norm(new_direction)
        
        return new_direction

    def _apply_boundary_reflection(self, bacterium):
        """
        应用边界反射条件。
        
        如果细菌位置超出[0, DOMAIN_SIZE]³范围，执行镜面反射：
        - 如果x < 0：x = -x，direction_x取反
        - 如果x > DOMAIN_SIZE：x = 2*DOMAIN_SIZE - x，direction_x取反
        """
        for i in range(3):
            pos = bacterium.position[i]
            
            if pos < 0:
                # 左边界反射
                bacterium.position[i] = -pos
                bacterium.direction[i] = -bacterium.direction[i]
            elif pos > self.domain_size:
                # 右边界反射
                bacterium.position[i] = 2 * self.domain_size - pos
                bacterium.direction[i] = -bacterium.direction[i]
