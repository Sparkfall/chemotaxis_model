# environment.py —— 模块A：环境场模块
"""
环境场模块（模块A）

模拟三维空间中底物颗粒释放引诱物后形成的浓度场。
包含扩散方程求解、对流项（可开关）以及浓度查询功能。

物理模型：
- 稳态情况下，单个球形底物颗粒周围的引诱物浓度场由Fick扩散方程的球对称解给出：
  C(r) = C₀ * R₀ / r  (r ≥ R₀)
  其中 C₀ 是颗粒表面浓度，R₀ 是颗粒半径，r 是到颗粒中心的距离。
  来源：Crank J (1975) "The Mathematics of Diffusion", Oxford University Press, Chapter 6.
"""

import numpy as np
from scipy.ndimage import gaussian_filter
import config
from utils import trilinear_interpolate, position_to_grid_index, get_concentration_at_position


class Environment:
    """
    环境浓度场类。
    
    管理三维空间中的引诱物浓度分布，包括：
    - 底物颗粒产生的稳态浓度场
    - 扩散更新
    - 对流效应（重力条件下）
    - 营养消耗反馈
    """

    def __init__(self, cfg=config):
        """
        用config.py中的参数初始化浓度场。
        
        参数：
            cfg: 配置模块（默认为config）
        """
        self.cfg = cfg
        self.domain_size = cfg.DOMAIN_SIZE
        self.grid_res = cfg.GRID_RESOLUTION
        self.grid_spacing = self.domain_size / (self.grid_res - 1)
        
        # 创建网格坐标
        self.x = np.linspace(0, self.domain_size, self.grid_res)
        self.y = np.linspace(0, self.domain_size, self.grid_res)
        self.z = np.linspace(0, self.domain_size, self.grid_res)
        self.X, self.Y, self.Z = np.meshgrid(self.x, self.y, self.z, indexing='ij')
        
        # 底物颗粒参数
        self.substrate_radius = cfg.SUBSTRATE_RADIUS
        self.substrate_positions = cfg.SUBSTRATE_POSITIONS
        self.surface_conc = cfg.SUBSTRATE_SURFACE_CONC
        
        # 物理参数
        self.diffusion_coeff = cfg.GLUCOSE_DIFFUSION_COEFF
        self.gravity_on = cfg.GRAVITY_ON
        self.convection_velocity = cfg.CONVECTION_VELOCITY
        
        # 初始化浓度场
        self.concentration_field = self._initialize_concentration_field()

    def _initialize_concentration_field(self):
        """
        初始化浓度场。
        
        基于球对称扩散方程的解析解，计算每个网格点的浓度。
        多个颗粒的贡献线性叠加。
        
        返回：
            初始浓度场（三维numpy数组）
        """
        conc = np.zeros((self.grid_res, self.grid_res, self.grid_res))
        
        for pos in self.substrate_positions:
            px, py, pz = pos
            
            # 计算到颗粒中心的距离
            r = np.sqrt((self.X - px)**2 + (self.Y - py)**2 + (self.Z - pz)**2)
            
            # 球对称扩散解：C(r) = C₀ * R₀ / r  (r ≥ R₀)
            # 在颗粒内部(r < R₀)，浓度等于表面浓度
            # 来源：Crank J (1975) "The Mathematics of Diffusion", Chapter 6
            mask_outside = r >= self.substrate_radius
            mask_inside = r < self.substrate_radius
            
            # 颗粒外部浓度
            conc_particle = np.zeros_like(conc)
            conc_particle[mask_outside] = self.surface_conc * self.substrate_radius / r[mask_outside]
            conc_particle[mask_inside] = self.surface_conc
            
            # 线性叠加多个颗粒的贡献
            conc += conc_particle
        
        return conc

    def get_concentration(self, position):
        """
        查询指定位置的浓度。
        
        使用最近邻或三线性插值从网格中查询浓度值。
        
        参数：
            position: [x, y, z]坐标（单位m）
        
        返回：
            该位置的引诱物浓度（单位M）
        """
        return get_concentration_at_position(
            self.concentration_field,
            position,
            self.domain_size
        )

    def update(self, dt, consumption_field):
        """
        推进浓度场一个时间步。
        
        包含：
        1. 减去营养消耗
        2. 扩散更新（有限差分法）
        3. 对流效应（重力条件下）
        4. 底物颗粒表面浓度保持恒定
        
        参数：
            dt: 时间步长（单位s）
            consumption_field: 各网格点的消耗量（单位mol），形状与浓度场相同
        """
        # 1. 减去营养消耗（转换为浓度变化）
        # 消耗量(mol) / 网格体积(m³) = 浓度变化(M)
        grid_volume = self.grid_spacing ** 3
        conc_change = consumption_field / grid_volume
        self.concentration_field -= conc_change
        
        # 确保浓度非负
        self.concentration_field = np.maximum(self.concentration_field, 0)
        
        # 2. 扩散更新（有限差分法）
        # 使用7点Laplacian模板
        # C_new = C + D * dt * laplacian(C)
        # 来源：数值PDE标准方法
        self._diffusion_step(dt)
        
        # 3. 对流效应（仅重力条件下）
        # 自然对流会在底物颗粒周围形成上升羽流
        # 简化处理：添加垂直方向的对流项
        if self.gravity_on:
            self._convection_step(dt)
        
        # 4. 底物颗粒表面浓度保持恒定（Dirichlet边界条件）
        self._apply_substrate_boundary()

    def _diffusion_step(self, dt):
        """
        执行一步扩散更新。
        
        使用有限差分法求解扩散方程：
        ∂C/∂t = D * ∇²C
        
        采用7点Laplacian模板：
        ∇²C ≈ (C[i+1,j,k] + C[i-1,j,k] + C[i,j+1,k] + C[i,j-1,k] + 
               C[i,j,k+1] + C[i,j,k-1] - 6*C[i,j,k]) / h²
        """
        C = self.concentration_field
        h = self.grid_spacing
        D = self.diffusion_coeff
        
        # 计算Laplacian（内部点）
        laplacian = np.zeros_like(C)
        
        # x方向二阶导数
        laplacian[1:-1, :, :] += (C[2:, :, :] - 2*C[1:-1, :, :] + C[:-2, :, :]) / h**2
        # y方向二阶导数
        laplacian[:, 1:-1, :] += (C[:, 2:, :] - 2*C[:, 1:-1, :] + C[:, :-2, :]) / h**2
        # z方向二阶导数
        laplacian[:, :, 1:-1] += (C[:, :, 2:] - 2*C[:, :, 1:-1] + C[:, :, :-2]) / h**2
        
        # 更新浓度场（显式欧拉法）
        # 稳定性条件：dt < h² / (6D)
        self.concentration_field += D * dt * laplacian

    def _convection_step(self, dt):
        """
        执行一步对流更新。
        
        在重力条件下，溶解的糖使液体密度增大，形成自然对流。
        简化处理：添加垂直方向（z轴）的对流项。
        
        对流速度随距离底物颗粒的距离衰减。
        """
        C = self.concentration_field
        
        # 计算到最近底物颗粒的距离
        distance_to_substrate = np.full_like(C, np.inf)
        for pos in self.substrate_positions:
            px, py, pz = pos
            r = np.sqrt((self.X - px)**2 + (self.Y - py)**2 + (self.Z - pz)**2)
            distance_to_substrate = np.minimum(distance_to_substrate, r)
        
        # 对流速度随距离衰减（在颗粒附近最强）
        # v_conv = v0 * exp(-r / R0)
        v_conv = self.convection_velocity * np.exp(-distance_to_substrate / self.substrate_radius)
        
        # 对流项：-v * ∂C/∂z（向上流动，z方向）
        dC_dz = np.zeros_like(C)
        dC_dz[:, :, 1:-1] = (C[:, :, 2:] - C[:, :, :-2]) / (2 * self.grid_spacing)
        # 边界处理（一阶差分）
        dC_dz[:, :, 0] = (C[:, :, 1] - C[:, :, 0]) / self.grid_spacing
        dC_dz[:, :, -1] = (C[:, :, -1] - C[:, :, -2]) / self.grid_spacing
        
        # 更新浓度场
        self.concentration_field -= v_conv * dt * dC_dz

    def _apply_substrate_boundary(self):
        """
        应用底物颗粒边界条件。
        
        Dirichlet边界条件：底物颗粒表面浓度保持恒定。
        """
        for pos in self.substrate_positions:
            px, py, pz = pos
            
            # 计算到颗粒中心的距离
            r = np.sqrt((self.X - px)**2 + (self.Y - py)**2 + (self.Z - pz)**2)
            
            # 在颗粒内部和表面，浓度保持恒定
            mask = r <= self.substrate_radius
            self.concentration_field[mask] = self.surface_conc
