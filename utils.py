# utils.py —— 公用工具函数
"""
公用工具函数模块

包含各模块共享的数学工具和辅助函数。
"""

import numpy as np


def rodrigues_rotation(v, k, theta):
    """
    使用Rodrigues旋转公式将向量v绕轴k旋转角度theta。
    
    参数：
        v: 要旋转的向量
        k: 旋转轴（单位向量）
        theta: 旋转角度（弧度）
    
    返回：
        旋转后的向量
    
    公式来源：Rodrigues' rotation formula
    v_rot = v*cos(theta) + (k×v)*sin(theta) + k*(k·v)*(1-cos(theta))
    """
    v = np.array(v, dtype=float)
    k = np.array(k, dtype=float)
    
    # 确保k是单位向量
    k = k / np.linalg.norm(k)
    
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    
    # Rodrigues公式
    v_rot = (v * cos_theta + 
             np.cross(k, v) * sin_theta + 
             k * np.dot(k, v) * (1 - cos_theta))
    
    return v_rot


def trilinear_interpolate(field, position, domain_size):
    """
    从三维网格场中进行三线性插值。
    
    参数：
        field: 三维numpy数组，浓度场
        position: [x, y, z]坐标（单位m）
        domain_size: 仿真域边长（单位m）
    
    返回：
        插值后的浓度值
    """
    grid_res = field.shape[0]
    grid_spacing = domain_size / (grid_res - 1)
    
    # 将物理坐标转换为网格坐标（包含边界）
    x, y, z = position
    gx = x / grid_spacing
    gy = y / grid_spacing
    gz = z / grid_spacing
    
    # 计算整数网格索引
    x0 = int(np.floor(gx))
    y0 = int(np.floor(gy))
    z0 = int(np.floor(gz))
    
    # 限制在有效范围内
    x0 = np.clip(x0, 0, grid_res - 2)
    y0 = np.clip(y0, 0, grid_res - 2)
    z0 = np.clip(z0, 0, grid_res - 2)
    
    x1 = x0 + 1
    y1 = y0 + 1
    z1 = z0 + 1
    
    # 计算插值权重
    xd = (gx - x0) / (x1 - x0) if x1 != x0 else 0
    yd = (gy - y0) / (y1 - y0) if y1 != y0 else 0
    zd = (gz - z0) / (z1 - z0) if z1 != z0 else 0
    
    # 三线性插值
    c000 = field[x0, y0, z0]
    c001 = field[x0, y0, z1]
    c010 = field[x0, y1, z0]
    c011 = field[x0, y1, z1]
    c100 = field[x1, y0, z0]
    c101 = field[x1, y0, z1]
    c110 = field[x1, y1, z0]
    c111 = field[x1, y1, z1]
    
    c00 = c000 * (1 - xd) + c100 * xd
    c01 = c001 * (1 - xd) + c101 * xd
    c10 = c010 * (1 - xd) + c110 * xd
    c11 = c011 * (1 - xd) + c111 * xd
    
    c0 = c00 * (1 - yd) + c10 * yd
    c1 = c01 * (1 - yd) + c11 * yd
    
    c = c0 * (1 - zd) + c1 * zd
    
    return c


def get_concentration_at_position(field, position, domain_size):
    """
    从三维网格场中获取指定位置的浓度值。
    
    优先使用最近邻插值（对于网格点），否则使用三线性插值。
    
    参数：
        field: 三维numpy数组，浓度场
        position: [x, y, z]坐标（单位m）
        domain_size: 仿真域边长（单位m）
    
    返回：
        该位置的浓度值
    """
    grid_res = field.shape[0]
    grid_spacing = domain_size / (grid_res - 1)
    
    # 将物理坐标转换为网格坐标
    x, y, z = position
    gx = x / grid_spacing
    gy = y / grid_spacing
    gz = z / grid_spacing
    
    # 检查是否正好在网格点上（允许小误差）
    ix = round(gx)
    iy = round(gy)
    iz = round(gz)
    
    if (abs(gx - ix) < 1e-10 and abs(gy - iy) < 1e-10 and abs(gz - iz) < 1e-10 and
        0 <= ix < grid_res and 0 <= iy < grid_res and 0 <= iz < grid_res):
        # 正好在网格点上，直接返回值
        return field[ix, iy, iz]
    
    # 否则使用三线性插值
    return trilinear_interpolate(field, position, domain_size)


def position_to_grid_index(position, domain_size, grid_resolution):
    """
    将物理坐标转换为网格索引。
    
    参数：
        position: [x, y, z]坐标（单位m）
        domain_size: 仿真域边长（单位m）
        grid_resolution: 网格分辨率
    
    返回：
        (i, j, k)网格索引
    """
    grid_spacing = domain_size / (grid_resolution - 1)
    x, y, z = position
    
    i = int(np.round(x / grid_spacing))
    j = int(np.round(y / grid_spacing))
    k = int(np.round(z / grid_spacing))
    
    # 限制在有效范围内
    i = np.clip(i, 0, grid_resolution - 1)
    j = np.clip(j, 0, grid_resolution - 1)
    k = np.clip(k, 0, grid_resolution - 1)
    
    return (i, j, k)
