# visualization.py —— 模块E（可视化部分）
"""
可视化模块

生成仿真结果的各种图表和可视化输出。
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import config


class Visualizer:
    """
    可视化类。
    
    生成以下可视化输出：
    1. 菌群空间分布的三维散点图（或二维投影）
    2. 底物颗粒周围菌密度随时间的变化曲线
    3. 累计营养消耗随时间的变化曲线
    4. 最终时刻的浓度场热力图（截面）
    """

    def __init__(self, cfg=config):
        """
        初始化可视化器。
        
        参数：
            cfg: 配置模块（默认为config）
        """
        self.cfg = cfg
        self.domain_size = cfg.DOMAIN_SIZE
        self.substrate_positions = cfg.SUBSTRATE_POSITIONS
        self.substrate_radius = cfg.SUBSTRATE_RADIUS
        self.proximity_radius = cfg.SUBSTRATE_PROXIMITY_RADIUS

    def plot_spatial_distribution(self, bacteria, time, save_path=None):
        """
        绘制菌群空间分布的三维散点图。
        
        参数：
            bacteria: 细菌列表
            time: 当前时间（单位s）
            save_path: 保存路径（如果为None则显示）
        """
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        
        # 提取活菌位置
        positions = np.array([b.position for b in bacteria if b.alive])
        
        if len(positions) > 0:
            ax.scatter(positions[:, 0] * 1e6,  # 转换为μm
                      positions[:, 1] * 1e6,
                      positions[:, 2] * 1e6,
                      c='blue', s=1, alpha=0.5)
        
        # 绘制底物颗粒
        for pos in self.substrate_positions:
            u = np.linspace(0, 2 * np.pi, 20)
            v = np.linspace(0, np.pi, 20)
            x = pos[0] * 1e6 + self.substrate_radius * 1e6 * np.outer(np.cos(u), np.sin(v))
            y = pos[1] * 1e6 + self.substrate_radius * 1e6 * np.outer(np.sin(u), np.sin(v))
            z = pos[2] * 1e6 + self.substrate_radius * 1e6 * np.outer(np.ones(np.size(u)), np.cos(v))
            ax.plot_surface(x, y, z, color='red', alpha=0.3)
        
        ax.set_xlabel('X (μm)')
        ax.set_ylabel('Y (μm)')
        ax.set_zlabel('Z (μm)')
        ax.set_title(f'Bacteria Spatial Distribution (t = {time:.0f} s)')
        
        # 设置坐标轴范围
        ax.set_xlim(0, self.domain_size * 1e6)
        ax.set_ylim(0, self.domain_size * 1e6)
        ax.set_zlim(0, self.domain_size * 1e6)
        
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            plt.close()
        else:
            plt.show()

    def plot_xy_projection(self, bacteria, time, save_path=None):
        """
        绘制菌群空间分布的XY平面投影。
        
        参数：
            bacteria: 细菌列表
            time: 当前时间（单位s）
            save_path: 保存路径（如果为None则显示）
        """
        fig, ax = plt.subplots(figsize=(8, 8))
        
        # 提取活菌位置
        positions = np.array([b.position for b in bacteria if b.alive])
        
        if len(positions) > 0:
            ax.scatter(positions[:, 0] * 1e6,  # 转换为μm
                      positions[:, 1] * 1e6,
                      c='blue', s=2, alpha=0.5)
        
        # 绘制底物颗粒投影（圆）
        for pos in self.substrate_positions:
            circle = plt.Circle((pos[0] * 1e6, pos[1] * 1e6),
                               self.substrate_radius * 1e6,
                               color='red', fill=False, linewidth=2)
            ax.add_patch(circle)
        
        ax.set_xlabel('X (μm)')
        ax.set_ylabel('Y (μm)')
        ax.set_title(f'Bacteria XY Projection (t = {time:.0f} s)')
        ax.set_aspect('equal')
        ax.set_xlim(0, self.domain_size * 1e6)
        ax.set_ylim(0, self.domain_size * 1e6)
        ax.grid(True, alpha=0.3)
        
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            plt.close()
        else:
            plt.show()

    def plot_density_near_substrate(self, time_points, density_data, labels, save_path=None):
        """
        绘制底物附近菌密度随时间的变化曲线（多组条件对比）。
        
        参数：
            time_points: 时间点数组
            density_data: 密度数据列表，每个元素对应一组条件
            labels: 条件标签列表
            save_path: 保存路径（如果为None则显示）
        """
        fig, ax = plt.subplots(figsize=(10, 6))
        
        colors = ['blue', 'red', 'green', 'cyan', 'magenta', 'orange']
        linestyles = ['-', '--', '-.', '-', '--', '-.']
        
        for i, (density, label) in enumerate(zip(density_data, labels)):
            color = colors[i % len(colors)]
            linestyle = linestyles[i % len(linestyles)]
            ax.plot(time_points, density, label=label, 
                   color=color, linestyle=linestyle, linewidth=2)
        
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Bacteria Density near Substrate (count/μm³)')
        ax.set_title('Bacteria Density near Substrate vs Time')
        ax.legend(loc='best')
        ax.grid(True, alpha=0.3)
        
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            plt.close()
        else:
            plt.show()

    def plot_total_consumption(self, time_points, consumption_data, labels, save_path=None):
        """
        绘制累计营养消耗随时间的变化曲线（多组条件对比）。
        
        参数：
            time_points: 时间点数组
            consumption_data: 消耗数据列表，每个元素对应一组条件
            labels: 条件标签列表
            save_path: 保存路径（如果为None则显示）
        """
        fig, ax = plt.subplots(figsize=(10, 6))
        
        colors = ['blue', 'red', 'green', 'cyan', 'magenta', 'orange']
        linestyles = ['-', '--', '-.', '-', '--', '-.']
        
        for i, (consumption, label) in enumerate(zip(consumption_data, labels)):
            color = colors[i % len(colors)]
            linestyle = linestyles[i % len(linestyles)]
            ax.plot(time_points, consumption, label=label,
                   color=color, linestyle=linestyle, linewidth=2)
        
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Total Glucose Consumed (mol)')
        ax.set_title('Total Glucose Consumption vs Time')
        ax.legend(loc='best')
        ax.grid(True, alpha=0.3)
        
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            plt.close()
        else:
            plt.show()

    def plot_concentration_slice(self, concentration_field, z_slice=None, save_path=None):
        """
        绘制浓度场的二维截面热力图。
        
        参数：
            concentration_field: 浓度场数组
            z_slice: z方向切片索引（如果为None则取中间切片）
            save_path: 保存路径（如果为None则显示）
        """
        if z_slice is None:
            z_slice = concentration_field.shape[2] // 2
        
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # 提取切片并转换为mM
        slice_data = concentration_field[:, :, z_slice] * 1e3
        
        # 绘制热力图
        im = ax.imshow(slice_data.T, origin='lower', cmap='viridis',
                      extent=[0, self.domain_size * 1e6, 0, self.domain_size * 1e6],
                      vmin=0, vmax=self.cfg.SUBSTRATE_SURFACE_CONC * 1e3)
        
        # 添加颜色条
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Glucose Concentration (mM)')
        
        # 标记底物颗粒位置
        for pos in self.substrate_positions:
            circle = plt.Circle((pos[0] * 1e6, pos[1] * 1e6),
                               self.substrate_radius * 1e6,
                               color='red', fill=False, linewidth=2)
            ax.add_patch(circle)
        
        ax.set_xlabel('X (μm)')
        ax.set_ylabel('Y (μm)')
        ax.set_title(f'Glucose Concentration Field (Z = {z_slice})')
        
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            plt.close()
        else:
            plt.show()

    def plot_results(self, time_points, density_data, consumption_data, labels, 
                    bacteria=None, concentration_field=None, output_dir=None):
        """
        生成完整的可视化结果。
        
        参数：
            time_points: 时间点数组
            density_data: 密度数据列表
            consumption_data: 消耗数据列表
            labels: 条件标签列表
            bacteria: 最终状态的细菌列表（可选）
            concentration_field: 最终浓度场（可选）
            output_dir: 输出目录（如果为None则显示）
        """
        # 1. 底物附近菌密度随时间变化
        if output_dir:
            save_path = f"{output_dir}/density_near_substrate.png"
        else:
            save_path = None
        self.plot_density_near_substrate(time_points, density_data, labels, save_path)
        
        # 2. 累计营养消耗随时间变化
        if output_dir:
            save_path = f"{output_dir}/total_consumption.png"
        else:
            save_path = None
        self.plot_total_consumption(time_points, consumption_data, labels, save_path)
        
        # 3. 浓度场热力图
        if concentration_field is not None:
            if output_dir:
                save_path = f"{output_dir}/concentration_field.png"
            else:
                save_path = None
            self.plot_concentration_slice(concentration_field, save_path=save_path)
        
        # 4. 最终菌群空间分布
        if bacteria is not None:
            if output_dir:
                save_path = f"{output_dir}/final_distribution_xy.png"
            else:
                save_path = None
            self.plot_xy_projection(bacteria, time_points[-1], save_path)
