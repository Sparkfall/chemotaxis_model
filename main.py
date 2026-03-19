# main.py —— 模块E（主控部分）
"""
主控模块

趋化性模拟主程序。
整合所有模块，执行仿真循环，记录数据并输出结果。

仿真循环（每个时间步）：
    1. 模块A：更新浓度场（扩散 + 对流 + 消耗反馈）
    2. 对每个细菌：
       a. 模块A：查询当前位置的浓度
       b. 模块B：根据浓度历史计算翻转概率
       c. 模块C：执行运动（run或tumble）
       d. 模块D：消耗当前位置的营养
    3. 模块D：更新细胞状态（生长/死亡）
    4. 模块E：记录数据，定期输出可视化
"""

import numpy as np
import os
import sys

# 导入配置
import config

# 导入各模块
from environment import Environment
from signaling import SignalingPathway
from movement import MovementEngine
from metabolism import MetabolismEngine
from bacteria import Bacterium
from visualization import Visualizer
from utils import position_to_grid_index


def run_simulation(chemotaxis_mode="wild_type", gravity_on=True, 
                   n_bacteria=None, total_time=None, output_dir=None):
    """
    运行一次仿真。
    
    参数：
        chemotaxis_mode: 趋化模式 ("wild_type", "no_chemotaxis", "enhanced")
        gravity_on: 是否开启重力 (True/False)
        n_bacteria: 细菌数量（默认使用config中的值）
        total_time: 总仿真时间（默认使用config中的值）
        output_dir: 输出目录
    
    返回：
        results: 包含仿真结果的字典
    """
    # 使用配置参数
    if n_bacteria is None:
        n_bacteria = config.N_BACTERIA
    if total_time is None:
        total_time = config.TOTAL_TIME
    
    dt = config.DT
    record_interval = config.RECORD_INTERVAL
    proximity_radius = config.SUBSTRATE_PROXIMITY_RADIUS
    
    print(f"\n{'='*60}")
    print(f"Starting simulation:")
    print(f"  Chemotaxis mode: {chemotaxis_mode}")
    print(f"  Gravity: {'ON' if gravity_on else 'OFF'}")
    print(f"  Bacteria count: {n_bacteria}")
    print(f"  Total time: {total_time} s")
    print(f"{'='*60}\n")
    
    # 1. 初始化所有模块
    # 临时修改config中的重力设置
    original_gravity = config.GRAVITY_ON
    config.GRAVITY_ON = gravity_on
    
    env = Environment(config)
    sig = SignalingPathway(config, chemotaxis_mode)
    mov = MovementEngine(config)
    met = MetabolismEngine(config)
    vis = Visualizer(config)
    
    # 恢复config
    config.GRAVITY_ON = original_gravity
    
    # 2. 创建细菌群体
    bacteria = []
    for i in range(n_bacteria):
        # 均匀随机分布
        pos = np.random.uniform(0, config.DOMAIN_SIZE, 3)
        bacteria.append(Bacterium(pos, i))
    
    print(f"Initialized {len(bacteria)} bacteria")
    
    # 3. 数据记录器
    time_points = []
    density_near_substrate = []  # 底物附近的菌密度
    total_consumption = []        # 累计营养消耗
    
    # 4. 主循环
    n_steps = int(total_time / dt)
    record_steps = int(record_interval / dt)
    
    print(f"Running {n_steps} steps...")
    
    for step in range(n_steps):
        t = step * dt
        
        # 创建消耗场（和浓度网格同形状）
        consumption_field = np.zeros((config.GRID_RESOLUTION,) * 3)
        
        # 处理每个细菌
        for bac in bacteria:
            if not bac.alive:
                continue
            
            # a. 查询浓度
            conc = env.get_concentration(bac.get_position())
            
            # b. 更新信号通路
            sig.update(bac, conc, dt)
            
            # c. 执行运动
            mov.step(bac, dt)
            
            # d. 营养消耗
            consumed = met.step(bac, conc, dt)
            
            # 将消耗映射到最近网格点
            grid_idx = position_to_grid_index(
                bac.get_position(), 
                config.DOMAIN_SIZE, 
                config.GRID_RESOLUTION
            )
            consumption_field[grid_idx] += consumed
        
        # 更新浓度场
        env.update(dt, consumption_field)
        
        # 记录数据（每record_interval秒记录一次）
        if step % record_steps == 0:
            # 计算底物附近菌密度
            substrate_pos = np.array(config.SUBSTRATE_POSITIONS[0])
            count_near = 0
            for bac in bacteria:
                if bac.alive:
                    dist = np.linalg.norm(bac.get_position() - substrate_pos)
                    if dist <= proximity_radius:
                        count_near += 1
            
            # 计算密度（count / volume）
            volume = (4/3) * np.pi * (proximity_radius ** 3)
            density = count_near / volume
            
            # 计算累计消耗
            total_consumed = sum(bac.nutrient_consumed for bac in bacteria)
            
            time_points.append(t)
            density_near_substrate.append(density)
            total_consumption.append(total_consumed)
            
            # 打印进度（每10个记录间隔打印一次，或至少每10秒打印一次）
            print_interval_steps = max(record_steps * 10, int(10.0 / dt))  # 至少每10秒打印一次
            if step % print_interval_steps == 0:
                print(f"  t = {t:.0f} s: density = {density:.2e} /μm³, "
                      f"consumed = {total_consumed:.2e} mol")
    
    print(f"\nSimulation completed!")
    print(f"  Final density near substrate: {density_near_substrate[-1]:.2e} /μm³")
    print(f"  Total glucose consumed: {total_consumption[-1]:.2e} mol")
    
    # 5. 返回结果
    results = {
        'time_points': np.array(time_points),
        'density_near_substrate': np.array(density_near_substrate),
        'total_consumption': np.array(total_consumption),
        'bacteria': bacteria,
        'concentration_field': env.concentration_field,
        'chemotaxis_mode': chemotaxis_mode,
        'gravity_on': gravity_on
    }
    
    return results


def run_all_conditions(output_dir="results"):
    """
    运行六组对比实验。
    
    条件组合：
    1. 有趋化 + 1g
    2. 无趋化 + 1g
    3. 增强趋化 + 1g
    4. 有趋化 + 微重力
    5. 无趋化 + 微重力
    6. 增强趋化 + 微重力
    """
    # 创建输出目录
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # 六组实验条件
    conditions = [
        ("wild_type", True, "WT + 1g"),
        ("no_chemotaxis", True, "No Chemo + 1g"),
        ("enhanced", True, "Enhanced + 1g"),
        ("wild_type", False, "WT + 0g"),
        ("no_chemotaxis", False, "No Chemo + 0g"),
        ("enhanced", False, "Enhanced + 0g"),
    ]
    
    all_results = []
    labels = []
    
    # 运行所有条件
    for chemotaxis_mode, gravity_on, label in conditions:
        results = run_simulation(
            chemotaxis_mode=chemotaxis_mode,
            gravity_on=gravity_on,
            output_dir=output_dir
        )
        all_results.append(results)
        labels.append(label)
    
    # 提取数据用于绘图
    time_points = all_results[0]['time_points']
    density_data = [r['density_near_substrate'] for r in all_results]
    consumption_data = [r['total_consumption'] for r in all_results]
    
    # 生成可视化
    print("\nGenerating visualizations...")
    vis = Visualizer(config)
    
    # 1. 底物附近菌密度对比
    vis.plot_density_near_substrate(
        time_points, density_data, labels,
        save_path=f"{output_dir}/density_comparison.png"
    )
    print(f"  Saved: {output_dir}/density_comparison.png")
    
    # 2. 累计消耗对比
    vis.plot_total_consumption(
        time_points, consumption_data, labels,
        save_path=f"{output_dir}/consumption_comparison.png"
    )
    print(f"  Saved: {output_dir}/consumption_comparison.png")
    
    # 3. 各条件的最终分布图
    for i, (results, label) in enumerate(zip(all_results, labels)):
        vis.plot_xy_projection(
            results['bacteria'], time_points[-1],
            save_path=f"{output_dir}/final_dist_{label.replace(' ', '_').replace('+', '')}.png"
        )
    
    # 4. 浓度场（使用第一个结果的浓度场）
    vis.plot_concentration_slice(
        all_results[0]['concentration_field'],
        save_path=f"{output_dir}/concentration_field.png"
    )
    print(f"  Saved: {output_dir}/concentration_field.png")
    
    print(f"\nAll results saved to: {output_dir}/")
    
    return all_results


def main():
    """
    主函数。
    """
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Bacterial Chemotaxis Simulation'
    )
    parser.add_argument(
        '--mode', type=str, default='all',
        choices=['all', 'wild_type', 'no_chemotaxis', 'enhanced'],
        help='Chemotaxis mode to simulate'
    )
    parser.add_argument(
        '--gravity', action='store_true',
        help='Enable gravity (1g)'
    )
    parser.add_argument(
        '--no-gravity', action='store_true',
        help='Disable gravity (microgravity)'
    )
    parser.add_argument(
        '--output', type=str, default='results',
        help='Output directory for results'
    )
    parser.add_argument(
        '--time', type=float, default=None,
        help='Total simulation time in seconds (default from config)'
    )
    parser.add_argument(
        '--bacteria', type=int, default=None,
        help='Number of bacteria (default from config)'
    )
    
    args = parser.parse_args()
    
    # 确定重力设置
    if args.gravity and args.no_gravity:
        print("Error: Cannot specify both --gravity and --no-gravity")
        sys.exit(1)
    elif args.gravity:
        gravity_on = True
    elif args.no_gravity:
        gravity_on = False
    else:
        gravity_on = config.GRAVITY_ON  # 使用config默认值
    
    # 运行仿真
    if args.mode == 'all':
        run_all_conditions(output_dir=args.output)
    else:
        results = run_simulation(
            chemotaxis_mode=args.mode,
            gravity_on=gravity_on,
            n_bacteria=args.bacteria,
            total_time=args.time,
            output_dir=args.output
        )
        
        # 生成可视化
        vis = Visualizer(config)
        label = f"{args.mode} + {'1g' if gravity_on else '0g'}"
        
        vis.plot_results(
            results['time_points'],
            [results['density_near_substrate']],
            [results['total_consumption']],
            [label],
            bacteria=results['bacteria'],
            concentration_field=results['concentration_field'],
            output_dir=args.output
        )
        
        print(f"\nResults saved to: {args.output}/")


if __name__ == "__main__":
    main()
