#!/usr/bin/env python3
"""
单元测试脚本

测试各模块的基本功能是否正确实现。
"""

import numpy as np
import sys
import os

# 添加父目录到路径
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import config
from bacteria import Bacterium, random_unit_vector
from environment import Environment
from signaling import SignalingPathway
from movement import MovementEngine
from metabolism import MetabolismEngine
from utils import rodrigues_rotation, trilinear_interpolate, position_to_grid_index


def test_random_unit_vector():
    """测试随机单位向量生成。"""
    print("Testing random_unit_vector()...")
    
    # 生成多个向量检查是否为单位长度
    for _ in range(100):
        v = random_unit_vector()
        assert np.abs(np.linalg.norm(v) - 1.0) < 1e-10, "Vector is not unit length"
    
    # 检查分布是否均匀（粗略检查）
    vectors = np.array([random_unit_vector() for _ in range(1000)])
    mean = np.mean(vectors, axis=0)
    assert np.all(np.abs(mean) < 0.1), "Distribution may not be uniform"
    
    print("  ✓ random_unit_vector() passed")


def test_rodrigues_rotation():
    """测试Rodrigues旋转公式。"""
    print("Testing rodrigues_rotation()...")
    
    # 测试1：绕自身旋转应该不变
    v = np.array([1, 0, 0])
    k = np.array([1, 0, 0])
    v_rot = rodrigues_rotation(v, k, np.pi/4)
    assert np.allclose(v, v_rot), "Rotation around self should not change vector"
    
    # 测试2：90度旋转
    v = np.array([1, 0, 0])
    k = np.array([0, 0, 1])
    v_rot = rodrigues_rotation(v, k, np.pi/2)
    assert np.allclose(v_rot, [0, 1, 0], atol=1e-10), "90 degree rotation failed"
    
    # 测试3：保持长度
    v = np.array([3, 4, 5])
    k = np.array([1, 1, 1])
    v_rot = rodrigues_rotation(v, k, np.pi/3)
    assert np.abs(np.linalg.norm(v) - np.linalg.norm(v_rot)) < 1e-10, "Length not preserved"
    
    print("  ✓ rodrigues_rotation() passed")


def test_trilinear_interpolate():
    """测试三线性插值。"""
    print("Testing trilinear_interpolate()...")
    
    # 创建简单测试场
    field = np.zeros((10, 10, 10))
    field[5, 5, 5] = 1.0
    
    # 在网格点应该返回精确值
    domain_size = 100e-6
    pos = np.array([50e-6, 50e-6, 50e-6])
    val = trilinear_interpolate(field, pos, domain_size)
    assert np.abs(val - 1.0) < 0.01, f"Interpolation at grid point failed: {val}"
    
    # 在零点应该返回0
    pos = np.array([0, 0, 0])
    val = trilinear_interpolate(field, pos, domain_size)
    assert val < 0.1, f"Interpolation at zero point failed: {val}"
    
    print("  ✓ trilinear_interpolate() passed")


def test_bacterium():
    """测试细菌类。"""
    print("Testing Bacterium class...")
    
    pos = [100e-6, 200e-6, 300e-6]
    bac = Bacterium(pos, 0)
    
    # 检查位置
    assert np.allclose(bac.get_position(), pos), "Position mismatch"
    
    # 检查方向是否为单位向量
    assert np.abs(np.linalg.norm(bac.direction) - 1.0) < 1e-10, "Direction not unit vector"
    
    # 检查初始状态
    assert bac.state == "run", "Initial state should be 'run'"
    assert bac.alive == True, "Initial alive should be True"
    
    # 测试翻转概率设置
    bac.set_tumble_probability(0.5)
    assert bac.tumble_probability == 0.5, "Tumble probability not set correctly"
    
    # 测试裁剪
    bac.set_tumble_probability(1.5)
    assert bac.tumble_probability == 1.0, "Tumble probability should be clipped to 1.0"
    bac.set_tumble_probability(-0.5)
    assert bac.tumble_probability == 0.0, "Tumble probability should be clipped to 0.0"
    
    print("  ✓ Bacterium class passed")


def test_environment():
    """测试环境模块。"""
    print("Testing Environment class...")
    
    # 使用较小的网格进行快速测试
    original_grid_res = config.GRID_RESOLUTION
    config.GRID_RESOLUTION = 20
    
    try:
        env = Environment(config)
        
        # 检查浓度场形状
        assert env.concentration_field.shape == (20, 20, 20), "Concentration field shape mismatch"
        
        # 检查底物表面浓度
        substrate_pos = np.array(config.SUBSTRATE_POSITIONS[0])
        conc_at_substrate = env.get_concentration(substrate_pos)
        assert np.abs(conc_at_substrate - config.SUBSTRATE_SURFACE_CONC) / config.SUBSTRATE_SURFACE_CONC < 0.1, \
            f"Concentration at substrate surface incorrect: {conc_at_substrate}"
        
        # 检查远处浓度趋近于0
        far_pos = np.array([100e-6, 100e-6, 100e-6])
        conc_far = env.get_concentration(far_pos)
        assert conc_far < config.SUBSTRATE_SURFACE_CONC * 0.5, "Concentration far from substrate should be lower"
        
        # 测试update方法
        consumption_field = np.zeros((20, 20, 20))
        env.update(0.1, consumption_field)
        
        print("  ✓ Environment class passed")
    finally:
        config.GRID_RESOLUTION = original_grid_res


def test_signaling():
    """测试信号通路模块。"""
    print("Testing SignalingPathway class...")
    
    # 测试无趋化模式
    sig_no_chemo = SignalingPathway(config, "no_chemotaxis")
    bac = Bacterium([0, 0, 0], 0)
    
    # 无论浓度如何，翻转概率应该恒定
    for conc in [0, 1e-6, 10e-6, 100e-6]:
        bac.receptor_activity = 0.5
        bac.methylation = 0.5
        sig_no_chemo.update(bac, conc, 0.1)
        assert np.abs(bac.tumble_probability - config.CW_BIAS_BASELINE) < 1e-10, \
            f"No chemotaxis mode: tumble probability should be constant at {config.CW_BIAS_BASELINE}"
    
    # 测试野生型模式
    sig_wt = SignalingPathway(config, "wild_type")
    bac_wt = Bacterium([0, 0, 0], 1)
    
    # 在恒定浓度下运行多步，受体活性应该回到基线附近（精确适应）
    constant_conc = 5e-6
    for _ in range(1000):
        sig_wt.update(bac_wt, constant_conc, 0.01)
    
    # 适应精度检查（允许5%误差）
    assert np.abs(bac_wt.receptor_activity - config.CW_BIAS_BASELINE) / config.CW_BIAS_BASELINE < 0.05, \
        f"Precise adaptation failed: activity = {bac_wt.receptor_activity}, expected ~{config.CW_BIAS_BASELINE}"
    
    print("  ✓ SignalingPathway class passed")


def test_movement():
    """测试运动模块。"""
    print("Testing MovementEngine class...")
    
    mov = MovementEngine(config)
    
    # 测试run阶段
    bac = Bacterium([1000e-6, 1000e-6, 1000e-6], 0)
    bac.direction = np.array([1, 0, 0])  # 沿x方向
    initial_pos = bac.get_position().copy()
    
    dt = 0.1
    mov.step(bac, dt)
    
    # 检查位移
    displacement = bac.get_position() - initial_pos
    expected_displacement = config.SWIM_SPEED * dt
    assert np.abs(displacement[0] - expected_displacement) < 1e-10, "Run displacement incorrect"
    assert np.abs(displacement[1]) < 1e-10, "Y displacement should be 0"
    assert np.abs(displacement[2]) < 1e-10, "Z displacement should be 0 (no gravity)"
    
    # 测试边界反射
    bac2 = Bacterium([-10e-6, 1000e-6, 1000e-6], 1)
    bac2.direction = np.array([-1, 0, 0])
    mov.step(bac2, dt)
    assert bac2.position[0] >= 0, "Boundary reflection failed: x < 0"
    assert bac2.direction[0] > 0, "Direction should be reversed after reflection"
    
    print("  ✓ MovementEngine class passed")


def test_metabolism():
    """测试代谢模块。"""
    print("Testing MetabolismEngine class...")
    
    met = MetabolismEngine(config)
    bac = Bacterium([0, 0, 0], 0)
    
    # 测试高浓度下的消耗
    high_conc = 100e-6  # 远高于Km
    dt = 1.0
    consumed1 = met.step(bac, high_conc, dt)
    expected_max = config.CONSUMPTION_RATE * dt
    assert np.abs(consumed1 - expected_max) / expected_max < 0.1, \
        f"High concentration consumption incorrect: {consumed1}, expected ~{expected_max}"
    
    # 测试低浓度下的消耗（应该减少）
    bac2 = Bacterium([0, 0, 0], 1)
    low_conc = 0.1e-6  # 远低于Km
    consumed2 = met.step(bac2, low_conc, dt)
    assert consumed2 < consumed1, "Low concentration should consume less"
    
    # 测试营养耗尽阈值
    bac3 = Bacterium([0, 0, 0], 2)
    very_low_conc = 1e-10  # 低于阈值
    consumed3 = met.step(bac3, very_low_conc, dt)
    assert consumed3 == 0.0, "Below threshold consumption should be 0"
    
    print("  ✓ MetabolismEngine class passed")


def run_all_tests():
    """运行所有测试。"""
    print("\n" + "="*60)
    print("Running Unit Tests")
    print("="*60 + "\n")
    
    tests = [
        test_random_unit_vector,
        test_rodrigues_rotation,
        test_trilinear_interpolate,
        test_bacterium,
        test_environment,
        test_signaling,
        test_movement,
        test_metabolism,
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            test()
            passed += 1
        except AssertionError as e:
            print(f"  ✗ {test.__name__} failed: {e}")
            failed += 1
        except Exception as e:
            print(f"  ✗ {test.__name__} error: {e}")
            failed += 1
    
    print("\n" + "="*60)
    print(f"Test Results: {passed} passed, {failed} failed")
    print("="*60 + "\n")
    
    return failed == 0


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
