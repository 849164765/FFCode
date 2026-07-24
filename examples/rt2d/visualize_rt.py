#!/usr/bin/env python3
"""生成RT不稳定性Bom=41.18 vs Bom=0的界面形态对比图
用法: python3 visualize_rt.py [step1 step2 ...]
不指定步数则自动找所有可用共同步
"""
import sys, os, glob, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_rt import collect_field, find_interface
import numpy as np

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False
    print('matplotlib not available')

RT2D_DIR = '/home/dlf/myCode/FerroRTinstability/examples/rt2d'
BOM0_DIR = '/home/dlf/myCode/FerroRTinstability/examples/rt2d_bom0'
T_STAR_SCALE = 36203

def get_common_steps():
    """获取两个案例共同的输出步数"""
    mag_vtms = set()
    for f in glob.glob(f'{RT2D_DIR}/vtkoutput/vtidata/*.vtm'):
        m = re.search(r'rt2d(\d+)\.vtm', f)
        if m: mag_vtms.add(int(m.group(1)))
    nomag_vtms = set()
    for f in glob.glob(f'{BOM0_DIR}/vtkoutput/vtidata/*.vtm'):
        m = re.search(r'rt2d(\d+)\.vtm', f)
        if m: nomag_vtms.add(int(m.group(1)))
    common = sorted(mag_vtms & nomag_vtms)
    return common

def visualize_steps(steps):
    if not HAS_MPL:
        print('matplotlib not available, cannot generate plots')
        return

    n = len(steps)
    if n == 0:
        print('No steps to visualize')
        return

    fig, axes = plt.subplots(2, n, figsize=(5*n, 12))
    if n == 1:
        axes = axes.reshape(2, 1)

    for col, step in enumerate(steps):
        step_int = int(step)
        t_star = step_int / T_STAR_SCALE

        # Bom=41.18
        phi_mag = collect_field(step_int, RT2D_DIR, 'PHI')
        # Bom=0
        phi_nomag = collect_field(step_int, BOM0_DIR, 'PHI')

        if phi_mag is None or phi_nomag is None:
            for row in range(2):
                axes[row, col].text(0.5, 0.5, f'T={step_int}\nNo data',
                    ha='center', va='center', transform=axes[row, col].transAxes)
            continue

        # 上行: Bom=41.18
        ax = axes[0, col]
        # 只显示界面附近区域 (y=300-700, 初始界面在512)
        y_lo, y_hi = 300, 750
        phi_show = phi_mag[y_lo:y_hi, :]
        ax.imshow(phi_show, cmap='RdBu_r', vmin=0, vmax=1, aspect='auto',
                  extent=[0, 256, y_lo, y_hi], origin='lower')
        ax.set_title(f'Bom=41.18  T={step_int} (t*={t_star:.2f})', fontsize=10)
        ax.set_xlabel('x')
        if col == 0:
            ax.set_ylabel('y (磁场)')
        # 画界面线
        iy = find_interface(phi_mag)
        valid = ~np.isnan(iy)
        ax.plot(np.where(valid)[0], iy[valid], 'k-', linewidth=0.5, alpha=0.5)

        # 下行: Bom=0
        ax = axes[1, col]
        phi_show = phi_nomag[y_lo:y_hi, :]
        ax.imshow(phi_show, cmap='RdBu_r', vmin=0, vmax=1, aspect='auto',
                  extent=[0, 256, y_lo, y_hi], origin='lower')
        ax.set_title(f'Bom=0  T={step_int} (t*={t_star:.2f})', fontsize=10)
        ax.set_xlabel('x')
        if col == 0:
            ax.set_ylabel('y (无磁场)')
        iy = find_interface(phi_nomag)
        valid = ~np.isnan(iy)
        ax.plot(np.where(valid)[0], iy[valid], 'k-', linewidth=0.5, alpha=0.5)

    plt.tight_layout()
    outpath = os.path.join(RT2D_DIR, 'rt_visual_comparison.png')
    plt.savefig(outpath, dpi=120, bbox_inches='tight')
    print(f'对比图已保存: {outpath}')
    plt.close()

    # 另外生成界面位置对比图
    fig, ax = plt.subplots(figsize=(10, 6))
    for step in steps:
        step_int = int(step)
        t_star = step_int / T_STAR_SCALE
        phi_mag = collect_field(step_int, RT2D_DIR, 'PHI')
        phi_nomag = collect_field(step_int, BOM0_DIR, 'PHI')
        if phi_mag is None or phi_nomag is None:
            continue
        iy_mag = find_interface(phi_mag)
        iy_nomag = find_interface(phi_nomag)
        x = np.arange(len(iy_mag))
        valid_m = ~np.isnan(iy_mag)
        valid_n = ~np.isnan(iy_nomag)
        ax.plot(x[valid_m], iy_mag[valid_m], '-', label=f'Bom=41.18 t*={t_star:.2f}', alpha=0.7)
        ax.plot(x[valid_n], iy_nomag[valid_n], '--', label=f'Bom=0 t*={t_star:.2f}', alpha=0.7)
    ax.set_xlabel('x')
    ax.set_ylabel('界面 y 位置')
    ax.set_title('界面位置对比 (实线=Bom=41.18, 虚线=Bom=0)')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(300, 750)
    outpath2 = os.path.join(RT2D_DIR, 'rt_interface_comparison.png')
    plt.savefig(outpath2, dpi=120, bbox_inches='tight')
    print(f'界面位置对比图: {outpath2}')
    plt.close()

if __name__ == '__main__':
    if len(sys.argv) > 1:
        steps = sys.argv[1:]
    else:
        steps = get_common_steps()
    print(f'可视化的时间步: {steps}')
    visualize_steps(steps)
