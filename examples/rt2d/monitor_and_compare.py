#!/usr/bin/env python3
"""监控RT不稳定性模拟进度并对比Bom=41.18和Bom=0的界面演化
用法:
  python3 monitor_and_compare.py           # 监控模式: 显示进度 + 可用对比
  python3 monitor_and_compare.py --full    # 完整对比: 所有可用时间步
  python3 monitor_and_compare.py --plot    # 生成对比图
"""
import base64, struct, glob, re, sys, os, subprocess
import numpy as np

# 复用 compare_rt.py 的函数
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_rt import collect_field, find_interface

RT2D_DIR = '/home/dlf/myCode/FerroRTinstability/examples/rt2d'
BOM0_DIR = '/home/dlf/myCode/FerroRTinstability/examples/rt2d_bom0'
T_STAR_SCALE = 36203  # t* = step / 36203

def get_latest_step(basedir):
    """获取最新输出的时间步"""
    vtms = sorted(glob.glob(f'{basedir}/vtkoutput/vtidata/*.vtm'),
                  key=lambda f: int(re.search(r'rt2d(\d+)', f).group(1)))
    if not vtms: return 0
    return int(re.search(r'rt2d(\d+)', vtms[-1]).group(1))

def get_sim_step(basedir, logfile):
    """从日志获取当前模拟步"""
    logpath = os.path.join(basedir, logfile)
    if not os.path.exists(logpath): return 0
    with open(logpath) as f:
        content = f.read()
    # 找最后一个 Diag T=XXXX
    matches = re.findall(r'Diag T=(\d+)', content)
    if matches: return int(matches[-1])
    return 0

def is_running(pattern='rt2d.exe'):
    """检查模拟是否在运行"""
    try:
        out = subprocess.check_output(['ps', '-ef'], text=True)
        return pattern in out
    except:
        return False

def analyze_step(step, verbose=True):
    """分析指定时间步的对比"""
    step_int = int(step)
    t_star = step_int / T_STAR_SCALE

    phi_mag = collect_field(step_int, RT2D_DIR, 'PHI')
    phi_nomag = collect_field(step_int, BOM0_DIR, 'PHI')

    if phi_mag is None or phi_nomag is None:
        return None

    iy_mag = find_interface(phi_mag)
    iy_nomag = find_interface(phi_nomag)

    valid_mag = iy_mag[~np.isnan(iy_mag)]
    valid_nomag = iy_nomag[~np.isnan(iy_nomag)]

    if len(valid_mag) == 0 or len(valid_nomag) == 0:
        return None

    result = {
        'step': step_int,
        't_star': t_star,
        'mag_spike': valid_mag.min(),       # 重相下探 (min y)
        'mag_bubble': valid_mag.max(),      # 轻相上涌 (max y)
        'mag_amp': valid_mag.max() - valid_mag.min(),
        'mag_mean': valid_mag.mean(),
        'nomag_spike': valid_nomag.min(),
        'nomag_bubble': valid_nomag.max(),
        'nomag_amp': valid_nomag.max() - valid_nomag.min(),
        'nomag_mean': valid_nomag.mean(),
    }
    # 磁场效果: spike差值 (负=磁场使spike更深=促进RT)
    result['spike_diff'] = result['mag_spike'] - result['nomag_spike']
    result['bubble_diff'] = result['mag_bubble'] - result['nomag_bubble']
    result['amp_diff'] = result['mag_amp'] - result['nomag_amp']

    if verbose:
        print(f'\n  T={step_int} (t*={t_star:.3f})')
        print(f'    Bom=41.18: spike={result["mag_spike"]:.1f}  bubble={result["mag_bubble"]:.1f}  amp={result["mag_amp"]:.1f}')
        print(f'    Bom=0:     spike={result["nomag_spike"]:.1f}  bubble={result["nomag_bubble"]:.1f}  amp={result["nomag_amp"]:.1f}')
        print(f'    差异:       spike Δ={result["spike_diff"]:+.1f}  bubble Δ={result["bubble_diff"]:+.1f}  amp Δ={result["amp_diff"]:+.1f}')
        if result['spike_diff'] < -5:
            print(f'    >>> 磁场使spike穿透更深 (促进RT不稳定性) <<<')
        elif result['spike_diff'] > 5:
            print(f'    >>> 磁场使spike穿透较浅 (抑制RT不稳定性) <<<')
        else:
            print(f'    >>> 磁场效果不明显 (差异<5格) <<<')

    return result

def monitor():
    """监控模式"""
    print('=' * 70)
    print('  RT不稳定性模拟监控')
    print('=' * 70)

    running = is_running()
    print(f'\n模拟进程: {"运行中" if running else "已停止"}')

    # Bom=41.18 进度
    mag_step = get_sim_step(RT2D_DIR, 'sim_log_80k_4p.txt')
    mag_output = get_latest_step(RT2D_DIR)
    print(f'\nBom=41.18 (磁场):')
    print(f'  当前步: T={mag_step} (t*={mag_step/T_STAR_SCALE:.3f})')
    print(f'  最新输出: T={mag_output}')
    print(f'  目标: T=80000 (t*=2.21)')
    print(f'  进度: {mag_step/80000*100:.1f}%')

    # Bom=0 进度
    nomag_step = get_sim_step(BOM0_DIR, 'sim_log_bom0_4p.txt')
    nomag_output = get_latest_step(BOM0_DIR)
    print(f'\nBom=0 (无磁场):')
    print(f'  当前步: T={nomag_step} (t*={nomag_step/T_STAR_SCALE:.3f})')
    print(f'  最新输出: T={nomag_output}')
    print(f'  目标: T=80000 (t*=2.21)')
    print(f'  进度: {nomag_step/80000*100:.1f}%')

    # H场诊断
    print(f'\nH场诊断 (Bom=41.18):')
    logpath = os.path.join(RT2D_DIR, 'sim_log_80k_4p.txt')
    if os.path.exists(logpath):
        with open(logpath) as f:
            content = f.read()
        matches = re.findall(r'Diag T=(\d+).*ratio=([\d.]+).*theory ([\d.]+)', content)
        if matches:
            last = matches[-1]
            print(f'  T={last[0]}: ratio={last[1]} (theory {last[2]})')

    # 可用对比
    print(f'\n可用对比时间步:')
    # 找共同的输出步 (Bom=41.18 OutputStep=2000, Bom=0 OutputStep=4000, 共同步=4000的倍数)
    common_steps = []
    for s in range(0, min(mag_output, nomag_output) + 1, 4000):
        if os.path.exists(f'{RT2D_DIR}/vtkoutput/vtidata/rt2d_T{s}_B0.vti') and \
           os.path.exists(f'{BOM0_DIR}/vtkoutput/vtidata/rt2d_T{s}_B0.vti'):
            common_steps.append(s)
    if common_steps:
        for s in common_steps:
            t = s / T_STAR_SCALE
            print(f'  T={s} (t*={t:.3f})')
    else:
        print('  (暂无共同输出)')

    return common_steps

def full_comparison(steps=None):
    """完整对比"""
    if steps is None:
        # 自动找所有可用共同步
        mag_out = get_latest_step(RT2D_DIR)
        nomag_out = get_latest_step(BOM0_DIR)
        max_step = min(mag_out, nomag_out)
        steps = list(range(0, max_step + 1, 4000))

    print('\n' + '=' * 70)
    print('  Bom=41.18 vs Bom=0 界面演化对比')
    print('=' * 70)
    print(f'\n  论文结论: 磁场加速RT不稳定性 (spike/bubble速度更快)')
    print(f'  初始界面: y=512, 扰动幅度=51.2 (0.1L*cos)')

    results = []
    for s in steps:
        r = analyze_step(s, verbose=True)
        if r:
            results.append(r)

    if len(results) >= 2:
        print('\n' + '=' * 70)
        print('  界面演化趋势分析')
        print('=' * 70)
        print(f'\n  {"T":>8} {"t*":>8} {"Bom=41.18 spike":>18} {"Bom=0 spike":>15} {"Δspike":>10} {"效果":>15}')
        for r in results:
            effect = '促进' if r['spike_diff'] < -5 else ('抑制' if r['spike_diff'] > 5 else '不明显')
            print(f'  {r["step"]:>8} {r["t_star"]:>8.3f} {r["mag_spike"]:>18.1f} {r["nomag_spike"]:>15.1f} {r["spike_diff"]:>+10.1f} {effect:>15}')

    return results

def make_plots(steps=None):
    """生成对比图"""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('matplotlib not available, skipping plots')
        return

    if steps is None:
        mag_out = get_latest_step(RT2D_DIR)
        nomag_out = get_latest_step(BOM0_DIR)
        max_step = min(mag_out, nomag_out)
        steps = list(range(0, max_step + 1, 4000))

    results = []
    for s in steps:
        r = analyze_step(s, verbose=False)
        if r:
            results.append(r)

    if len(results) < 2:
        print('Not enough data for plots')
        return

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # Plot 1: Spike penetration (重相下探深度)
    ax = axes[0]
    t_stars = [r['t_star'] for r in results]
    mag_spikes = [512 - r['mag_spike'] for r in results]  # 穿透深度 = 512 - spike_y
    nomag_spikes = [512 - r['nomag_spike'] for r in results]
    ax.plot(t_stars, mag_spikes, 'ro-', label='Bom=41.18 (磁场)', linewidth=2)
    ax.plot(t_stars, nomag_spikes, 'bs--', label='Bom=0 (无磁场)', linewidth=2)
    ax.set_xlabel('t* = t/√(L/(g·At))')
    ax.set_ylabel('Spike 穿透深度 (格)')
    ax.set_title('重相下探深度 (Spike)')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 2: Bubble rise (轻相上涌高度)
    ax = axes[1]
    mag_bubbles = [r['mag_bubble'] - 512 for r in results]
    nomag_bubbles = [r['nomag_bubble'] - 512 for r in results]
    ax.plot(t_stars, mag_bubbles, 'ro-', label='Bom=41.18 (磁场)', linewidth=2)
    ax.plot(t_stars, nomag_bubbles, 'bs--', label='Bom=0 (无磁场)', linewidth=2)
    ax.set_xlabel('t* = t/√(L/(g·At))')
    ax.set_ylabel('Bubble 上涌高度 (格)')
    ax.set_title('轻相上涌高度 (Bubble)')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 3: Amplitude
    ax = axes[2]
    mag_amps = [r['mag_amp'] for r in results]
    nomag_amps = [r['nomag_amp'] for r in results]
    ax.plot(t_stars, mag_amps, 'ro-', label='Bom=41.18 (磁场)', linewidth=2)
    ax.plot(t_stars, nomag_amps, 'bs--', label='Bom=0 (无磁场)', linewidth=2)
    ax.set_xlabel('t* = t/√(L/(g·At))')
    ax.set_ylabel('扰动幅度 (格)')
    ax.set_title('界面扰动幅度')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plotpath = os.path.join(RT2D_DIR, 'rt_comparison.png')
    plt.savefig(plotpath, dpi=150, bbox_inches='tight')
    print(f'\n对比图已保存: {plotpath}')
    plt.close()

if __name__ == '__main__':
    if '--full' in sys.argv:
        full_comparison()
    elif '--plot' in sys.argv:
        make_plots()
    else:
        common = monitor()
        if common and len(common) > 1:
            print('\n' + '=' * 70)
            print('  当前可用对比')
            print('=' * 70)
            full_comparison(common)
