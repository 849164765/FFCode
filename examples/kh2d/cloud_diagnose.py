#!/usr/bin/env python3
"""
云端诊断脚本：快速判断云端运行的是 Dirichlet BC (修复版) 还是 Neumann BC (旧版).

用法:
    cd /path/to/FerroKHinstability/examples/kh2d
    python3 cloud_diagnose.py [vti输出目录]

如果不指定目录，默认查找 ./vtkoutput/vtidata

诊断逻辑:
    1. 检查 git commit 版本
    2. 检查 kh2d.ini 配置
    3. 分析 VTK 输出中的 H 场统计
    4. 根据 |H|/H0 比值判断 BC 类型:
       - > 90%: Dirichlet BC (修复版, 正确)
       - 40-60%: Neumann BC (旧版, 有bug)
       - < 5%: MF 未演化或更早版本
"""
import os
import sys
import re
import math
import glob
import subprocess
import xml.etree.ElementTree as ET
import base64
import struct


def check_git_version():
    """检查 git commit 版本."""
    print("=" * 70)
    print("[1/4] Git 代码版本检查")
    print("=" * 70)
    # 找到项目根目录
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # 向上查找 .git
    proj_root = script_dir
    for _ in range(5):
        if os.path.isdir(os.path.join(proj_root, '.git')):
            break
        proj_root = os.path.dirname(proj_root)

    try:
        # 当前 HEAD
        head = subprocess.check_output(
            ['git', '-C', proj_root, 'log', '--oneline', '-1'],
            stderr=subprocess.DEVNULL
        ).decode().strip()
        print(f"  当前 HEAD: {head}")

        # 检查是否包含 Dirichlet BC 修复
        fix_commit = "3d84ac1"
        result = subprocess.run(
            ['git', '-C', proj_root, 'log', '--oneline', '--all'],
            capture_output=True, text=True
        )
        if fix_commit in result.stdout:
            # 检查当前 HEAD 是否包含该修复
            contains = subprocess.run(
                ['git', '-C', proj_root, 'branch', '--contains', fix_commit],
                capture_output=True, text=True
            )
            # 检查当前分支
            branch = subprocess.check_output(
                ['git', '-C', proj_root, 'rev-parse', '--abbrev-ref', 'HEAD'],
                stderr=subprocess.DEVNULL
            ).decode().strip()
            print(f"  当前分支: {branch}")

            if fix_commit in head:
                print(f"  ✓ 已包含 Dirichlet BC 修复 ({fix_commit})")
            else:
                # 检查 HEAD 是否在 fix_commit 之后
                is_ancestor = subprocess.run(
                    ['git', '-C', proj_root, 'merge-base', '--is-ancestor', fix_commit, 'HEAD'],
                    capture_output=True
                )
                if is_ancestor.returncode == 0:
                    print(f"  ✓ HEAD 在 Dirichlet BC 修复 ({fix_commit}) 之后")
                else:
                    print(f"  ✗ HEAD 不包含 Dirichlet BC 修复 ({fix_commit})!")
                    print(f"    >>> 这是问题根因! 请执行: git pull origin {branch}")
                    return False
        else:
            print(f"  ✗ 仓库中未找到 Dirichlet BC 修复 ({fix_commit})")
            return False
    except Exception as e:
        print(f"  ⚠ 无法检查 git 版本: {e}")
    print()
    return True


def check_config():
    """检查 kh2d.ini 配置."""
    print("=" * 70)
    print("[2/4] kh2d.ini 配置检查")
    print("=" * 70)
    ini_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'kh2d.ini')
    if not os.path.exists(ini_path):
        print(f"  ⚠ 未找到 {ini_path}")
        return

    expected = {
        'Ni': ('512', '网格数X (应为512, L=256半域)'),
        'Nj': ('512', '网格数Y (应为512)'),
        'Bom': ('636', '磁Bond数'),
        'TotalStep': ('80000', '总步数'),
        'OutputStep': ('2000', '输出间隔'),
    }

    with open(ini_path) as f:
        content = f.read()

    for key, (expected_val, desc) in expected.items():
        m = re.search(rf'{key}\s*=\s*(\S+)', content)
        if m:
            actual = m.group(1).split(';')[0].split('#')[0].strip()
            ok = actual == expected_val
            status = '✓' if ok else '⚠'
            print(f"  {status} {key} = {actual}  (期望: {expected_val})  {desc}")
        else:
            print(f"  ✗ {key} 未找到  {desc}")
    print()


def parse_vti(path):
    """解析 VTI 文件, 返回 {fieldname: values}."""
    tree = ET.parse(path)
    piece = tree.getroot().find('.//Piece')
    if piece is None:
        return None
    ext = list(map(int, piece.get('Extent').split()))
    nx = ext[1] - ext[0] + 1
    ny = ext[3] - ext[2] + 1
    pds = piece.find('PointData')
    if pds is None:
        return None
    result = {}
    for da in pds.findall('DataArray'):
        name = da.get('Name')
        ncomp = int(da.get('NumberOfComponents', '1'))
        text = ''.join(da.text.split())
        data = base64.b64decode(text[8:])
        nvals = nx * ny * ncomp
        if len(data) >= nvals * 8:
            vals = struct.unpack(f'<{nvals}d', data[:nvals * 8])
        else:
            raw = base64.b64decode(text)
            vals = struct.unpack(f'<{nvals}d', raw[4:4 + nvals * 8])
        result[name] = vals
    return result, nx, ny


def analyze_hfield(vtidir):
    """分析 H 场统计."""
    print("=" * 70)
    print("[3/4] VTK 输出 H 场分析")
    print("=" * 70)

    if not os.path.isdir(vtidir):
        print(f"  ✗ 目录不存在: {vtidir}")
        return None

    vtms = sorted(glob.glob(os.path.join(vtidir, 'kh2d*.vtm')),
                  key=lambda x: int(re.search(r'kh2d(\d+)', os.path.basename(x)).group(1)))

    if not vtms:
        print(f"  ✗ 未找到 VTM 文件 in {vtidir}")
        return None

    print(f"  找到 {len(vtms)} 个 VTM 文件")
    print(f"  分析最后3个时间步...\n")

    # 分析最后3个 VTM (最稳定的统计)
    results = []
    for vtm in vtms[-3:]:
        step = int(re.search(r'kh2d(\d+)', os.path.basename(vtm)).group(1))

        tree = ET.parse(vtm)
        datasets = tree.getroot().findall('.//DataSet')
        base_dir = os.path.dirname(vtm)

        hmag_all = []
        hmag_interface = []
        for ds in datasets:
            fp = ds.get('file')
            if not fp:
                continue
            full = os.path.join(base_dir, fp) if not os.path.isabs(fp) else fp
            if not os.path.exists(full):
                continue
            result = parse_vti(full)
            if result is None:
                continue
            fields, _, _ = result
            if 'PHI' not in fields or 'HX' not in fields:
                continue
            phi = fields['PHI']
            hx = fields['HX']
            hy = fields.get('HY', (0.0,) * len(hx))
            for idx in range(len(phi)):
                hmag = math.sqrt(hx[idx]**2 + hy[idx]**2)
                hmag_all.append(hmag)
                if 0.1 < phi[idx] < 0.9:
                    hmag_interface.append(hmag)

        if not hmag_all:
            continue

        hmag_mean = sum(hmag_all) / len(hmag_all)
        hmag_max = max(hmag_all)
        hmag_min = min(hmag_all)

        results.append({
            'step': step,
            'hmag_mean': hmag_mean,
            'hmag_max': hmag_max,
            'hmag_min': hmag_min,
            'n_interface': len(hmag_interface),
            'hmag_interface_mean': sum(hmag_interface) / len(hmag_interface) if hmag_interface else 0,
        })

        print(f"  step={step}: |H| mean={hmag_mean:.6e} max={hmag_max:.6e}  "
              f"界面处mean={sum(hmag_interface)/len(hmag_interface) if hmag_interface else 0:.6e}")

    return results


def diagnose(results, git_ok):
    """根据分析结果诊断 BC 类型."""
    print("=" * 70)
    print("[4/4] 诊断结论")
    print("=" * 70)

    # 理论 H0 计算 (根据 kh2d.ini 参数)
    ini_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'kh2d.ini')
    params = {'Ni': 512, 'Bom': 636.0, 'U': 0.02, 'We': 10000.0,
              'rho_h': 1.0, 'Fr': 1.0, 'Cell_Len': 1.0}
    if os.path.exists(ini_path):
        with open(ini_path) as f:
            content = f.read()
        for key in params:
            m = re.search(rf'^\s*{key}\s*=\s*([0-9.eE+-]+)', content, re.MULTILINE)
            if m:
                try:
                    val = m.group(1)
                    if key == 'Ni':
                        params[key] = int(val)
                    else:
                        params[key] = float(val)
                except ValueError:
                    pass

    Ni = params['Ni']; Bom = params['Bom']; U = params['U']
    We = params['We']; rho_h = params['rho_h']; Fr = params['Fr']
    Cell_Len = params['Cell_Len']

    L = Ni * Cell_Len * 0.5
    sigma = rho_h * L * U * U / We
    gravity = U * U / (Fr * L)
    H0_theory = math.sqrt(2.0 * sigma * Bom / L)

    print(f"\n  理论参数: Ni={Ni}, L={L}, Bom={Bom}, H0={H0_theory:.6e}")
    print(f"  t* = step × {math.sqrt(gravity/L):.6e}  (t*=6 需 {int(6.0/math.sqrt(gravity/L))} 步)\n")

    if not results:
        print("  ⚠ 无可用数据")
        return

    # 取最后一个时间步的统计
    last = results[-1]
    ratio = last['hmag_mean'] / H0_theory if H0_theory > 0 else 0
    step = last['step']
    t_star = step * math.sqrt(gravity / L)

    print(f"  最后时间步: step={step}, t*={t_star:.4f}")
    print(f"  |H|/H0 (mean) = {ratio:.4f} ({ratio*100:.1f}%)")
    print(f"  |H|/H0 (max)  = {last['hmag_max']/H0_theory:.4f}")
    print()

    if ratio > 0.90:
        print("  ✓✓✓ 诊断: Dirichlet BC (修复版, 正确)")
        print("       H 场维持良好 (>90% of H0)")
        print("       磁场力应有效稳定 KH 不稳定性")
        if t_star < 2.0:
            print(f"       当前 t*={t_star:.2f} 应处于早期生长阶段 (正常)")
        else:
            print(f"       当前 t*={t_star:.2f}")
    elif 0.40 < ratio < 0.70:
        print("  ✗✗✗ 诊断: Neumann BC (旧版, 有bug!)")
        print(f"       H 场衰减至 {ratio*100:.1f}% of H0 (应 >90%)")
        print("       磁场力不足, KH演化加速约3.84倍")
        print(f"       step={step} 外观 ≈ 修复版 step={int(step*3.84)} (t*={step*3.84*math.sqrt(gravity/L):.2f})")
        print()
        print("  >>> 修复方法:")
        print("      git pull origin FerroKHinstability")
        print("      确认 git log -1 显示 commit 3d84ac1")
        print("      重新编译并运行仿真")
    elif ratio < 0.05:
        print("  ✗✗✗ 诊断: H 场几乎为零!")
        print("       MF 格子可能未演化 (Bom<=0 或更早版本代码)")
        print("       完全无磁场稳定化, KH演化最快")
    else:
        print(f"  ⚠ 诊断: H 场异常 ({ratio*100:.1f}% of H0)")
        print("       请检查代码版本和配置")

    print()
    print("=" * 70)
    if not git_ok:
        print("  ⚠ Git 版本检查发现问题, 请先解决代码版本问题")
    print("  诊断完成")
    print("=" * 70)


def main():
    print()
    print("╔" + "═" * 68 + "╗")
    print("║  云端 KH 不稳定性仿真诊断脚本" + " " * 36 + "║")
    print("║  检查: 代码版本 / 配置 / H场状态 / BC类型判断" + " " * 22 + "║")
    print("╚" + "═" * 68 + "╝")
    print()

    # 1. Git 版本
    git_ok = check_git_version()

    # 2. 配置
    check_config()

    # 3. H 场分析
    vtidir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        os.path.dirname(os.path.abspath(__file__)), 'vtkoutput', 'vtidata')
    results = analyze_hfield(vtidir)

    # 4. 诊断
    diagnose(results, git_ok)


if __name__ == '__main__':
    main()
