#!/usr/bin/env python3
"""
比较 Bom=636 (修复后) 与 Bom=0 的界面振幅增长率.

使用 VTI 文件的 Origin 属性确定块的全局位置, 正确组装全域 phi.
"""
import os
import sys
import struct
import base64
import xml.etree.ElementTree as ET
import math
import glob
import numpy as np


def parse_vti_block(path):
    """解析 VTI 文件, 返回 (origin_x, origin_y, nx, ny, phi_values)."""
    tree = ET.parse(path)
    root = tree.getroot()
    piece = root.find('.//Piece')
    if piece is None:
        return None
    ext = list(map(int, piece.get('Extent').split()))
    nx = ext[1] - ext[0] + 1
    ny = ext[3] - ext[2] + 1

    # 获取 Origin
    img = root.find('ImageData')
    if img is not None:
        origin = list(map(float, img.get('Origin', '0 0 0').split()))
        ox, oy = origin[0], origin[1]
    else:
        ox, oy = 0.0, 0.0

    pds = piece.find('PointData')
    if pds is None:
        return None
    for da in pds.findall('DataArray'):
        if da.get('Name') != 'PHI':
            continue
        text = ''.join(da.text.split())
        hdr_b64 = text[:8]
        hdr = base64.b64decode(hdr_b64)
        data_b64 = text[8:]
        data = base64.b64decode(data_b64)
        nvals = nx * ny
        if len(data) >= nvals * 8:
            vals = struct.unpack(f'<{nvals}d', data[:nvals * 8])
        else:
            raw = base64.b64decode(text)
            vals = struct.unpack(f'<{nvals}d', raw[4:4 + nvals * 8])
        return (ox, oy, nx, ny, list(vals))
    return None


def load_phi_global(vtm_path):
    """加载 VTM, 使用 Origin 组装全域 phi 2D 数组."""
    tree = ET.parse(vtm_path)
    root = tree.getroot()
    datasets = root.findall('.//DataSet')
    base_dir = os.path.dirname(vtm_path)

    blocks = []
    for ds in datasets:
        fp = ds.get('file')
        if not fp:
            continue
        full = os.path.join(base_dir, fp) if not os.path.isabs(fp) else fp
        if not os.path.exists(full):
            continue
        result = parse_vti_block(full)
        if result is not None:
            blocks.append(result)

    if not blocks:
        return None

    # 确定全域大小
    # Origin 是块左下角(含ghost), 内部格点从 Origin+1 开始
    max_x = 0
    max_y = 0
    for (ox, oy, nx, ny, _) in blocks:
        max_x = max(max_x, int(round(ox)) + nx)
        max_y = max(max_y, int(round(oy)) + ny)

    phi = np.zeros((max_y, max_x))

    for (ox, oy, nx, ny, vals) in blocks:
        arr = np.array(vals).reshape(ny, nx)
        # Origin 是 ghost cell 的坐标, 内部格点从 Origin+1 开始
        ix0 = int(round(ox))
        iy0 = int(round(oy))
        phi[iy0:iy0+ny, ix0:ix0+nx] = arr

    # 提取内部区域 (去掉 ghost cells)
    # ghost cells 通常在边界, 内部是 [1:max-1]
    # 但不同块的 ghost 位置不同, 简单处理: 取 [1:255, 1:255]
    return phi[1:max_y-1, 1:max_x-1]


def extract_interface_amplitude(phi):
    """提取界面振幅.

    界面位置: phi=0.5 等值线
    对每个 x 列, 找到 phi=0.5 的 y 位置 (线性插值)
    返回界面位置的 y 坐标数组和振幅
    """
    ny, nx = phi.shape
    interface_y = np.full(nx, np.nan)

    for i in range(nx):
        col = phi[:, i]
        for j in range(ny - 1):
            if (col[j] - 0.5) * (col[j+1] - 0.5) <= 0:
                if abs(col[j+1] - col[j]) > 1e-10:
                    frac = (0.5 - col[j]) / (col[j+1] - col[j])
                    interface_y[i] = j + frac
                else:
                    interface_y[i] = j
                break

    # 去除 NaN
    valid = ~np.isnan(interface_y)
    if np.sum(valid) < 2:
        return np.array([]), 0, 0

    iy = interface_y[valid]
    amplitude_std = np.std(iy)
    amplitude_range = np.max(iy) - np.min(iy)

    return iy, amplitude_std, amplitude_range


def analyze_growth(vtm_dir, label):
    """分析一个目录下所有 VTM 的界面振幅增长率."""
    vtms = sorted(glob.glob(os.path.join(vtm_dir, 'kh2d*.vtm')),
                  key=lambda x: int(os.path.basename(x).replace('kh2d','').replace('.vtm','') or 0))

    print(f"\n{'='*60}")
    print(f"  {label}")
    print(f"{'='*60}")

    steps = []
    amplitudes = []

    for vtm in vtms:
        step = int(os.path.basename(vtm).replace('kh2d','').replace('.vtm','') or 0)
        phi = load_phi_global(vtm)
        if phi is None:
            continue

        iy, amp_std, amp_range = extract_interface_amplitude(phi)

        steps.append(step)
        amplitudes.append(amp_std)
        print(f"  Step {step:5d}: amplitude(std)={amp_std:.4f}  range={amp_range:.4f}  n_valid={len(iy)}")

    if len(steps) < 2:
        print("  数据不足, 无法计算增长率")
        return steps, amplitudes

    # 计算增长率 (线性拟合 log(amplitude) vs step)
    steps_arr = np.array(steps, dtype=float)
    amps_arr = np.array(amplitudes)

    # 过滤掉振幅为0的点
    mask = amps_arr > 1e-6
    if np.sum(mask) < 2:
        print("  有效数据不足")
        return steps, amplitudes

    log_amps = np.log(amps_arr[mask])
    fit_steps = steps_arr[mask]

    # 线性拟合: log(A) = gamma * step + C
    coeffs = np.polyfit(fit_steps, log_amps, 1)
    gamma = coeffs[0]

    # 理论增长率 (无磁场): gamma_theory = k * DeltaU / 2
    L = 128
    U = 0.02
    k = 2 * math.pi / L
    gamma_theory = k * U / 2

    print(f"\n  增长率分析:")
    print(f"    测量 gamma = {gamma:.6e} /step")
    print(f"    理论 gamma = {gamma_theory:.6e} /step (无磁场, k*DeltaU/2)")
    print(f"    比值 gamma_meas/gamma_theory = {gamma/gamma_theory:.4f}")

    return steps, amplitudes, gamma


def main():
    base_dir = os.path.dirname(os.path.abspath(__file__))

    # Bom=636 修复后
    bom636_dir = os.path.join(base_dir, 'vtkoutput_bom636_fixed', 'vtidata')
    if os.path.exists(bom636_dir):
        result636 = analyze_growth(bom636_dir, "Bom=636 (Dirichlet BC 修复后)")
    else:
        print("Bom=636 修复后结果不存在")
        result636 = ([], [], 0)

    # Bom=0
    bom0_dir = os.path.join(base_dir, 'vtkoutput', 'vtidata')
    if os.path.exists(bom0_dir):
        result0 = analyze_growth(bom0_dir, "Bom=0 (无磁场)")
    else:
        print("Bom=0 结果不存在")
        result0 = ([], [], 0)

    # 对比
    if len(result636[0]) > 1 and len(result0[0]) > 1:
        print(f"\n{'='*60}")
        print(f"  增长率对比")
        print(f"{'='*60}")
        gamma636 = result636[2]
        gamma0 = result0[2]
        print(f"  Bom=636: gamma = {gamma636:.6e} /step")
        print(f"  Bom=0:   gamma = {gamma0:.6e} /step")
        if gamma0 > 0:
            print(f"  比值 gamma_636/gamma_0 = {gamma636/gamma0:.4f}")
            print(f"  磁场使增长率降低: {(1-gamma636/gamma0)*100:.1f}%")

        # 振幅对比
        print(f"\n  振幅对比:")
        s636, a636 = result636[0], result636[1]
        s0, a0 = result0[0], result0[1]
        for i in range(min(len(s636), len(s0))):
            if s636[i] == s0[i] and a0[i] > 1e-6:
                ratio = a636[i] / a0[i]
                print(f"    Step {s636[i]:5d}: Bom=636={a636[i]:.4f}  Bom=0={a0[i]:.4f}  ratio={ratio:.4f}")


if __name__ == '__main__':
    main()
