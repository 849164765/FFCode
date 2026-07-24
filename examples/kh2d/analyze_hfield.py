#!/usr/bin/env python3
"""
分析磁场 H 分布, 验证磁场力是否生效.

直接扫描所有 VTI block, 找到界面格点(0.1<phi<0.9),
检查 |H| 在界面附近是否有梯度 (∇|H|≠0 才能产生磁场力 F=χ|H|∇|H|).
无需重建全域网格, 避免块坐标问题.
"""
import os
import sys
import struct
import base64
import xml.etree.ElementTree as ET
import math
import glob


def parse_vti_fields(path):
    """解析 VTI 文件, 返回 {fieldname: ((ox,oy,nx,ny), values)} 字典."""
    tree = ET.parse(path)
    root = tree.getroot()
    piece = root.find('.//Piece')
    if piece is None:
        return None
    ext = list(map(int, piece.get('Extent').split()))
    nx = ext[1] - ext[0] + 1
    ny = ext[3] - ext[2] + 1

    result = {}
    pds = piece.find('PointData')
    if pds is None:
        return result
    for da in pds.findall('DataArray'):
        name = da.get('Name')
        ncomp = int(da.get('NumberOfComponents', '1'))
        text = ''.join(da.text.split())
        # VTK binary+base64: 前 8 字符为 header(uint32), 其余为 data
        hdr_b64 = text[:8]
        hdr = base64.b64decode(hdr_b64)
        data_b64 = text[8:]
        data = base64.b64decode(data_b64)
        nvals = nx * ny * ncomp
        if len(data) >= nvals * 8:
            vals = struct.unpack(f'<{nvals}d', data[:nvals * 8])
        else:
            raw = base64.b64decode(text)
            vals = struct.unpack(f'<{nvals}d', raw[4:4 + nvals * 8])
        result[name] = ((ext[0], ext[2], nx, ny), vals)
    return result


def analyze_step(vtm_path):
    """分析单个时间步, 扫描所有 block 找界面处 H 场统计."""
    tree = ET.parse(vtm_path)
    root = tree.getroot()
    datasets = root.findall('.//DataSet')
    base_dir = os.path.dirname(vtm_path)

    # 全域统计
    hmag_all = []         # 所有格点 |H|
    hmag_interface = []    # 界面处 |H| (0.1<phi<0.9)
    hmag_ferro = []        # 铁磁流体内部 |H| (phi>0.99)
    hmag_upper = []        # 上层 |H| (phi<0.01)

    # 界面处 ∇|H| (用块内中心差分)
    grad_hmag_interface = []

    # H 方向统计 (界面处)
    hx_interface = []
    hy_interface = []

    block_count = 0
    for ds in datasets:
        fp = ds.get('file')
        if not fp:
            continue
        full = os.path.join(base_dir, fp) if not os.path.isabs(fp) else fp
        if not os.path.exists(full):
            continue
        fields = parse_vti_fields(full)
        if fields is None or 'PHI' not in fields or 'HX' not in fields:
            continue
        block_count += 1

        (ox, oy, nx, ny), phi = fields['PHI']
        _, hx = fields['HX']
        _, hy = fields['HY']

        # 计算 |H|
        hmag = [0.0] * (nx * ny)
        for idx in range(nx * ny):
            hmag[idx] = math.sqrt(hx[idx] ** 2 + hy[idx] ** 2)
            hmag_all.append(hmag[idx])
            p = phi[idx]
            if 0.1 < p < 0.9:
                hmag_interface.append(hmag[idx])
                hx_interface.append(hx[idx])
                hy_interface.append(hy[idx])
            elif p > 0.99:
                hmag_ferro.append(hmag[idx])
            elif p < 0.01:
                hmag_upper.append(hmag[idx])

        # 界面处 ∇|H| (块内中心差分, 避开边界)
        for j in range(1, ny - 1):
            for i in range(1, nx - 1):
                idx = j * nx + i
                p = phi[idx]
                if 0.1 < p < 0.9:
                    dhdx = (hmag[j * nx + (i + 1)] - hmag[j * nx + (i - 1)]) * 0.5
                    dhdy = (hmag[(j + 1) * nx + i] - hmag[(j - 1) * nx + i]) * 0.5
                    grad_mag = math.sqrt(dhdx ** 2 + dhdy ** 2)
                    grad_hmag_interface.append(grad_mag)

    if not hmag_all:
        print(f"  无法加载 {vtm_path}")
        return

    # 统计
    hmag_min = min(hmag_all)
    hmag_max = max(hmag_all)
    hmag_mean = sum(hmag_all) / len(hmag_all)
    hmag_range = hmag_max - hmag_min

    # 估算 H0 (理论值)
    # H0 = sqrt(2*sigma*Bom/L)
    # sigma = rho_h*L*U^2/We, L=128 (256x256域), U=0.02, We=10000, Bom=636
    L = 128
    U = 0.02
    We = 10000.0
    rho_h = 1.0
    Bom = 636.0
    Fr = 1.0
    sigma = rho_h * L * U * U / We
    gravity = U * U / (Fr * L)
    H0_theory = math.sqrt(2.0 * sigma * Bom / L)

    print(f"\n  === {os.path.basename(vtm_path)} ===")
    print(f"  加载块数: {block_count}")
    print(f"  总格点数: {len(hmag_all)}")
    print(f"  理论 H0 = {H0_theory:.6e}")

    print(f"  |H| 全域: min={hmag_min:.6e}  max={hmag_max:.6e}  mean={hmag_mean:.6e}")
    print(f"  |H| range = {hmag_range:.6e}  (相对值: {hmag_range/hmag_mean*100 if hmag_mean>0 else 0:.2f}%)")
    print(f"  |H|/H0 (mean) = {hmag_mean/H0_theory:.4f}")

    if hmag_ferro:
        fm = sum(hmag_ferro) / len(hmag_ferro)
        fmin, fmax = min(hmag_ferro), max(hmag_ferro)
        print(f"  |H| 铁磁流体内部(phi>0.99): mean={fm:.6e}  range={fmax-fmin:.6e}  (n={len(hmag_ferro)})")
    if hmag_upper:
        um = sum(hmag_upper) / len(hmag_upper)
        umin, umax = min(hmag_upper), max(hmag_upper)
        print(f"  |H| 上层内部(phi<0.01): mean={um:.6e}  range={umax-umin:.6e}  (n={len(hmag_upper)})")
    if hmag_interface:
        im = sum(hmag_interface) / len(hmag_interface)
        imin, imax = min(hmag_interface), max(hmag_interface)
        print(f"  |H| 界面处(0.1<phi<0.9): mean={im:.6e}  range={imax-imin:.6e}  (n={len(hmag_interface)})")

    if grad_hmag_interface:
        gmax = max(grad_hmag_interface)
        gmean = sum(grad_hmag_interface) / len(grad_hmag_interface)
        print(f"  ∇|H| 界面处: max={gmax:.6e}  mean={gmean:.6e}  (n={len(grad_hmag_interface)})")

        # 估算磁场力: F = chi * |H| * grad|H|
        chi = 1.0
        if hmag_interface:
            im = sum(hmag_interface) / len(hmag_interface)
            fmag_typical = chi * im * gmean
            fmag_max = chi * imax * gmax
        else:
            fmag_typical = 0
            fmag_max = 0

        print(f"  磁场力 F_mag: 典型值≈{fmag_typical:.6e}  最大值≈{fmag_max:.6e}")
        print(f"  对比: 重力 g = {gravity:.6e}")
        print(f"  F_mag/g 比值: 典型={fmag_typical/gravity if gravity > 0 else 0:.4f}  最大={fmag_max/gravity if gravity > 0 else 0:.4f}")

    if hx_interface:
        hxa = [abs(x) for x in hx_interface]
        hya = [abs(y) for y in hy_interface]
        print(f"  H 方向(界面处): |HX|mean={sum(hxa)/len(hxa):.6e}  |HY|mean={sum(hya)/len(hya):.6e}  比值={sum(hya)/sum(hxa) if sum(hxa)>0 else 0:.4f}")

    # 关键判断
    print(f"\n  >>> 诊断:")
    if hmag_range / hmag_mean < 0.01 if hmag_mean > 0 else True:
        print(f"  ⚠ |H| 几乎均匀 (range/mean < 1%), 磁场力 ≈ 0!")
        print(f"  ⚠ 这解释了 Bom=636 与 Bom=0 增长率几乎相同!")
    elif grad_hmag_interface:
        gmean = sum(grad_hmag_interface) / len(grad_hmag_interface)
        if hmag_interface:
            im = sum(hmag_interface) / len(hmag_interface)
            fmag = 1.0 * im * gmean
            ratio = fmag / gravity if gravity > 0 else 0
            if ratio < 0.1:
                print(f"  ⚠ 磁场力远小于重力 (F_mag/g={ratio:.4f} < 0.1), 稳定化效果微弱!")
            else:
                print(f"  ✓ 磁场力可观 (F_mag/g={ratio:.4f}), 稳定化应有效果")


def main():
    vtidir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'vtkoutput', 'vtidata')
    if len(sys.argv) > 1:
        vtms = [sys.argv[1]]
    else:
        vtms = sorted(glob.glob(os.path.join(vtidir, 'kh2d*.vtm')),
                      key=lambda x: int(''.join(c for c in os.path.basename(x) if c.isdigit()) or 0))

    print(f"找到 {len(vtms)} 个 VTM 文件")
    for vtm in vtms:
        analyze_step(vtm)


if __name__ == '__main__':
    main()
