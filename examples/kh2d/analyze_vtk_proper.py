#!/usr/bin/env python3
"""使用标准库正确解析 VTI 文件 (base64+binary inline 格式) 并分析 KH 界面演化.

用法:
  python3 analyze_vtk_proper.py [vti_dir] [L]
默认:
  vti_dir = vtkoutput/vtidata
  L = Ni*Cell_Len/2 (从 kh2d.ini 推断)
"""
import os
import sys
import re
import base64
import struct
import math
import xml.etree.ElementTree as ET
import numpy as np


def parse_ini():
    """从 kh2d.ini 读取 Ni, Cell_Len, U, Fr, Bom."""
    cfg = {'Ni':256,'Cell_Len':1.0,'U':0.02,'Fr':1.0,'Bom':636}
    try:
        with open('kh2d.ini','r') as f:
            section=None
            for ln in f:
                s=ln.strip()
                if s.startswith('[') and s.endswith(']'):
                    section=s[1:-1]; continue
                if '=' in s and not s.startswith(';'):
                    k,v=s.split('=',1)
                    k=k.strip(); v=v.split(';',1)[0].split('#',1)[0].strip()
                    if k in cfg:
                        try:
                            if isinstance(cfg[k],int): cfg[k]=int(v)
                            else: cfg[k]=float(v)
                        except: pass
    except FileNotFoundError:
        pass
    return cfg


def decode_base64_binary(b64_bytes, dtype='<f8'):
    """解码 VTK inline base64 binary 数据.
    VTK inline binary base64 由多个独立 base64 块组成, 每块自带 padding.
    格式: [header块(uint32 size)] + [data块(raw bytes)]
    每块以 '=' padding 结尾, 块之间无分隔符.
    """
    # 切分为多个独立 base64 块 (按 '=' padding 分割)
    chunks = []
    i = 0
    n = len(b64_bytes)
    while i < n:
        end = i
        while end < n and b64_bytes[end:end+1] != b'=':
            end += 1
        # 包含所有 padding
        while end < n and b64_bytes[end:end+1] == b'=':
            end += 1
        if end > i:
            chunks.append(b64_bytes[i:end])
            i = end
        else:
            break
    if not chunks:
        return np.array([], dtype=dtype)

    # 解码并拼接
    all_raw = b''
    for c in chunks:
        try:
            all_raw += base64.b64decode(c)
        except Exception:
            return np.array([], dtype=dtype)

    # 第一块是 uint32 size header (4 字节)
    if len(all_raw) < 4:
        return np.array([], dtype=dtype)
    hdr_size = struct.unpack('<I', all_raw[:4])[0]
    data_bytes = all_raw[4:4+hdr_size]
    if len(data_bytes) == 0:
        # 可能没有 header, 直接是数据
        data_bytes = all_raw
    return np.frombuffer(data_bytes, dtype=dtype)


def read_vti_phi(filepath):
    """读取 VTI 文件中的 PHI 字段, 返回 (phi_array, local_ext, origin).
    直接用正则提取 base64 内容, 避免 XML parser 截断.
    返回:
      phi: (ny, nx) numpy 数组
      local_ext: [xmin,xmax,ymin,ymax,0,0] 块内局部坐标
      origin: (Ox, Oy) 全局块原点位置 (cell 0 左下角的全局坐标)
    """
    with open(filepath, 'rb') as f:
        data = f.read()

    # 用正则提取 WholeExtent
    m = re.search(rb'WholeExtent="([^"]+)"', data)
    if not m:
        return None
    ext = list(map(int, m.group(1).split()))
    nx = ext[1]-ext[0]+1
    ny = ext[3]-ext[2]+1

    # 提取 Origin (块的全局起点位置, 左下角坐标)
    m_org = re.search(rb'Origin="([^"]+)"', data)
    if not m_org:
        return None
    org_parts = list(map(float, m_org.group(1).split()))
    Ox, Oy = org_parts[0], org_parts[1]

    # 提取 PHI DataArray 内容
    m = re.search(rb'<DataArray[^>]*Name="PHI"[^>]*>([^<]+)</DataArray>', data)
    if not m:
        return None
    b64_content = m.group(1).strip()

    # 检测类型
    m_type = re.search(rb'<DataArray[^>]*Name="PHI"[^>]*type="([^"]+)"', data)
    dtype_str = m_type.group(1).decode() if m_type else 'Float64'
    dtype = '<f8' if dtype_str == 'Float64' else '<f4'

    try:
        values = decode_base64_binary(b64_content, dtype)
    except Exception:
        return None

    if len(values) < nx*ny:
        return None
    phi = values[:nx*ny].reshape(ny, nx)
    return phi, ext, (Ox, Oy)


def find_interface(phi, ext):
    """找 phi=0.5 等高线位置."""
    ny, nx = phi.shape
    interface_y = []
    for i in range(nx):
        col = phi[:, i]
        crossings = []
        for j in range(ny-1):
            if (col[j]-0.5)*(col[j+1]-0.5) < 0:
                frac = (0.5-col[j])/(col[j+1]-col[j])
                y_pos = (ext[2]+j) + frac
                crossings.append(y_pos)
        if crossings:
            interface_y.append(crossings[len(crossings)//2])
        else:
            interface_y.append(float('nan'))
    return np.array(interface_y)


def read_vtm_all_phi(vtm_path, vti_dir):
    """读取 VTM 文件引用的所有 VTI, 合并成全局 phi.
    使用每个 VTI 的 Origin 字段计算全局位置:
      global_gx = round(Ox + 0.5) + local_i
      global_gy = round(Oy + 0.5) + local_j
    (因为 VTK Origin 是 cell 0 左下角坐标, Spacing=1, cell i 中心在 Origin+(i+0.5))
    """
    try:
        tree = ET.parse(vtm_path)
    except:
        return None
    # 收集所有 VTI 文件 (去重)
    vti_files = []
    seen = set()
    for ds in tree.iter('DataSet'):
        f = ds.get('file')
        if f and f.endswith('.vti') and f not in seen:
            vti_files.append(f)
            seen.add(f)
    if not vti_files:
        return None

    # 读取所有 VTI
    all_phi = []
    for vf in vti_files:
        vti_path = os.path.join(vti_dir, vf)
        if not os.path.exists(vti_path):
            continue
        result = read_vti_phi(vti_path)
        if result is None:
            continue
        phi, ext, (Ox, Oy) = result
        # 全局 cell 索引 (左下角)
        gx0 = int(round(Ox + 0.5))
        gy0 = int(round(Oy + 0.5))
        all_phi.append((phi, gx0, gy0))

    if not all_phi:
        return None

    # 找全局域范围
    min_gx = min(gx for _, gx, _ in all_phi)
    max_gx = max(gx + phi.shape[1] - 1 for phi, gx, _ in all_phi)
    min_gy = min(gy for _, _, gy in all_phi)
    max_gy = max(gy + phi.shape[0] - 1 for phi, _, gy in all_phi)

    gnx = max_gx - min_gx + 1
    gny = max_gy - min_gy + 1
    global_phi = np.full((gny, gnx), np.nan)
    for phi, gx, gy in all_phi:
        ny, nx = phi.shape
        x0 = gx - min_gx
        y0 = gy - min_gy
        global_phi[y0:y0+ny, x0:x0+nx] = phi

    return global_phi, [min_gx, max_gx, min_gy, max_gy]


def main():
    vti_dir = sys.argv[1] if len(sys.argv) > 1 else 'vtkoutput/vtidata'
    L_override = float(sys.argv[2]) if len(sys.argv) > 2 else None

    cfg = parse_ini()
    Ni = cfg['Ni']
    Cell_Len = cfg['Cell_Len']
    L = L_override if L_override else Ni*Cell_Len*0.5
    U = cfg['U']
    Fr = cfg['Fr']
    Bom = cfg['Bom']

    g = U**2/(Fr*L)
    t_star_factor = math.sqrt(g/L)

    # 理论增长率 (无磁场)
    wavelength = L
    k = 2*math.pi/wavelength
    deltaU = U
    gamma_theory = k*deltaU/2

    print("="*70)
    print(f"KH 不稳定性增长率分析")
    print("="*70)
    print(f"VTI 目录: {vti_dir}")
    print(f"网格: {Ni}x{Ni}, L={L}")
    print(f"U={U}, Fr={Fr}, g={g:.6e}, Bom={Bom}")
    print(f"理论增长率 (无磁场): γ={gamma_theory:.6e} /步  (τ={1/gamma_theory:.0f}步)")
    print(f"t* 因子 = {t_star_factor:.6e}, t*=6 对应 ~{6/t_star_factor:.0f} 步")
    print("-"*70)
    print(f"{'步数':>6} | {'t*':>7} | {'振幅':>8} | {'峰位':>8} | {'谷位':>8} | {'|H|max估计':>10}")
    print("-"*70)

    # 找所有 VTM 文件
    vtms = sorted([f for f in os.listdir(vti_dir) if f.endswith('.vtm')],
                  key=lambda x: int(re.search(r'kh2d(\d+)\.vtm', x).group(1)))

    results = []
    for vtm_name in vtms:
        m = re.search(r'kh2d(\d+)\.vtm', vtm_name)
        if not m: continue
        step = int(m.group(1))
        vtm_path = os.path.join(vti_dir, vtm_name)
        result = read_vtm_all_phi(vtm_path, vti_dir)
        if result is None:
            print(f"{step:6d} | 读取失败")
            continue
        global_phi, ext = result
        interface_y = find_interface(global_phi, ext)
        valid_y = interface_y[~np.isnan(interface_y)]
        if len(valid_y) > 0:
            amplitude = (np.max(valid_y)-np.min(valid_y))/2
            peak = np.max(valid_y)
            trough = np.min(valid_y)
        else:
            amplitude = peak = trough = float('nan')
        t_star = step*t_star_factor
        # 估计 |H|max: 从 sim log 知道 max 约 5.6e-2 (Bom=636, 256x256)
        # 此处仅打印占位
        hmax_str = "-"
        print(f"{step:6d} | {t_star:7.4f} | {amplitude:8.2f} | {peak:8.2f} | {trough:8.2f} | {hmax_str:>10}")
        results.append({'step':step, 't_star':t_star, 'amplitude':amplitude,
                       'peak':peak, 'trough':trough})

    print("-"*70)
    if len(results) >= 2:
        # 计算实际增长率 (取前几个点, 线性区)
        early = [r for r in results if r['step'] <= 3000 and not math.isnan(r['amplitude']) and r['amplitude'] > 1e-3]
        if len(early) >= 2:
            steps = np.array([r['step'] for r in early], dtype=float)
            amps = np.array([r['amplitude'] for r in early])
            log_amps = np.log(amps)
            fit = np.polyfit(steps, log_amps, 1)
            gamma_actual = fit[0]
            print(f"\n实际增长率 (线性区): γ_actual = {gamma_actual:.6e} /步")
            print(f"理论增长率 (无磁场): γ_theory  = {gamma_theory:.6e} /步")
            print(f"比值: γ_actual/γ_theory = {gamma_actual/gamma_theory:.3f}x")
            print(f"稳定化效应: {(1-gamma_actual/gamma_theory)*100:.1f}% 减速")

            # 推算 t*=6 所需步数
            # 振幅从 0.1*L 开始, t*=6 时振幅约 L (达到饱和)
            amp_0 = 0.1*L
            amp_target = 0.5*L  # 振幅达到 0.5L 时认为接近 t*=6
            steps_to_t6 = math.log(amp_target/amp_0)/gamma_actual
            print(f"\n推算 t*=6 (振幅达 0.5L) 所需步数: ~{steps_to_t6:.0f} 步")
            print(f"当前 20000 步对应 t* = {20000*t_star_factor:.3f}")
            amp_at_20k = amp_0 * math.exp(gamma_actual*20000)
            print(f"推算 20000 步时振幅: {amp_at_20k:.2f} = {amp_at_20k/L*100:.1f}% L")


if __name__ == '__main__':
    main()
