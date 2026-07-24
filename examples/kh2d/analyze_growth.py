#!/usr/bin/env python3
"""分析KH不稳定性的界面增长率，与理论值比较"""
import os
import struct
import math
import re
import numpy as np

def read_vti_phi(filepath):
    """读取VTI文件中的PHI字段"""
    with open(filepath, 'rb') as f:
        data = f.read()

    # 解析XML头部
    text = data[:5000].decode('ascii', errors='ignore')

    # 提取WholeExtent
    m = re.search(r'WholeExtent="([^"]+)"', text)
    if not m:
        return None
    ext = list(map(int, m.group(1).split()))
    nx = ext[1] - ext[0] + 1
    ny = ext[3] - ext[2] + 1

    # 找到PHI scalar的offset和size
    # 寻找 appended data 标记
    m_offset = re.search(r'offset="(\d+)"', text)

    # 寻找 Float64 标记
    if 'Float64' in text:
        dtype_size = 8
        dtype_fmt = 'd'
    else:
        dtype_size = 4
        dtype_fmt = 'f'

    # 找到数据起始位置（appended data after XML header）
    # 寻找第一个 <AppendedData encoding="raw"> 标记
    m_appended = re.search(r'<AppendedData encoding="raw">', text)
    if m_appended:
        data_start = m_appended.end()
        # 跳过可能的换行符
        while data_start < len(data) and data[data_start:data_start+1] in [b'\n', b'\r']:
            data_start += 1
    else:
        # 尝试内联数据
        m_data = re.search(r'>\s*_', text)
        if m_data:
            data_start = text.find('_', m_data.start()) + 1
        else:
            return None

    # 读取数据 - 先读header(4 bytes for size in appended mode)
    # VTK appended data format: 4-byte size header + data
    if m_offset:
        offset = int(m_offset.group(1))
        # 在appended data区域中，offset是相对于appended data开始的
        actual_offset = data_start + offset
        # 读取4字节size header
        if actual_offset + 4 > len(data):
            return None
        size = struct.unpack('i', data[actual_offset:actual_offset+4])[0]
        actual_data_start = actual_offset + 4
    else:
        actual_data_start = data_start
        size = nx * ny * dtype_size

    n_vals = size // dtype_size
    if n_vals != nx * ny:
        # 可能包含overlap
        n_vals = min(n_vals, nx * ny)

    values = struct.unpack(f'{n_vals}{dtype_fmt}', data[actual_data_start:actual_data_start+n_vals*dtype_size])
    phi = np.array(values[:nx*ny]).reshape(ny, nx)
    return phi, ext

def find_interface(phi, ext):
    """找到phi=0.5的等高线位置(界面位置)"""
    ny, nx = phi.shape
    # 对于每一列x，找到phi最接近0.5的y位置
    interface_y = []
    for i in range(nx):
        col = phi[:, i]
        # 找到phi从>0.5变为<0.5的位置(或反之)
        crossings = []
        for j in range(ny-1):
            if (col[j] - 0.5) * (col[j+1] - 0.5) < 0:
                # 线性插值找到精确位置
                frac = (0.5 - col[j]) / (col[j+1] - col[j])
                y_pos = (ext[2] + j) + frac
                crossings.append(y_pos)
        if crossings:
            # 取中间的crossing（界面）
            interface_y.append(crossings[len(crossings)//2])
        else:
            interface_y.append(float('nan'))
    return np.array(interface_y)

def main():
    vtidir = "vtkoutput_bom636_fixed/vtidata"
    if not os.path.exists(vtidir):
        print(f"Error: {vtidir} not found")
        return

    # 找到所有VTM文件
    vtms = sorted([f for f in os.listdir(vtidir) if f.endswith('.vtm')])
    if not vtms:
        print("No VTM files found")
        return

    # 参数
    U = 0.02
    Fr = 1.0
    Ni = 256
    Cell_Len = 1.0
    L = Ni * Cell_Len * 0.5  # 半域 = 128

    g = U**2 / (Fr * L)
    t_star_factor = math.sqrt(g / L)

    # 理论增长率 (无磁场)
    wavelength = L  # 扰动波长 = L
    k = 2 * math.pi / wavelength
    deltaU = U  # 剪切速度 = U
    gamma = k * deltaU / 2  # KH增长率

    print(f"=== KH Instability Growth Rate Analysis ===")
    print(f"Domain: {Ni}x{Ni}, L={L}")
    print(f"U={U}, Fr={Fr}, g={g:.6e}")
    print(f"Wavelength={wavelength}, k={k:.6f}")
    print(f"Shear velocity ΔU={deltaU}")
    print(f"Theoretical growth rate γ={gamma:.6e} per step")
    print(f"E-folding time τ={1/gamma:.1f} steps")
    print(f"t* factor = {t_star_factor:.6e}")
    print()

    # 分析每个时间步
    results = []
    for vtm_name in vtms:
        # 提取步数
        m = re.search(r'kh2d(\d+)\.vtm', vtm_name)
        if not m:
            continue
        step = int(m.group(1))

        # 读取VTM文件获取block列表
        vtm_path = os.path.join(vtidir, vtm_name)
        with open(vtm_path, 'r') as f:
            vtm_content = f.read()

        # 找到所有VTI文件引用
        vti_refs = re.findall(r'kh2d_T\d+_B\d+\.vti', vtm_content)

        # 合并所有block的phi数据
        all_phi = {}
        for vti_name in vti_refs:
            vti_path = os.path.join(vtidir, vti_name)
            if not os.path.exists(vti_path):
                continue
            result = read_vti_phi(vti_path)
            if result is None:
                continue
            phi, ext = result

            # 提取block id
            m_b = re.search(r'_B(\d+)\.vti', vti_name)
            if m_b:
                block_id = int(m_b.group(1))
                all_phi[block_id] = (phi, ext)

        if not all_phi:
            continue

        # 找到全局域范围
        min_x = min(ext[0] for _, ext in all_phi.values())
        max_x = max(ext[1] for _, ext in all_phi.values())
        min_y = min(ext[2] for _, ext in all_phi.values())
        max_y = max(ext[3] for _, ext in all_phi.values())

        global_nx = max_x - min_x + 1
        global_ny = max_y - min_y + 1

        # 合并到全局数组
        global_phi = np.full((global_ny, global_nx), np.nan)
        for phi, ext in all_phi.values():
            nx = ext[1] - ext[0] + 1
            ny = ext[3] - ext[2] + 1
            x0 = ext[0] - min_x
            y0 = ext[2] - min_y
            global_phi[y0:y0+ny, x0:x0+nx] = phi

        # 找界面
        interface_y = find_interface(global_phi, [min_x, max_x, min_y, max_y])

        # 计算界面振幅 (max - min) / 2
        valid_y = interface_y[~np.isnan(interface_y)]
        if len(valid_y) > 0:
            amplitude = (np.max(valid_y) - np.min(valid_y)) / 2
            mean_pos = np.mean(valid_y)
            peak_pos = np.max(valid_y)
            trough_pos = np.min(valid_y)
        else:
            amplitude = float('nan')
            mean_pos = float('nan')
            peak_pos = float('nan')
            trough_pos = float('nan')

        t_star = step * t_star_factor
        theoretical_amp = 0.1 * L * math.exp(gamma * step) if step > 0 else 0.1 * L

        results.append({
            'step': step,
            't_star': t_star,
            'amplitude': amplitude,
            'peak': peak_pos,
            'trough': trough_pos,
            'mean': mean_pos,
            'theory_amp': theoretical_amp,
        })

        print(f"Step {step:6d} | t*={t_star:7.4f} | amp={amplitude:8.2f} | "
              f"peak={peak_pos:8.2f} | trough={trough_pos:8.2f} | "
              f"theory_amp={theoretical_amp:8.2f}")

    print()
    if len(results) >= 2:
        # 计算实际增长率
        r0 = results[0]
        r1 = results[-1]
        if r0['amplitude'] > 0 and r1['amplitude'] > 0:
            dt = r1['step'] - r0['step']
            if dt > 0 and r0['amplitude'] > 0:
                growth_ratio = r1['amplitude'] / r0['amplitude']
                actual_gamma = math.log(growth_ratio) / dt if growth_ratio > 0 else 0
                print(f"=== Growth Rate Comparison ===")
                print(f"Steps: {r0['step']} → {r1['step']} (Δt={dt})")
                print(f"Amplitude: {r0['amplitude']:.2f} → {r1['amplitude']:.2f} (ratio={growth_ratio:.4f})")
                print(f"Actual growth rate: γ_actual={actual_gamma:.6e} per step")
                print(f"Theory growth rate: γ_theory ={gamma:.6e} per step")
                print(f"Ratio: γ_actual/γ_theory = {actual_gamma/gamma:.2f}x")
                if actual_gamma > 2 * gamma:
                    print("⚠️  WARNING: Growth rate is MUCH higher than theory!")
                elif actual_gamma < 0.5 * gamma:
                    print("⚠️  WARNING: Growth rate is MUCH lower than theory!")
                else:
                    print("✓ Growth rate matches theory (within 2x)")

if __name__ == '__main__':
    main()
