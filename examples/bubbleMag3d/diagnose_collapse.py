#!/usr/bin/env python3
# diagnose_collapse.py — 分析 bubbleMag3d 气泡崩溃的诊断脚本
#
# 用法:
#   python3 diagnose_collapse.py [vtkdir] [--steps 0,20,25,30]
#
# 功能:
#   1. 读取 vtm/vti (base64 二进制格式), 组装 128 个 block 为完整域
#   2. 每步检查: NaN/Inf, PHI 越界, 气泡体积/质心/球形度, |H| 尖峰,
#      PSI 远场线性度, 速度/Ma, 密度范围
#   3. 周期性接缝检查: x/y 边界 ghost 平面 vs 对侧物理平面 (128进程跨rank同步是否失效)
#   4. 块间连续性: 相邻 block 共享面的场跳变
#   5. 高频振荡 (棋盘格) 指标: HMAG 的逐格跳变
import re, base64, glob, sys, argparse
import numpy as np

CS2 = 1.0/3.0
CS  = np.sqrt(CS2)

def read_vti(path):
    """解析 vti: 4字节长度前缀(base64, 恒为8字符含'==') + payload(base64)"""
    txt = open(path, 'rb').read()
    ext = re.findall(rb'Extent="(\d+) (\d+) (\d+) (\d+) (\d+) (\d+)"', txt)[0]
    shp = tuple(int(b)-int(a)+1 for a, b in zip(ext[::2], ext[1::2]))
    origin = [float(v) for v in re.findall(rb'Origin="([-\d.]+) ([-\d.]+) ([-\d.]+)"', txt)[0]]
    fields = {}
    for m in re.findall(rb'<DataArray\s+([^>]*)>\s*([A-Za-z0-9+/=\r\n]+)\s*</DataArray>', txt):
        attrs, b64 = m
        nm = re.search(rb'Name="([^"]+)"', attrs).group(1).decode()
        ncomp = int(re.search(rb'NumberOfComponents="(\d+)"', attrs).group(1))
        s = b64.decode('ascii').strip()
        size = int.from_bytes(base64.b64decode(s[:8]), 'little')
        payload = base64.b64decode(s[8:] + '=' * (-len(s[8:]) % 4))
        arr = np.frombuffer(payload[:size], dtype='<f8')
        if ncomp == 1:
            fields[nm] = arr.reshape(shp)
        else:
            fields[nm] = arr.reshape(shp + (ncomp,))
    return fields, origin

def assemble(vtkdir, step, Ni, Nj, Nk):
    """组装所有 block 的物理单元 (去掉 1 层 ghost) 到全局数组; 同时保留 ghost 平面用于接缝检查"""
    g = {}
    gfull = {}
    layout = []  # (block_id, gx0, gy0, gz0, nx, ny, nz)
    for p in sorted(glob.glob(f'{vtkdir}/bubbleMag3d_T{step}_B*.vti')):
        bid = int(re.search(r'_B(\d+)\.vti', p).group(1))
        f, origin = read_vti(p)
        gx0, gy0, gz0 = (int(o + 0.5) for o in origin)
        shp = f['PHI'].shape
        for nm, arr in f.items():
            if nm not in g:
                g[nm] = np.zeros((Nk, Nj, Ni) if arr.ndim == 3 else (Nk, Nj, Ni, arr.shape[-1]))
            g[nm][gz0:gz0+shp[2]-2, gy0:gy0+shp[1]-2, gx0:gx0+shp[0]-2] = arr[1:-1, 1:-1, 1:-1]
            if nm in ('PHI','PSI','HX','HY','HZ','HMAG'):
                gfull.setdefault(nm, {})[bid] = arr
        layout.append((bid, gx0, gy0, gz0, shp[0], shp[1], shp[2]))
    return g, gfull, layout

def report_stats(g, step, H0):
    print(f"==================== step = {step} ====================")
    # --- NaN / Inf ---
    bad = 0
    for nm, arr in g.items():
        n = int(np.isnan(arr).sum() + np.isinf(arr).sum())
        bad += n
        if n: print(f"  [!!] {nm}: {n} NaN/Inf")
    print(f"  NaN/Inf total: {bad}")
    # --- PHI ---
    phi = g['PHI']
    print(f"  PHI    : min={phi.min():.4f} max={phi.max():.4f}  "
          f"overshoot(phi<0|phi>1)={np.sum((phi<0)|(phi>1))}  "
          f"mean={phi.mean():.4f}")
    vol = np.sum(phi < 0.5)
    # 质心 (气泡 = phi<0.5)
    idx = np.argwhere(phi < 0.5)
    if len(idx) > 0:
        com = idx.mean(axis=0)
        # 球形度: 到质心的距离分布
        d = np.linalg.norm(idx - com, axis=1)
        Rmean, Rstd = d.mean(), d.std()
        print(f"  Bubble : vol={vol} cells  COM=({com[2]:.1f},{com[1]:.1f},{com[0]:.1f})  "
              f"Rmean={Rmean:.2f} Rstd={Rstd:.3f} (Rstd/Rmean={Rstd/Rmean:.3f})")
        # 沿三个方向的延伸半径 (collapse/拉长检测)
        for ax, name in [(0,'z'), (1,'y'), (2,'x')]:
            span = idx[:,ax].max() - idx[:,ax].min() + 1
            print(f"      span_{name} = {span}")
    else:
        print("  Bubble : 不存在 (已完全崩溃!)")
    # --- PSI 远场线性度 (取远离气泡的角落柱) ---
    psi = g['PSI']
    z = np.arange(psi.shape[0])
    # 在 (x,y)=(4,4) 角落的柱: 应满足 psi ≈ -H0*(z+0.5)
    col = psi[2:40, 4, 4]
    zc = z[2:40] + 0.5
    if abs(H0) > 0:
        dev = np.abs(col + H0*zc).max()
        print(f"  PSI    : min={psi.min():.5f} max={psi.max():.5f}  "
              f"far-field linearity dev={dev:.5f} ({dev/max(abs(H0),1e-12)*100:.2f}% of H0)")
    # --- H 场 ---
    hm = g['HMAG']
    imax = np.unravel_index(np.argmax(hm), hm.shape)
    print(f"  HMAG   : max={hm.max():.5f} ({hm.max()/max(H0,1e-30):.2f}xH0) at (z,y,x)={imax}  "
          f"mean={hm.mean():.5f}")
    # 逐格最大跳变 (棋盘格指标)
    dj = max(np.abs(np.diff(hm, axis=0)).max(),
             np.abs(np.diff(hm, axis=1)).max(),
             np.abs(np.diff(hm, axis=2)).max())
    print(f"  HMAG   : 最大逐格跳变 = {dj:.5f} ({dj/max(H0,1e-30):.2f}xH0)")
    # --- 速度 ---
    vel = g['Velocity']
    umag = np.sqrt((vel**2).sum(-1))
    print(f"  |u|    : max={umag.max():.5f}  Ma={umag.max()/CS:.4f}  mean={umag.mean():.6f}")
    # --- 密度 ---
    rho = g['Density']
    print(f"  rho    : min={rho.min():.4f} max={rho.max():.4f}")
    return

def seam_check(gfull, layout, step, H0):
    """周期性接缝检查: x/y 边界 block 的 ghost 平面 vs 对侧物理平面
    若 128 进程跨 rank 时 SyncMFPeriodicGhosts 失效, 这里会暴露 PSI/H ghost 陈旧"""
    print(f"----- 周期性接缝检查 (step={step}) -----")
    # 建立 block 索引: (gx0//32, gy0//32, gz0//32) -> bid
    bmap = {}
    for bid, gx0, gy0, gz0, nx, ny, nz in layout:
        bmap[(gx0//32, gy0//32, gz0//32)] = bid
    NX, NY, NZ = 4, 4, 8
    for nm in ['PSI', 'HX', 'HY', 'HZ', 'HMAG']:
        errs = []
        # x 方向: (0,y,z) 的 ghost i=0 vs (3,y,z) 的物理 i=32 ; (3) 的 ghost i=33 vs (0) 的物理 i=1
        for zz in range(NZ):
            for yy in range(NY):
                bL = bmap[(0, yy, zz)]; bR = bmap[(NX-1, yy, zz)]
                aL, aR = gfull[nm][bL], gfull[nm][bR]
                errs.append(np.abs(aL[1:-1, 1:-1, 0] - aR[1:-1, 1:-1, 32]).max())
                errs.append(np.abs(aR[1:-1, 1:-1, 33] - aL[1:-1, 1:-1, 1]).max())
        # y 方向: (x,0,z) 的 ghost j=0 vs (x,3,z) 的物理 j=32
        for zz in range(NZ):
            for xx in range(NX):
                bF = bmap[(xx, 0, zz)]; bB = bmap[(xx, NY-1, zz)]
                aF, aB = gfull[nm][bF], gfull[nm][bB]
                errs.append(np.abs(aF[1:-1, 0, 1:-1] - aB[1:-1, 32, 1:-1]).max())
                errs.append(np.abs(aB[1:-1, 33, 1:-1] - aF[1:-1, 1, 1:-1]).max())
        e = max(errs)
        ref = max(abs(H0), np.abs(gfull[nm][list(gfull[nm].keys())[0]]).max())
        print(f"  {nm:5s}: 周期性接缝最大偏差 = {e:.5f}  ({e/ref*100:.2f}% of ref={ref:.5f})")
    return

def block_continuity(gfull, layout, step):
    """块间连续性: 相邻 block 重叠 ghost 与物理单元的最大差值 (MPI ghost 同步是否正常)"""
    print(f"----- 块间连续性 (step={step}) -----")
    bmap = {}
    for bid, gx0, gy0, gz0, nx, ny, nz in layout:
        bmap[(gx0//32, gy0//32, gz0//32)] = bid
    for nm in ['PHI', 'PSI', 'HMAG']:
        worst = 0.0
        for (ix, iy, iz), bid in bmap.items():
            a = gfull[nm][bid]
            # +x 邻居: 本块 ghost i=33 vs 邻居物理 i=1
            if (ix+1, iy, iz) in bmap:
                b = gfull[nm][bmap[(ix+1, iy, iz)]]
                worst = max(worst, np.abs(a[1:-1, 1:-1, 33] - b[1:-1, 1:-1, 1]).max())
            if (ix, iy+1, iz) in bmap:
                b = gfull[nm][bmap[(ix, iy+1, iz)]]
                worst = max(worst, np.abs(a[1:-1, 33, 1:-1] - b[1:-1, 1, 1:-1]).max())
            if (ix, iy, iz+1) in bmap:
                b = gfull[nm][bmap[(ix, iy, iz+1)]]
                worst = max(worst, np.abs(a[33, 1:-1, 1:-1] - b[1, 1:-1, 1:-1]).max())
        print(f"  {nm:5s}: 相邻块共享面最大跳变 = {worst:.5e}")
    return

def interface_slice(phi, step, H0):
    """中心切片: 气泡轮廓提取 + 界面厚度估计 + 轮廓径向均匀性"""
    zc, yc, xc = phi.shape[0]//2, phi.shape[1]//2, phi.shape[2]//2
    sl = phi[:, yc, :]  # x-z 平面
    print(f"----- 界面分析 (step={step}, y={yc} 中心切片 x-z) -----")
    # 界面厚度: |∇φ| 峰值区宽度
    gx = np.abs(np.diff(sl, axis=1)); gz = np.abs(np.diff(sl, axis=0))
    gmag = np.zeros_like(sl)
    gmag[:, 1:] += gx; gmag[:, :-1] += gx
    gmag[1:, :] += gz; gmag[:-1, :] += gz
    gmag *= 0.5
    print(f"  |∇φ|max={gmag.max():.4f}  (理想tanh界面 W=2: max|∇φ|≈1/W=0.5)")
    # 界面上的 phi 分布
    iface = (phi > 0.05) & (phi < 0.95)
    print(f"  界面单元数={iface.sum()}  界面内 phi 范围=[{phi[iface].min():.3f},{phi[iface].max():.3f}]")
    # 径向剖面: 以 (xc, zc) 为中心
    zz, xx = np.mgrid[0:phi.shape[0], 0:phi.shape[2]]
    r = np.sqrt(((zz - zc)/1.0)**2 + ((xx - xc)/1.0)**2)
    rbin = np.linspace(0, r.max(), 60)
    prof = np.zeros(len(rbin)-1); profn = np.zeros(len(rbin)-1)
    ri = np.digitize(r.ravel(), rbin)
    phif = phi[:, yc, :].ravel()
    for i in range(1, len(rbin)):
        m = ri == i
        if m.sum():
            prof[i-1] = phif[m].mean()
            profn[i-1] = m.sum()
    # 找 phi=0.5 处半径
    cross = np.where(np.diff(prof > 0.5))[0]
    print(f"  径向剖面 phi=0.5 交叉点(距中心): "
          f"{[round((rbin[i]+rbin[i+1])/2,1) for i in cross][:6]}")
    return

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('vtkdir', nargs='?', default='./vtkoutput/vtidata')
    ap.add_argument('--steps', default='0,20,25,30')
    ap.add_argument('--mesh', default='128,128,256')
    args = ap.parse_args()
    Ni, Nj, Nk = (int(v) for v in args.mesh.split(','))
    H0 = 0.0200  # 由 Eo=20,Re=40,Bom=1.94,R=25 计算 (sqrt(2*Bom*sigma/D))
    steps = [int(s) for s in args.steps.split(',')]
    allg = {}
    for s in steps:
        g, gfull, layout = assemble(args.vtkdir, s, Ni, Nj, Nk)
        allg[s] = (g, gfull, layout)
        report_stats(g, s, H0)
        seam_check(gfull, layout, s, H0)
        block_continuity(gfull, layout, s)
        interface_slice(g['PHI'], s, H0)
    # 气泡演化汇总
    print("==================== 气泡演化汇总 ====================")
    for s in steps:
        phi = allg[s][0]['PHI']
        vol = np.sum(phi < 0.5)
        idx = np.argwhere(phi < 0.5)
        com = idx.mean(axis=0) if len(idx) else np.array([-1,-1,-1])
        print(f"  step={s:4d}  vol={vol:8d}  COM_z={com[0]:7.2f}  "
              f"phi_max={phi.max():.4f} phi_min={phi.min():.4f}")

if __name__ == '__main__':
    main()
