#!/usr/bin/env python3
# validate_field.py — 磁场解析解验证
# 对比 field_validate 输出的 ψ/H 与球体匀强磁场解析解:
#   μ_in=1 (气泡), μ_out=9 (铁磁流体), 均匀场 H0 沿 +z
#   球内: H = 3μ_out/(μ_in+2μ_out)·H0 ẑ (均匀)
#   球外: ψ = -H0·z + C·H0·R³·z/r³, C=(μ_in-μ_out)/(μ_in+2μ_out)
#   极轴: H_z = H0·(1 - 2|C|(R/r)³), 赤道: H_r = H0·(1 + |C|(R/r)³)
# 用法: python3 validate_field.py [vtkdir] [--mesh 64,64,128] [--radius 12] [--center 32,32,48] [--h0 0.0417]
import re, base64, glob, sys, argparse
import numpy as np

def read_vti(path):
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
        fields[nm] = arr.reshape(shp) if ncomp == 1 else arr.reshape(shp + (ncomp,))
    return fields, origin

def assemble(vtkdir, step, Ni, Nj, Nk):
    g = {}
    for p in sorted(glob.glob(f'{vtkdir}/fieldcheck_T{step}_B*.vti')):
        f, origin = read_vti(p)
        gx0, gy0, gz0 = (int(o + 0.5) for o in origin)
        shp = f['PHI'].shape
        for nm, arr in f.items():
            if nm not in g:
                g[nm] = np.zeros((Nk, Nj, Ni) if arr.ndim == 3 else (Nk, Nj, Ni, arr.shape[-1]))
            g[nm][gz0:gz0+shp[2]-2, gy0:gy0+shp[1]-2, gx0:gx0+shp[0]-2] = arr[1:-1, 1:-1, 1:-1]
    return g

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('vtkdir', nargs='?', default='./vtkoutput/vtidata')
    ap.add_argument('--mesh', default='64,64,128')
    ap.add_argument('--radius', type=float, default=12.0)
    ap.add_argument('--center', default='32,32,48')
    ap.add_argument('--h0', type=float, default=0.0)
    ap.add_argument('--muin', type=float, default=1.0)
    ap.add_argument('--muout', type=float, default=9.0)
    ap.add_argument('--step', type=int, default=0)
    args = ap.parse_args()
    Ni, Nj, Nk = (int(v) for v in args.mesh.split(','))
    xc, yc, zc = (float(v) for v in args.center.split(','))
    R = args.radius

    g = assemble(args.vtkdir, args.step, Ni, Nj, Nk)
    phi = g['PHI']; psi = g['PSI']; Hx = g['HX']; Hy = g['HY']; Hz = g['HZ']

    # 球心位置 (最近格点)
    iz, iy, ix = int(round(zc)), int(round(yc)), int(round(xc))
    zz, yy, xx = np.mgrid[0:Nk, 0:Nj, 0:Ni]
    r = np.sqrt((zz-zc)**2 + (yy-yc)**2 + (xx-xc)**2)

    # H0: 用远场实测 (角落均值) 或参数
    if args.h0 > 0:
        H0 = args.h0
    else:
        H0 = np.abs(Hz[2:8, 2:8, 2:8]).mean()
    print(f"H0 = {H0:.6f}")

    C = (args.muin - args.muout) / (args.muin + 2*args.muout)
    H_in = 3*args.muout/(args.muin + 2*args.muout) * H0
    print(f"C = {C:.5f}, 解析 H_in = {H_in:.6f} = {H_in/H0:.5f}·H0")

    # ---- 1. 球内均匀性 (气泡内部 φ<0.05, r<R-2) ----
    inside = (r < R - 2.0) & (phi < 0.05)  # 气泡内部 (μ_in=1, φ=0)
    if inside.sum():
        Hz_in = Hz[inside]
        err_in = np.abs(Hz_in - H_in).max() / H0
        print(f"\n[1] 球内 (r<R-2, φ<0.05, {inside.sum()}格): "
              f"H_z 实测 mean={Hz_in.mean():.5f} (±{Hz_in.std():.5f})  "
              f"解析 {H_in:.5f}  相对误差 max={err_in*100:.2f}%")
        Hmag_in = np.sqrt(Hx[inside]**2+Hy[inside]**2+Hz[inside]**2)
        print(f"    |H| 实测 {Hmag_in.mean():.5f}  横向 H_x/H_y max="
              f"{np.abs(Hx[inside]).max():.2e}/{np.abs(Hy[inside]).max():.2e} (应≈0)")

    # ---- 2. 极轴 (x=y=0 柱): H_z vs 解析 ----
    axis = (np.abs(xx-xc) < 0.6) & (np.abs(yy-yc) < 0.6) & (np.abs(zz-zc) > R)
    ra = np.abs(zz[axis] - zc)
    Hz_a = Hz[axis]; phia = phi[axis]
    sel = (ra > R + 2.0) & (ra < 40)  # 远离界面与壁面
    if sel.sum():
        ra_s = ra[sel]; Hz_s = Hz_a[sel]
        ana = H0 * (1 - 2*abs(C) * (R/ra_s)**3)
        err_axis = np.abs(Hz_s - ana).max() / H0
        print(f"\n[2] 极轴 (r={ra_s.min():.1f}..{ra_s.max():.1f}): "
              f"H_z 实测 vs 解析 最大偏差={err_axis*100:.2f}%H0")
        for i in range(0, len(ra_s), max(1, len(ra_s)//6)):
            print(f"    r={ra_s[i]:6.1f}: 实测 {Hz_s[i]:.5f}  解析 {ana[i]:.5f}")

    # ---- 3. 赤道 (z=zc 平面): H_r vs 解析 ----
    eq = (np.abs(zz-zc) < 0.6) & (r > R)
    re_ = r[eq]; phie = phi[eq]
    Hr_e = np.sqrt(Hx[eq]**2 + Hy[eq]**2 + Hz[eq]**2)  # 赤道处 H 近似径向
    sel2 = (re_ > R + 2.0) & (re_ < 40)
    if sel2.sum():
        re_s = re_[sel2]; Hr_s = Hr_e[sel2]
        ana2 = H0 * (1 + abs(C) * (R/re_s)**3)
        err_eq = np.abs(Hr_s - ana2).max() / H0
        print(f"\n[3] 赤道 (r={re_s.min():.1f}..{re_s.max():.1f}): "
              f"|H| 实测 vs 解析 最大偏差={err_eq*100:.2f}%H0")
        for i in range(0, len(re_s), max(1, len(re_s)//6)):
            print(f"    r={re_s[i]:6.1f}: 实测 {Hr_s[i]:.5f}  解析 {ana2[i]:.5f}")

    # ---- 4. 界面边界条件: H 切向连续 + μH 法向连续 ----
    print("\n[4] 界面边界条件 (r 在 R±1.5 的壳层):")
    shell = (np.abs(r - R) < 1.5)
    if shell.sum():
        # 法向 n = r̂; 切向 = H - (H·n)n
        nx_ = (xx-xc)/r; ny_ = (yy-yc)/r; nz_ = (zz-zc)/r
        Hn = Hx*nx_ + Hy*ny_ + Hz*nz_
        mu = args.muin + phi*(args.muout - args.muin)
        # 取界面两侧各一层
        in_shell = shell & (r < R)
        out_shell = shell & (r > R)
        muHn_in = (mu*Hn)[in_shell]; muHn_out = (mu*Hn)[out_shell]
        Ht_in = np.sqrt((Hx**2+Hy**2+Hz**2 - Hn**2))[in_shell]
        Ht_out = np.sqrt((Hx**2+Hy**2+Hz**2 - Hn**2))[out_shell]
        print(f"    μH_n: 内侧 mean={muHn_in.mean():.5f}  外侧 mean={muHn_out.mean():.5f}"
              f"  (应相等)")
        print(f"    H_t : 内侧 mean={Ht_in.mean():.5f}  外侧 mean={Ht_out.mean():.5f}"
              f"  (应相等)")

    # ---- 5. ψ 远场线性度 (绝对坐标: ψ = -H0·z + C·H0·R³·(z-zc)/r³) ----
    far = (r > 3*R) & (r < 40)
    if far.sum():
        ana_psi = -H0*zz[far] + C*H0*R**3*(zz[far]-zc)/r[far]**3
        ref = np.abs(ana_psi).max()
        err_psi = np.abs(psi[far] - ana_psi).max() / ref
        print(f"\n[5] ψ 远场 (r>3R, {far.sum()}格): 与解析最大偏差 = {err_psi*100:.2f}% "
              f"(相对远场 ψ 量级 {ref:.4f})")

if __name__ == '__main__':
    main()
