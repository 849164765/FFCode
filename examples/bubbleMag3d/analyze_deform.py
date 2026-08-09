#!/usr/bin/env python3
# analyze_deform.py — 气泡变形/崩溃分析 (bubbleMag3d 本地复现用)
#
# 用法:
#   python3 analyze_deform.py [vtkdir] [--mesh 64,64,128] [--block 32] [--steps 0,50,...]
#                            [--radius 16] [--h0 0.0289]
#
# 输出每个时间步:
#   - 气泡 x/y/z 跨度、赤道半径 vs 极向半径 (左右被压/上下拉伸检测)
#   - 体积、质心、球形度
#   - |u|max/Ma、|F|max、|H|max、NaN 统计
#   - 赤道界面厚度
#   - 周期性接缝偏差 (跨 rank 同步正确性)
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
    gfull = {}
    layout = []
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

def seam_check(gfull, layout, Ni, Nj, Nk, block):
    bmap = {}
    for bid, gx0, gy0, gz0, nx, ny, nz in layout:
        bmap[(gx0//block, gy0//block, gz0//block)] = bid
    NX, NY, NZ = Ni//block, Nj//block, Nk//block
    res = {}
    for nm in ['PSI', 'HX', 'HY', 'HZ', 'HMAG']:
        errs = []
        for zz in range(NZ):
            for yy in range(NY):
                if (0, yy, zz) in bmap and (NX-1, yy, zz) in bmap:
                    bL, bR = bmap[(0, yy, zz)], bmap[(NX-1, yy, zz)]
                    aL, aR = gfull[nm][bL], gfull[nm][bR]
                    errs.append(np.abs(aL[1:-1,1:-1,0]-aR[1:-1,1:-1,block]).max())
                    errs.append(np.abs(aR[1:-1,1:-1,block+1]-aL[1:-1,1:-1,1]).max())
            for xx in range(NX):
                if (xx, 0, zz) in bmap and (xx, NY-1, zz) in bmap:
                    bF, bB = bmap[(xx, 0, zz)], bmap[(xx, NY-1, zz)]
                    aF, aB = gfull[nm][bF], gfull[nm][bB]
                    errs.append(np.abs(aF[1:-1,0,1:-1]-aB[1:-1,block,1:-1]).max())
                    errs.append(np.abs(aB[1:-1,block+1,1:-1]-aF[1:-1,1,1:-1]).max())
        res[nm] = max(errs) if errs else float('nan')
    return res

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('vtkdir', nargs='?', default='./vtkoutput/vtidata')
    ap.add_argument('--mesh', default='64,64,128')
    ap.add_argument('--block', type=int, default=32)
    ap.add_argument('--steps', default=None)
    ap.add_argument('--radius', type=float, default=16.0)
    ap.add_argument('--h0', type=float, default=0.0)
    args = ap.parse_args()
    Ni, Nj, Nk = (int(v) for v in args.mesh.split(','))
    CS = np.sqrt(1/3)
    if args.steps:
        steps = [int(s) for s in args.steps.split(',')]
    else:
        steps = sorted(int(re.search(r'_T(\d+)_', p).group(1))
                       for p in glob.glob(f'{args.vtkdir}/bubbleMag3d_T*_B0.vti'))
        steps = sorted(set(steps))
    if not steps:
        print("no timestep files found"); return

    print(f"=== 变形/崩溃分析: mesh {Ni}x{Nj}x{Nk}, R={args.radius}, 步: {steps} ===")
    print(f"{'step':>5} {'vol':>8} {'span_x':>7} {'span_y':>7} {'span_z':>7} "
          f"{'R_eq':>6} {'R_pol':>6} {'R_eq/R_pol':>9} {'COM_z':>6} "
          f"{'|u|max':>7} {'Ma':>6} {'|F|max':>9} {'|H|max':>7} {'NaN':>5}")
    prev = None
    for s in steps:
        g, gfull, layout = assemble(args.vtkdir, s, Ni, Nj, Nk)
        phi = g['PHI']
        nanu = int(np.isnan(g['Velocity']).sum() + np.isnan(g['Force']).sum())
        idx = np.argwhere(phi < 0.5)
        if len(idx) == 0:
            print(f"{s:5d} 气泡消失 (完全崩溃)")
            continue
        com = idx.mean(axis=0)
        spans = idx.max(axis=0) - idx.min(axis=0) + 1
        # 赤道半径: 过质心水平面 (x-y) 内 phi=0.5 平均半径
        zc = int(round(com[0]))
        eq = phi[zc]
        yy, xx = np.mgrid[0:Nj, 0:Ni]
        r = np.sqrt((xx-com[2])**2 + (yy-com[1])**2)
        rbins = np.arange(0, r.max(), 0.5)
        prof = np.array([eq[(r>=a)&(r<a+0.5)].mean() for a in rbins[:-1]])
        ok = ~np.isnan(prof)
        if ok.sum() > 2:
            Req = rbins[:-1][ok][np.argmin(np.abs(prof[ok]-0.5))]
        else:
            Req = float('nan')
        # 极向半径: 过质心垂直轴 phi=0.5
        col = phi[:, int(round(com[1])), int(round(com[2]))]
        zz = np.arange(Nk)
        cross = np.where(np.diff((col > 0.5).astype(int)))[0]
        Rpol = float('nan')
        if len(cross) >= 2:
            Rpol = (cross[-1] - cross[0]) / 2.0
        umag = np.sqrt((g['Velocity']**2).sum(-1))
        Fmag = np.sqrt((g['Force']**2).sum(-1))
        hm = g['HMAG']
        ratio = Req/Rpol if Rpol and not np.isnan(Req) else float('nan')
        h0ref = args.h0 if args.h0 > 0 else 1.0
        print(f"{s:5d} {len(idx):8d} {spans[2]:7d} {spans[1]:7d} {spans[0]:7d} "
              f"{Req:6.2f} {Rpol:6.2f} {ratio:9.3f} {com[0]:6.1f} "
              f"{umag.max():7.4f} {umag.max()/CS:6.3f} {Fmag.max():9.2e} "
              f"{hm.max()/h0ref:7.2f}xH0 {nanu:5d}")
        # 崩溃签名
        if nanu > 0:
            nz_, ny_, nx_ = np.where(np.isnan(g['Velocity']).any(-1))
            if len(nz_):
                print(f"    !! 速度 NaN 区域: z={nz_.min()}..{nz_.max()} y={ny_.min()}..{ny_.max()} x={nx_.min()}..{nx_.max()}")
            seam = seam_check(gfull, layout, Ni, Nj, Nk, args.block)
            worst = max(seam.values())
            print(f"    接缝偏差: PSI={seam['PSI']:.2e} HX={seam['HX']:.2e} HMAG={seam['HMAG']:.2e}")
        prev = (phi, com)

if __name__ == '__main__':
    main()
