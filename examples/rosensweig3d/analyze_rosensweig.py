#!/usr/bin/env python3
# analyze_rosensweig.py — Rosensweig instability 3D peak analysis
#
# Reads the rosensweig3d VTK output, extracts the phi=0.5 interface height
# field z_if(x,y), and reports peak statistics (peak count, peak height,
# valley depth, mean spacing / wavelength) at each output step.
#
# Usage:
#   python3 analyze_rosensweig.py [vtkdir] [--mesh 128,128,192] [--lambda 25]
#                                 [--iface 64] [--plot]
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
        # VTK stores point data x-fastest; reshape to (z,y,x[,c]) so the first
        # axis is the vertical (interface scan) direction. (For a cubic block
        # reshape((nx,ny,nz)) accidentally gives the same result, which is why
        # this only surfaced on non-cubic blocks.)
        if ncomp == 1:
            fields[nm] = arr.reshape((shp[2], shp[1], shp[0]))
        else:
            fields[nm] = arr.reshape((shp[2], shp[1], shp[0], ncomp))
    return fields, origin

def assemble(vtkdir, step, Ni, Nj, Nk):
    g = {}
    layout = []
    for p in sorted(glob.glob(f'{vtkdir}/rosensweig3d_T{step}_B*.vti')):
        bid = int(re.search(r'_B(\d+)\.vti', p).group(1))
        f, origin = read_vti(p)
        gx0, gy0, gz0 = (int(o + 0.5) for o in origin)
        shp = f['PHI'].shape  # (nz, ny, nx) after the (z,y,x) reshape
        for nm, arr in f.items():
            if nm not in g:
                g[nm] = np.zeros((Nk, Nj, Ni) if arr.ndim == 3 else (Nk, Nj, Ni, arr.shape[-1]))
            g[nm][gz0:gz0+shp[0]-2, gy0:gy0+shp[1]-2, gx0:gx0+shp[2]-2] = arr[1:-1, 1:-1, 1:-1]
        layout.append((bid, gx0, gy0, gz0, shp[0], shp[1], shp[2]))
    # sanity: every global cell covered?
    covered = (g['PHI'] != 0).any(axis=0) if 'PHI' in g else None
    return g, layout

def interface_height(phi, thresh=0.5):
    """Return z_if(x,y): first phi>=thresh cell scanning from the top (crest)."""
    Nk, Nj, Ni = phi.shape
    # mask: phi >= thresh (ferrofluid)
    above = phi >= thresh
    # for each column, highest k with above
    zidx = np.argmax(above, axis=0)  # first True scanning upward (lowest ferrofluid cell)
    has = above.any(axis=0)
    # refine with linear interp at the crossing just above the highest ferrofluid cell
    zc = np.full((Nj, Ni), np.nan)
    for j in range(Nj):
        for i in range(Ni):
            if not has[j, i]:
                continue
            k = zidx[j, i]
            # scan up from k until phi < thresh
            kk = k
            while kk < Nk-1 and phi[kk, j, i] >= thresh:
                kk += 1
            if kk == Nk-1:
                zc[j, i] = Nk-1
            else:
                p0, p1 = phi[kk-1, j, i], phi[kk, j, i]
                w = (thresh - p0) / (p1 - p0) if p1 != p0 else 0.5
                zc[j, i] = (kk-1) + w
    return zc

def count_peaks(z, height_thresh, Ni, Nj):
    """Count peak *clusters* (connected components above the height threshold)
    and their mean nearest-neighbour spacing. Connected-component labelling is
    done with a two-pass union-find-free scan (no scipy needed)."""
    zc = np.nan_to_num(z, nan=np.nanmin(z))
    mask = zc > height_thresh
    # flood-fill labels (8-connectivity)
    labels = np.zeros_like(mask, int)
    label = 0
    stack = []
    for j in range(Nj):
        for i in range(Ni):
            if mask[j, i] and labels[j, i] == 0:
                label += 1
                stack.append((j, i))
                labels[j, i] = label
                while stack:
                    y, x = stack.pop()
                    for dy in (-1, 0, 1):
                        yy = (y + dy) % Nj
                        for dx in (-1, 0, 1):
                            if dy == 0 and dx == 0:
                                continue
                            xx = (x + dx) % Ni
                            if mask[yy, xx] and labels[yy, xx] == 0:
                                labels[yy, xx] = label
                                stack.append((yy, xx))
    if label == 0:
        return 0, float('nan')
    # centroid of each cluster
    cents = []
    for lab in range(1, label + 1):
        ys, xs = np.nonzero(labels == lab)
        cents.append((np.mean(ys), np.mean(xs)))
    # mean nearest-neighbour spacing with periodic wrap
    sp = []
    for p in cents:
        dmin = 1e9
        for q in cents:
            if q is p:
                continue
            dx = abs(p[1]-q[1]); dy = abs(p[0]-q[0])
            dx = min(dx, Ni-dx); dy = min(dy, Nj-dy)
            d = np.hypot(dx, dy)
            if d < dmin:
                dmin = d
        if dmin < 1e8:
            sp.append(dmin)
    spacing = float(np.mean(sp)) if sp else float('nan')
    return label, spacing

def dominant_wavelength(z, Ni, Nj):
    """Dominant wavelength of the interface pattern via 2D FFT (periodic box)."""
    zc = np.nan_to_num(z, nan=np.nanmean(z))
    zc = zc - zc.mean()
    F = np.fft.fft2(zc)
    kys, kxs = np.mgrid[0:Nj, 0:Ni]
    kx = kxs - (kxs > Ni//2) * Ni
    ky = kys - (kys > Nj//2) * Nj
    kr = np.sqrt(kx**2 + ky**2)
    kmax = min(Ni, Nj) // 2
    kvals = np.arange(1, kmax)
    p = np.array([np.abs(F[(kr > k-0.5) & (kr < k+0.5)]).sum() for k in kvals])
    if p.size == 0 or p.max() == 0:
        return float('nan')
    km = kvals[np.argmax(p)]
    return Ni / km

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('vtkdir', nargs='?', default='./vtkoutput/vtidata')
    ap.add_argument('--mesh', default='128,128,192')
    ap.add_argument('--block', type=int, default=32)
    ap.add_argument('--lambda', dest='lam', type=float, default=25.0)
    ap.add_argument('--iface', type=float, default=64.0)
    ap.add_argument('--steps', default=None)
    ap.add_argument('--plot', action='store_true')
    args = ap.parse_args()
    Ni, Nj, Nk = (int(v) for v in args.mesh.split(','))
    if args.steps:
        steps = [int(s) for s in args.steps.split(',')]
    else:
        steps = sorted(int(re.search(r'_T(\d+)_', p).group(1))
                       for p in glob.glob(f'{args.vtkdir}/rosensweig3d_T*_B0.vti'))
        steps = sorted(set(steps))
    if not steps:
        print("no timestep files found"); return

    print(f"=== Rosensweig peak analysis: mesh {Ni}x{Nj}x{Nk}, lambda_c={args.lam}, "
          f"iface={args.iface} ===")
    print(f"{'step':>6} {'zmax':>7} {'zmin':>7} {'amp':>7} {'zmean':>7} "
          f"{'nclust':>7} {'spacing':>8} {'waveleng':>8} {'|u|max':>8}")
    for s in steps:
        g, _ = assemble(args.vtkdir, s, Ni, Nj, Nk)
        phi = g['PHI']
        zc = interface_height(phi)
        # peak height above the mean interface
        zmean = np.nanmean(zc)
        zmax = np.nanmax(zc)
        zmin = np.nanmin(zc)
        amp = zmax - zmin
        # count peak clusters: connected regions above zmean + 0.25*amp
        thresh = zmean + 0.25 * amp
        npeak, spacing = count_peaks(zc, thresh, Ni, Nj)
        wl = dominant_wavelength(zc, Ni, Nj)
        # |u|max from Velocity (components), if present
        umax = float('nan')
        if 'Velocity' in g:
            umax = np.sqrt((g['Velocity']**2).sum(-1)).max()
        print(f"{s:6d} {zmax:7.2f} {zmin:7.2f} {amp:7.2f} {zmean:7.2f} "
              f"{npeak:7d} {spacing:8.1f} {wl:8.1f} {umax:8.4f}")
        if args.plot:
            try:
                import matplotlib
                matplotlib.use('Agg')
                import matplotlib.pyplot as plt
                from mpl_toolkits.mplot3d import Axes3D  # noqa
                x = np.arange(Ni); y = np.arange(Nj)
                X, Y = np.meshgrid(x, y)
                fig = plt.figure(figsize=(10, 4))
                ax = fig.add_subplot(1, 2, 1, projection='3d')
                ax.plot_surface(X, Y, zc, cmap='viridis', linewidth=0)
                ax.set_title(f'interface height step {s}')
                ax2 = fig.add_subplot(1, 2, 2)
                im = ax2.imshow(zc, origin='lower', extent=[0, Ni, 0, Nj], cmap='viridis')
                ax2.set_title(f'peak map (npeaks={npeak})')
                plt.colorbar(im, ax=ax2)
                plt.savefig(f'rosensweig_peaks_{s}.png', dpi=110)
                print(f"  -> rosensweig_peaks_{s}.png")
            except ImportError:
                print("  (matplotlib not available, skip plot)")

if __name__ == '__main__':
    main()
