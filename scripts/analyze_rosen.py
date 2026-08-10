#!/usr/bin/env python3
"""analyze_rosen.py — Rosensweig 3D interface analysis from rosenMag3d vti output

Usage: python3 analyze_rosen.py <workdir> [z0] [--png dir]

Reads PHI/HMAG from block-structured 3D VTI files of one output step,
assembles the global domain, computes the interface height map h(x,y)
(phi=0.5 crossing) and reports:
  - h_max / h_min / peak amplitude (hmax-hmin)/2
  - dominant wavelength (2D FFT of h)
  - bulk |H| uniformity check (mean |H| in a solvent bulk slab)
"""

import sys, os, glob, re, base64, argparse
import numpy as np

FIELD_TYPES = {
    "Float32": np.float32,
    "Float64": np.float64,
}


def decode_payload(payload):
    """Writer emits two concatenated base64 blobs: 4-byte size prefix + data."""
    p = re.sub(r"\s+", "", payload)
    out = b""
    while p:
        i = p.find("==")
        if i < 0:
            pad = "=" * ((-len(p)) % 4)
            out += base64.b64decode(p + pad)
            break
        out += base64.b64decode(p[: i + 2])
        p = p[i + 2:]
    return out


def parse_vti(path):
    with open(path, "rb") as f:
        data = f.read()
    text = data.decode("ascii", errors="ignore")
    origin = [float(v) for v in re.search(r'Origin="([^"]*)"', text).group(1).split()]
    spacing = [float(v) for v in re.search(r'Spacing="([^"]*)"', text).group(1).split()]
    extent = [int(v) for v in re.search(r'Extent="([^"]*)"', text).group(1).split()]
    arrays = {}
    for m in re.finditer(r'<DataArray([^>]*)>(.*?)</DataArray>', text, re.S):
        attrs, payload = m.group(1), m.group(2)
        name = re.search(r'Name="([^"]*)"', attrs).group(1)
        dtype = FIELD_TYPES[re.search(r'type="([^"]*)"', attrs).group(1)]
        comps = int(re.search(r'NumberOfComponents="([^"]*)"', attrs).group(1))
        raw = decode_payload(payload)
        nbytes = np.frombuffer(raw[:4], dtype="<i4")[0]
        vals = np.frombuffer(raw[4:4 + nbytes], dtype="<f8").astype(np.float64)
        if comps > 1:
            vals = vals.reshape(-1, comps)
        arrays[name] = vals
    return {"origin": origin, "spacing": spacing, "extent": extent, "arrays": arrays}


def load_step(step, vtkdir):
    files = sorted(glob.glob(f"{vtkdir}/vtidata/rosenMag3d_T{step}_B*.vti"))
    if not files:
        return None
    blocks = [parse_vti(p) for p in files]
    mn = np.array([1e30] * 3)
    mx = np.array([-1e30] * 3)
    bmin, bmax = [], []
    for b in blocks:
        o = np.array(b["origin"]) + 0.5 * np.array(b["spacing"])
        ex = b["extent"]
        counts = np.array([ex[1] - ex[0] + 1, ex[3] - ex[2] + 1, ex[5] - ex[4] + 1])
        # local 0 = halo row; physical cells are local 1..counts-2
        bmin.append(o.astype(int))
        bmax.append((o + counts - 2).astype(int))
        mn = np.minimum(mn, bmin[-1]); mx = np.maximum(mx, bmax[-1])
    shape = tuple((mx - mn).astype(int))
    mn = mn.astype(int)
    phi = np.zeros(shape); hmag = np.zeros(shape)
    for b, b0, b1 in zip(blocks, bmin, bmax):
        a = b["arrays"]
        ex = b["extent"]
        # block memory layout id = k*nx*ny + j*nx + i  (z slowest, x fastest)
        ph = a["PHI"].reshape(ex[5] + 1, ex[3] + 1, ex[1] + 1).transpose(2, 1, 0)
        hm = a["HMAG"].reshape(ex[5] + 1, ex[3] + 1, ex[1] + 1).transpose(2, 1, 0)
        # drop the 1-row halo on each side (physical = local 1..counts-2)
        ph = ph[1:-1, 1:-1, 1:-1]
        hm = hm[1:-1, 1:-1, 1:-1]
        sl = tuple(slice(x0 - mn[d], x1 - mn[d]) for d, (x0, x1) in enumerate(zip(b0, b1)))
        phi[sl] = ph
        hmag[sl] = hm
    return {"origin": mn, "phi": phi, "hmag": hmag, "blocks": len(files)}


def interface_height(phi, z0_guess):
    nz = phi.shape[2]
    above = phi >= 0.5
    cross = np.argmax(~above, axis=2)
    h = np.zeros(phi.shape[:2])
    for iy in range(phi.shape[1]):
        for ix in range(phi.shape[0]):
            k = int(cross[iy, ix])
            if k == 0 or k >= nz - 1:
                h[iy, ix] = z0_guess
                continue
            p0, p1 = phi[iy, ix, k - 1], phi[iy, ix, k]
            if p1 - p0 == 0:
                h[iy, ix] = k
            else:
                h[iy, ix] = (k - 1) + (0.5 - p0) / (p1 - p0)
    return h


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("workdir", help="simulation working dir containing vtkoutput/")
    ap.add_argument("z0", type=float, nargs="?", default=42.67, help="initial interface z")
    ap.add_argument("--png", metavar="DIR", default=None, help="write h(x,y) PNG snapshots")
    args = ap.parse_args()

    vtkdir = os.path.join(args.workdir, "vtkoutput")
    steps = sorted(set(int(re.search(r"_T(\d+)_", p).group(1))
                       for p in glob.glob(f"{vtkdir}/vtidata/rosenMag3d_T*_B0.vti")))
    if not steps:
        print("no vti files found"); return
    print(f"steps found: {steps}")

    for s in steps:
        d = load_step(s, vtkdir)
        if d is None:
            continue
        phi = d["phi"]
        ix = slice(1, -1) if phi.shape[0] > 2 else slice(None)
        iy = slice(1, -1) if phi.shape[1] > 2 else slice(None)
        h = interface_height(phi[ix, iy], args.z0)
        amp = 0.5 * (h.max() - h.min())
        hc = h - h.mean()
        F = np.abs(np.fft.fftshift(np.fft.fft2(hc)))
        F[F.shape[0] // 2, F.shape[1] // 2] = 0
        ky, kx = np.unravel_index(np.argmax(F), F.shape)
        ny, nx = F.shape
        lam_x = nx / abs(kx - nx // 2) if kx != nx // 2 else float("nan")
        lam_y = ny / abs(ky - ny // 2) if ky != ny // 2 else float("nan")
        bulk_hmag = d["hmag"][:, :, 5:15].mean()
        print(f"t={s:5d}  hmax={h.max():8.3f} hmin={h.min():8.3f} amp={amp:7.3f} "
              f"lambda_x={lam_x:6.2f} lambda_y={lam_y:6.2f} bulk|H|={bulk_hmag:.4f}")
        if args.png:
            try:
                import matplotlib
                matplotlib.use("Agg")
                import matplotlib.pyplot as plt
                fig, ax = plt.subplots(figsize=(6, 6))
                im = ax.imshow(h, origin="lower", cmap="viridis")
                plt.colorbar(im, ax=ax, label="interface z")
                ax.set_title(f"t={s}  amp={amp:.2f}")
                ax.set_xlabel("x"); ax.set_ylabel("y")
                fig.savefig(os.path.join(args.png, f"h_t{s:05d}.png"), dpi=100)
                plt.close(fig)
            except Exception as e:
                print(f"  (png skip: {e})")


if __name__ == "__main__":
    main()
