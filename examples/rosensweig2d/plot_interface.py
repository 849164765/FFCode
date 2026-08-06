#!/usr/bin/env python3
# plot_interface.py — view the Rosensweig interface from the VTK output.
#
# Usage:
#   python3 plot_interface.py <case_dir> [step]
#   e.g. python3 plot_interface.py case_H0_10 9000
#        (step omitted -> last output step)
#
# Coordinate convention: ORIGIN AT THE BOTTOM-LEFT CORNER of the domain (0,0).
#   x in [0, Lx] mm  (Lx = 21 mm),  y in [0, Ly] mm  (Ly = 7 mm)
# (NOT at the left-boundary midpoint; the ghost-column offset of the vti
#  writer is stripped: only the physical cells are mapped.)
#
# Outputs:
#   interface_<step>.csv : x_mm, y_iface_mm  (phi=0.5 crossing per column)
#   interface_<step>.png : PHI field image + interface profile (mm axes)
import base64, re, glob, sys, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

NI = 512
NJ = 171
LX_MM = 21.0          # physical domain length (mm)
LY_MM = 7.0           # physical domain height (mm)
DX_MM = LX_MM / NI    # mm per cell (x)
DY_MM = LY_MM / NJ    # mm per cell (y)

def decode_vti_payload(s):
    s = re.sub(rb'\s', b'', s)
    head = s[:8]                      # size-prefix chunk (padded)
    rest = s[8:]
    rest += b'=' * (-len(rest) % 4)
    return base64.b64decode(head) + base64.b64decode(rest)

def read_vti_phi(fn):
    data = open(fn, 'rb').read()
    m = re.search(rb'Origin="([-\d.eE]+) ([-\d.eE]+)', data)
    ox, oy = float(m.group(1)), float(m.group(2))
    m = re.search(rb'Extent="(\d+) (\d+) (\d+) (\d+)', data)
    x0, x1, y0, y1 = map(int, m.groups())
    pat = rb'<DataArray type="Float64" Name="(\w+)" format="binary" encoding="base64" NumberOfComponents="\d+">(.*?)</DataArray>'
    for mm in re.finditer(pat, data, re.S):
        if mm.group(1).decode() != 'PHI':
            continue
        raw = decode_vti_payload(mm.group(2))
        n = int.from_bytes(raw[:4], 'little')
        phi = np.frombuffer(raw[4:4 + n], dtype='<f8')
        return ox, oy, phi.reshape(y1 - y0 + 1, x1 - x0 + 1)
    return ox, oy, None

def assemble_phi(vtidata):
    """Assemble the full PHI field (physical cells only)."""
    phi = np.full((NJ, NI), np.nan)
    for f in sorted(glob.glob(os.path.join(vtidata, '*.vti'))):
        ox, oy, p = read_vti_phi(f)
        if p is None:
            continue
        ny, nx = p.shape
        for j in range(ny):
            y = oy + j
            if not (0.5 <= y <= NJ - 0.5):
                continue
            for i in range(nx):
                x = ox + i
                if 0.5 <= x <= NI - 0.5:
                    phi[int(y), int(x)] = p[j, i]
    return phi

def interface_profile(phi):
    """phi=0.5 crossing per x column, in mm from the bottom-left corner."""
    y_if = np.full(NI, np.nan)
    for x in range(NI):
        col = phi[:, x]
        for j in range(NJ - 1):
            a, b = col[j], col[j + 1]
            if not np.isnan(a) and not np.isnan(b) and (a - 0.5) * (b - 0.5) <= 0 and a != b:
                # physical y (mm) from the bottom edge: cell row j spans
                # [j*dy,(j+1)*dy], crossing at j+0.5+frac cells
                y_if[x] = (j + 0.5 + (0.5 - a) / (b - a)) * DY_MM
                break
    return y_if

def last_step(vtidata):
    files = sorted(glob.glob(os.path.join(vtidata, 'rosensweig2d_T*_B0.vti')))
    if not files:
        return None
    return int(re.search(r'_T(\d+)_', os.path.basename(files[-1])).group(1))

def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    case_dir = sys.argv[1]
    vtidata = os.path.join(case_dir, 'vtkoutput', 'vtidata')
    step = int(sys.argv[2]) if len(sys.argv) > 2 else last_step(vtidata)
    if step is None:
        print(f'[error] no VTK output found in {case_dir}')
        sys.exit(1)

    phi = assemble_phi(vtidata)
    y_if = interface_profile(phi)
    x_mm = (np.arange(NI) + 0.5) * DX_MM

    with open(f'interface_{step}.csv', 'w') as fo:
        fo.write('# x_mm  y_iface_mm  (origin: bottom-left corner of the domain)\n')
        for x in range(NI):
            fo.write(f'{x_mm[x]:.6f} {y_if[x]:.6f}\n')

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 8), gridspec_kw={'height_ratios': [2, 1]})
    im = ax1.imshow(phi, origin='lower', aspect='auto', cmap='jet', vmin=0, vmax=1,
                    extent=[0, LX_MM, 0, LY_MM], interpolation='nearest')
    plt.colorbar(im, ax=ax1, fraction=0.03, label='phi')
    ax1.set_xlabel('x (mm)'); ax1.set_ylabel('y (mm)')
    ax1.set_title(f'PHI field, step {step} (origin: bottom-left corner)')

    ax2.plot(x_mm, y_if, 'k-', lw=1.2)
    ax2.axhline(57.0 * DY_MM, color='gray', ls='--', lw=0.8)
    ax2.set_xlabel('x (mm)'); ax2.set_ylabel('y (mm)')
    ax2.set_title('interface y(x), phi=0.5')
    ax2.set_xlim(0, LX_MM); ax2.set_ylim(0, LY_MM)
    plt.tight_layout()
    plt.savefig(f'interface_{step}.png', dpi=150)
    print(f'step {step}: interface y range {np.nanmin(y_if):.3f}..{np.nanmax(y_if):.3f} mm')
    print(f'peak height = {np.nanmax(y_if)-np.nanmean(y_if):.4f} mm, '
          f'valley depth = {np.nanmean(y_if)-np.nanmin(y_if):.4f} mm')
    print('saved interface_%d.csv / interface_%d.png' % (step, step))

if __name__ == '__main__':
    main()
