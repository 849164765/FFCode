#!/usr/bin/env python3
# plot_peak_vs_H0.py — peak height & valley depth vs H0 (paper Fig. 23 style).
#
# Usage:
#   python3 plot_peak_vs_H0.py [case_dir_prefix] [h0_list]
#   e.g. python3 plot_peak_vs_H0.py case_H0_ "0 2 4 6 8 10 12 14 15"
#
# For each case dir (case_H0_<h0>/), finds the LAST output step, reads the PHI
# field from the .vti files, extracts the interface y_if(x) at phi=0.5, and
# computes:
#   peak height  = max(y_if) - mean(y_if)
#   valley depth = mean(y_if) - min(y_if)
# (both in mm; dx = Lx_phys/Ni = 21mm/512 = 41.02um)
# Writes peak_vs_H0.csv and peak_vs_H0.png.
import base64, re, glob, sys, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

NI = 512
NJ = 171
LX_MM = 21.0          # physical domain length (mm)
DX_MM = LX_MM / NI    # mm per cell

def decode_vti_payload(s):
    s = re.sub(rb'\s', b'', s)
    head = s[:8]                      # size-prefix chunk (padded)
    rest = s[8:]
    rest += b'=' * (-len(rest) % 4)
    return base64.b64decode(head) + base64.b64decode(rest)

def read_phi_vti(fn):
    data = open(fn, 'rb').read()
    m = re.search(rb'Origin="([-\d.eE]+) ([-\d.eE]+)', data)
    ox, oy = float(m.group(1)), float(m.group(2))
    m = re.search(rb'Extent="(\d+) (\d+) (\d+) (\d+)', data)
    x0, x1, y0, y1 = map(int, m.groups())
    pat = rb'<DataArray type="Float64" Name="(\w+)" format="binary" encoding="base64" NumberOfComponents="\d+">(.*?)</DataArray>'
    arrs = {}
    for mm in re.finditer(pat, data, re.S):
        if mm.group(1).decode() != 'PHI':
            continue
        raw = decode_vti_payload(mm.group(2))
        n = int.from_bytes(raw[:4], 'little')
        arrs['PHI'] = np.frombuffer(raw[4:4 + n], dtype='<f8')
    nx = x1 - x0 + 1
    ny = y1 - y0 + 1
    return ox, oy, arrs['PHI'].reshape(ny, nx)

def interface_profile(stepdir):
    """Assemble the full PHI field of one output step and return y_if(x)."""
    phi = np.full((NJ, NI), np.nan)
    for f in sorted(glob.glob(os.path.join(stepdir, '*.vti'))):
        ox, oy, p = read_phi_vti(f)
        ny, nx = p.shape
        for j in range(ny):
            y = oy + j
            if not (0.5 <= y <= NJ - 0.5):
                continue
            for i in range(nx):
                x = ox + i
                if 0.5 <= x <= NI - 0.5:
                    phi[int(y), int(x)] = p[j, i]
    y_if = np.full(NI, np.nan)
    for x in range(NI):
        col = phi[:, x]
        for j in range(NJ - 1):
            a, b = col[j], col[j + 1]
            if not np.isnan(a) and not np.isnan(b) and (a - 0.5) * (b - 0.5) <= 0 and a != b:
                y_if[x] = j + (0.5 - a) / (b - a)
                break
    return y_if

def last_step_dir(case_dir):
    """Return the vtidata directory of the last output step in case_dir."""
    files = sorted(glob.glob(os.path.join(case_dir, 'vtkoutput', 'vtidata',
                                          'rosensweig2d_T*_B0.vti')))
    if not files:
        return None
    last = files[-1]
    step = int(re.search(r'_T(\d+)_', os.path.basename(last)).group(1))
    return os.path.join(case_dir, 'vtkoutput', 'vtidata'), step

def peak_valley(case_dir):
    vtidata, step = last_step_dir(case_dir)
    if vtidata is None:
        print(f'  [skip] no vtk output in {case_dir}')
        return None, None, None
    y_if = interface_profile(vtidata)
    m = np.nanmean(y_if)
    ph = np.nanmax(y_if) - m          # peak height (cells)
    vd = m - np.nanmin(y_if)          # valley depth (cells)
    return step, ph * DX_MM, vd * DX_MM

def main():
    prefix = sys.argv[1] if len(sys.argv) > 1 else 'case_H0_'
    h0_list = (sys.argv[2].split() if len(sys.argv) > 2
               else ['0', '2', '4', '6', '8', '10', '12', '14', '15'])
    rows = []
    for h0 in h0_list:
        case = f'{prefix}{h0}'
        print(f'H0={h0} kA/m: {case}')
        step, ph, vd = peak_valley(case)
        if step is None:
            continue
        print(f'  step {step}: peak={ph:.4f} mm  valley={vd:.4f} mm')
        rows.append((float(h0), step, ph, vd))

    if not rows:
        print('no data found')
        return
    rows = np.array(rows)
    with open('peak_vs_H0.csv', 'w') as fo:
        fo.write('H0_kAm,step,peak_height_mm,valley_depth_mm\n')
        for r in rows:
            fo.write(f'{r[0]:g},{int(r[1])},{r[2]:.6f},{r[3]:.6f}\n')

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(rows[:, 0], rows[:, 2], 'o-', color='tab:red', label='peak height')
    ax.plot(rows[:, 0], rows[:, 3], 's-', color='tab:blue', label='valley depth')
    ax.axvline(4.7, color='k', ls='--', lw=0.8)
    ax.text(4.75, ax.get_ylim()[1] * 0.9 if False else rows[:, 2].max() * 0.95,
            '$H_c$=4.7', fontsize=9)
    ax.set_xlabel(r'$H_0$ (kA/m)')
    ax.set_ylabel('height (mm)')
    ax.set_title('Rosensweig instability: peak height & valley depth vs $H_0$')
    ax.grid(alpha=0.3)
    ax.legend()
    plt.tight_layout()
    plt.savefig('peak_vs_H0.png', dpi=150)
    print('saved peak_vs_H0.csv / peak_vs_H0.png')

if __name__ == '__main__':
    main()
