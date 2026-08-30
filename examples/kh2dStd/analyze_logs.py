#!/usr/bin/env python3
"""
analyze_logs.py — parse kh2dStd out{}.log files and answer two questions:

1) Time-scale calibration
   Old logs print  t*_plot = k * step*sqrt(g/L)  with k=3.
   For a visual anchor "step n1..n2 matches paper t* = T", the implied k is
       k = T / (n * sqrt(g/L)).
   New logs also print t*_eq, so k is simply t*_plot/t*_eq.

2) Magnetic-force magnitude
   Old logs print [DEBUG] Hmag / Fmag (Kelvin force only, local to rank 0).
   New logs print [KH2D_MAG] Hmag / F_kelvin / F_press / F_total (global MPI
   reductions), where F_press = -0.5*Hx^2*grad(mu).

Usage:
   python3 analyze_logs.py [dir] [--anchor T n1 n2] [--k SCALE]
Examples:
   python3 analyze_logs.py .                    # k=4.0 default
   python3 analyze_logs.py . --anchor 8 25600 25600
   python3 analyze_logs.py . --k 4.0
"""
import sys, os, re, glob, argparse
import numpy as np

OLD_DBG = re.compile(
    r'\[KH2D_DBG\] step=(\d+) t\*=([\d.eE+-]+)\s+hmin=([\d.eE+-]+) hmax=([\d.eE+-]+) amp=([\d.eE+-]+)\s+ncross=(\d+)')
NEW_DBG = re.compile(
    r'\[KH2D_DBG\] step=(\d+) t\*_eq=([\d.eE+-]+) t\*_plot=([\d.eE+-]+)\s+hmin=([\d.eE+-]+) hmax=([\d.eE+-]+) amp=([\d.eE+-]+)\s+ncross=(\d+)')
MAG = re.compile(
    r'\[KH2D_MAG\] step=(\d+) Hmag: min=([\d.eE+-]+) max=([\d.eE+-]+) mean=([\d.eE+-]+) range=([\d.eE+-]+) \(H0=([\d.eE+-]+)\)\s*')
FORCE = re.compile(
    r'\[KH2D_MAG\] step=(\d+) F_kelvin: max=([\d.eE+-]+) mean=([\d.eE+-]+) \| F_press: max=([\d.eE+-]+) mean=([\d.eE+-]+) \| F_total: max=([\d.eE+-]+) mean=([\d.eE+-]+)')
OLD_FORCE = re.compile(r'\[DEBUG step=(\d+)\] Fmag: max=([\d.eE+-]+) mean=([\d.eE+-]+)')
OLD_HMAG = re.compile(
    r'\[DEBUG step=(\d+)\] Hmag: min=([\d.eE+-]+) max=([\d.eE+-]+) mean=([\d.eE+-]+)\s+range=([\d.eE+-]+) \(H0=([\d.eE+-]+)\)')


def load(path):
    rows, mag = {}, {}
    with open(path) as f:
        for line in f:
            m = NEW_DBG.search(line)
            if m:
                s = int(m.group(1))
                rows[s] = dict(t_eq=float(m.group(2)), t_plot=float(m.group(3)),
                               hmin=float(m.group(4)), hmax=float(m.group(5)),
                               amp=float(m.group(6)), ncross=int(m.group(7)))
                continue
            m = OLD_DBG.search(line)
            if m:
                s = int(m.group(1))
                rows[s] = dict(t_plot=float(m.group(2)), hmin=float(m.group(3)),
                               hmax=float(m.group(4)), amp=float(m.group(5)),
                               ncross=int(m.group(6)))
                continue
            m = FORCE.search(line)
            if m:
                mag[int(m.group(1))] = dict(Fk_max=float(m.group(2)), Fk_mean=float(m.group(3)),
                                            Fp_max=float(m.group(4)), Fp_mean=float(m.group(5)),
                                            Ft_max=float(m.group(6)), Ft_mean=float(m.group(7)))
                continue
            m = MAG.search(line)
            if m:
                mag.setdefault(int(m.group(1)), {}).update(
                    dict(Hmin=float(m.group(2)), Hmax=float(m.group(3)), Hmean=float(m.group(4)),
                         Hrange=float(m.group(5)), H0=float(m.group(6))))
                continue
            m = OLD_FORCE.search(line)
            if m:
                mag.setdefault(int(m.group(1)), {}).update(
                    dict(Fk_max=float(m.group(2)), Fk_mean=float(m.group(3))))
            m = OLD_HMAG.search(line)
            if m:
                mag.setdefault(int(m.group(1)), {}).update(
                    dict(Hmin=float(m.group(2)), Hmax=float(m.group(3)), Hmean=float(m.group(4)),
                         Hrange=float(m.group(5)), H0=float(m.group(6))))
    return rows, mag


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('dir', nargs='?', default='.')
    ap.add_argument('--anchor', nargs=3, type=float, default=[8.0, 25600.0, 25600.0],
                    help='paper_tstar step_min step_max (default: 8 25600 25600)')
    ap.add_argument('--k', type=float, default=4.0)
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.dir, 'out*.log')))
    if not files:
        print('no out*.log found'); return
    data = {}
    for path in files:
        m = re.search(r'out(\d+)\.log', os.path.basename(path))
        if not m:
            continue
        bom = int(m.group(1))
        rows, mag = load(path)
        data[bom] = (rows, mag)

    g_guess = 1.5625e-6
    if args.anchor:
        Tp, n1, n2 = args.anchor
        for n in (n1, n2):
            k = Tp / (n * np.sqrt(g_guess / 256.0))
            print(f'k for paper t*={Tp:g} at step {n:.0f}: {k:.4f}')
    print(f'\nUsing k = {args.k} (t*_plot = k*step*sqrt(g/L)):')
    print(f'{"step":>6s} {"t*":>6s}', end='')
    for b in sorted(data):
        print(f' | {"B"+str(b):>6s} amp', end='')
    print()
    steps = sorted(data[0][0].keys())
    for s in steps[::2]:
        if s % 2000:
            pass
        t = s * args.k * np.sqrt(g_guess / 256.0)
        line = f'{s:6d} {t:6.3f}'
        for b in sorted(data):
            r = data[b][0].get(s)
            line += f' | {r["amp"] if r else float("nan"):10.3f}'
        print(line)

    print('\nAmplitude at selected paper t* (linear interp in step):')
    print(f'{"t*":>5s}', end='')
    for b in sorted(data):
        print(f' | {"B"+str(b):>6s} amp  ncross', end='')
    print()
    for tt in [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6]:
        qstep = tt / (args.k * np.sqrt(g_guess / 256.0))
        line = f'{tt:5.2f}'
        for b in sorted(data):
            rows = data[b][0]
            xs = np.array(sorted(rows))
            amp = np.interp(qstep, xs, [rows[int(x)]['amp'] for x in xs])
            ncr = np.interp(qstep, xs, [rows[int(x)]['ncross'] for x in xs])
            line += f' | {amp:8.2f} {ncr:6.0f}'
        print(line)

    # magnetic force summary
    if any(data[b][1] for b in data):
        print('\nMagnetic force diagnostics (latest available step per Bom):')
        for b in sorted(data):
            mag = data[b][1]
            if not mag:
                continue
            s = max(mag)
            m = mag[s]
            print(f'B{b:>5d} step={s}:', '  '.join(f'{k}={v:.3e}' for k, v in m.items()))
            if 'Fk_mean' in m and 'Fp_mean' in m and m['Fk_mean'] > 0:
                print(f'          F_press/F_kelvin: max={m.get("Fp_max",0)/m["Fk_max"] if m.get("Fk_max") else 0:.3f} '
                      f'mean={m["Fp_mean"]/m["Fk_mean"]:.3f}')


if __name__ == '__main__':
    main()
