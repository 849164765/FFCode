#!/usr/bin/env python3
"""analyze_magnetic.py — Debug magnetic fields from VTK output"""
import sys, os, glob, argparse
import numpy as np
from plot_mag_bubble import assemble_field

def analyze(vtm_path):
    vtm_dir = vtm_path if os.path.isdir(vtm_path) else os.path.dirname(vtm_path)
    if os.path.isdir(vtm_dir):
        vtms = sorted(glob.glob(os.path.join(vtm_dir, '*.vtm')))
        steps = []
        for v in vtms:
            bn = os.path.basename(v)
            try:
                s = int(bn.split('.')[0].split('_T')[-1] if '_T' in bn else bn.split('.')[0].replace('bubbleMag2d',''))
                steps.append((s, v))
            except: pass
        if steps:
            steps.sort(key=lambda x: x[0])
            vtm_path = steps[-1][1]
            print(f"Using step {steps[-1][0]}: {vtm_path}")

    fields_to_check = ['PHI', 'PSI', 'HX', 'HY', 'Force']
    available = {}
    for fname in fields_to_check:
        data, extent, dx = assemble_field(vtm_path, fname)
        if data is not None:
            available[fname] = data
            print(f"  Loaded {fname}: shape={data.shape}, range=[{np.nanmin(data):.6f}, {np.nanmax(data):.6f}]")

    if 'HY' in available:
        hy = available['HY']
        h0_expected = 0.013  # from the verification runs
        print(f"\n  HY stats: mean={np.nanmean(hy):.6f}, median={np.nanmedian(hy):.6f}")
        print(f"  Expected H0 = {h0_expected}")
        deviation = np.abs(np.nanmean(hy) - h0_expected) / h0_expected
        print(f"  Deviation from H0: {deviation*100:.1f}%")

    if 'PSI' in available:
        psi = available['PSI']
        y_min = 0
        dy = 1.0  # assuming dx=1
        gy = psi.shape[0]
        # Check linearity: psi should be -H0*y
        y_vals = np.arange(gy) * dy
        # Take a column at the center
        cx = psi.shape[1] // 2
        psi_col = psi[:, cx]
        valid = ~np.isnan(psi_col)
        if valid.sum() > 0:
            y_fit = y_vals[valid]
            psi_fit = psi_col[valid]
            # Linear fit: psi = a*y + b
            coeffs = np.polyfit(y_fit, psi_fit, 1)
            print(f"\n  PSI linear fit (center column): psi = {coeffs[0]:.6f}*y + {coeffs[1]:.6f}")
            print(f"  Expected: psi = {-h0_expected}*y + 0")
            h_actual = -coeffs[0]
            print(f"  Effective H (from slope): {h_actual:.6f} (expected {h0_expected})")

    if 'Force' in available:
        force = available['Force']
        fx = force[:, :, 0] if force.ndim == 3 else force
        fy = force[:, :, 1] if force.ndim == 3 else None
        print(f"\n  Force X: range=[{np.nanmin(fx):.2e}, {np.nanmax(fx):.2e}], mean={np.nanmean(np.abs(fx)):.2e}")
        if fy is not None:
            print(f"  Force Y: range=[{np.nanmin(fy):.2e}, {np.nanmax(fy):.2e}], mean={np.nanmean(np.abs(fy)):.2e}")

    if 'PHI' in available:
        phi = available['PHI']
        gy, gx = phi.shape
        y_vals = np.arange(gy)
        # Center of mass in Y
        bubble_mask = phi < 0.5
        if bubble_mask.sum() > 0:
            cy = np.sum(y_vals[:, None] * bubble_mask) / np.sum(bubble_mask)
            cx = np.sum(np.arange(gx)[None, :] * bubble_mask) / np.sum(bubble_mask)
            bubble_h = np.max(np.where(bubble_mask.any(axis=1))[0]) - np.min(np.where(bubble_mask.any(axis=1))[0])
            bubble_w = np.max(np.where(bubble_mask.any(axis=0))[0]) - np.min(np.where(bubble_mask.any(axis=0))[0])
            print(f"\n  Bubble: center=({cx:.1f}, {cy:.1f}), size={bubble_w}x{bubble_h}")
            print(f"  Aspect ratio (W/H): {bubble_w/bubble_h:.3f}")

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('input', help='VTM file or vtkoutput directory')
    args = p.parse_args()
    analyze(args.input)
