#!/usr/bin/env python3
"""Comprehensive analysis: peak height, phi leakage, and artifact check."""
import sys, os, glob
sys.path.insert(0, os.path.dirname(__file__))
import plot_vti
import numpy as np

def load_field(step, vti_dir=None):
    """Load field from a specific vti_dir."""
    if vti_dir is None:
        vti_dir = plot_vti.VTI_DIR
    pattern = os.path.join(vti_dir, f'rosensweig2d_T{step}_B*.vti')
    files = sorted(glob.glob(pattern), key=lambda f: int(f.split('_B')[1].split('.')[0]))
    if not files:
        return None
    all_data = []
    max_x = max_y = 0
    for f in files:
        ox, oy, nx, ny, data = plot_vti.parse_vti(f)
        ix0, iy0 = int(round(ox)), int(round(oy))
        all_data.append((ix0, iy0, nx, ny, data))
        max_x = max(max_x, ix0 + nx)
        max_y = max(max_y, iy0 + ny)
    field = np.full((max_y, max_x), np.nan)
    for ix0, iy0, nx, ny, data in all_data:
        field[iy0:iy0+ny, ix0:ix0+nx] = data.reshape(ny, nx)
    return field

def measure_peaks(field):
    """Measure peak amplitude: interface at phi=0.5, find max deviation from y0=28."""
    interior = field[1:-1, 1:-1]
    ny_int, nx_int = interior.shape
    y0 = 28
    amplitudes = []
    for i in range(nx_int):
        col = interior[:, i]
        for j in range(ny_int - 1):
            if (col[j] - 0.5) * (col[j+1] - 0.5) <= 0:
                denom = col[j+1] - col[j]
                frac = (0.5 - col[j]) / denom if abs(denom) > 1e-12 else 0.5
                amplitudes.append(j + frac - y0)
                break
    if amplitudes:
        return max(amplitudes), min(amplitudes), np.mean(np.array(amplitudes))
    return 0, 0, 0

def check_leakage(field):
    """Check phi leakage at bottom rows (should be ~1.0 for heavy ferrofluid)."""
    results = {}
    for j in range(6):
        row = field[j, 1:-1]
        results[j] = (float(row.min()), float(row.mean()))
    return results

def check_artifacts(field):
    """Check for deep trough artifacts (interface >10 cells below y0=28)."""
    interior = field[1:-1, 1:-1]
    ny_int, nx_int = interior.shape
    y0 = 28
    abnormal = []
    for i in range(nx_int):
        col = interior[:, i]
        for j in range(ny_int - 1):
            if (col[j] - 0.5) * (col[j+1] - 0.5) <= 0:
                denom = col[j+1] - col[j]
                frac = (0.5 - col[j]) / denom if abs(denom) > 1e-12 else 0.5
                y_cross = j + frac
                if y_cross - y0 < -10:
                    abnormal.append((int(i+1), float(y_cross)))
                break
    return abnormal

steps = [0, 1000, 5000, 10000, 20000, 30000, 50000, 70000, 100000]
mm_per_cell = 21.0 / 252

print("=" * 80)
print("COMPREHENSIVE ANALYSIS: Peak Height + Leakage + Artifacts")
print("=" * 80)
print(f"{'Step':>6}  {'Peak_mm':>8}  {'Trough_mm':>10}  {'Row4_min':>9}  {'Artifacts':>10}")
print("-" * 80)

for step in steps:
    field = load_field(step)
    if field is None:
        continue
    peak, trough, mean = measure_peaks(field)
    leakage = check_leakage(field)
    artifacts = check_artifacts(field)
    row4_min = leakage[4][0] if 4 in leakage else 0
    n_art = len(artifacts)
    print(f"{step:>6}  {peak*mm_per_cell:>+8.4f}  {trough*mm_per_cell:>+10.4f}  {row4_min:>9.4f}  {n_art:>10}")

print("\n" + "=" * 80)
print("DETAILED LEAKAGE CHECK (bottom rows, phi should be ~1.0)")
print("=" * 80)
for step in [1000, 10000, 50000, 100000]:
    field = load_field(step)
    if field is None:
        continue
    leakage = check_leakage(field)
    print(f"\nStep {step}:")
    for j in range(6):
        minv, meanv = leakage[j]
        print(f"  row {j}: min={minv:.4f} mean={meanv:.4f}")

print("\n" + "=" * 80)
print("ARTIFACT CHECK (columns with interface >10 cells below y0=28)")
print("=" * 80)
for step in [50000, 70000, 100000]:
    field = load_field(step)
    if field is None:
        continue
    artifacts = check_artifacts(field)
    if artifacts:
        print(f"Step {step}: {len(artifacts)} artifacts at x={[(x,round(yc,1)) for x,yc in artifacts]}")
    else:
        print(f"Step {step}: NO artifacts (clean!)")
