#!/usr/bin/env python3
"""Analyze peak height trajectory from VTK output."""
import sys
import os
sys.path.insert(0, os.path.dirname(__file__))
from plot_vti import assemble_field
import numpy as np

def measure_peaks(field):
    """Measure peak amplitude: interface at phi=0.5, find max deviation from y0=28."""
    ny, nx = field.shape
    interior = field[1:-1, 1:-1]
    ny_int, nx_int = interior.shape
    y0 = 28  # interface position (84/3)

    amplitudes = []
    for i in range(nx_int):
        col = interior[:, i]
        for j in range(ny_int - 1):
            if (col[j] - 0.5) * (col[j+1] - 0.5) <= 0:
                denom = col[j+1] - col[j]
                if abs(denom) > 1e-12:
                    frac = (0.5 - col[j]) / denom
                else:
                    frac = 0.5
                y_cross = j + frac
                amplitudes.append(y_cross - y0)
                break

    if amplitudes:
        return max(amplitudes), min(amplitudes), np.mean(np.array(amplitudes))
    return 0, 0, 0

steps = [0, 1000, 2000, 3000, 4000, 5000, 7000, 10000, 15000, 20000, 30000, 50000, 70000, 100000]
mm_per_cell = 21.0 / 252
print(f"{'Step':>6}  {'Peak':>8}  {'Trough':>8}  {'Mean':>8}  {'Peak_mm':>8}")
print("-" * 52)
for step in steps:
    field = assemble_field(step)
    if field is None:
        print(f"  Step {step}: no data")
        continue
    peak, trough, mean = measure_peaks(field)
    print(f"{step:>6}  {peak:>+8.3f}  {trough:>+8.3f}  {mean:>+8.4f}  {peak*mm_per_cell:>8.4f}")
