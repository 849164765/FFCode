#!/usr/bin/env python3
"""Parse VTI files and plot the PHI field for the Rosensweig2D example."""
import os
import sys
import base64
import struct
import glob
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

VTI_DIR = os.path.join(os.path.dirname(__file__), 'vtkoutput', 'vtidata')

def parse_vti(filepath):
    """Parse a single VTI file, return (origin_x, origin_y, nx, ny, data_1d)."""
    tree = ET.parse(filepath)
    root = tree.getroot()
    img = root.find('.//ImageData')
    origin = list(map(float, img.get('Origin').split()))
    piece = root.find('.//Piece')
    extent = list(map(int, piece.get('Extent').split()))
    nx = extent[1] - extent[0] + 1
    ny = extent[3] - extent[2] + 1
    da = root.find('.//DataArray')
    text = da.text.strip()
    # VTK "binary"+"base64" format: header and data are SEPARATE base64 blocks.
    # Header is always 4 bytes (int32) → 8 base64 chars (with '==' padding).
    header_b64 = text[:8]
    data_b64 = text[8:]
    header_bytes = base64.b64decode(header_b64)
    nbytes = struct.unpack('<i', header_bytes)[0]
    data_bytes = base64.b64decode(data_b64)
    data = np.frombuffer(data_bytes[:nbytes], dtype='<f8')
    return origin[0], origin[1], nx, ny, data

def assemble_field(step):
    """Assemble the full 2D field from all block VTI files for a given step."""
    pattern = os.path.join(VTI_DIR, f'rosensweig2d_T{step}_B*.vti')
    files = sorted(glob.glob(pattern), key=lambda f: int(f.split('_B')[1].split('.')[0]))
    if not files:
        print(f"No VTI files found for step {step}")
        return None
    all_data = []
    max_x = 0
    max_y = 0
    for f in files:
        ox, oy, nx, ny, data = parse_vti(f)
        ix0 = int(round(ox))      # cell-center position → 0-based index
        iy0 = int(round(oy))
        all_data.append((ix0, iy0, nx, ny, data))
        max_x = max(max_x, ix0 + nx)
        max_y = max(max_y, iy0 + ny)
    field = np.full((max_y, max_x), np.nan)
    for ix0, iy0, nx, ny, data in all_data:
        block = data.reshape(ny, nx)
        field[iy0:iy0+ny, ix0:ix0+nx] = block
    return field

def plot_step(step, outpath):
    field = assemble_field(step)
    if field is None:
        return
    fig, ax = plt.subplots(figsize=(14, 4))
    im = ax.imshow(field, origin='lower', aspect='equal', cmap='RdBu_r',
                   vmin=-1, vmax=1, interpolation='nearest')
    ax.set_title(f'PHI field — Step {step}')
    ax.set_xlabel('X (lattice)')
    ax.set_ylabel('Y (lattice)')
    plt.colorbar(im, ax=ax, label='PHI')
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()
    print(f'Saved {outpath}')

if __name__ == '__main__':
    steps = [0, 200] if len(sys.argv) < 2 else [int(s) for s in sys.argv[1:]]
    for step in steps:
        outpath = os.path.join(os.path.dirname(__file__), f'phi_step{step}.png')
        plot_step(step, outpath)
