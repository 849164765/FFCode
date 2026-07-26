#!/usr/bin/env python3
"""Analyze bubble deformation from VTK output files."""
import os
import sys
import base64
import struct
import numpy as np
import xml.etree.ElementTree as ET
from glob import glob

def parse_vti_block(vti_path):
    """Parse a VTI file and extract PHI field data."""
    tree = ET.parse(vti_path)
    root = tree.getroot()

    # Get extent and origin
    image_data = root.find('ImageData')
    whole_extent = list(map(int, image_data.get('WholeExtent').split()))
    origin = list(map(float, image_data.get('Origin').split()))
    spacing = list(map(float, image_data.get('Spacing').split()))

    piece = image_data.find('Piece')
    piece_extent = list(map(int, piece.get('Extent').split()))

    # Get PHI data array
    point_data = piece.find('PointData')
    if point_data is None:
        return None, None

    phi_array = None
    for data_array in point_data.findall('DataArray'):
        if data_array.get('Name') == 'PHI':
            phi_array = data_array
            break

    if phi_array is None:
        return None, None

    # Decode base64 data
    encoding = phi_array.get('encoding', 'base64')
    data_format = phi_array.get('format', 'binary')
    raw_text = phi_array.text.strip()

    if encoding != 'base64':
        return None, None

    # VTK base64 binary format: header and data are SEPARATE base64 blocks
    # concatenated in the XML text.
    # Header block: base64-encoded uint32 (4 bytes -> 8 base64 chars with ==)
    # Data block: base64-encoded array (N bytes -> ceil(N/3)*4 base64 chars)
    # Must decode them SEPARATELY due to padding.

    if data_format == 'binary':
        # Decode header (first 8 base64 chars = 4 bytes uint32)
        header_b64 = raw_text[:8]
        header_bytes = base64.b64decode(header_b64)
        header_size = struct.unpack('<I', header_bytes)[0]

        # Decode data (remaining base64 chars)
        data_b64 = raw_text[8:]
        data_bytes = base64.b64decode(data_b64)

        if len(data_bytes) < header_size:
            return None, None

        phi = np.frombuffer(data_bytes[:header_size], dtype=np.float64)
    else:
        return None, None

    # Reshape to 2D grid
    x0, x1, y0, y1 = piece_extent[0], piece_extent[1], piece_extent[2], piece_extent[3]
    nx = x1 - x0 + 1
    ny = y1 - y0 + 1

    if len(phi) != nx * ny:
        # Try without header (appended format)
        phi = np.frombuffer(all_data, dtype=np.float64)
        if len(phi) != nx * ny:
            return None, None

    phi_grid = phi.reshape(ny, nx)

    return (piece_extent, origin, spacing), phi_grid


def load_timestep(timestep, vtidir='./vtkoutput/vtidata'):
    """Load all blocks for a given timestep and combine into global grid."""
    pattern = os.path.join(vtidir, f'bubbleMag2d_T{timestep}_B*.vti')
    files = sorted(glob(pattern), key=lambda x: int(x.split('_B')[1].split('.')[0]))

    if not files:
        print(f"No files found for timestep {timestep}")
        return None, None

    blocks = []
    for f in files:
        result = parse_vti_block(f)
        if result[0] is not None:
            blocks.append(result)

    if not blocks:
        return None, None

    # Compute global cell positions using Origin field
    # Origin gives the position of the bottom-left corner of the block
    # Cell center at local index (i,j) is at global (origin_x+0.5+i, origin_y+0.5+j)
    # So global cell index = int(origin + 0.5) + local_index
    all_x0 = []
    all_y0 = []
    for (extent, origin, spacing), phi in blocks:
        gx0 = int(round(origin[0] + 0.5))
        gy0 = int(round(origin[1] + 0.5))
        all_x0.append(gx0)
        all_y0.append(gy0)

    min_x = min(all_x0)
    min_y = min(all_y0)
    max_x = max(x0 + phi.shape[1] - 1 for (x0, phi) in zip(all_x0, [b[1] for b in blocks]))
    max_y = max(y0 + phi.shape[0] - 1 for (y0, phi) in zip(all_y0, [b[1] for b in blocks]))

    global_nx = max_x - min_x + 1
    global_ny = max_y - min_y + 1

    phi_global = np.full((global_ny, global_nx), np.nan)

    for (extent, origin, spacing), phi in blocks:
        gx0 = int(round(origin[0] + 0.5))
        gy0 = int(round(origin[1] + 0.5))
        ny, nx = phi.shape
        phi_global[gy0-min_y:gy0-min_y+ny, gx0-min_x:gx0-min_x+nx] = phi

    return phi_global, (min_x, max_x, min_y, max_y)


def measure_bubble(phi, threshold=0.5):
    """Measure bubble dimensions from phi field.

    phi=0 inside bubble, phi=1 outside (based on initialization).
    Bubble is where phi < threshold.
    """
    # Find bubble region (phi < threshold)
    bubble_mask = phi < threshold

    if not bubble_mask.any():
        return None

    # Find bounding box
    rows = np.any(bubble_mask, axis=1)
    cols = np.any(bubble_mask, axis=0)
    y_indices = np.where(rows)[0]
    x_indices = np.where(cols)[0]

    if len(y_indices) == 0 or len(x_indices) == 0:
        return None

    # Bubble dimensions
    Dx = x_indices[-1] - x_indices[0] + 1  # horizontal diameter
    Dy = y_indices[-1] - y_indices[0] + 1  # vertical diameter

    # Center
    cx = (x_indices[0] + x_indices[-1]) / 2.0
    cy = (y_indices[0] + y_indices[-1]) / 2.0

    # Elongation (Dy/Dx > 1 means elongated in y/field direction)
    elong = Dy / Dx if Dx > 0 else 0

    return {
        'Dx': Dx,
        'Dy': Dy,
        'elong': elong,
        'cx': cx,
        'cy': cy,
        'area': bubble_mask.sum()
    }


def main():
    vtidir = './vtkoutput/vtidata'
    if len(sys.argv) > 1:
        vtidir = sys.argv[1]

    # Find all available timesteps
    all_files = glob(os.path.join(vtidir, 'bubbleMag2d_T*_B0.vti'))
    timesteps = sorted(set(int(f.split('_T')[1].split('_B')[0]) for f in all_files))

    print(f"Found timesteps: {timesteps}")
    print(f"{'Step':>6} | {'Dx':>6} | {'Dy':>6} | {'elong':>8} | {'cx':>7} | {'cy':>7} | {'area':>8}")
    print("-" * 65)

    for t in timesteps:
        phi, extent = load_timestep(t, vtidir)
        if phi is None:
            continue

        result = measure_bubble(phi)
        if result:
            print(f"{t:6d} | {result['Dx']:6d} | {result['Dy']:6d} | {result['elong']:8.4f} | {result['cx']:7.1f} | {result['cy']:7.1f} | {result['area']:8d}")
        else:
            print(f"{t:6d} | No bubble found")


if __name__ == '__main__':
    main()
