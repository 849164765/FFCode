#!/usr/bin/env python3
"""plot_mag_bubble.py — VTK VTM/VTI → PNG bubble shape visualization
Usage: python3 plot_mag_bubble.py <vtm_file_or_dir> [--step N] [--field PHI]
Reads PHI field from block-structured VTI files, plots phi=0.5 contour"""

import sys, os, glob, xml.etree.ElementTree as ET, base64, struct, argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def parse_vti(fpath):
    """Parse binary base64 ImageData VTI, return {name: np.array} dict"""
    tree = ET.parse(fpath)
    root = tree.getroot()
    image = root.find('ImageData')
    piece = image.find('Piece')
    ext = [int(x) for x in (piece.get('Extent') or image.get('WholeExtent')).split()]
    nx, ny = ext[1]-ext[0], ext[3]-ext[2]  # +1 for inclusive? VTI ext is inclusive
    origin = [float(x) for x in (image.get('Origin') or '0 0 0').split()]
    spacing = [float(x) for x in (image.get('Spacing') or '1 1 1').split()]

    fields = {}
    for da in (piece.find('PointData') or piece.find('CellData') or []):
        name = da.get('Name')
        ncomp = int(da.get('NumberOfComponents', 1))
        fmt = da.get('format', 'ascii')
        raw = da.text.strip() if da.text else ''
        if fmt == 'binary':
            # VTI "appended" format: first 8 chars = base64-encoded Int32 header
            # e.g. '4BwAAA==' → 4 bytes → byte count of payload
            clean = raw.replace('\n','').strip()
            header_b64 = clean[:8]
            header_bytes = base64.b64decode(header_b64)
            if len(header_bytes) < 4:
                header_bytes = header_bytes + b'\x00' * (4 - len(header_bytes))
            header_size = struct.unpack('<I', header_bytes)[0]
            data_b64 = clean[8:]
            bin_data = base64.b64decode(data_b64)
            vti_type = da.get('type', 'Float64')
            ftype = 'f' if vti_type=='Float32' else 'd'
            fsize = struct.calcsize(ftype)
            actual_n = len(bin_data) // fsize
            arr = np.array(struct.unpack(f"<{actual_n}{ftype}", bin_data))
            if ncomp > 1:
                arr = arr.reshape((-1, ncomp))
            fields[name] = arr, nx, ny, origin, spacing
        else:
            vals = [float(x) for x in raw.split()]
            arr = np.array(vals)
            if ncomp > 1:
                arr = arr.reshape((-1, ncomp))
            fields[name] = arr, nx, ny, origin, spacing
    return fields

def assemble_field(vtm_path, field_name='PHI'):
    """Read all blocks from VTM, assemble field array"""
    vtm_dir = os.path.dirname(os.path.abspath(vtm_path))
    tree = ET.parse(vtm_path)
    root = tree.getroot()
    mbd = root.find('vtkMultiBlockDataSet')

    # First pass: read metadata
    blocks = []
    for block in mbd.findall('Block'):
        ds = block.find('DataSet')
        fname = ds.get('file')
        fidx = ds.get('index', '0')
        vti_path = os.path.join(vtm_dir, fname)
        if not os.path.exists(vti_path):
            continue
        tree2 = ET.parse(vti_path)
        img = tree2.find('ImageData') or tree2.getroot()
        piece = img.find('Piece')
        ext = [int(x) for x in (piece.get('Extent') or img.get('WholeExtent')).split()]
        origin = [float(x) for x in (img.get('Origin') or '0 0 0').split()]
        spacing = [float(x) for x in (img.get('Spacing') or '1 1 1').split()]
        blocks.append({'file': vti_path, 'extent': ext, 'origin': origin, 'spacing': spacing})

    if not blocks:
        return None, None, None

    # Determine global grid
    x_min = min(b['origin'][0] for b in blocks)
    y_min = min(b['origin'][1] for b in blocks)
    x_max = max(b['origin'][0] + b['extent'][1]*b['spacing'][0] for b in blocks)
    y_max = max(b['origin'][1] + b['extent'][3]*b['spacing'][1] for b in blocks)
    dx = blocks[0]['spacing'][0]
    dy = blocks[0]['spacing'][1]

    gx = int((x_max - x_min) / dx) + 1
    gy = int((y_max - y_min) / dy) + 1

    field = np.zeros((gy, gx))
    mask = np.zeros((gy, gx), dtype=bool)

    for b in blocks:
        parsed = parse_vti(b['file'])
        if field_name not in parsed:
            continue
        arr, nx, ny, ox, sp = parsed[field_name]
        # arr is flat, shape (n_points, ncomp) or (n_points,)
        if arr.ndim > 1 and arr.shape[1] > 1:
            arr = arr[:, 0]  # take first component
        arr = arr.flatten()
        npts = (nx+1)*(ny+1)
        if len(arr) > npts:
            arr = arr[:npts]
        elif len(arr) < npts:
            arr2 = np.zeros(npts)
            arr2[:len(arr)] = arr
            arr = arr2
        arr2d = arr.reshape((ny+1, nx+1)).T  # VTI point order: x fastest

        ix = int(round((ox[0] - x_min) / dx))
        iy = int(round((ox[1] - y_min) / dy))
        ex = min(ix + nx + 1, gx)
        ey = min(iy + ny + 1, gy)

        bx = ex - ix
        by = ey - iy
        field[iy:ey, ix:ex] = arr2d[:bx, :by].T  # transpose for (y,x) order
        mask[iy:ey, ix:ex] = True

    field[~mask] = np.nan
    return field, (x_min, x_max, y_min, y_max), dx

def plot_bubble(vtm_path, out_path=None, field_name='PHI', contour_val=0.5):
    field, extent, dx = assemble_field(vtm_path, field_name)
    if field is None:
        print(f"ERROR: field '{field_name}' not found in {vtm_path}")
        return

    x_min, x_max, y_min, y_max = extent
    h = y_max - y_min

    fig, ax = plt.subplots(1, 1, figsize=(8, h/x_max*8 if x_max > 0 else 8))
    x = np.linspace(x_min, x_max + dx, field.shape[1])
    y = np.linspace(y_min, y_max + dx, field.shape[0])

    # Fill background
    im = ax.pcolormesh(x, y, field, shading='auto', cmap='RdBu', vmin=0, vmax=1)
    # phi=0.5 contour (bubble interface)
    ax.contour(x, y, field, levels=[contour_val], colors='black', linewidths=2)
    ax.set_aspect('equal')
    ax.set_title(f'{field_name} field (contour at {contour_val})')
    plt.colorbar(im, ax=ax, label=field_name)

    if out_path is None:
        out_dir = os.path.dirname(os.path.abspath(vtm_path))
        out_path = os.path.join(out_dir, f'bubble_shape_{field_name}.png')
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f"Saved: {out_path}")
    plt.close()

def main():
    p = argparse.ArgumentParser(description='Plot bubble shape from VTK VTM/VTI files')
    p.add_argument('input', help='VTM file path or vtkoutput directory')
    p.add_argument('--step', type=int, default=None, help='time step to plot (default: latest)')
    p.add_argument('--field', default='PHI', help='field name (default: PHI)')
    p.add_argument('--output', default=None, help='output PNG path')
    args = p.parse_args()

    vtm = args.input
    if os.path.isdir(vtm):
        vtms = sorted(glob.glob(os.path.join(vtm, '*.vtm')))
        if args.step is not None:
            # filter to files matching the requested step
            step_vtms = [v for v in vtms if f'_T{args.step}_' in v or v.endswith(f'{args.step}.vtm')]
            if step_vtms:
                vtm = step_vtms[0]
                print(f"Using step {args.step}: {vtm}")
            else:
                print(f"ERROR: No VTM found for step {args.step}, available: {[os.path.basename(v) for v in vtms]}")
                return
        else:
            # Find max step
            steps = []
            for v in vtms:
                bn = os.path.basename(v)
                try:
                    s = int(bn.split('.')[0].split('_T')[-1] if '_T' in bn else bn.split('.')[0].replace('bubbleMag2d',''))
                    steps.append((s, v))
                except: pass
            if steps:
                steps.sort(key=lambda x: x[0])
                vtm = steps[-1][1]
                print(f"Using step {steps[-1][0]}: {vtm}")
            elif vtms:
                vtm = vtms[-1]
    plot_bubble(vtm, args.output, args.field)

if __name__ == '__main__':
    main()
