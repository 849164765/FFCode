#!/usr/bin/env python3
"""Quick KH instability growth rate analyzer.
Parses VTI files, extracts phi=0.5 interface, measures amplitude growth.
"""
import base64, struct, math, os, sys, glob, re

def parse_vti(path):
    """Parse a single VTI block, return (extent, phi_array)."""
    import xml.etree.ElementTree as ET
    tree = ET.parse(path)
    root = tree.getroot()
    img = root.find('ImageData')
    piece = img.find('Piece')
    ext = list(map(int, piece.get('Extent').split()))
    nx = ext[1]-ext[0]+1
    ny = ext[3]-ext[2]+1
    da = piece.find('PointData').find('DataArray')
    name = da.get('Name')
    if name != 'PHI':
        return None  # skip non-PHI blocks
    # VTK binary+base64: header and data may be separate base64 blocks
    # Remove ALL whitespace from the text
    raw_b64 = ''.join(da.text.split())
    raw = base64.b64decode(raw_b64)
    # VTK binary: first 4 bytes = data size (uint32), then data
    n = struct.unpack('<I', raw[:4])[0]
    actual = len(raw) - 4
    if actual < n:
        # Header and data are separate base64 blocks; re-decode data portion
        # Try splitting at the first '==' or base64 boundary
        # Actually VTK encodes header and data as separate base64 sequences
        # joined without separator. Decode them separately.
        return None  # fallback
    vals = struct.unpack(f'<{n//8}d', raw[4:4+n])
    return (ext[0], ext[2], nx, ny), list(vals)

def parse_vti_robust(path):
    """Robust VTI parser handling VTK's separate base64 header/data encoding."""
    import xml.etree.ElementTree as ET
    tree = ET.parse(path)
    root = tree.getroot()
    img = root.find('ImageData')
    piece = img.find('Piece')
    ext = list(map(int, piece.get('Extent').split()))
    nx = ext[1]-ext[0]+1
    ny = ext[3]-ext[2]+1
    da = piece.find('PointData').find('DataArray')
    name = da.get('Name')
    if name != 'PHI':
        return None
    text = ''.join(da.text.split())  # remove whitespace
    # VTK binary+base64: first base64 block encodes the header (uint32),
    # second block encodes the data. They are concatenated without separator.
    # The header is always 4 bytes -> encoded as ceil(4/3)*4 = 8 base64 chars
    # (with padding if needed). Decode first 8 chars as header.
    hdr_b64 = text[:8]
    hdr = base64.b64decode(hdr_b64)
    n = struct.unpack('<I', hdr[:4])[0]
    # Remaining base64 is the data
    data_b64 = text[8:]
    data = base64.b64decode(data_b64)
    nvals = nx * ny
    if len(data) >= nvals * 8:
        vals = struct.unpack(f'<{nvals}d', data[:nvals*8])
    else:
        # Fallback: try decoding all at once
        raw = base64.b64decode(text)
        vals = struct.unpack(f'<{nvals}d', raw[4:4+nvals*8])
    return (ext[0], ext[2], nx, ny), list(vals)

def assemble_phi(vtk_dir, step):
    """Assemble all blocks for a timestep into a 2D phi array."""
    # Find all VTI files for this step
    pattern = os.path.join(vtk_dir, 'vtidata', f'kh2d_T{step}_B*.vti')
    files = sorted(glob.glob(pattern), key=lambda p: int(re.search(r'_B(\d+)', p).group(1)))
    if not files:
        return None
    blocks = []
    max_x, max_y = 0, 0
    min_x, min_y = 1e9, 1e9
    for f in files:
        result = parse_vti_robust(f)
        if result is None:
            continue
        (ox, oy, nx, ny), vals = result
        blocks.append((ox, oy, nx, ny, vals))
        max_x = max(max_x, ox+nx)
        max_y = max(max_y, oy+ny)
        min_x = min(min_x, ox)
        min_y = min(min_y, oy)
    # Assemble
    width = max_x - min_x
    height = max_y - min_y
    phi = [[0.0]*width for _ in range(height)]
    for ox, oy, nx, ny, vals in blocks:
        for j in range(ny):
            for i in range(nx):
                phi[oy+j-min_y][ox+i-min_x] = vals[j*nx+i]
    return phi, width, height

def extract_interface(phi, width, height):
    """Extract interface y-position at each x column (phi=0.5 contour)."""
    interface_y = []
    for i in range(width):
        prev_phi = None
        y_cross = None
        for j in range(height):
            p = phi[j][i]
            if prev_phi is not None:
                # Check if 0.5 is crossed
                if (prev_phi <= 0.5 <= p) or (prev_phi >= 0.5 >= p):
                    # Linear interpolation
                    t = (0.5 - prev_phi) / (p - prev_phi) if abs(p-prev_phi) > 1e-12 else 0.5
                    y_cross = (j-1) + t
                    break
            prev_phi = p
        interface_y.append(y_cross if y_cross is not None else height/2)
    return interface_y

def main():
    vtk_dir = '/home/dlf/myCode/FerroKHinstability/examples/kh2d/vtkoutput'
    steps = [0, 500, 1000]
    print("=" * 60)
    print("KH Instability Growth Rate Analysis (Bom=0, 256x256, L=128)")
    print("=" * 60)

    amplitudes = {}
    for step in steps:
        result = assemble_phi(vtk_dir, step)
        if result is None:
            print(f"Step {step}: no data")
            continue
        phi, width, height = result
        interface = extract_interface(phi, width, height)
        # Amplitude = max - min of interface
        y_max = max(interface)
        y_min = min(interface)
        amp = (y_max - y_min) / 2.0
        # Also compute mean and std
        y_mean = sum(interface) / len(interface)
        y_std = math.sqrt(sum((y - y_mean)**2 for y in interface) / len(interface))
        amplitudes[step] = (amp, y_mean, y_std, y_max, y_min)
        print(f"Step {step:5d}: amp={amp:8.3f}  mean_y={y_mean:7.2f}  std={y_std:7.3f}  max={y_max:7.2f}  min={y_min:7.2f}")

    # Theoretical growth rate
    L = 128  # half-domain for 256x256
    U = 0.02
    k = 2 * math.pi / L  # perturbation wavelength = L
    deltaU = U
    gamma_theory = k * deltaU / 2  # pure KH growth rate
    print(f"\n--- Theoretical (pure KH, no magnetic) ---")
    print(f"L={L}, U={U}, k={k:.6f}, ΔU={deltaU}")
    print(f"γ_theory = k·ΔU/2 = {gamma_theory:.6e} per step")
    print(f"e-folding time = {1/gamma_theory:.1f} steps")

    # Initial perturbation amplitude
    A0 = 0.1 * L  # 0.1*L from Eq.78
    print(f"\nInitial perturbation amplitude A0 = 0.1·L = {A0:.1f}")

    # Compare measured vs theoretical
    if 0 in amplitudes and 500 in amplitudes:
        amp0 = amplitudes[0][0]
        amp500 = amplitudes[500][0]
        if amp0 > 0.1:
            ratio = amp500 / amp0
            gamma_meas = math.log(ratio) / 500
            print(f"\n--- Measured growth rate (step 0→500) ---")
            print(f"amp(0)={amp0:.3f}, amp(500)={amp500:.3f}")
            print(f"ratio={ratio:.4f}, γ_meas = ln(ratio)/500 = {gamma_meas:.6e}")
            print(f"γ_meas/γ_theory = {gamma_meas/gamma_theory:.3f}")

    if 500 in amplitudes and 1000 in amplitudes:
        amp500 = amplitudes[500][0]
        amp1000 = amplitudes[1000][0]
        if amp500 > 0.1:
            ratio = amp1000 / amp500
            gamma_meas = math.log(ratio) / 500
            print(f"\n--- Measured growth rate (step 500→1000) ---")
            print(f"amp(500)={amp500:.3f}, amp(1000)={amp1000:.3f}")
            print(f"ratio={ratio:.4f}, γ_meas = ln(ratio)/500 = {gamma_meas:.6e}")
            print(f"γ_meas/γ_theory = {gamma_meas/gamma_theory:.3f}")

if __name__ == '__main__':
    main()
