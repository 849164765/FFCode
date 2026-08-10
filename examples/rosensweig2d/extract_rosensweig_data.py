#!/usr/bin/env python3
"""Extract time-series diagnostics from Rosensweig VTI files on the cluster.

Only the requested numeric output steps are read. The original VTI files are
never copied or rewritten. The result is one long-format CSV with one row per
case and time step.

Examples:
    # Zero-argument run with built-in defaults (steps 0:end:1000, 128 blocks,
    # PHI+HMAG, strict mode, output rosensweig_diagnostics.csv):
    python3 extract_rosensweig_data.py

    # All available output steps through T20000
    python3 extract_rosensweig_data.py --steps all

    # Only selected steps
    python3 extract_rosensweig_data.py --steps 0,5000,10000,15000,20000

    # Every fifth available output step
    python3 extract_rosensweig_data.py --steps all --stride 5

The CSV is intentionally compact. It stores field statistics and regional
diagnostics, not every grid point. Use it to locate the first bad time step;
the original VTI can then be inspected only for that step if necessary.
"""

import argparse
import base64
import csv
import math
import re
import sys
from pathlib import Path

import numpy as np


NI = 512
NJ = 171
LX_MM = 21.0
LY_MM = 7.0
DX_MM = LX_MM / NI
DY_MM = LY_MM / NJ
WALL_BAND_CELLS = 12
DEFAULT_FIELDS = ("PHI", "HMAG")

# Defaults for a zero-argument run from the rosensweig2d example directory.
# Override any of them via the matching command-line option.
DEFAULTS = {
    "root": ".",
    "steps": "0:end:1000",
    "stride": 1,
    "case_glob": "case_H0_*",
    "expected_blocks": 128,
    "out": "rosensweig_diagnostics.csv",
    "fields": ",".join(DEFAULT_FIELDS),
    "strict": True,
}

VTI_STEP_RE = re.compile(r"_T(\d+)_")
VTI_BLOCK_RE = re.compile(r"_B(\d+)\.vti$")
ATTRIBUTE_RE = re.compile(rb'([A-Za-z_:][\w:.-]*)="([^"]*)"')
DATA_ARRAY_RE = re.compile(rb"<DataArray\b([^>]*)>(.*?)</DataArray>", re.S)


def numeric_step(path):
    match = VTI_STEP_RE.search(path.name)
    return int(match.group(1)) if match else None


def numeric_block(path):
    match = VTI_BLOCK_RE.search(path.name)
    return int(match.group(1)) if match else None


def decode_vti_payload(payload):
    """Decode the base64 payload emitted by FreeLB's binary XML writer."""
    compact = re.sub(rb"\s", b"", payload)
    if len(compact) < 8:
        raise ValueError("invalid or empty VTI base64 payload")
    head = compact[:8]
    rest = compact[8:]
    rest += b"=" * (-len(rest) % 4)
    return base64.b64decode(head) + base64.b64decode(rest)


def parse_attributes(blob):
    return {
        key.decode("ascii"): value.decode("utf-8")
        for key, value in ATTRIBUTE_RE.findall(blob)
    }


def decode_array(payload, type_name, expected_bytes):
    raw = decode_vti_payload(payload)
    if type_name == "Float64":
        dtype = np.dtype("<f8")
    elif type_name == "Float32":
        dtype = np.dtype("<f4")
    else:
        raise ValueError(f"unsupported VTI data type: {type_name}")

    # FreeLB currently uses a four-byte length prefix. Accept an eight-byte
    # prefix too, for VTK configurations using UInt64 headers.
    for prefix_size in (4, 8):
        if len(raw) < prefix_size:
            continue
        payload_size = int.from_bytes(raw[:prefix_size], "little")
        if payload_size != expected_bytes:
            continue
        end = prefix_size + payload_size
        if end <= len(raw):
            return np.frombuffer(raw[prefix_size:end], dtype=dtype).copy()

    raise ValueError(
        f"binary payload length mismatch: expected {expected_bytes} bytes, "
        f"decoded {len(raw)} bytes"
    )


def read_vti(path, wanted_fields):
    data = path.read_bytes()

    image_match = re.search(
        rb'<ImageData\b[^>]*\bOrigin="([^\"]+)"', data
    )
    if image_match is None:
        raise ValueError(f"missing ImageData Origin in {path}")
    origin = [float(value) for value in image_match.group(1).split()[:2]]

    piece_match = re.search(rb'<Piece\b[^>]*\bExtent="([^\"]+)"', data)
    if piece_match is None:
        piece_match = re.search(
            rb'<ImageData\b[^>]*\bWholeExtent="([^\"]+)"', data
        )
    if piece_match is None:
        raise ValueError(f"missing VTI extent in {path}")
    extent = [int(value) for value in piece_match.group(1).split()]
    if len(extent) < 4:
        raise ValueError(f"invalid VTI extent in {path}: {extent}")
    x0, x1, y0, y1 = extent[:4]
    nx = x1 - x0 + 1
    ny = y1 - y0 + 1

    arrays = {}
    wanted = set(wanted_fields)
    for match in DATA_ARRAY_RE.finditer(data):
        attrs = parse_attributes(match.group(1))
        name = attrs.get("Name")
        if name not in wanted:
            continue
        if attrs.get("format") != "binary" or attrs.get("encoding") != "base64":
            raise ValueError(f"unsupported {name} format in {path}")
        components = int(attrs.get("NumberOfComponents", "1"))
        type_name = attrs.get("type", "Float64")
        item_size = 8 if type_name == "Float64" else 4
        expected_bytes = nx * ny * components * item_size
        values = decode_array(match.group(2), type_name, expected_bytes)
        shape = (ny, nx, components) if components > 1 else (ny, nx)
        arrays[name] = (values.reshape(shape), components)

    return origin[0], origin[1], (x0, x1, y0, y1), arrays


def assemble_step(files, wanted_fields):
    """Assemble requested fields from one numeric output step."""
    fields = {}

    for path in files:
        ox, oy, extent, arrays = read_vti(path, wanted_fields)
        x0, x1, y0, y1 = extent
        nx = x1 - x0 + 1
        ny = y1 - y0 + 1

        x_coordinates = ox + x0 + np.arange(nx, dtype=float)
        y_coordinates = oy + y0 + np.arange(ny, dtype=float)
        x_mask = (x_coordinates >= 0.5) & (x_coordinates <= NI - 0.5)
        y_mask = (y_coordinates >= 0.5) & (y_coordinates <= NJ - 0.5)
        x_indices = np.floor(x_coordinates[x_mask]).astype(int)
        y_indices = np.floor(y_coordinates[y_mask]).astype(int)
        local_x = np.flatnonzero(x_mask)
        local_y = np.flatnonzero(y_mask)

        for name, (values, components) in arrays.items():
            if name not in fields:
                shape = (NJ, NI, components) if components > 1 else (NJ, NI)
                fields[name] = np.full(shape, np.nan, dtype=float)
            target = fields[name]
            for local_row, global_row in zip(local_y, y_indices):
                if components == 1:
                    target[global_row, x_indices] = values[local_row, local_x]
                else:
                    target[global_row, x_indices, :] = values[local_row, local_x, :]

    return fields


def extract_interface(phi):
    """Return the phi=0.5 interface in lattice y coordinates."""
    y_interface = np.full(NI, np.nan, dtype=float)
    for x in range(NI):
        column = phi[:, x]
        for row in range(NJ - 1):
            lower = column[row]
            upper = column[row + 1]
            if not (np.isfinite(lower) and np.isfinite(upper)):
                continue
            if (lower - 0.5) * (upper - 0.5) <= 0 and lower != upper:
                y_interface[x] = (
                    row + 0.5 + (0.5 - lower) / (upper - lower)
                )
                break
    return y_interface


def finite_values(values, mask=None):
    if mask is not None:
        values = values[mask]
    return values[np.isfinite(values)]


def add_stats(row, prefix, values, mask=None):
    finite = finite_values(values, mask)
    if finite.size == 0:
        row[f"{prefix}_min"] = math.nan
        row[f"{prefix}_max"] = math.nan
        row[f"{prefix}_mean"] = math.nan
        return
    row[f"{prefix}_min"] = float(np.min(finite))
    row[f"{prefix}_max"] = float(np.max(finite))
    row[f"{prefix}_mean"] = float(np.mean(finite))


def add_max(row, name, values, mask=None):
    finite = finite_values(values, mask)
    row[name] = float(np.max(finite)) if finite.size else math.nan


def case_h0(case):
    match = re.match(r"case_H0_([-+0-9.eE]+)$", case.name)
    return float(match.group(1)) if match else math.nan


def available_steps(vtidata):
    steps = {numeric_step(path) for path in vtidata.glob("rosensweig2d_T*_B*.vti")}
    return sorted(step for step in steps if step is not None)


def parse_step_spec(spec, available, stride):
    if stride < 1:
        raise ValueError("--stride must be at least 1")
    if spec.lower() == "all":
        return available[::stride]

    selected = set()
    for token in (part.strip() for part in spec.split(",")):
        if not token:
            continue
        if ":" not in token:
            selected.add(int(token))
            continue
        parts = token.split(":")
        if len(parts) not in (2, 3):
            raise ValueError(f"invalid step range: {token}")
        start = int(parts[0]) if parts[0] else available[0]
        raw_stop = parts[1]
        if raw_stop.lower() in ("end", "max"):
            stop = available[-1] if available else -1
        else:
            stop = int(raw_stop)
        range_stride = int(parts[2]) if len(parts) == 3 and parts[2] else 1
        if range_stride < 1:
            raise ValueError(f"invalid step range stride: {token}")
        if available and stop < 0:
            continue
        selected.update(range(start, stop + 1, range_stride))
    return sorted(selected)[::stride]


def region_masks():
    wall = np.zeros((NJ, NI), dtype=bool)
    wall[:WALL_BAND_CELLS, :] = True
    wall[-WALL_BAND_CELLS:, :] = True
    interior = ~wall

    seam = np.zeros((NJ, NI), dtype=bool)
    seam[:, 0] = True
    seam[:, -1] = True
    return wall, interior, seam


def max_abs_jump(values):
    if values is None or values.shape[1] < 2:
        return math.nan
    difference = np.abs(values[:, 0] - values[:, -1])
    finite = difference[np.isfinite(difference)]
    return float(np.max(finite)) if finite.size else math.nan


def extract_case_step(case, step, wanted_fields, expected_blocks):
    vtidata = case / "vtkoutput" / "vtidata"
    files = sorted(
        path
        for path in vtidata.glob(f"rosensweig2d_T{step}_B*.vti")
        if numeric_step(path) == step
    )
    row = {
        "case": case.name,
        "H0_kAm": case_h0(case),
        "step": step,
        "status": "ok",
        "vti_files": len(files),
        "expected_vti_files": expected_blocks,
    }
    if not files:
        row["status"] = "missing_step"
        row["available_steps"] = ";".join(map(str, available_steps(vtidata)))
        return row

    block_ids = {numeric_block(path) for path in files}
    expected_ids = set(range(expected_blocks)) if expected_blocks else None
    if expected_ids is not None and (
        len(files) != expected_blocks or block_ids != expected_ids
    ):
        row["status"] = "block_count_mismatch"
        row["block_ids"] = ";".join(
            str(block_id) for block_id in sorted(block_ids) if block_id is not None
        )
        return row

    fields = assemble_step(files, wanted_fields)
    phi = fields.get("PHI")
    if phi is None:
        raise ValueError(f"PHI is missing from {case}")

    wall, interior, seam = region_masks()
    phi_finite = np.isfinite(phi)
    row["phi_cells"] = int(phi_finite.sum())
    row["phi_nan_cells"] = int(phi.size - phi_finite.sum())
    add_stats(row, "phi", phi)
    row["phi_floor_fraction"] = (
        float(np.mean(phi[phi_finite] <= 0.010001)) if phi_finite.any() else math.nan
    )
    row["phi_cap_fraction"] = (
        float(np.mean(phi[phi_finite] >= 0.998999)) if phi_finite.any() else math.nan
    )
    add_stats(row, "phi_wall", phi, wall)
    add_stats(row, "phi_interior", phi, interior)
    add_stats(row, "phi_seam", phi, seam)
    add_stats(row, "phi_bottom", phi[:WALL_BAND_CELLS, :])
    add_stats(row, "phi_top", phi[-WALL_BAND_CELLS:, :])
    add_stats(row, "phi_bottom_deviation", np.abs(phi[:WALL_BAND_CELLS, :] - 1.0))
    add_stats(row, "phi_top_deviation", np.abs(phi[-WALL_BAND_CELLS:, :]))
    row["phi_seam_jump_max"] = max_abs_jump(phi)

    y_interface = extract_interface(phi)
    valid_interface = y_interface[np.isfinite(y_interface)]
    row["interface_columns"] = int(valid_interface.size)
    if valid_interface.size:
        interface_mean = float(np.mean(valid_interface))
        row["interface_mean_mm"] = interface_mean * DY_MM
        row["peak_height_mm"] = (
            float(np.max(valid_interface)) - interface_mean
        ) * DY_MM
        row["valley_depth_mm"] = (
            interface_mean - float(np.min(valid_interface))
        ) * DY_MM
    else:
        row["interface_mean_mm"] = math.nan
        row["peak_height_mm"] = math.nan
        row["valley_depth_mm"] = math.nan

    for name in ("PSI", "HX", "HY", "Density"):
        if name in fields:
            add_stats(row, name.lower(), fields[name])

    hmag = fields.get("HMAG")
    if hmag is None and "HX" in fields and "HY" in fields:
        hmag = np.sqrt(fields["HX"] ** 2 + fields["HY"] ** 2)
    if hmag is not None:
        add_stats(row, "hmag", hmag)
        add_stats(row, "hmag_wall", hmag, wall)
        add_stats(row, "hmag_interior", hmag, interior)
        add_stats(row, "hmag_seam", hmag, seam)
        row["hmag_seam_jump_max"] = max_abs_jump(hmag)

    force_mag = None
    if "Force" in fields and fields["Force"].ndim == 3:
        force = fields["Force"]
        force_mag = np.sqrt(force[:, :, 0] ** 2 + force[:, :, 1] ** 2)
        add_stats(row, "force_mag", force_mag)
        add_max(row, "force_wall_max", force_mag, wall)
        add_max(row, "force_interior_max", force_mag, interior)
        add_max(row, "force_seam_max", force_mag, seam)

    if "Velocity" in fields and fields["Velocity"].ndim == 3:
        velocity = fields["Velocity"]
        velocity_mag = np.sqrt(velocity[:, :, 0] ** 2 + velocity[:, :, 1] ** 2)
        add_stats(row, "velocity_mag", velocity_mag)
        add_max(row, "velocity_wall_max", velocity_mag, wall)
        add_max(row, "velocity_interior_max", velocity_mag, interior)

    return row


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(DEFAULTS["root"]))
    steps = parser.add_mutually_exclusive_group()
    steps.add_argument("--step", type=int, help="one exact output step")
    steps.add_argument(
        "--steps",
        default=DEFAULTS["steps"],
        help="all, comma-separated steps, or inclusive ranges such as 0:end:1000",
    )
    parser.add_argument("--stride", type=int, default=DEFAULTS["stride"])
    parser.add_argument("--case-glob", default=DEFAULTS["case_glob"])
    parser.add_argument(
        "--expected-blocks",
        type=int,
        default=DEFAULTS["expected_blocks"],
        help="expected B*.vti block count; use 0 to disable the check",
    )
    parser.add_argument("--out", type=Path, default=Path(DEFAULTS["out"]))
    parser.add_argument(
        "--fields",
        default=DEFAULTS["fields"],
        help="comma-separated VTI fields to load",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        default=DEFAULTS["strict"],
        help="return nonzero if any requested case/step is missing or fails",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    if args.step is not None:
        step_spec = str(args.step)
    else:
        step_spec = args.steps
    wanted_fields = tuple(
        field.strip() for field in args.fields.split(",") if field.strip()
    )
    root = args.root.resolve()
    cases = sorted(
        (path for path in root.glob(args.case_glob) if path.is_dir()),
        key=lambda path: (math.isnan(case_h0(path)), case_h0(path), path.name),
    )
    if not cases:
        print(f"[error] no case directories matching {args.case_glob!r} under {root}")
        return 2

    rows = []
    failures = 0
    for case in cases:
        vtidata = case / "vtkoutput" / "vtidata"
        available = available_steps(vtidata)
        try:
            selected_steps = parse_step_spec(step_spec, available, args.stride)
        except ValueError as error:
            print(f"[error] {error}", file=sys.stderr)
            return 2
        if not selected_steps:
            failures += 1
            rows.append(
                {
                    "case": case.name,
                    "H0_kAm": case_h0(case),
                    "step": "",
                    "status": "no_steps",
                    "available_steps": ";".join(map(str, available)),
                }
            )
            print(f"[skip] {case.name}: no selected steps")
            continue

        for step in selected_steps:
            try:
                row = extract_case_step(
                    case, step, wanted_fields, args.expected_blocks or None
                )
            except Exception as error:  # report one bad case/step and continue
                failures += 1
                row = {
                    "case": case.name,
                    "H0_kAm": case_h0(case),
                    "step": step,
                    "status": f"error:{type(error).__name__}",
                    "vti_files": 0,
                    "expected_vti_files": args.expected_blocks,
                }
                print(f"[error] {case.name} T{step}: {error}", file=sys.stderr)
            rows.append(row)
            if row["status"] == "ok":
                print(
                    f"[ok] {case.name} T{step}: {row['vti_files']} VTI, "
                    f"peak={row.get('peak_height_mm', math.nan):.6f} mm, "
                    f"valley={row.get('valley_depth_mm', math.nan):.6f} mm, "
                    f"Hmax={row.get('hmag_max', math.nan):.6g}"
                )
            else:
                if row["status"] in ("missing_step", "block_count_mismatch"):
                    failures += 1
                print(
                    f"[skip] {case.name} T{step}: {row['status']} "
                    f"({row.get('available_steps', 'no available-step list')})"
                )

        output = args.out or root / "rosensweig_diagnostics.csv"
        output.parent.mkdir(parents=True, exist_ok=True)
    preferred = ["case", "H0_kAm", "step", "status", "vti_files"]
    all_keys = set().union(*(row.keys() for row in rows))
    fieldnames = preferred + sorted(all_keys - set(preferred))
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    print(f"saved unified diagnostic CSV: {output}")
    return 2 if args.strict and failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
