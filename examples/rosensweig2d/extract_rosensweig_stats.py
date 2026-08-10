#!/usr/bin/env python3
"""Extract time-series statistics from Rosensweig VTI files on the cluster.

Successor of extract_rosensweig_data.py with two fixes:

1. The expected block count is auto-detected per case from the VTI files
   themselves, so sweeps run with any MPI rank count extract cleanly. The old
   script hard-coded 128 blocks, and the current sweep (256 MPI ranks)
   produced a CSV where every step was flagged block_count_mismatch and no
   statistics were computed. Override with --expected-blocks to force an
   explicit count instead of auto-detection.
2. Extraction runs in parallel over (case, step) pairs with a process pool
   (--workers). A full sweep produces tens of thousands of VTI files
   (e.g. 11 cases x 21 steps x 256 blocks = 59136 files).

Only the requested numeric output steps are read. The original VTI files are
never copied or rewritten. The result is one long-format CSV with one row per
case and time step. Column names are identical to extract_rosensweig_data.py,
so existing plotting pipelines keep working.

Examples:
    # Zero-argument run from the rosensweig2d example directory:
    # auto-detect block count, steps 0:end:1000, PHI+HMAG:
    python3 extract_rosensweig_stats.py

    # All available output steps through T20000, every fifth step:
    python3 extract_rosensweig_stats.py --steps all --stride 5

    # 64 parallel workers on the cluster:
    python3 extract_rosensweig_stats.py --workers 64

    # Force an explicit block count instead of auto-detection:
    python3 extract_rosensweig_stats.py --expected-blocks 256

    # Return nonzero if any case/step is missing or fails (slurm/CI):
    python3 extract_rosensweig_stats.py --strict
"""

import argparse
import base64
import csv
import math
import os
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor
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
    "expected_blocks": 0,  # 0 = auto-detect per case from the VTI files
    "workers": min(32, os.cpu_count() or 1),
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
    """Return the phi=0.5 interface in lattice y coordinates.

    Primary method: linear interpolation of the phi=0.5 crossing between
    adjacent rows. Fallback for columns with no 0.5 crossing (different phi
    convention, e.g. 0..255 or step functions): the steepest-gradient row.
    """
    y_interface = np.full(NI, np.nan, dtype=float)
    for x in range(NI):
        column = phi[:, x]
        found = None
        for row in range(NJ - 1):
            lower = column[row]
            upper = column[row + 1]
            if not (np.isfinite(lower) and np.isfinite(upper)):
                continue
            if (lower - 0.5) * (upper - 0.5) <= 0 and lower != upper:
                found = row + 0.5 + (0.5 - lower) / (upper - lower)
                break
        if found is None:
            best_row, best_slope = 0, -1.0
            for row in range(NJ - 1):
                lower = column[row]
                upper = column[row + 1]
                if np.isfinite(lower) and np.isfinite(upper):
                    slope = abs(upper - lower)
                    if slope > best_slope:
                        best_slope, best_row = slope, row
            if best_slope > 0:
                found = best_row + 0.5
        y_interface[x] = found
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


def case_block_ids(vtidata):
    """One glob per case: all output steps and the union of block ids."""
    steps = set()
    blocks = set()
    for path in vtidata.glob("rosensweig2d_T*_B*.vti"):
        step = numeric_step(path)
        if step is None:
            continue
        steps.add(step)
        block = numeric_block(path)
        if block is not None:
            blocks.add(block)
    return sorted(steps), blocks


def probe_layout(vtidata, step):
    """Read only the XML header of one block file: Origin, Extent, field names.

    Cheap layout diagnostics for a case whose extracted interface looks wrong
    (e.g. all-zero peak/valley): confirms whether blocks carry global extents
    and which arrays are actually present.
    """
    files = sorted(vtidata.glob(f"rosensweig2d_T{step}_B*.vti"))
    if not files:
        return None
    head = files[0].read_bytes()
    origin = re.search(rb'Origin="([^"]+)"', head)
    extent = re.search(rb'Extent="([^"]+)"', head)
    names = re.findall(rb'Name="([A-Za-z0-9_]+)"', head)
    return (
        files[0].name,
        origin.group(1).decode("ascii") if origin else "?",
        extent.group(1).decode("ascii") if extent else "?",
        sorted(set(name.decode("ascii") for name in names)),
    )


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


def extract_case_step(case, step, wanted_fields, expected_blocks,
                      reference_ids, reference_steps):
    vtidata = case / "vtkoutput" / "vtidata"
    files = sorted(
        path
        for path in vtidata.glob(f"rosensweig2d_T{step}_B*.vti")
        if numeric_step(path) == step
    )
    expected_count = expected_blocks or len(reference_ids)
    row = {
        "case": case.name,
        "H0_kAm": case_h0(case),
        "step": step,
        "status": "ok",
        "vti_files": len(files),
        "expected_vti_files": expected_count,
    }
    if not files:
        row["status"] = "missing_step"
        row["available_steps"] = ";".join(map(str, reference_steps))
        return row

    found_ids = {numeric_block(path) for path in files}
    row["block_ids"] = ";".join(
        str(block) for block in sorted(found_ids) if block is not None
    )

    if expected_blocks:
        expected_ids = set(range(expected_blocks))
        if found_ids != expected_ids:
            row["status"] = "block_count_mismatch"
            missing = expected_ids - found_ids
            if missing:
                row["missing_ids"] = ";".join(map(str, sorted(missing)))
            return row
    elif found_ids != reference_ids:
        missing = reference_ids - found_ids
        row["status"] = (
            "partial" if found_ids < reference_ids else "block_count_mismatch"
        )
        if missing:
            row["missing_ids"] = ";".join(map(str, sorted(missing)))
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

    # Interface extraction assumes a 0..1 order parameter. If the stored phi
    # uses another convention (e.g. 0..255, or diverged values), normalize it
    # so the 0.5 crossing lands at the midpoint of the transition ramp.
    phi_interface = phi
    if phi_finite.any():
        phi_values = phi[phi_finite]
        phi_span_values = float(np.max(phi_values) - np.min(phi_values))
        if phi_span_values > 1.5:
            phi_min = float(np.min(phi_values))
            phi_interface = (phi - phi_min) / phi_span_values
    row["phi_floor_fraction"] = (
        float(np.mean(phi[phi_finite] <= 0.010001)) if phi_finite.any() else math.nan
    )
    row["phi_cap_fraction"] = (
        float(np.mean(phi[phi_finite] >= 0.998999)) if phi_finite.any() else math.nan
    )
    if phi_finite.any():
        phi_values = phi[phi_finite]
        row["phi_span"] = float(np.max(phi_values) - np.min(phi_values))
    else:
        row["phi_span"] = math.nan
    add_stats(row, "phi_wall", phi, wall)
    add_stats(row, "phi_interior", phi, interior)
    add_stats(row, "phi_seam", phi, seam)
    add_stats(row, "phi_bottom", phi[:WALL_BAND_CELLS, :])
    add_stats(row, "phi_top", phi[-WALL_BAND_CELLS:, :])
    add_stats(row, "phi_bottom_deviation", np.abs(phi[:WALL_BAND_CELLS, :] - 1.0))
    add_stats(row, "phi_top_deviation", np.abs(phi[-WALL_BAND_CELLS:, :]))
    row["phi_seam_jump_max"] = max_abs_jump(phi)

    y_interface = extract_interface(phi_interface)
    valid_interface = y_interface[np.isfinite(y_interface)]
    row["interface_columns"] = int(valid_interface.size)
    if valid_interface.size:
        interface_mean = float(np.mean(valid_interface))
        row["interface_mean_mm"] = interface_mean * DY_MM
        row["interface_min_mm"] = float(np.min(valid_interface)) * DY_MM
        row["interface_max_mm"] = float(np.max(valid_interface)) * DY_MM
        row["peak_height_mm"] = (
            float(np.max(valid_interface)) - interface_mean
        ) * DY_MM
        row["valley_depth_mm"] = (
            interface_mean - float(np.min(valid_interface))
        ) * DY_MM
    else:
        row["interface_mean_mm"] = math.nan
        row["interface_min_mm"] = math.nan
        row["interface_max_mm"] = math.nan
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


def run_task(task):
    case, step, wanted_fields, expected_blocks, reference_ids, reference_steps = task
    try:
        return extract_case_step(case, step, wanted_fields, expected_blocks,
                                 reference_ids, reference_steps)
    except Exception as error:
        return {
            "case": case.name,
            "H0_kAm": case_h0(case),
            "step": step,
            "status": f"error:{type(error).__name__}",
            "vti_files": 0,
            "expected_vti_files": expected_blocks or len(reference_ids),
            "error": str(error)[:200],
        }


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
        help="expected B*.vti block count; 0 (default) = auto-detect per case",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=DEFAULTS["workers"],
        help="parallel extraction workers (default: min(32, ncpu))",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="print per-case block layout (Origin/Extent/array names) of one file",
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
    tasks = []
    for case in cases:
        vtidata = case / "vtkoutput" / "vtidata"
        reference_steps, reference_ids = case_block_ids(vtidata)
        if not reference_steps:
            rows.append(
                {
                    "case": case.name,
                    "H0_kAm": case_h0(case),
                    "step": "",
                    "status": "no_steps",
                    "available_steps": "",
                }
            )
            print(f"[skip] {case.name}: no output steps found")
            continue
        try:
            selected_steps = parse_step_spec(step_spec, reference_steps, args.stride)
        except ValueError as error:
            print(f"[error] {error}", file=sys.stderr)
            return 2
        if not selected_steps:
            rows.append(
                {
                    "case": case.name,
                    "H0_kAm": case_h0(case),
                    "step": "",
                    "status": "no_steps",
                    "available_steps": ";".join(map(str, reference_steps)),
                }
            )
            print(f"[skip] {case.name}: no selected steps")
            continue

        expected = args.expected_blocks or len(reference_ids)
        print(
            f"[info] {case.name}: {len(reference_ids)} blocks, "
            f"{len(reference_steps)} steps -> {len(selected_steps)} selected"
        )
        if args.debug:
            layout = probe_layout(vtidata, reference_steps[0])
            if layout is not None:
                name, origin, extent, arrays = layout
                print(
                    f"[layout] {case.name}: {name} Origin={origin} "
                    f"Extent={extent} arrays={','.join(arrays)}"
                )
            else:
                print(f"[layout] {case.name}: no VTI files at T{reference_steps[0]}")
        for step in selected_steps:
            tasks.append(
                (case, step, wanted_fields, args.expected_blocks,
                 reference_ids, reference_steps)
            )

    started = time.time()
    completed = 0
    if tasks:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            for row in pool.map(run_task, tasks):
                rows.append(row)
                completed += 1
                if row["status"] == "ok":
                    coverage = row.get("interface_columns", 0) / NI
                    warning = ""
                    if coverage < 0.8:
                        warning = (
                            f" [warn: interface in only {row['interface_columns']}"
                            f"/{NI} columns]"
                        )
                    print(
                        f"[ok] {row['case']} T{row['step']}: "
                        f"{row['vti_files']} VTI, "
                        f"peak={row.get('peak_height_mm', math.nan):.6f} mm, "
                        f"valley={row.get('valley_depth_mm', math.nan):.6f} mm, "
                        f"Hmax={row.get('hmag_max', math.nan):.6g}"
                        f"{warning}"
                    )
                else:
                    print(
                        f"[skip] {row['case']} T{row['step']}: {row['status']} "
                        f"({row.get('available_steps', '')})"
                    )

    failures = sum(1 for row in rows if row["status"] != "ok")
    output = args.out if args.out.is_absolute() else root / args.out
    output.parent.mkdir(parents=True, exist_ok=True)
    preferred = ["case", "H0_kAm", "step", "status", "vti_files"]
    all_keys = set().union(*(row.keys() for row in rows))
    fieldnames = preferred + sorted(all_keys - set(preferred))
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    print(f"saved unified diagnostic CSV: {output}")
    print(
        f"done in {time.time() - started:.1f}s: {len(rows)} rows, "
        f"{len(rows) - failures} ok, {failures} failed"
    )
    return 2 if args.strict and failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
