#!/usr/bin/env python3
"""
Build PHASE.DAT for the antenna profile from an ANTEX .atx file.

Uses PCO (Phase Center Offset) + PCV (Phase Center Variation) per
lib/snr/data/ant/gemini-phase.md:
  total_correction_mm = PCO_projection(azimuth, zenith) + PCV(azimuth, zenith)
  phase_deg = 360 * (total_correction_mm / wavelength_mm)

For a single principal-plane cut (profile basic format), NOAZI (azimuth-
independent) values are used and phase is written as boresight_deg, phase_deg.

Usage:
  python atx_phase_to_profile.py TRM598.atx TRM59800.00 NONE L1
  python atx_phase_to_profile.py TRM598.atx TRM59800.00 NONE L1 -o lib/snr/data/ant/profile

Output: TRM59800.00__NONE__L1__RHCP__PHASE.DAT (and optionally LHCP) in profile dir.
"""

import argparse
import math
import os
import re
import sys

# GPS L1 wavelength in mm (c / 1575.42e6 * 1000)
WAVELENGTH_MM = {
    'L1': 299792458.0 / 1575.42e6 * 1000,
    'L2': 299792458.0 / 1227.60e6 * 1000,
    'L5': 299792458.0 / 1176.45e6 * 1000,
}


def parse_atx_frequency_block(lines, start_idx):
    """
    Parse one frequency block (from START OF FREQUENCY to END OF FREQUENCY).
    Returns (pco_north_mm, pco_east_mm, pco_up_mm, zen_deg_list, azi_deg_list, pcv_mm_grid)
    where pcv_mm_grid[azi_idx][zen_idx] = PCV in mm.
    """
    i = start_idx
    pco_north = pco_east = pco_up = None
    zen_deg = []
    azi_deg = []
    grid = []  # list of rows: each row is [pcv at zen 0, 5, ..., 90]

    while i < len(lines):
        line = lines[i]
        if 'END OF FREQUENCY' in line:
            break
        # NORTH / EAST / UP line: "      0.61      1.87     86.60                              NORTH / EAST / UP"
        if 'NORTH / EAST / UP' in line:
            parts = line.split()
            if len(parts) >= 3:
                pco_north = float(parts[0])
                pco_east = float(parts[1])
                pco_up = float(parts[2])
            i += 1
            continue
        # Data row: "   NOAZI    0.00   -0.38  ..." or "     0.0    0.00   -0.37  ..."
        if pco_north is not None and line.strip():
            parts = line.split()
            if len(parts) >= 2:
                try:
                    azi_label = parts[0]
                    values = [float(x) for x in parts[1:]]
                    if azi_label == 'NOAZI':
                        zen_deg = list(range(0, min(91, len(values) * 5), 5))[:len(values)]
                        if not zen_deg and values:
                            zen_deg = list(range(0, 91, 5))[:len(values)]
                    azi_val = 0.0 if azi_label == 'NOAZI' else float(azi_label)
                    if values:
                        azi_deg.append(azi_val)
                        grid.append(values)
                except (ValueError, IndexError):
                    pass
        i += 1

    return pco_north, pco_east, pco_up, zen_deg, azi_deg, grid


def find_antenna_and_freq(path, antenna_type, freq_code):
    """
    Find antenna block and frequency block in ATX file.
    antenna_type e.g. 'TRM59800.00', freq_code e.g. 'G01' for L1.
    Returns (pco_north, pco_east, pco_up, zen_deg, azi_deg, pcv_grid) or None.
    """
    with open(path, 'r') as f:
        lines = f.readlines()

    in_antenna = False
    ant_match = re.compile(r'^\s*' + re.escape(antenna_type.strip()) + r'\s', re.IGNORECASE)

    for i, line in enumerate(lines):
        if 'START OF ANTENNA' in line:
            in_antenna = False
        if ant_match.search(line):
            in_antenna = True
        if in_antenna and 'START OF FREQUENCY' in line:
            # Frequency code is on same line: "   G01     START OF FREQUENCY"
            if freq_code.upper() in line.split():
                result = parse_atx_frequency_block(lines, i + 1)
                if result[0] is not None:
                    return result
    return None


def freq_name_to_atx_code(freq_name):
    """L1 -> G01, L2 -> G02, L5 -> G05."""
    m = {'L1': 'G01', 'L2': 'G02', 'L5': 'G05', 'R1': 'R01', 'R2': 'R02'}
    return m.get(freq_name.upper(), 'G01')


def pco_projection_mm(north_mm, east_mm, up_mm, azimuth_deg, zenith_deg):
    """
    Project PCO onto signal direction (unit vector toward satellite).
    azimuth_deg: clockwise from North; zenith_deg: from vertical (0=up, 90=horizon).
    Per gemini-phase.md: pco_proj = north*cos(az)*sin(zen) + east*sin(az)*sin(zen) + up*cos(zen).
    """
    az_rad = math.radians(azimuth_deg)
    zen_rad = math.radians(zenith_deg)
    return (north_mm * math.cos(az_rad) * math.sin(zen_rad) +
            east_mm * math.sin(az_rad) * math.sin(zen_rad) +
            up_mm * math.cos(zen_rad))


def interpolate_pcv(zen_deg_list, azi_deg_list, grid, azimuth_deg, zenith_deg):
    """Bilinear interpolation of PCV grid. Use NOAZI row if only one row."""
    if not grid or not zen_deg_list:
        return 0.0
    # If we only have NOAZI (one row), use that row for all azimuths
    if len(grid) == 1:
        row = grid[0]
    else:
        # Row order matches azi_deg_list; find bracketing azimuth indices
        idx_hi = None
        for k, a in enumerate(azi_deg_list):
            if a >= azimuth_deg:
                idx_hi = k
                break
        if idx_hi is None:
            row = grid[-1]
        elif idx_hi == 0:
            row = grid[0]
        else:
            w = (azimuth_deg - azi_deg_list[idx_hi - 1]) / (azi_deg_list[idx_hi] - azi_deg_list[idx_hi - 1])
            row = [(1 - w) * grid[idx_hi - 1][j] + w * grid[idx_hi][j] for j in range(len(grid[0]))]

    # Interpolate in zenith
    nzen = len(zen_deg_list)
    if zenith_deg <= zen_deg_list[0]:
        return row[0]
    if zenith_deg >= zen_deg_list[-1]:
        return row[-1]
    for j in range(nzen - 1):
        if zen_deg_list[j + 1] >= zenith_deg:
            w = (zenith_deg - zen_deg_list[j]) / (zen_deg_list[j + 1] - zen_deg_list[j])
            return (1 - w) * row[j] + w * row[j + 1]
    return row[-1]


def total_phase_correction_mm(pco_north, pco_east, pco_up, zen_list, azi_list, grid, azimuth_deg, zenith_deg):
    """PCO projection + interpolated PCV (gemini-phase.md)."""
    pco_proj = pco_projection_mm(pco_north, pco_east, pco_up, azimuth_deg, zenith_deg)
    pcv = interpolate_pcv(zen_list, azi_list, grid, azimuth_deg, zenith_deg)
    return pco_proj + pcv


def mm_to_phase_deg(total_mm, wavelength_mm):
    """Convert path-length correction (mm) to phase in degrees."""
    return 360.0 * (total_mm / wavelength_mm)


def write_phase_dat(filepath, boresight_deg, phase_deg):
    """
    Write profile-format PHASE.DAT: first line NaN and delta (0), then boresight_deg, phase_deg.
    """
    with open(filepath, 'w') as f:
        f.write('NaN\t0\n')
        for b, p in zip(boresight_deg, phase_deg):
            f.write(f'{b}\t{p}\n')


def main():
    parser = argparse.ArgumentParser(
        description='Build PHASE.DAT from ANTEX .atx (PCO + PCV) for antenna profile'
    )
    parser.add_argument('atx_path', help='Path to .atx file (e.g. TRM598.atx)')
    parser.add_argument('model', help='Antenna model (e.g. TRM59800.00)')
    parser.add_argument('radome', default='NONE', nargs='?', help='Radome code')
    parser.add_argument('freq', default='L1', nargs='?', help='Frequency (L1, L2, L5)')
    parser.add_argument('-o', '--output-dir', default=None,
                        help='Profile directory; default: lib/snr/data/ant/profile')
    parser.add_argument('--polar', default='RHCP', choices=['RHCP', 'LHCP'],
                        help='Polarization suffix for filename (default RHCP)')
    args = parser.parse_args()

    atx_path = os.path.abspath(args.atx_path)
    if not os.path.isfile(atx_path):
        print(f'Error: ATX file not found: {atx_path}', file=sys.stderr)
        sys.exit(1)

    freq_code = freq_name_to_atx_code(args.freq)
    result = find_antenna_and_freq(atx_path, args.model.strip(), freq_code)
    if result is None:
        print(f'Error: Antenna {args.model} / frequency {freq_code} not found in {atx_path}', file=sys.stderr)
        sys.exit(1)

    pco_north, pco_east, pco_up, zen_list, azi_list, grid = result
    if not zen_list:
        zen_list = list(range(0, 91, 5))

    wavelength_mm = WAVELENGTH_MM.get(args.freq.upper(), WAVELENGTH_MM['L1'])

    # Principal-plane cut: azimuth = 0 (North). Boresight 0 = zenith, 180 = nadir.
    # For boresight > 90 we mirror zenith: zenith = 180 - boresight (nadir = back direction).
    boresight_points = list(range(0, 181, 5))
    phase_deg_list = []
    for boresight in boresight_points:
        zenith = boresight if boresight <= 90 else (180 - boresight)
        total_mm = total_phase_correction_mm(
            pco_north, pco_east, pco_up, zen_list, azi_list, grid,
            azimuth_deg=0.0, zenith_deg=zenith
        )
        phase_deg_list.append(mm_to_phase_deg(total_mm, wavelength_mm))

    if args.output_dir is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        args.output_dir = os.path.join(script_dir, '..', 'lib', 'snr', 'data', 'ant', 'profile')
    out_dir = os.path.normpath(args.output_dir)
    os.makedirs(out_dir, exist_ok=True)

    phase_filename = f'{args.model}__{args.radome}__{args.freq}__{args.polar}__PHASE.DAT'
    phase_path = os.path.join(out_dir, phase_filename)
    write_phase_dat(phase_path, boresight_points, phase_deg_list)
    print(f'Wrote: {phase_path}')
    if args.polar == 'RHCP':
        # Optionally write LHCP same file (phase pattern is same; tool applies ±90°)
        lhcp_path = os.path.join(out_dir, f'{args.model}__{args.radome}__{args.freq}__LHCP__PHASE.DAT')
        write_phase_dat(lhcp_path, boresight_points, phase_deg_list)
        print(f'Wrote: {lhcp_path}')


if __name__ == '__main__':
    main()
