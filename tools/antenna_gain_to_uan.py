#!/usr/bin/env python3
"""
Convert antenna gain pattern from SNR simulator GAIN.DAT/DATX format
to Remcom Wireless InSite .uan (User-Defined Antenna) format.

Usage:
    python antenna_gain_to_uan.py TRM55971.00 NONE L1 RHCP
    python antenna_gain_to_uan.py TRM29659.00 NONE L2 RHCP -o Trimble_L2.uan

The antenna gain data files are in lib/snr/data/ant/profile/
Format: MODEL__RADOME__FREQ__POLAR__GAIN.DAT or .DATX
"""

import argparse
import os
import sys

# GPS/GNSS frequencies in Hz
FREQUENCIES = {
    'L1': 154 * 10.23e6,   # 1575.42 MHz
    'L2': 120 * 10.23e6,   # 1227.60 MHz
    'L5': 115 * 10.23e6,   # 1176.45 MHz
    'R1': 120 * 10.23e6,   # GLONASS
    'R2': 120 * 10.23e6,
}


def load_gain_profile(filepath):
    """
    Load antenna gain from .DAT or .DATX file.
    Returns (theta_deg, phi_deg, gain_db, delta).
    theta = angle from zenith (0-180), phi = azimuth (0-360)
    """
    import numpy as np

    data = np.loadtxt(filepath)
    if data.ndim == 1:
        data = data.reshape(1, -1)

    ncols = data.shape[1]
    if ncols == 2:
        # Basic format: boresight_angle, gain
        data = np.column_stack([data[:, 0], np.zeros(len(data)), data[:, 1]])

    # First row: NaN, NaN/azi, delta
    delta = float(data[0, 2]) if not np.isnan(data[0, 2]) else 0.0
    data = data[1:]

    ang = data[:, 0]  # boresight angle (0=zenith, 180=nadir, 360=wrap)
    azi = data[:, 1]
    val = data[:, 2]  # pseudo-gain
    gain_db = val - delta  # true gain

    # Convert to theta (0-180) and phi (0-360) per angle_boresight_positive
    theta_deg = np.where(ang < 0, np.abs(ang), ang)
    theta_deg = np.where(theta_deg > 180, 360 - theta_deg, theta_deg)
    phi_deg = np.where((ang < 0) | (ang > 180), azi + 180, azi)
    phi_deg = np.mod(phi_deg, 360)

    return theta_deg, phi_deg, gain_db, delta


def write_uan(filepath, theta, phi, gain_db, phase_deg, freq_hz, polarization='RHCP'):
    """
    Write Remcom Wireless InSite UAN format.
    For circular polarization (RHCP/LHCP), split gain -3dB between theta and phi.
    Phase: use phase_deg if provided (from PHASE.DAT), else derive from polarization.
    For Remcom/IEEE: RHCP phase_phi - phase_theta = -90°; LHCP = +90°.
    """
    import numpy as np

    # For circular pol: equal power in theta and phi components
    gain_theta_db = gain_db - 3.0
    gain_phi_db = gain_db - 3.0

    # Phase: theta-phi relationship for circular polarization.
    # IEEE/Remcom: RHCP = E_theta - j*E_phi => phase_phi lags theta by 90° (phase_phi - phase_theta = -90°).
    # So RHCP: phase_phi = phase_theta - 90°; LHCP: phase_phi = phase_theta + 90°.
    if phase_deg is not None:
        # Use measured phase for both components (phase pattern from file)
        phase_theta = np.asarray(phase_deg, dtype=float)
        phase_phi = phase_theta + (-90.0 if polarization.upper() == 'RHCP' else 90.0)
    else:
        # No phase data: derive from polarization convention
        phase_theta = np.zeros_like(gain_db)
        phase_phi = np.full_like(gain_db, -90.0 if polarization.upper() == 'RHCP' else 90.0)

    theta_vals = np.unique(theta)
    phi_vals = np.unique(phi)
    theta_inc = np.min(np.diff(np.sort(theta_vals))) if len(theta_vals) > 1 else 5.0
    phi_inc = np.min(np.diff(np.sort(phi_vals))) if len(phi_vals) > 1 else 5.0
    if np.isnan(theta_inc) or theta_inc <= 0:
        theta_inc = 5.0
    if np.isnan(phi_inc) or phi_inc <= 0:
        phi_inc = 5.0

    # maximum_gain: peak gain in dBi (max of theta/phi components)
    maximum_gain = float(np.nanmax(np.maximum(gain_theta_db, gain_phi_db)))

    # phi_max: 355 when phi_inc=5 (0,5,...,355) to match Wireless InSite convention
    phi_max = 360 - phi_inc if phi_inc > 0 else 360

    with open(filepath, 'w') as f:
        f.write('begin_<parameters> \n')
        f.write('format free\n')
        f.write(f'phi_min 0.000000\n')
        f.write(f'phi_max {phi_max:.6f}\n')
        f.write(f'phi_inc {phi_inc:.6f}\n')
        f.write('theta_min 0.000000\n')
        f.write('theta_max 180.0000\n')
        f.write(f'theta_inc {theta_inc:.6f}\n')
        f.write('complex\n')
        f.write('mag_phase\n')
        f.write('pattern gain\n')
        f.write('magnitude dB\n')
        f.write(f'maximum_gain {maximum_gain:.6f}\n')
        f.write('phase degrees\n')
        f.write('direction degrees\n')
        f.write('polarization theta_phi\n')
        f.write('end_<parameters>\n')
        for i in range(len(theta)):
            f.write(f'{theta[i]:.6f}  {phi[i]:.6f}  {gain_theta_db[i]:.6f}  '
                    f'{gain_phi_db[i]:.6f}  {phase_theta[i]:.6f}  {phase_phi[i]:.6f}\n')


def expand_azimuth_symmetric(theta_deg, gain_db, theta_step=5, phi_step=5):
    """Expand 1D (azimuth-symmetric) pattern to full theta-phi grid.
    Order: phi varies first (0,5,...,355) then theta, to match UAN convention.
    """
    import numpy as np

    theta_grid = np.arange(0, 181, theta_step)
    # phi 0,5,...,355 to match Wireless InSite UAN phi_max convention
    phi_grid = np.arange(0, 360 - phi_step + 1e-9, phi_step)
    phi_mesh, theta_mesh = np.meshgrid(phi_grid, theta_grid)
    theta_flat = theta_mesh.ravel()
    phi_flat = phi_mesh.ravel()

    # Interpolate gain (constant in phi for azimuth-symmetric)
    gain_flat = np.interp(theta_flat, theta_deg, gain_db)
    return theta_flat, phi_flat, gain_flat


def main():
    import numpy as np

    parser = argparse.ArgumentParser(
        description='Convert antenna gain to Remcom Wireless InSite .uan format'
    )
    parser.add_argument('model', help='Antenna model (e.g., TRM55971.00)')
    parser.add_argument('radome', help='Radome code (e.g., NONE)')
    parser.add_argument('freq', help='Frequency band (L1, L2, L5)')
    parser.add_argument('polar', choices=['RHCP', 'LHCP'], help='Polarization')
    parser.add_argument('-o', '--output', help='Output .uan file path')
    parser.add_argument('--data-dir', default=None,
                        help='Path to ant/profile directory')
    args = parser.parse_args()

    # Find data directory
    if args.data_dir:
        data_dir = args.data_dir
    else:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        data_dir = os.path.join(script_dir, '..', 'lib', 'snr', 'data', 'ant', 'profile')
    data_dir = os.path.normpath(data_dir)

    gain_filename = f'{args.model}__{args.radome}__{args.freq}__{args.polar}__GAIN'
    gain_filepath = os.path.join(data_dir, gain_filename + '.DAT')
    if not os.path.exists(gain_filepath):
        gain_filepath = os.path.join(data_dir, gain_filename + '.DATX')
    if not os.path.exists(gain_filepath):
        print(f'Error: Antenna file not found: {gain_filename}.DAT or .DATX', file=sys.stderr)
        print(f'  Looked in: {data_dir}', file=sys.stderr)
        sys.exit(1)

    theta_deg, phi_deg, gain_db, delta = load_gain_profile(gain_filepath)

    # Load phase if available (optional - most antennas have gain only)
    phase_deg = None
    phase_filename = f'{args.model}__{args.radome}__{args.freq}__{args.polar}__PHASE'
    phase_filepath = os.path.join(data_dir, phase_filename + '.DAT')
    if not os.path.exists(phase_filepath):
        phase_filepath = os.path.join(data_dir, phase_filename + '.DATX')
    if os.path.exists(phase_filepath):
        theta_phase, _, phase_vals, _ = load_gain_profile(phase_filepath)  # same format as gain
        phase_deg = phase_vals

    # Check if azimuth-symmetric (single phi value or all same)
    phi_unique = len(set(phi_deg.round(2)))
    if phi_unique < 3:
        # Deduplicate: for symmetric pattern, take unique theta and mean gain
        theta_unique, idx = np.unique(theta_deg, return_inverse=True)
        gain_unique = np.array([np.mean(gain_db[idx == i]) for i in range(len(theta_unique))])
        theta_deg, phi_deg, gain_db = expand_azimuth_symmetric(theta_unique, gain_unique)
        if phase_deg is not None:
            # Phase may be on a different theta grid (e.g. 0:5:180); interpolate onto theta_unique
            theta_phase_u = np.unique(theta_phase)
            phase_at_unique = np.array([np.mean(phase_deg[theta_phase == t]) for t in theta_phase_u])
            phase_unique = np.interp(theta_unique, theta_phase_u, phase_at_unique)
            _, _, phase_deg = expand_azimuth_symmetric(theta_unique, phase_unique)

    freq_hz = FREQUENCIES.get(args.freq.upper(), FREQUENCIES['L1'])

    if args.output:
        outpath = args.output
    else:
        outpath = f'{args.model}_{args.freq}_{args.polar}.uan'
    if not outpath.endswith('.uan'):
        outpath += '.uan'

    write_uan(outpath, theta_deg, phi_deg, gain_db, phase_deg, freq_hz, args.polar)
    print(f'Wrote: {outpath}')


if __name__ == '__main__':
    main()
