function antenna_gain_to_uan (model, radome, freq_name, polar, varargin)
%ANTENNA_GAIN_TO_UAN  Convert antenna gain pattern to Remcom Wireless InSite .uan format
%
% antenna_gain_to_uan (model, radome, freq_name, polar)
% antenna_gain_to_uan (model, radome, freq_name, polar, output_path)
%
% Converts antenna gain pattern from the SNR simulator's GAIN.DAT/DATX format
% to Remcom Wireless InSite's .uan (User-Defined Antenna) file format.
%
% Inputs:
%   model     - Antenna model code (e.g., 'TRM55971.00', 'TRM29659.00')
%   radome    - Radome code (e.g., 'NONE')
%   freq_name - Frequency band ('L1', 'L2', 'L5', 'R1', 'R2')
%   polar     - Polarization ('RHCP' or 'LHCP')
%   output_path - (optional) Output .uan file path; default: <model>_<freq>_<polar>.uan
%
% Example:
%   antenna_gain_to_uan('TRM55971.00', 'NONE', 'L1', 'RHCP')
%   antenna_gain_to_uan('TRM29659.00', 'NONE', 'L2', 'RHCP', 'Trimble_choke_ring_L2.uan')
%
% The .uan file can be imported into Remcom Wireless InSite for propagation
% modeling with the actual antenna radiation pattern.
%
% See also: snr_setup_ant_profile_load, snr_setup_ant_comp

    if (nargin < 5) || isempty(varargin{1})
        output_path = [];
    else
        output_path = varargin{1};
    end

    % Load antenna profile data
    filename = snr_setup_ant_filename('gain', model, radome, freq_name, polar);
    data_dir = snr_setup_ant_path();
    filepath = fullfile(data_dir, 'profile', filename);

    % Try extended format first
    if exist([filepath 'X'], 'file')
        filepath = [filepath 'X'];
    end
    if ~exist(filepath, 'file')
        error('snr:ant:uan:fileNotFound', ...
            'Antenna gain file not found: %s', filepath);
    end

    profile = snr_setup_ant_profile_load(filepath, 'gain');

    % Get frequency in Hz for UAN header
    freq_hz = get_frequency_hz(freq_name);

    % Build theta, phi, gain arrays for UAN format
    % SNR uses: elev = 90 - ang (elevation from horizon), azim = azimuth
    % UAN uses: theta = angle from zenith (0-180), phi = azimuth (0-360)
    % So: theta = 90 - elev = ang (boresight angle)
    theta_deg = profile.ang;
    phi_deg = profile.azim;
    gain_db = profile.gain_db;

    % For 2-column (azimuth-symmetric) data, expand to full theta-phi grid
    % Order: phi varies first (0,5,...,355) then theta, to match UAN convention
    if isscalar(unique(phi_deg)) || numel(unique(phi_deg)) < 3
        theta_domain = linspace(0, 180, 37);
        phi_domain = linspace(0, 355, 72);
        [phi_grid, theta_grid] = meshgrid(phi_domain, theta_domain);
        theta_grid = theta_grid(:);
        phi_grid = phi_grid(:);
        % Interpolate gain (azimuth-independent)
        gain_grid = interp1(theta_deg, gain_db, theta_grid, 'linear', 'extrap');
    else
        % Full 3D data from DATX - need to map to theta-phi
        [theta_grid, phi_grid, gain_grid] = expand_profile_to_grid(profile);
    end

    % For circular polarization (RHCP/LHCP), split gain between theta and phi
    % components. Equal power split: -3 dB each.
    gain_theta_db = gain_grid - 3;
    gain_phi_db = gain_grid - 3;
    phase_theta_deg = zeros(size(gain_grid));
    phase_phi_deg = zeros(size(gain_grid));

    % Determine angle increments for UAN header
    theta_vals = unique(theta_grid);
    phi_vals = unique(phi_grid);
    theta_inc = min(diff(sort(theta_vals)));
    if isempty(theta_inc) || theta_inc == 0, theta_inc = 5; end
    phi_inc = min(diff(sort(phi_vals)));
    if isempty(phi_inc) || phi_inc == 0, phi_inc = 5; end

    % Write UAN file
    if isempty(output_path)
        output_path = sprintf('%s_%s_%s.uan', model, freq_name, polar);
    end
    if isempty(regexpi(output_path, '\.uan$'))
        output_path = [output_path '.uan'];
    end

    write_uan_file(output_path, theta_grid, phi_grid, ...
        gain_theta_db, gain_phi_db, phase_theta_deg, phase_phi_deg, ...
        theta_inc, phi_inc);

    fprintf('Wrote: %s\n', output_path);
end

function freq_hz = get_frequency_hz(freq_name)
    switch upper(freq_name)
        case 'L1'
            freq_hz = 154 * 10.23e6;  % 1575.42 MHz
        case 'L2'
            freq_hz = 120 * 10.23e6;  % 1227.60 MHz
        case 'L5'
            freq_hz = 115 * 10.23e6;  % 1176.45 MHz
        case {'R1', 'R2'}
            freq_hz = 120 * 10.23e6;  % GLONASS ~L2
        otherwise
            freq_hz = 1575.42e6;     % Default L1
    end
end

function [theta_grid, phi_grid, gain_grid] = expand_profile_to_grid(profile)
    % Map profile (elev, azim) to UAN (theta, phi)
    % theta = 90 - elev, phi = azim (wrap to 0-360)
    theta_deg = 90 - profile.elev;
    phi_deg = mod(profile.azim, 360);

    theta_domain = linspace(0, 180, 37);
    phi_domain = linspace(0, 355, 72);
    [phi_grid, theta_grid] = meshgrid(phi_domain, theta_domain);
    theta_grid = theta_grid(:);
    phi_grid = phi_grid(:);

    F = scatteredInterpolant(theta_deg, phi_deg, profile.gain_db, ...
        'linear', 'nearest');
    gain_grid = F(theta_grid, phi_grid);
end

function write_uan_file(filepath, theta, phi, g_theta, g_phi, ...
    phase_theta, phase_phi, theta_inc, phi_inc)
    % Write Remcom Wireless InSite UAN format
    fid = fopen(filepath, 'w');
    if fid < 0
        error('snr:ant:uan:writeFailed', 'Cannot open file: %s', filepath);
    end

    maximum_gain = max(nanmax(g_theta), nanmax(g_phi));
    phi_max = 360 - phi_inc;

    % UAN header (Wireless InSite format)
    fprintf(fid, 'begin_<parameters> \n');
    fprintf(fid, 'format free\n');
    fprintf(fid, 'phi_min 0.000000\n');
    fprintf(fid, 'phi_max %.6f\n', phi_max);
    fprintf(fid, 'phi_inc %.6f\n', phi_inc);
    fprintf(fid, 'theta_min 0.000000\n');
    fprintf(fid, 'theta_max 180.0000\n');
    fprintf(fid, 'theta_inc %.6f\n', theta_inc);
    fprintf(fid, 'complex\n');
    fprintf(fid, 'mag_phase\n');
    fprintf(fid, 'pattern gain\n');
    fprintf(fid, 'magnitude dB\n');
    fprintf(fid, 'maximum_gain %.6f\n', maximum_gain);
    fprintf(fid, 'phase degrees\n');
    fprintf(fid, 'direction degrees\n');
    fprintf(fid, 'polarization theta_phi\n');
    fprintf(fid, 'end_<parameters>\n');

    % Data rows: theta, phi, G_theta, G_phi, phase_theta, phase_phi
    for i = 1:numel(theta)
        fprintf(fid, '%.6f  %.6f  %.6f  %.6f  %.6f  %.6f\n', ...
            theta(i), phi(i), g_theta(i), g_phi(i), ...
            phase_theta(i), phase_phi(i));
    end

    fclose(fid);
end
