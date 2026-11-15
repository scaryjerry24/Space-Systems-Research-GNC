close all; clear; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ORBIT FOV CRATERS SCRIPT
% Builds on Alec Delgado's least_squares_main.m and Alan Figueroa's orbit
% (Orbit_and_ErrorAnalysis.m). Put this file in Alec's folder and run it.
%
% WHAT IT DOES
%  At every true anomaly in nu_list:
%   1. Spacecraft position from the orbit, and the sun angle at that spot
%      (Sun direction fixed; the sun angle changes as the spacecraft moves).
%   2. Alec's angular_error_calc.m for that distance and sun angle -> which
%      craters are detected and the error stats for each.
%   3. Alec's crater_pos.m -> places those craters inside the view and
%      gives sigma_el, sigma_az and noisy az/el for each.
%   4. Drops craters on the night side. LS / WLS / BCWLS position fix.
%  So the craters CHANGE at every position: each view gets its own set,
%  straight from Alec's code. Crater IDs are unique across the whole run
%  ("p3_c12" = position 3, crater 12).
%
% FIELD OF VIEW
%  The view is the area where Alec's table has detections (craters out to
%  the largest bearing bin with data). Optional: cam_half_angle limits it
%  further to a camera cone pointing straight down.
%
% OUTPUTS
%  - Figure per position: Moon, orbit, view (cyan), LOS lines to that
%    position's craters, colored by sigma_bearing.
%  - Overview figure: every crater from every position on one Moon,
%    colored by position, so you can see the craters change with the view.
%  - fov_craters_all.csv     : one row per crater per position
%  - fov_summary.csv         : one row per position
%  - fov_craters_results.xlsx: both tables as sheets
%
% USES (Alec's files, unchanged): angular_error_calc.m,
% discretized_results.mat, crater_pos.m, LS_estimate.m, WLS_estimate.m,
% BCWLS_estimate.m
%
% Author: Jeremiah Symington
% Date: 2026-11-08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INPUT
mu_moon   = 4.9048695e12;          % (m^3/s^2) (Alan)
orbit_a   = 1.74e6 + 8260e3;       % (m) 10,000 km from center (about 8,260 km altitude)
orbit_e   = 0;
orbit_i   = deg2rad(30);           % (rad)
orbit_LAN = deg2rad(20);           % (rad)
orbit_w   = deg2rad(0);            % (rad)

nu_list   = 5:20:105;              % (deg) true anomalies (1 deg = about 4.2 min at this orbit)
sun_dir_I = [1; 0; 0];             % Sun direction from the Moon, fixed
cam_half_angle = [];               % (deg) optional camera cone, straight down. [] = use Alec's table only
make_figures   = true;
rng(1);                            % same results every run

%% Constants
radius_Moon = 1.7374e6;            % (m)
bin_w   = 15;                      % (deg) bearing bin width in Alec's table
sun_hat = sun_dir_I/norm(sun_dir_I);
n_pos   = length(nu_list);
load('discretized_results.mat', 'discretize_moon_angle', 'discretize_range')   % [A] bin edges

all_rows = {};  summary_rows = [];  res = cell(1,n_pos);

fprintf("------------------------- FOV Craters Test -------------------------\n")

%% Loop over true anomaly
for p = 1:n_pos

    % --- 1. Position and sun angle ---
    r_sc_I = orbit_state(orbit_a, orbit_e, orbit_i, orbit_LAN, orbit_w, deg2rad(nu_list(p)), mu_moon);
    d_km   = norm(r_sc_I)/1000;
    sun_angle = acosd(dot(r_sc_I/norm(r_sc_I), sun_hat));          % [O] Alec's moon_angle input
    sun_bin   = discretize_moon_angle(find(discretize_moon_angle <= sun_angle, 1, 'last'));
    range_bin = discretize_range(find(discretize_range <= d_km, 1, 'last'));

    % --- 2. Alec's crater table for this distance and sun angle ---
    [~,repeat_matrix_detections] = angular_error_calc(d_km, min(sun_angle,179.9));   % [A]
    moon_ang_area  = 2*rad2deg(asin(1737.4/d_km));                                   % [A] LSM L90
    sc_bearings    = repeat_matrix_detections(:,3);                                  % (deg) bin edges
    std_dev_norm   = repeat_matrix_detections(:,2);
    crater_radius  = repeat_matrix_detections(:,4);
    std_dev_unnorm = deg2rad(std_dev_norm*moon_ang_area/100);                         % [A] LSM L91
    keep = ~isnan(std_dev_unnorm) & std_dev_unnorm > 0;                               % [A] LSM L94-105
    sc_bearings = sc_bearings(keep);  crater_radius = crater_radius(keep);
    std_dev_unnorm = std_dev_unnorm(keep);  std_dev_norm = std_dev_norm(keep);
    bearing_bin = sc_bearings;
    sc_bearings = min(sc_bearings + bin_w*rand(size(sc_bearings)), 89.9);             % spread inside 15 deg bin
    bearing_max = max(bearing_bin) + bin_w;                                           % edge of the view

    % --- 3. Alec's crater_pos: places craters, sigmas, noisy az/el ---
    if isempty(sc_bearings)
        idx = [];
    else
        [std_el, std_az, az_meas, el_meas, r_c_LOS, R_LOS_I] = ...
            crater_pos(r_sc_I, sc_bearings, std_dev_unnorm);                           % [A]
        std_el = std_el(:);  std_az = std_az(:);  az_meas = az_meas(:)';  el_meas = el_meas(:)';
        r_crater_I = r_sc_I + R_LOS_I*r_c_LOS;                                         % (m) on the Moon

        % --- 4. Keep lit craters (and inside the camera cone, if set) ---
        n_c       = r_crater_I./vecnorm(r_crater_I);
        solar_inc = acosd(sun_hat'*n_c)';
        off_nadir = acosd((r_c_LOS'*[0;1;0])./vecnorm(r_c_LOS)');                      % angle from straight down (+y = Moon center)
        idx = find(solar_inc < 90);
        if ~isempty(cam_half_angle), idx = idx(off_nadir(idx) <= cam_half_angle); end
    end
    n_det = length(idx);

    % --- Angular errors and position fix ---
    errs = NaN(1,3);
    if n_det > 0
        az_nom = atan2(r_c_LOS(2,idx), r_c_LOS(1,idx));
        el_nom = atan2(r_c_LOS(3,idx), vecnorm(r_c_LOS(1:2,idx)));
        az_err = rad2deg(angle(exp(1i*(az_meas(idx) - az_nom))));                    % measured - true
        el_err = rad2deg(el_meas(idx) - el_nom);
        if n_det >= 3
            [pos_LS,errs(1),A_mat,z_mat] = LS_estimate(az_meas(idx), el_meas(idx), r_c_LOS(:,idx));
            [pos_WLS,errs(2),W_block,d_o] = WLS_estimate(std_az(idx), std_el(idx), az_meas(idx), el_meas(idx), ...
                                                         r_c_LOS(:,idx), pos_LS, A_mat, z_mat);
            [~,errs(3)] = BCWLS_estimate(std_az(idx), std_el(idx), az_meas(idx), el_meas(idx), pos_WLS, A_mat, W_block, d_o);
        end
        lat = asind(r_crater_I(3,idx)'/radius_Moon);  lon = atan2d(r_crater_I(2,idx)', r_crater_I(1,idx)');
        for k = 1:n_det
            j = idx(k);
            all_rows(end+1,:) = {nu_list(p), round(sun_angle,2), sun_bin, range_bin, ...
                sprintf('p%d_c%d', p, k), round(lat(k),3), round(lon(k),3), crater_radius(j), ...
                round(sc_bearings(j),2), bearing_bin(j), round(norm(r_c_LOS(:,j))/1000,1), ...
                round(std_dev_norm(j),5), round(rad2deg(std_dev_unnorm(j)),6), ...
                round(rad2deg(std_el(j)),6), round(rad2deg(std_az(j)),6), ...
                round(az_err(k),6), round(el_err(k),6), round(solar_inc(j),2)}; %#ok<SAGROW>
        end
        res{p} = struct('r_sc_I',r_sc_I, 'r_crater_I',r_crater_I(:,idx), ...
            'sigma_b',rad2deg(std_dev_unnorm(idx)), 'sun_angle',sun_angle, 'bearing_max',bearing_max);
    else
        res{p} = struct('r_sc_I',r_sc_I, 'r_crater_I',zeros(3,0), 'sigma_b',[], ...
            'sun_angle',sun_angle, 'bearing_max',bearing_max);
    end

    summary_rows = [summary_rows; nu_list(p), round(sun_angle,1), sun_bin, n_det, ...
        round(median(rad2deg(std_dev_unnorm(idx))),5), round(errs/1000,3)]; %#ok<AGROW>
    % Trace: shows Alec's code is rerun and which slice of the .mat it used
    fprintf("nu = %5.1f | sun angle %5.1f -> sun bin %3d, range bin %5d | table rows %3d, bearings %2d-%2d deg | median sigma_b %.5f deg | %3d lit craters | WLS %.2f km\n", ...
        nu_list(p), sun_angle, sun_bin, range_bin, length(bearing_bin), min([bearing_bin; NaN]), ...
        max([bearing_bin; NaN]) + bin_w, median(rad2deg(std_dev_unnorm)), n_det, errs(2)/1000)
end

%% Tables
All_craters = cell2table(all_rows, 'VariableNames', {'nu_deg','sun_angle_deg','sun_angle_bin_deg', ...
    'range_bin_km','crater_id','lat_deg','lon_deg','radius_bin','bearing_deg','bearing_bin_deg', ...
    'range_km','std_norm','sigma_bearing_deg','sigma_el_deg','sigma_az_deg', ...
    'az_error_deg','el_error_deg','solar_inc_deg'});
Summary = array2table(summary_rows, 'VariableNames', {'nu_deg','sun_angle_deg','sun_angle_bin_deg', ...
    'craters_detected','median_sigma_bearing_deg','LS_error_km','WLS_error_km','BCWLS_error_km'});

fprintf('\n=== Summary (one row per position) ===\n');  disp(Summary)

writetable(All_craters, 'fov_craters_all.csv');
writetable(Summary,     'fov_summary.csv');
out_file = 'fov_craters_results.xlsx';
if exist(out_file,'file'), delete(out_file); end
writetable(Summary,     out_file, 'Sheet', 'Summary');
writetable(All_craters, out_file, 'Sheet', 'All_craters');
fprintf('Tables saved to fov_craters_all.csv, fov_summary.csv and %s\n', out_file);

%% Figures
if make_figures
    km = 1/1000;
    nu_arc = linspace(0, 2*pi, 360);
    r_arc = zeros(3, length(nu_arc));
    for k = 1:length(nu_arc)
        r_arc(:,k) = orbit_state(orbit_a, orbit_e, orbit_i, orbit_LAN, orbit_w, nu_arc(k), mu_moon);
    end
    r_pos = cell2mat(cellfun(@(s) s.r_sc_I, res, 'UniformOutput', false));
    % Same camera as orbit_positions_main.m: looking down on the point
    % under the spacecraft at nu = 60 deg, from the side and slightly behind
    [r_v, v_v] = orbit_state(orbit_a, orbit_e, orbit_i, orbit_LAN, orbit_w, deg2rad(60), mu_moon);
    up    = r_v/norm(r_v);
    along = v_v - (v_v'*up)*up;  along = along/norm(along);
    side  = cross(up, along);
    cam   = 0.75*up + 0.55*side - 0.35*along;  cam = cam/norm(cam);
    all_sig = cell2mat(cellfun(@(s) s.sigma_b(:), res, 'UniformOutput', false)');
    if isempty(all_sig), c_lim = [0 1];
    else, c_lim = [min(all_sig) prctile(all_sig,95)+1e-9]; end                % top 5% don't wash out colors
    [xs, ys, zs] = sphere(120);

    % --- One figure per position ---
    for p = 1:n_pos
        figure('Color','k', 'Position', [100 100 900 750]); hold on
        moon_shaded(xs, ys, zs, radius_Moon*km, sun_hat)
        plot3(r_arc(1,:)*km, r_arc(2,:)*km, r_arc(3,:)*km, 'w-')
        plot3(r_pos(1,:)*km, r_pos(2,:)*km, r_pos(3,:)*km, 'w.', 'MarkerSize', 12)
        rs = res{p}.r_sc_I*km;
        ring = footprint_ring(res{p}.r_sc_I, res{p}.bearing_max, radius_Moon)*km;
        if ~isempty(ring), fill3(ring(1,:), ring(2,:), ring(3,:), 'c', 'FaceAlpha', 0.25, 'EdgeColor', 'c'); end
        rc = res{p}.r_crater_I*1.002*km;
        for j = 1:size(rc,2)
            plot3([rs(1) rc(1,j)], [rs(2) rc(2,j)], [rs(3) rc(3,j)], 'g-', 'LineWidth', 0.5)
        end
        if ~isempty(rc)
            scatter3(rc(1,:), rc(2,:), rc(3,:), 20, res{p}.sigma_b, 'filled', 'MarkerEdgeColor', 'k')
            text(rc(1,:), rc(2,:), rc(3,:), arrayfun(@(k) sprintf(' %d', k), 1:size(rc,2), ...
                'UniformOutput', false), 'Color', 'y', 'FontSize', 7)   % k = crater_id p<p>_c<k> in the CSV
        end
        plot3(rs(1), rs(2), rs(3), 'p', 'Color', 'b', 'MarkerFaceColor', 'c', 'MarkerSize', 14)
        caxis(c_lim);  axis equal off
        camproj('perspective'); ctr = 0.5*mean(r_pos,2)*km;
        camtarget(ctr'); campos((ctr + 40000*cam)'); camup(up'); camva(30)
        title(sprintf('(%d)  \\nu = %g°, sun angle %.0f°, %d detected', p, nu_list(p), ...
            res{p}.sun_angle, size(rc,2)), 'Color', 'w', 'FontWeight', 'normal', ...
            'FontName', 'Times New Roman', 'FontSize', 12)
        cb = colorbar; cb.Label.String = '\sigma_{bearing} (deg)'; cb.Color = 'w';
    end

    % --- Overview: every position's craters on one Moon ---
    figure('Color','k', 'Position', [100 100 900 750]); hold on
    moon_shaded(xs, ys, zs, radius_Moon*km, sun_hat)
    plot3(r_arc(1,:)*km, r_arc(2,:)*km, r_arc(3,:)*km, 'w-')
    colors = lines(n_pos);  h = gobjects(1, n_pos);
    for p = 1:n_pos
        rs = res{p}.r_sc_I*km;  rc = res{p}.r_crater_I*1.002*km;
        h(p) = plot3(rs(1), rs(2), rs(3), 'p', 'Color', colors(p,:), 'MarkerFaceColor', colors(p,:), 'MarkerSize', 12);
        if ~isempty(rc)
            plot3(rc(1,:), rc(2,:), rc(3,:), '.', 'Color', colors(p,:), 'MarkerSize', 10)
        end
    end
    axis equal off
    camproj('perspective'); ctr = 0.5*mean(r_pos,2)*km;
    camtarget(ctr'); campos((ctr + 40000*cam)'); camup(up'); camva(30)
    legend(h, arrayfun(@(v) sprintf('\\nu = %g°', v), nu_list, 'UniformOutput', false), ...
        'TextColor', 'w', 'Color', 'k', 'Location', 'eastoutside')
    title('Craters seen from each position (color = position)', 'Color', 'w', ...
        'FontWeight', 'normal', 'FontName', 'Times New Roman', 'FontSize', 12)
end

%% Local functions

function [r_I, v_I] = orbit_state(a, e, inc, LAN, w, nu, mu)
% Position/velocity at true anomaly nu (same rotation as Alan's code)
p_semi = a*(1 - e^2);
r_mag  = p_semi/(1 + e*cos(nu));
R1 = @(th)[cos(th) sin(th) 0; -sin(th) cos(th) 0; 0 0 1];
R2 = @(th)[1 0 0; 0 cos(th) sin(th); 0 -sin(th) cos(th)];
R_total = R1(-LAN)*R2(-inc)*R1(-w);
r_I = R_total*[r_mag*cos(nu); r_mag*sin(nu); 0];
v_I = R_total*sqrt(mu/p_semi)*[-sin(nu); e + cos(nu); 0];
end

function ring = footprint_ring(r_sc, bearing_max, R)
% Edge of the view on the surface (3xN, m), crater_pos.m geometry
ring = [];
if isempty(bearing_max) || isnan(bearing_max), return, end
d = norm(r_sc);  b = deg2rad(min(bearing_max, 90));
gamma = b - asin(R*sin(b)/d);                          % central angle to the edge
u  = r_sc/d;
e1 = cross([0;0;1], u);  e1 = e1/norm(e1);  e2 = cross(u, e1);
t  = linspace(0, 2*pi, 120);
ring = 1.002*R*(cos(gamma)*u + sin(gamma)*(e1*cos(t) + e2*sin(t)));
end

function moon_shaded(xs, ys, zs, R, sun_hat)
% Moon shaded straight from the Sun direction (no MATLAB lights, so the
% lit side always faces the Sun). Brightness = 0.15 + 0.75*cos(sun angle).
shade = 0.15 + 0.75*max(0, xs*sun_hat(1) + ys*sun_hat(2) + zs*sun_hat(3));
C = repmat(shade, 1, 1, 3);
surf(R*xs, R*ys, R*zs, C, 'EdgeColor', 'none', 'FaceColor', 'interp', 'FaceLighting', 'none')
end