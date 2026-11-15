close all; clear; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ORBIT POSITIONS SCRIPT
% Builds on Alec Delgado's least_squares_main.m and Alan Figueroa's orbit
% (Orbit_and_ErrorAnalysis.m). Put this file in Alec's folder and run it.
%
% WHAT IT DOES
%  1. Runs Alec's code (angular_error_calc + crater_pos) once at nu_ref.
%     The lit craters it places become one fixed, numbered crater group.
%  2. Moves the spacecraft along Alan's orbit to each nu in nu_list. At
%     each spot, every crater gets its sigma from Alec's table (distance,
%     sun angle, bearing, radius) and noisy az/el like crater_pos.m.
%     A crater is "detected" if it is lit and Alec's table has data for
%     its bearing. LS / WLS / BCWLS estimates use the detected craters.
%
% OUTPUTS
%  - Figure: 3-D Moon close-up at each nu. Cyan = visible area, star =
%    spacecraft, green = LOS lines, dots = craters colored by sigma_el
%    (white circles = not detected). Numbers match the tables.
%  - Tables in the Command Window and in orbit_positions_results.xlsx:
%      Summary       : one row per position
%      Sigma_bearing : total sigma of each crater at each nu (NaN = not detected)
%      Bearing       : bearing angle of each crater at each nu
%      Sigma_el : sigma_el of each crater at each nu (blank = not detected)
%      Sigma_az : same for azimuth
%      All_data : one row per crater per position, with Alec's table
%                 stats behind each sigma (same as the CSV below)
%  - sigma_bearing_by_position.csv, bearing_by_position.csv : the two
%    tables above as CSV
%  - crater_sigmas_by_position.csv : sigma and angular-error stats for
%    every crater at every position (fixed Sun, sun angle changes as the
%    spacecraft moves)
%
% WHERE EACH NUMBER COMES FROM  ([A] = Alec's code, [O] = this script)
%  sun angle     [O] acosd(r_sc_hat . sun_hat) at each position. This is
%                    Alec's "moon_angle" / "solar_angle" input
%                    (least_squares_main.m L65 and Solar Phase Experiment
%                    L56-59, L497, L540).
%  table bins    [A] angular_error_calc.m L22-30: loads
%                    discretized_results.mat and slices discretized_data at
%                    the sun-angle bin (discretize_moon_angle, at or below)
%                    and range bin (discretize_range, at or below).
%  mean/std/count[A] angular_error_calc.m L33-44: count, normalized mean
%                    and std for each (bearing bin, radius bin).
%  sigma_bearing [A] least_squares_main.m L90-92: std * Moon angular
%                    size / 100 -> radians ("un-normalize").
%  sigma_el/az   [A] crater_pos.m L146-147: sigma_el = sigma/sqrt((4-pi)/2),
%                    sigma_az = sigma_el / cos(elevation).
%  angular error [A] crater_pos.m L99-100: Rayleigh draw with
%                    b = sigma/sqrt((4-pi)/2), random direction.
%  crater spots  [A] crater_pos.m (random theta, L58) at nu_ref only;
%                    after that the craters stay fixed [O].
%  estimates     [A] LS_estimate.m, WLS_estimate.m, BCWLS_estimate.m
%
% USES (Alec's files, unchanged): angular_error_calc.m,
% discretized_results.mat, crater_pos.m, LS_estimate.m, WLS_estimate.m,
% BCWLS_estimate.m
%
% NOTES
%  - Sun direction is fixed; the sun angle changes because the spacecraft
%    moves. Moon rotation is ignored.
%  - Orbit is 10,000 km from center (Alan used 200 km altitude). Alec's
%    table has far more detections here (craters out to 60 deg bearing).
%
% Author: Jeremiah Symington
% Date: 2026-10-24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INPUT
mu_moon   = 4.9048695e12;          % (m^3/s^2) (Alan)
orbit_a   = 1.74e6 + 8260e3;      % (m) 10,000 km from center (about 8,260 km altitude)
orbit_e   = 0;
orbit_i   = deg2rad(30);           % (rad)
orbit_LAN = deg2rad(20);           % (rad)
orbit_w   = deg2rad(0);            % (rad)

nu_ref    = 60;                    % (deg) where the crater group is built
nu_list   = 5:20:105;              % (deg) 6 positions (1 deg = about 4.2 min at this orbit)
sun_dir_I = [1; 0; 0];             % Sun direction from the Moon, fixed
rng(1);                            % same craters every run

% Optional real Moon image (NOT needed). If this file is in the folder it
% is wrapped on the Moon; otherwise the Moon is plain gray. Free download:
% https://svs.gsfc.nasa.gov/4720 -> lroc_color_poles_4k.tif
moon_texture_file = 'lroc_color_poles_4k.tif';

%% Constants
radius_Moon = 1.7374e6;            % (m)
bin_w   = 15;                      % (deg) bearing bin width in Alec's table
sun_hat = sun_dir_I/norm(sun_dir_I);
n_pos   = length(nu_list);

load('discretized_results.mat', 'discretize_moon_angle', 'discretize_range')   % [A] table bin edges

%% 1. Crater group from Alec's code at nu_ref
r_ref = orbit_state(orbit_a, orbit_e, orbit_i, orbit_LAN, orbit_w, deg2rad(nu_ref), mu_moon);
sun_angle_ref = acosd(dot(r_ref/norm(r_ref), sun_hat));
[~,repeat_matrix_detections] = angular_error_calc(norm(r_ref)/1000, sun_angle_ref);   % [A] table rows
moon_ang_area  = 2*rad2deg(asin(1737.4/(norm(r_ref)/1000)));                      % [A] LSM L90
std_dev_unnorm = deg2rad(repeat_matrix_detections(:,2)*moon_ang_area/100);        % [A] LSM L91
keep = ~isnan(std_dev_unnorm) & std_dev_unnorm > 0;
sc_bearings   = repeat_matrix_detections(keep,3) + bin_w*rand(sum(keep),1);   % spread inside bin
crater_radius = repeat_matrix_detections(keep,4);

[~,~,~,~,r_crater_LOS_ref,R_LOS_I_ref] = crater_pos(r_ref, min(sc_bearings,89.9), std_dev_unnorm(keep));  % [A] places craters
field_I = r_ref + R_LOS_I_ref*r_crater_LOS_ref;            % (m) crater positions, inertial
lit     = (sun_hat'*(field_I./vecnorm(field_I)))' > 0;
field_I = field_I(:,lit);  crater_radius = crater_radius(lit);
n_field = size(field_I,2);
crater_lat = asind(field_I(3,:)'/radius_Moon);  crater_lon = atan2d(field_I(2,:)', field_I(1,:)');

%% 2. Same craters seen from each position
sigma_el_grid = NaN(n_field, n_pos);   sigma_az_grid = NaN(n_field, n_pos);
sigma_b_grid  = NaN(n_field, n_pos);   bearing_grid  = NaN(n_field, n_pos);
summary_rows = [];  res = cell(1,n_pos);
csv_rows = [];                          % crater_sigmas_by_position.csv

for p = 1:n_pos
    r_sc_I = orbit_state(orbit_a, orbit_e, orbit_i, orbit_LAN, orbit_w, deg2rad(nu_list(p)), mu_moon);
    d_km   = norm(r_sc_I)/1000;
    sun_angle = acosd(dot(r_sc_I/norm(r_sc_I), sun_hat));           % [O] Alec's moon_angle input

    % Alec's table at this distance and sun angle
    [~,repeat_matrix_detections] = angular_error_calc(d_km, sun_angle);   % [A] AEC L22-44
    sun_bin   = discretize_moon_angle(find(discretize_moon_angle <= sun_angle, 1, 'last'));  % [A] AEC L26
    range_bin = discretize_range(find(discretize_range <= d_km, 1, 'last'));                % [A] AEC L25
    moon_ang_area = 2*rad2deg(asin(1737.4/d_km));                     % [A] LSM L90
    has_data = repeat_matrix_detections(:,2) > 0;
    bearing_max = min(max(repeat_matrix_detections(has_data,3)) + bin_w, 90);   % edge of visible area

    % Each crater: bearing, lighting, sigma from the table
    n_c       = field_I./vecnorm(field_I);
    c_to_sc   = r_sc_I - field_I;
    bearing   = acosd(sum(n_c.*c_to_sc,1)./vecnorm(c_to_sc))';
    solar_inc = acosd(sun_hat'*n_c)';
    sigma = NaN(n_field,1);  stats = NaN(n_field,4);      % [count, mean_norm, std_norm, mean_deg]
    for j = 1:n_field
        if bearing(j) < 90 && solar_inc(j) < 90
            [sigma(j), stats(j,:)] = table_sigma(repeat_matrix_detections, bearing(j), ...
                                                 crater_radius(j), moon_ang_area, bin_w);
        end
    end
    idx = find(~isnan(sigma));

    % Noisy measurements (crater_pos.m error model) and Alec's estimators
    [std_el, std_az, az_meas, el_meas, r_c_LOS] = fixed_crater_meas(r_sc_I, field_I(:,idx), sigma(idx));  % [A] crater_pos model
    az_nom = atan2(r_c_LOS(2,:), r_c_LOS(1,:));
    el_nom = atan2(r_c_LOS(3,:), vecnorm(r_c_LOS(1:2,:)));
    az_err = NaN(n_field,1);  el_err = NaN(n_field,1);
    az_err(idx) = rad2deg(angle(exp(1i*(az_meas - az_nom))));      % measured - true
    el_err(idx) = rad2deg(el_meas - el_nom);
    errs = NaN(1,3);
    if length(idx) >= 3
        [pos_LS,errs(1),A_mat,z_mat] = LS_estimate(az_meas, el_meas, r_c_LOS);
        [pos_WLS,errs(2),W_block,d_o] = WLS_estimate(std_az, std_el, az_meas, el_meas, r_c_LOS, pos_LS, A_mat, z_mat);
        [~,errs(3)] = BCWLS_estimate(std_az, std_el, az_meas, el_meas, pos_WLS, A_mat, W_block, d_o);
    end

    sigma_el_grid(idx,p) = rad2deg(std_el(:));
    sigma_az_grid(idx,p) = rad2deg(std_az(:));
    sigma_b_grid(:,p)    = rad2deg(sigma);       % total (bearing) sigma, before the az/el split
    bearing_grid(:,p)    = bearing;              % (deg) bearing angle of each crater from here
    summary_rows = [summary_rows; nu_list(p), round(sun_angle,1), length(idx), ...
        round(median(sigma_el_grid(idx,p)),4), round(errs/1000,3)]; %#ok<AGROW>
    csv_rows = [csv_rows; repmat([nu_list(p), round(sun_angle,2), sun_bin, range_bin], n_field, 1), ...
        (1:n_field)', round([crater_lat, crater_lon],3), crater_radius, round(bearing,2), ...
        floor(min(bearing,89.99)/bin_w)*bin_w, ~isnan(sigma), stats(:,1), round(stats(:,2:4),5), ...
        round(rad2deg(sigma),5), round([sigma_el_grid(:,p), sigma_az_grid(:,p), az_err, el_err],5)]; %#ok<AGROW>
    res{p} = struct('r_sc_I',r_sc_I, 'idx',idx, 'sun_angle',sun_angle, 'bearing_max',bearing_max);
end

%% 3. Tables
nu_names = arrayfun(@(v) strrep(sprintf('nu_%g', v), '.', 'p'), nu_list, 'UniformOutput', false);
Summary  = array2table(summary_rows, 'VariableNames', {'nu_deg','sun_angle_deg','craters_detected', ...
               'median_sigma_el_deg','LS_error_km','WLS_error_km','BCWLS_error_km'});
Sigma_el = [table((1:n_field)', crater_radius, 'VariableNames', {'crater','radius'}), ...
            array2table(round(sigma_el_grid,4), 'VariableNames', nu_names)];
Sigma_az = [table((1:n_field)', crater_radius, 'VariableNames', {'crater','radius'}), ...
            array2table(round(sigma_az_grid,4), 'VariableNames', nu_names)];
All_data = array2table(csv_rows, 'VariableNames', {'nu_deg','sun_angle_deg','sun_angle_bin_deg', ...
    'range_bin_km','crater','lat_deg','lon_deg','radius_bin','bearing_deg','bearing_bin_deg','detected', ...
    'table_count','table_mean_norm','table_std_norm','mean_error_deg','sigma_bearing_deg', ...
    'sigma_el_deg','sigma_az_deg','az_error_deg','el_error_deg'});
Sigma_bearing = [table((1:n_field)', crater_radius, 'VariableNames', {'crater','radius'}), ...
                 array2table(round(sigma_b_grid,5), 'VariableNames', nu_names)];
Bearing = [table((1:n_field)', crater_radius, 'VariableNames', {'crater','radius'}), ...
           array2table(round(bearing_grid,2), 'VariableNames', nu_names)];

fprintf('\n=== Summary (one row per position) ===\n');           disp(Summary)
fprintf('=== sigma_bearing (deg): total angular-error sigma of each crater at each nu (NaN = not detected) ===\n')
disp(Sigma_bearing)
fprintf('=== bearing angle (deg) of each crater at each nu ===\n'); disp(Bearing)

% CSV files
writetable(All_data,      'crater_sigmas_by_position.csv');
writetable(Sigma_bearing, 'sigma_bearing_by_position.csv');
writetable(Bearing,       'bearing_by_position.csv');

% Excel file (deleted and rebuilt each run, so define it first)
out_file = 'orbit_positions_results.xlsx';
if exist(out_file,'file'), delete(out_file); end
writetable(Summary,       out_file, 'Sheet', 'Summary');
writetable(Sigma_bearing, out_file, 'Sheet', 'Sigma_bearing');
writetable(Bearing,       out_file, 'Sheet', 'Bearing');
writetable(Sigma_el,      out_file, 'Sheet', 'Sigma_el');
writetable(Sigma_az,      out_file, 'Sheet', 'Sigma_az');
writetable(All_data,      out_file, 'Sheet', 'All_data');
fprintf('Tables saved to %s and 3 CSV files\n', out_file);

%% 4. Figure: 3-D Moon close-up at each position
km = 1/1000;
moon_img = [];
if exist(moon_texture_file,'file')
    moon_img = im2double(imread(moon_texture_file));
    if size(moon_img,3) == 1, moon_img = repmat(moon_img,1,1,3); end
end

% Orbit arc over the craters, and local directions at the crater group
nu_arc = linspace(0, 2*pi, 360);
r_arc = zeros(3,length(nu_arc));
for k = 1:length(nu_arc)
    r_arc(:,k) = orbit_state(orbit_a, orbit_e, orbit_i, orbit_LAN, orbit_w, nu_arc(k), mu_moon);
end
r_pos = cell2mat(cellfun(@(s) s.r_sc_I, res, 'UniformOutput', false));
up    = mean(field_I,2);  up = up/norm(up);
[~,v_ref] = orbit_state(orbit_a, orbit_e, orbit_i, orbit_LAN, orbit_w, deg2rad(nu_ref), mu_moon);
along = v_ref - (v_ref'*up)*up;  along = along/norm(along);
side  = cross(up, along);
fc    = radius_Moon*up*km;                              % center of the crater group (km)
cam   = 0.75*up + 0.55*side - 0.35*along;  cam = cam/norm(cam);
lat0  = asind(up(3));  lon0 = atan2d(up(2), up(1));
fl    = field_I*1.001*km;                               % craters, just above the surface
c_lim = [min(sigma_el_grid(:)) max(sigma_el_grid(:))+1e-9];



for p = 1:n_pos
    figure('Color','k', 'Position', [100 100 900 750]); hold on
    moon_surface([-90 90], [-180 180], 0.5, radius_Moon*km, moon_img)
    light('Position', sun_hat', 'Style', 'infinite')

    ring = footprint_ring(res{p}.r_sc_I, res{p}.bearing_max, radius_Moon)*km;
    fill3(ring(1,:), ring(2,:), ring(3,:), 'c', 'FaceAlpha', 0.3, 'EdgeColor', 'c')
    plot3(r_arc(1,:)*km, r_arc(2,:)*km, r_arc(3,:)*km, 'w-')
    plot3(r_pos(1,:)*km, r_pos(2,:)*km, r_pos(3,:)*km, 'w.', 'MarkerSize', 10)

    rs = res{p}.r_sc_I*km;  below = radius_Moon*km*rs/norm(rs);
    plot3([rs(1) below(1)], [rs(2) below(2)], [rs(3) below(3)], 'w--')
    idx = res{p}.idx;  nd = setdiff(1:n_field, idx);
    for j = idx(:)'
        plot3([rs(1) fl(1,j)], [rs(2) fl(2,j)], [rs(3) fl(3,j)], 'g-')
    end
    if ~isempty(nd), plot3(fl(1,nd), fl(2,nd), fl(3,nd), 'wo', 'MarkerSize', 4); end
    if ~isempty(idx)
        scatter3(fl(1,idx), fl(2,idx), fl(3,idx), 28, sigma_el_grid(idx,p), 'filled', 'MarkerEdgeColor', 'k')
    end
    text(fl(1,:), fl(2,:), fl(3,:), arrayfun(@(j) sprintf(' %d',j), 1:n_field, 'UniformOutput', false), ...
        'Color', 'y', 'FontSize', 7)
    plot3(rs(1), rs(2), rs(3), 'p', 'Color', 'b', 'MarkerFaceColor', 'c', 'MarkerSize', 14)

    caxis(c_lim); axis equal off
    camproj('perspective'); ctr = 0.5*mean(r_pos,2)*km;
    camtarget(ctr'); campos((ctr + 40000*cam)'); camup(up'); camva(30)
    title(sprintf('(%d)  \\nu = %g°, sun angle %.0f°, %d detected', p, nu_list(p), res{p}.sun_angle, length(idx)), ...
        'Color', 'w', 'FontWeight', 'normal', 'FontName', 'Times New Roman', 'FontSize', 12)
    cb = colorbar; cb.Label.String = '\sigma_{el} (deg)'; cb.Color = 'w';
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

function [s, stats] = table_sigma(T, bearing, radius, moon_ang_area, bin_w)
% Sigma (rad) from Alec's table (rows of repeat_matrix_detections:
% [mean_norm std_norm bearing_bin radius_bin], one row per detected crater)
% for this crater's bearing bin, closest radius bin. NaN if no detections.
% stats = [count, mean_norm, std_norm, mean_deg] of that table bin.
s = NaN;  stats = NaN(1,4);
rows = find(T(:,3) == floor(min(bearing,89.99)/bin_w)*bin_w & T(:,2) > 0);
if isempty(rows), return, end
[~,k] = min(abs(T(rows,4) - radius));
r = T(rows(k),:);
s = deg2rad(r(2)*moon_ang_area/100);                   % [A] LSM L91 un-normalize
stats = [sum(T(:,3) == r(3) & T(:,4) == r(4)), r(1), r(2), r(1)*moon_ang_area/100];   % [A] LSM L92 mean
end

function [std_el, std_az, az_meas, el_meas, r_c_LOS] = fixed_crater_meas(r_sc_I, r_crater_I, sigma)
% crater_pos.m frames and error model, for craters at known positions
x_aux = r_sc_I/norm(r_sc_I);
y_aux = cross([0;0;1], x_aux);  y_aux = y_aux/norm(y_aux);
z_aux = cross(x_aux, y_aux);
R_LOS_I = [y_aux, -x_aux, z_aux];                     % LOS -> inertial
r_c_LOS = R_LOS_I'*(r_crater_I - r_sc_I);             % (m) spacecraft -> crater
sigma   = sigma(:);
err     = raylrnd(sigma/sqrt((4-pi)/2));              % Rayleigh angle error
th      = 2*pi*rand(length(sigma),1);
r_meas  = zeros(size(r_c_LOS));
for i = 1:length(sigma)
    x_e = r_c_LOS(:,i)/norm(r_c_LOS(:,i));
    y_e = cross([0;0;1], x_e);  y_e = y_e/norm(y_e);
    z_e = cross(x_e, y_e);
    r_meas(:,i) = [x_e y_e z_e]*norm(r_c_LOS(:,i))*[cos(err(i)); sin(err(i))*cos(th(i)); sin(err(i))*sin(th(i))];
end
el_nom  = atan2(r_c_LOS(3,:), vecnorm(r_c_LOS(1:2,:)));
az_meas = atan2(r_meas(2,:), r_meas(1,:));
el_meas = atan2(r_meas(3,:), vecnorm(r_meas(1:2,:)));
std_el  = sigma/sqrt((4-pi)/2);                        % (rad) crater_pos.m
std_az  = std_el./cos(el_nom(:));
end

function ring = footprint_ring(r_sc, bearing_max, R)
% Edge of the visible area on the surface (3xN, m), crater_pos.m geometry
d = norm(r_sc);  b = deg2rad(bearing_max);
gamma = b - asin(R*sin(b)/d);                          % central angle to the edge
u  = r_sc/d;
e1 = cross([0;0;1], u);  e1 = e1/norm(e1);  e2 = cross(u, e1);
t  = linspace(0, 2*pi, 120);
ring = 1.002*R*(cos(gamma)*u + sin(gamma)*(e1*cos(t) + e2*sin(t)));
end

function moon_surface(lat_lim, lon_lim, step, R, img)
% Patch of the Moon (lat/lon box, deg). img = optional equirectangular
% Moon image; empty -> plain gray.
[LON, LAT] = meshgrid(lon_lim(1):step:lon_lim(2), lat_lim(1):step:lat_lim(2));
X = R*cosd(LAT).*cosd(LON);  Y = R*cosd(LAT).*sind(LON);  Z = R*sind(LAT);
if isempty(img)
    C = 0.6*ones([size(X) 3]);
else
    [H, W, ~] = size(img);
    row = min(max(round((90 - LAT)/180*(H-1)) + 1, 1), H);
    col = mod(round(mod(LON + 180, 360)/360*(W-1)), W) + 1;
    C = zeros([size(X) 3]);
    for k = 1:3, band = img(:,:,k); C(:,:,k) = band(sub2ind([H W], row, col)); end
end
surf(X, Y, Z, C, 'EdgeColor', 'none', 'FaceColor', 'interp', 'FaceLighting', 'gouraud', ...
    'AmbientStrength', 0.45, 'DiffuseStrength', 0.7, 'SpecularStrength', 0)
end