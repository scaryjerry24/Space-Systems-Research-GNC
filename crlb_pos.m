function [cov_p, sig_p] = crlb_pos(r_crater_LOS_nom, std_azimuth, std_elevation)
% Position-only bound (3x3), LOS frame. Same model as LS/WLS_estimate,
% evaluated at the true crater angles.
    N = size(r_crater_LOS_nom, 2);
    A = zeros(2*N, 3);
    W = zeros(2*N);
    for i = 1:N
        r  = r_crater_LOS_nom(:,i);
        az = atan2(r(2), r(1));
        el = atan2(r(3), hypot(r(1), r(2)));
        u_th  = [ sin(az); -cos(az); 0];
        u_psi = [-cos(az)*sin(el); -sin(az)*sin(el); cos(el)];
        A(2*i-1:2*i,:) = [u_th'; u_psi'];
        d  = norm(r);
        Do = diag([-d*cos(el), d]);
        R_eta = Do * diag([std_azimuth(i)^2, std_elevation(i)^2]) * Do';
        W(2*i-1:2*i, 2*i-1:2*i) = inv(R_eta);
    end
    cov_p = inv(A' * W * A);
    sig_p = sqrt(diag(cov_p));     % 1-sigma: [x (HCA); y (LOS); z (VCA)]
end