function [T_K, p_Pa] = cloud_base_T_p(t_win, z_msl_m, P)
% CLOUD_BASE_T_P  INTERPSONDE temperature and pressure at a given height, for a list of 30-min windows.
%
%   [T_K, p_Pa] = cloud_base_T_p(t_win, z_msl_m, P)
%
%   For each window t_win(i), takes the 30-min mean INTERPSONDE profile of the same
%   window (P from interpsonde_30min_profiles) and interpolates it linearly in
%   height to z_msl_m(i), the cloud base (Dedrick et al. 2026, p. 1065: cloud-base
%   thermodynamics from INTERPSONDE at the ceilometer cloud-base height). Pressure
%   is interpolated in ln(p), which is linear in height for an isothermal layer
%   (the hypsometric equation, e.g. Wallace & Hobbs 2006, Ch. 3; I have not
%   checked the equation number); over 20-m levels the difference from linear
%   interpolation of p is negligible (my estimate).
%
%   Inputs
%     t_win     window starts (datetime column)
%     z_msl_m   cloud-base height for each window, m on INTERPSONDE's height scale
%               (NaN gives NaN)
%     P         struct from interpsonde_30min_profiles
%   Outputs
%     T_K, p_Pa cloud-base temperature (K) and pressure (Pa); NaN where the window
%               has no INTERPSONDE profile or z is outside the profile

T_K = NaN(numel(t_win), 1);
p_Pa = NaN(numel(t_win), 1);
[has, loc] = ismember(t_win, P.time);
for i = find(has & isfinite(z_msl_m))'
    T_prof = P.T_K(loc(i), :);
    p_prof = P.p_Pa(loc(i), :);
    ok = isfinite(T_prof) & isfinite(p_prof);
    if nnz(ok) < 2, continue, end
    T_K(i) = interp1(P.height_m(ok), T_prof(ok), z_msl_m(i), 'linear', NaN);
    p_Pa(i) = exp(interp1(P.height_m(ok), log(p_prof(ok)), z_msl_m(i), 'linear', NaN));
end

end
