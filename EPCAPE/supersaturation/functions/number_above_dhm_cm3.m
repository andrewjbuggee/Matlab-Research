function n_above_cm3 = number_above_dhm_cm3(N_cm3, Dg_um, sigma_g, d_hm_um)
% NUMBER_ABOVE_DHM_CM3  Number of particles larger than the Hoppel minimum, from fitted lognormal modes.
%
%   n_above_cm3 = number_above_dhm_cm3(N_cm3, Dg_um, sigma_g, d_hm_um)
%
%   For a lognormal mode with number N, median diameter Dg and geometric
%   standard deviation sigma_g, the number of particles larger than D is
%
%       N(D' > D) = N/2 * erfc( ln(D / Dg) / (sqrt(2) ln sigma_g) )
%
%   Summing over the Aitken, accumulation and coarse modes at D = D_HM gives the
%   number of particles that activated in recent clouds, the target of the
%   s_CCN retrieval (Dedrick et al. 2026, p. 1065; D_HM defined on p. 1062).
%   This is the cumulative size distribution (Rogers & Yau 1989, p. 91, eq. 6.11)
%   counted from D_HM upward; for a lognormal mode (a normal distribution in ln D)
%   the integral is the erfc expression above (a mathematical identity).
%
%   Inputs (one row per time; columns = Aitken, accumulation, coarse)
%     N_cm3     mode number concentrations, cm-3 (K x 3)
%     Dg_um     mode median diameters, um (K x 3)
%     sigma_g   mode geometric standard deviations (K x 3)
%     d_hm_um   Hoppel minimum diameter, um (K x 1)
%   Output
%     n_above_cm3  (K x 1), cm-3. NaN where D_HM is missing. A mode with a missing
%                  parameter contributes 0 (as in the Python notebook).

n_above_cm3 = zeros(size(d_hm_um));
for m = 1:3
    part = 0.5 .* N_cm3(:, m) .* erfc(log(d_hm_um ./ Dg_um(:, m)) ./ (sqrt(2) .* log(sigma_g(:, m))));
    part(~isfinite(part)) = 0;   % a mode with a missing parameter adds nothing
    n_above_cm3 = n_above_cm3 + part;
end
n_above_cm3(~isfinite(d_hm_um)) = NaN;

end
