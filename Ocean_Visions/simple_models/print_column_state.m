function print_column_state(R, label)
% PRINT_COLUMN_STATE  One-screen summary of the final state of a run.
P = R.P;
fprintf('  %s:\n', label);
fprintf('    T_sk %.1f K | T_atm [%s] K | T_ml - T_f %.3f K | h %.3f m\n', R.T_sk_K(end), ...
    strtrim(sprintf('%.1f ', R.T_atm_K(:,end))), R.T_ml_K(end) - P.ocean.T_f_K, R.h_m(end));
fprintf('    LWD %.1f | LWU %.1f | SH_dn %+.1f | LH_dn %+.1f | cond up to skin %.1f | cond at base %.1f | ocean %.1f W/m2 | growth %.2f cm/day\n', ...
    R.LWD_sfc_W_m2(end), R.LWU_sfc_W_m2(end), R.SH_down_W_m2(end), R.LH_down_W_m2(end), ...
    R.F_cond_top_up_W_m2(end), R.F_cond_bot_up_W_m2(end), R.F_ocean_up_W_m2(end), 100*86400*R.dhdt_bot_m_s(end));
h_eff = R.h_m(end) + P.ice.h_snow_m * P.ice.k_W_m_K / P.ice.k_snow_W_m_K;
fprintf('    damping: lambda_LW %.2f, lambda_SH %.2f, lambda_LH %.2f W/m2/K; steady lambda_G = k/h_eff = %.2f', ...
    R.lambda_LW_W_m2_K(end), R.lambda_SH_W_m2_K(end), R.lambda_LH_W_m2_K(end), P.ice.k_W_m_K/h_eff);
if isfinite(R.C_H(end)), fprintf(' | C_H %.2e, Ri_b %.3f, zeta %.2f', R.C_H(end), R.Ri_b(end), R.zeta(end)); end
fprintf('\n');
end
