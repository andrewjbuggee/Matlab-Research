function P = model_config(level, varargin)
% MODEL_CONFIG  Build the parameter/option struct for one complexity level.
%
%   P = model_config(3)                                   the 22 Oct 2026 model
%   P = model_config(1, 'lw_emission', 'linearized')      simplest structure, linear radiation
%   P = model_config(4, 'turb_scheme', 'most', 'lapse_kind', 'moist', 'h_snow_m', 0.2)
%   P = model_config(..., 'quiet', true)                  no summary printout
%
% Thin front door to seb_column_params: 'complexity' selects the layer
% structure (seb_column_finalize.m builds P.atm.layers), every other name/value
% pair overrides a leaf field anywhere in P (physics options in P.opt, physical
% parameters in P.ice / P.atm / P.ocean / P.cloud / P.num).  A summary of the
% configuration is printed unless 'quiet' is true.  To edit a single layer, do
% it on P.atm.layers(k) and call P = seb_column_finalize(P).

quiet = false;
keep = true(size(varargin));
for i = 1:2:numel(varargin)
    if strcmpi(varargin{i}, 'quiet'), quiet = varargin{i+1}; keep(i:i+1) = false; end
end
varargin = varargin(keep);
P = seb_column_params('complexity', level, varargin{:});
if ~quiet, model_summary(P); end
end

function model_summary(P)
o = P.opt;  a = P.atm;
fprintf('\n--- model configuration: complexity level %d ---\n', o.complexity);
fprintf('  atmosphere: %d layer(s), profile = %s', a.n_layers, o.atm_profile);
if strcmpi(o.atm_profile, 'lapse')
    if isnumeric(o.lapse_kind), fprintf(' (%g K/km)', 1000*o.lapse_kind); else, fprintf(' (%s)', o.lapse_kind); end
end
fprintf('\n');
for k = a.n_layers:-1:1
    L = a.layers(k);
    fprintf('    [%d] %-24s dp = %4g hPa (%4.0f m)  eps_gas = %.2f  T_bg = %5.1f K  mix-up = %4.1f%s\n', ...
        k, L.name, L.dp_hPa, L.dz_m, L.emis_gas, L.T_bg_K, L.gamma_mix_above_W_m2_K, ...
        tern(L.has_cloud, sprintf('  + cloud (%s %g g/m2, eps %.3f -> slab %.3f)', P.cloud.phase, P.cloud.WP_g_m2, P.cloud.emis, P.cloud.emis_slab), ''));
end
fprintf('  longwave  : grey slabs, emission = %s\n', o.lw_emission);
fprintf('  turbulence: %s', o.turb_scheme);
switch lower(o.turb_scheme)
    case 'bulk_const', fprintf(' (C_H = %.1e, U = %g m/s)', a.C_H, a.U_m_s);
    case 'most',       fprintf(' (z_ref = %g m, z0m = %.1e, z0h = %.1e, U = %g m/s)', o.z_ref_m, o.z0m_m, o.z0h_m, a.U_m_s);
    case 'linear',     fprintf(' (lambda_SH = %.2f W/m2/K)', o.lambda_SH_lin_W_m2_K);
end
fprintf(', latent = %s\n', o.latent_scheme);
if strcmpi(o.ice_model, 'slab')
    fprintf('  ice       : single slab, h = %.2f m, growth = %s\n', P.ice.h_init_m, o.ice_growth);
else
    fprintf('  ice       : skin + %d layer(s) [%s], h = %.2f m, snow = %.2f m, growth = %s%s\n', P.ice.N_layers, ...
        strtrim(sprintf('%g ', P.ice.layer_fractions)), P.ice.h_init_m, P.ice.h_snow_m, o.ice_growth, tern(P.ice.use_brine, ', brine corrections', ''));
end
fprintf('  ocean     : %s (ML %g m; deep %s)\n', o.ocean_model, P.ocean.h_ml_m, tern(P.ocean.deep_prognostic, 'prognostic', 'fixed'));
fprintf('  advection : tau_adv = %.2g d; atmosphere mode = %s\n', a.tau_adv_s/86400, a.mode);
end

function s = tern(c, a, b)
if c, s = a; else, s = b; end
end
