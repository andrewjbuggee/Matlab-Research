function [fig, png_path] = draw_model_diagram(P, out_dir)
% DRAW_MODEL_DIAGRAM  Render the column as configured in P, in the style of the
%                     26 Oct 2026 cartoons: boxes for the layers that exist,
%                     arrows for the fluxes that connect them, and the physics
%                     options written on the right.
%
%   [fig, png] = draw_model_diagram(P)            figure only (invisible)
%   [fig, png] = draw_model_diagram(P, out_dir)   also writes model_diagram_levelN.png
%
% What is drawn is read from P, so the picture always matches the configuration
% that will run: the number and names of atmospheric layers, whether the cloud is
% a layer of its own or an emissivity of a thicker slab, the skin, snow, the ice
% as one slab or N layers of the configured fractions, the ocean slabs, and
% every flux path -- longwave both ways between neighbours, SH/LH between the
% lowest air and the surface, conduction through the ice, the ocean heat flux,
% freezing at the base, mixing between coupled air layers, and advection into
% each air layer.

if nargin < 2, out_dir = ''; end
a = P.atm;  o = P.opt;  ice = P.ice;
esc = @(str) strrep(str, '_', '\_');            % option names carry underscores; tex would subscript them
slab_ice = strcmpi(o.ice_model, 'slab');
zero_layer = strcmpi(o.ice_model, 'zero_layer');
has_snow = ice.h_snow_m > 0 && ~slab_ice;
deep = P.ocean.deep_prognostic;

% ---- vertical layout (arbitrary units, top to bottom) ---------------------
rows = {};                                  % {label, height, colour, kind}
rows(end+1,:) = {'space', 0.35, [0.05 0.05 0.15], 'space'};
for k = a.n_layers:-1:1
    L = a.layers(k);
    lab = sprintf('%s   (%g hPa, \\epsilon_{gas} = %.2f)', L.name, L.dp_hPa, L.emis_gas);
    col = [0.80 0.88 0.97];
    if L.has_cloud
        lab = sprintf('%s   (%g hPa, \\epsilon_{gas} = %.2f; cloud: %s %g g m^{-2}, \\epsilon_{slab} = %.2f)', ...
            L.name, L.dp_hPa, L.emis_gas, P.cloud.phase, P.cloud.WP_g_m2, P.cloud.emis_slab);
        col = [0.86 0.86 0.92];
    end
    rows(end+1,:) = {lab, 0.9, col, sprintf('atm%d', k)};
end
if ~slab_ice
    rows(end+1,:) = {'skin  (no heat capacity: fluxes sum to zero)', 0.12, [0.55 0.55 0.55], 'skin'};
    if has_snow, rows(end+1,:) = {sprintf('snow  %.2f m, k = %.2f W m^{-1} K^{-1}', ice.h_snow_m, ice.k_snow_W_m_K), 0.35, [0.95 0.95 0.95], 'snow'}; end
    fr = ice.layer_fractions;
    if zero_layer
        rows(end+1,:) = {sprintf('ice  (linear profile, no heat capacity; Semtner 0-layer / Thorndike 1992; h = %.1f m)', ice.h_init_m), 1.4, [0.78 0.90 0.95], 'ice1'};
    else
        for l = 1:ice.N_layers
            rows(end+1,:) = {sprintf('ice layer %d  (%.0f %% of h)', l, 100*fr(l)), max(0.25, 2.0*fr(l)), [0.78 0.90 0.95], sprintf('ice%d', l)};
        end
    end
else
    rows(end+1,:) = {sprintf('ice slab  (one temperature, h = %.1f m)', ice.h_init_m), 1.4, [0.78 0.90 0.95], 'ice1'};
end
rows(end+1,:) = {sprintf('ocean mixed layer  (%g m)', P.ocean.h_ml_m), 0.6, [0.55 0.70 0.90], 'ml'};
if deep, rows(end+1,:) = {sprintf('deep ocean  (%g m reservoir)', P.ocean.h_deep_m), 0.6, [0.35 0.50 0.80], 'deep'}; end

nR = size(rows,1);
y_top = zeros(nR,1); y_bot = zeros(nR,1);  y = 0;
for i = 1:nR, y_top(i) = y; y_bot(i) = y - rows{i,2}; y = y_bot(i) - 0.06; end
kinds = rows(:,4);
idx = @(k) find(strcmp(kinds, k), 1);
ymid = @(i) 0.5*(y_top(i) + y_bot(i));

fig = new_figure_clean('Position', [60 60 1500 900]);
ax = axes(fig); hold(ax, 'on'); axis(ax, 'off');
x0 = 3.0; x1 = 7.6;                          % box horizontal extent
for i = 1:nR
    rectangle(ax, 'Position', [x0 y_bot(i) x1-x0 rows{i,2}], 'FaceColor', rows{i,3}, 'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 1.0);
    tc = 'k'; if strcmp(kinds{i}, 'space') || strcmp(kinds{i}, 'deep'), tc = 'w'; end
    text(ax, 0.5*(x0+x1), ymid(i), rows{i,1}, 'HorizontalAlignment', 'center', 'FontSize', 9.5, 'Color', tc, 'Interpreter', 'tex');
end

% ---- fluxes ----------------------------------------------------------------
cLW = [0.70 0.09 0.17]; cT = [0.30 0.45 0.69]; cC = [0.51 0.45 0.70]; cO = [0.20 0.50 0.80]; cA = [0.95 0.60 0.10];
% longwave both ways between every pair of neighbouring radiative layers
rad = [idx('space'); find(startsWith(kinds, 'atm')); idx(tern(slab_ice, 'ice1', 'skin'))];
for i = 1:numel(rad)-1
    u = rad(i); d = rad(i+1);
    arrow(ax, x0+0.55, y_top(d)+0.02, x0+0.55, y_bot(u)-0.02, cLW);         % upward LW
    if ~strcmp(kinds{u}, 'space')
        arrow(ax, x0+0.85, y_bot(u)-0.02, x0+0.85, y_top(d)+0.02, cLW);     % downward LW
    end
end
text(ax, x0+0.70, y_bot(1)-0.02, 'LW', 'Color', cLW, 'FontSize', 8.5, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
% SH, LH between the lowest air layer and the surface
lo = idx('atm1'); sf = idx(tern(slab_ice, 'ice1', 'skin'));
arrow(ax, x0+1.9, y_bot(lo)-0.02, x0+1.9, y_top(sf)+0.02, cT);  text(ax, x0+1.9, y_bot(lo)+0.08, 'SH', 'Color', cT, 'FontSize', 8.5, 'HorizontalAlignment', 'center');
if ~strcmpi(o.latent_scheme, 'off')
    arrow(ax, x0+2.3, y_bot(lo)-0.02, x0+2.3, y_top(sf)+0.02, [0.33 0.66 0.41]);  text(ax, x0+2.3, y_bot(lo)+0.08, 'LH', 'Color', [0.33 0.66 0.41], 'FontSize', 8.5, 'HorizontalAlignment', 'center');
end
% conduction down through snow and ice to the base, and the ocean heat flux up
cond = find(startsWith(kinds, 'skin') | startsWith(kinds, 'snow') | startsWith(kinds, 'ice'));
for i = 1:numel(cond)-1
    arrow(ax, x0+3.2, y_bot(cond(i))-0.02, x0+3.2, y_top(cond(i+1))+0.02, cC);
end
ib = cond(end); ml = idx('ml');
arrow(ax, x0+3.2, y_bot(ib)-0.02, x0+3.2, y_top(ml)+0.02, cC);  text(ax, x0+3.45, y_bot(ib)-0.08, 'k dT/dz', 'Color', cC, 'FontSize', 8.5);
arrow(ax, x0+3.8, y_top(ml)+0.02, x0+3.8, y_bot(ib)-0.02, cO);  text(ax, x0+3.95, y_top(ml)+0.10, 'F_{ocean}', 'Color', cO, 'FontSize', 8.5);
text(ax, x1+0.1, y_bot(ib), sprintf('T = T_f = %.2f K: freezing, \\rho L dh/dt = k dT/dz - F_{ocean}%s', P.ocean.T_f_K, tern(strcmpi(o.ice_growth,'fixed'), '  [growth OFF: fixed h]', '')), ...
    'FontSize', 8.5, 'Color', [0.1 0.1 0.1], 'VerticalAlignment', 'middle');
if deep
    dp = idx('deep');
    arrow(ax, x0+3.5, y_top(dp)+0.02, x0+3.5, y_bot(ml)-0.02, cO);  text(ax, x0+3.7, ymid(dp)+0.2, sprintf('\\gamma_{deep} = %g W m^{-2} K^{-1}', P.ocean.gamma_deep_W_m2_K), 'Color', cO, 'FontSize', 8.5);
end
% advection into each air layer; mixing between coupled layers
for k = 1:a.n_layers
    i = idx(sprintf('atm%d', k));
    arrow(ax, x0-0.9, ymid(i), x0-0.05, ymid(i), cA);
    if strcmpi(o.atm_forcing, 'heating')
        adv_txt = sprintf('advection: heating %.1f W m^{-2} (D = %g)  ', a.D_adv_W_m2*a.D_weight_vec(k), a.D_adv_W_m2);
    else
        adv_txt = sprintf('advection: relax to %.0f K, \\tau = %.1f d  ', a.layers(k).T_bg_K, a.tau_adv_s/86400);
    end
    text(ax, x0-0.95, ymid(i), adv_txt, 'Color', cA, 'FontSize', 8.5, 'HorizontalAlignment', 'right');
    if a.layers(k).gamma_mix_above_W_m2_K > 0 && k < a.n_layers
        j = idx(sprintf('atm%d', k+1));
        arrow(ax, x1-0.5, y_top(i)+0.02, x1-0.5, y_bot(j)-0.02, [0.4 0.4 0.4]);
        arrow(ax, x1-0.3, y_bot(j)-0.02, x1-0.3, y_top(i)+0.02, [0.4 0.4 0.4]);
        text(ax, x1-0.4, y_bot(j)+0.1, sprintf('mixing \\gamma = %g', a.layers(k).gamma_mix_above_W_m2_K), 'Color', [0.3 0.3 0.3], 'FontSize', 8, 'HorizontalAlignment', 'center');
    end
end
text(ax, x0+0.55, y_top(1)+0.12, 'OLR', 'Color', cLW, 'FontSize', 8.5, 'HorizontalAlignment', 'center');

% ---- the physics options, on the right ---------------------------------------
opt_lines = {
    sprintf('\\bf complexity level %d', o.complexity)
    sprintf('longwave: grey slabs, emission = %s', esc(o.lw_emission))
    sprintf('   surface \\epsilon = %.3f; cloud \\epsilon = 1 - exp(-k_{abs} WP), k_{abs} = %.3f m^2 g^{-1}', ice.emis_sfc, P.cloud.k_abs_m2_g)
    sprintf('turbulent fluxes: %s%s', esc(o.turb_scheme), tern(strcmpi(o.turb_scheme,'bulk_const'), sprintf(' (C_H = %.1e)', a.C_H), tern(strcmpi(o.turb_scheme,'most'), sprintf(' (z_{ref} = %g m)', o.z_ref_m), tern(strcmpi(o.turb_scheme,'linear'), sprintf(' (\\lambda_{SH} = %.1f)', o.lambda_SH_lin_W_m2_K), ''))))
    sprintf('   latent: %s; U = %g m s^{-1}; RH_{ice} = %.2f', esc(o.latent_scheme), a.U_m_s, a.RH_ice_la)
    sprintf('air profile: %s%s; advection: %s', esc(o.atm_profile), tern(strcmpi(o.atm_profile,'lapse'), sprintf(' (%s)', tern(isnumeric(o.lapse_kind), sprintf('%g K/km', 1000*o.lapse_kind), o.lapse_kind)), ''), esc(o.atm_forcing))
    sprintf('ice: %s%s', tern(slab_ice, 'single slab', tern(zero_layer, 'zero-layer (no heat capacity)', sprintf('skin + %d layers [%s]', ice.N_layers, strtrim(sprintf('%g ', ice.layer_fractions))))), tern(ice.use_brine, ', brine corrections', ''))
    sprintf('   k = %.2f W m^{-1} K^{-1}, \\rho c = %.2g J m^{-3} K^{-1}, L = %.3g J kg^{-1}', ice.k_W_m_K, ice.rho_c_J_m3_K, P.const.L_fus_J_kg)
    sprintf('   growth: %s', esc(o.ice_growth))
    sprintf('ocean: %s; F_{ocean} = \\rho c c_h u_* (T_{ml} - T_f), coeff %.0f W m^{-2} K^{-1}', esc(o.ocean_model), P.ocean.coeff_iceocean_W_m2_K)
    sprintf('time step %g s; no solar (polar night)', P.num.dt_s)};
text(ax, x1+0.15, y_top(1), opt_lines, 'FontSize', 9, 'VerticalAlignment', 'top', 'Interpreter', 'tex');

xlim(ax, [-1.6 13.5]); ylim(ax, [y_bot(end)-0.3, y_top(1)+0.45]);
title(ax, sprintf('Column model, complexity level %d: energy moves only up and down', o.complexity), 'FontSize', 12);

png_path = '';
if ~isempty(out_dir)
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end
    png_path = fullfile(out_dir, sprintf('model_diagram_level%d.png', o.complexity));
    exportgraphics(fig, png_path, 'Resolution', 150);
end
end

function arrow(ax, x1, y1, x2, y2, col)
% a straight arrow in data units with a small filled head
plot(ax, [x1 x2], [y1 y2], '-', 'Color', col, 'LineWidth', 1.4);
L = hypot(x2-x1, y2-y1); if L == 0, return; end
ux = (x2-x1)/L; uy = (y2-y1)/L;  hl = 0.09; hw = 0.05;
px = x2 - hl*ux; py = y2 - hl*uy;
patch(ax, [x2, px + hw*uy, px - hw*uy], [y2, py - hw*ux, py + hw*ux], col, 'EdgeColor', col);
end

function s = tern(c, a, b)
if c, s = a; else, s = b; end
end
