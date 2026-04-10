%% Plot the Droplet Size Distribution as an imagesc plot versus Altitude
%
%   Creates a 2D color image where the x-axis is droplet radius (microns),
%   the y-axis is altitude (m), and the color represents the number
%   concentration Nc (cm^-3) at each size bin and altitude level.
%   A black solid line marks the weighted mean radius at each altitude, and
%   dashed black lines mark the +/- 1 weighted standard deviation bounds.
%
%   Inputs:
%     vert_profile     - single VOCALS-REx profile struct (e.g. ensemble_profiles{1,1})
%     instrument       - 'all'  : combined CDP + 2DC distribution
%                        'CDP'  : CDP size bins only (radius <= 25 um)
%     three_panel_flag - false : single imagesc panel only
%                        true  : three panels: imagesc | re vs altitude |
%                                effective variance vs altitude
%
% By Andrew John Buggee
%%

function [] = plot_droplet_size_dist_vs_altitude(vert_profile, instrument, three_panel_flag)

if nargin ~= 3
    error([newline, 'Wrong number of inputs. Need 3: vert_profile, instrument, three_panel_flag', newline])
end

if ~isstruct(vert_profile) || ~isscalar(vert_profile)
    error([newline, 'vert_profile must be a single (scalar) struct.', newline])
end

if ~ismember(instrument, {'all', 'CDP'})
    error([newline, 'instrument must be ''all'' or ''CDP''.', newline])
end

% -------------------------------------------------------------------------
% Font sizes
% -------------------------------------------------------------------------
ttl_fnt = 18;
ax_fnt  = 18;
cb_fnt  = 16;
ln_sz   = 2;

% -------------------------------------------------------------------------
% Extract data
% -------------------------------------------------------------------------
r_bins_all = double(vert_profile.drop_radius_bin_center);   % 1x91, microns
altitude   = vert_profile.altitude;                          % 1x28, meters
Nc_full    = double(vert_profile.Nc);                        % 91x28, cm^-3

% -------------------------------------------------------------------------
% Select instrument range
% -------------------------------------------------------------------------
% CDP typically measures droplet radii up to ~25 um (50 um diameter).
% Bins beyond this are measured by the 2DC instrument.
CDP_r_max_um = 25;

use_2DC = vert_profile.flag_2DC_data_is_conforming && strcmp(instrument, 'all');

if strcmp(instrument, 'CDP') || ~vert_profile.flag_2DC_data_is_conforming
    cdp_mask     = r_bins_all <= CDP_r_max_um;
    r_bins_plot  = r_bins_all(cdp_mask);
    Nc_plot      = Nc_full(cdp_mask, :);
    re_plot      = vert_profile.re_CDP;        % precomputed effective radius
    instr_label  = 'CDP';
    if ~vert_profile.flag_2DC_data_is_conforming && strcmp(instrument, 'all')
        instr_label = 'CDP (2DC non-conforming)';
    end
else
    r_bins_plot  = r_bins_all;
    Nc_plot      = Nc_full;
    re_plot      = vert_profile.re;            % precomputed effective radius
    instr_label  = 'CDP + 2DC';
end

n_bins = length(r_bins_plot);
n_alt  = length(altitude);

% -------------------------------------------------------------------------
% Compute weighted mean radius and +/- 1 weighted std dev at each altitude
%
% Note: the distributions at each level may not be Gaussian (normFit,
% logNormFit, and gammaFit in the struct capture the best-fit shapes), but
% the weighted mean and std dev are still a useful descriptor of the
% distribution width regardless of shape.
% -------------------------------------------------------------------------
r_col  = r_bins_plot(:);
r_mean = NaN(1, n_alt);
r_std  = NaN(1, n_alt);

for kk = 1:n_alt
    nc_col = Nc_plot(:, kk);
    nc_col(nc_col < 0) = 0;
    total = sum(nc_col);
    if total > 0
        mu         = sum(r_col .* nc_col) / total;
        sigma      = sqrt(sum((r_col - mu).^2 .* nc_col) / total);
        r_mean(kk) = mu;
        r_std(kk)  = sigma;
    end
end

r_lo = max(r_mean - r_std, min(r_bins_plot));
r_hi = min(r_mean + r_std, max(r_bins_plot));

% -------------------------------------------------------------------------
% Prepare log10 of Nc for display (zeros/NaNs -> NaN so they show as background)
% -------------------------------------------------------------------------
Nc_log = log10(Nc_plot);
Nc_log(~isfinite(Nc_log)) = NaN;

% imagesc expects the matrix as [n_y x n_x] = [n_alt x n_bins]
Nc_img = Nc_log';    % n_alt x n_bins

% -------------------------------------------------------------------------
% Compute effective variance at each altitude for three-panel mode
%   v_e(z) = sum[(r - r_e)^2 * r^2 * n(r)] / (r_e^2 * sum[r^2 * n(r)])
% -------------------------------------------------------------------------
v_eff = NaN(1, n_alt);
if three_panel_flag
    for kk = 1:n_alt
        nc_col = Nc_plot(:, kk);
        nc_col(nc_col < 0) = 0;
        r_e_k  = re_plot(kk);
        denom  = sum(r_col.^2 .* nc_col);
        if denom > 0 && isfinite(r_e_k) && r_e_k > 0
            v_eff(kk) = sum((r_col - r_e_k).^2 .* r_col.^2 .* nc_col) ...
                        / (r_e_k^2 * denom);
        end
    end
end

% -------------------------------------------------------------------------
% Figure layout
% -------------------------------------------------------------------------
if three_panel_flag
    figure;
    set(gcf, 'Position', [0 0 1600 650])
else
    figure;
    set(gcf, 'Position', [0 0 750 600])
end

% =========================================================================
% Panel 1 — imagesc of droplet size distribution
% =========================================================================
if three_panel_flag
    ax1 = subplot(1, 3, 1);
else
    ax1 = axes; %#ok<LAXES>
end

imagesc(ax1, r_bins_plot, altitude, Nc_img);
set(ax1, 'YDir', 'normal');     % altitude increases upward
colormap(ax1, 'turbo');

cb1 = colorbar(ax1);
cb1.Label.String     = '$\log_{10}(N_c \ [\mathrm{cm}^{-3}])$';
cb1.Label.Interpreter = 'latex';
cb1.Label.FontSize   = cb_fnt;

xlabel(ax1, 'Droplet radius ($\mu$m)', 'Interpreter', 'latex', 'FontSize', ax_fnt)
ylabel(ax1, 'Altitude (m)',            'Interpreter', 'latex', 'FontSize', ax_fnt)
title(ax1, ['$N_c$ distribution — ', instr_label], ...
    'Interpreter', 'latex', 'FontSize', ttl_fnt)

% overlay weighted mean and +/- 1 sigma lines
hold(ax1, 'on')
plot(ax1, r_mean, altitude, 'k-',  'LineWidth', ln_sz,   'DisplayName', 'Mean radius')
plot(ax1, r_lo,   altitude, 'k--', 'LineWidth', ln_sz-0.5, 'DisplayName', '$\mu \pm \sigma$')
plot(ax1, r_hi,   altitude, 'k--', 'LineWidth', ln_sz-0.5, 'HandleVisibility', 'off')
hold(ax1, 'off')

legend(ax1, 'Interpreter', 'latex', 'Location', 'best', ...
    'Color', [0.15 0.15 0.15], 'TextColor', 'w', 'FontSize', 13)

% =========================================================================
% Panels 2 & 3 (only when three_panel_flag is true)
% =========================================================================
if three_panel_flag

    % --- Panel 2: effective radius vs altitude ----------------------------
    ax2 = subplot(1, 3, 2);

    plot(ax2, re_plot, altitude, '.-', ...
        'Color', [0.2 0.4 0.8], 'LineWidth', ln_sz, 'MarkerSize', 16)
    grid(ax2, 'on'); grid(ax2, 'minor');

    if use_2DC
        re_xlabel = '$r_e$ ($\mu$m)';
    else
        re_xlabel = '$r_e$ ($\mu$m) — CDP only';
    end
    xlabel(ax2, re_xlabel,        'Interpreter', 'latex', 'FontSize', ax_fnt)
    ylabel(ax2, 'Altitude (m)',   'Interpreter', 'latex', 'FontSize', ax_fnt)
    title(ax2,  'Effective radius', 'Interpreter', 'latex', 'FontSize', ttl_fnt)

    % --- Panel 3: effective variance vs altitude --------------------------
    ax3 = subplot(1, 3, 3);

    plot(ax3, v_eff, altitude, '.-', ...
        'Color', [0.8 0.2 0.2], 'LineWidth', ln_sz, 'MarkerSize', 16)
    grid(ax3, 'on'); grid(ax3, 'minor');

    xlabel(ax3, 'Effective variance $v_e$', 'Interpreter', 'latex', 'FontSize', ax_fnt)
    ylabel(ax3, 'Altitude (m)',             'Interpreter', 'latex', 'FontSize', ax_fnt)
    title(ax3,  'Effective variance',       'Interpreter', 'latex', 'FontSize', ttl_fnt)

    % super-title with thresholds
    if isfield(vert_profile, 'LWC_threshold') && isfield(vert_profile, 'Nc_threshold')
        sgtitle(['$LWC \geq$ ', num2str(vert_profile.LWC_threshold), ' $\mathrm{g\ m}^{-3}$', ...
                 '   $N_c \geq$ ', num2str(vert_profile.Nc_threshold), ' $\mathrm{cm}^{-3}$'], ...
            'Interpreter', 'latex', 'FontSize', ttl_fnt)
    end

    % link all three y-axes so zooming is synchronized
    linkaxes([ax1 ax2 ax3], 'y')

end

end
