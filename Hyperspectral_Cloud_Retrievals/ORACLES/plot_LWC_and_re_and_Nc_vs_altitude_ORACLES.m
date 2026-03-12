%% Plot LWC, effective radius, and droplet number concentration vs altitude
%  for ORACLES vertical profiles

% INPUTS:
% -------
%   vert_profiles      - struct array OR cell array of ORACLES vertical profiles
%                        (output of find_verticalProfiles_ORACLES, or ensemble_profiles)
%
%   indices            - vector of indices into vert_profiles to plot
%
%   normalize_altitude - logical; if true, normalize each profile's altitude
%                        to [0, 1] so that cloud bottom = 0 and cloud top = 1

% By Andrew John Buggee

%%

function [] = plot_LWC_and_re_and_Nc_vs_altitude_ORACLES(vert_profiles, indices, normalize_altitude)


if nargin ~= 3
    error([newline, 'Wrong number of inputs. Need 3: vertical profiles, indices to plot, ', ...
        'normalization flag', newline])
end


N_curves = length(indices);
legend_str = cell(1, N_curves);

% Font sizes
ttl_fnt = 20;
ax_fnt  = 30;
lgnd_fnt = 18;

% Starting color index (cycles through mySavedColors palette)
clr_start = 62;


% -------------------------------------------------------------------------
% Helper: extract a single profile regardless of cell or struct array
% -------------------------------------------------------------------------
if iscell(vert_profiles)
    getProf = @(k) vert_profiles{k};
else
    getProf = @(k) vert_profiles(k);
end


figure;

for nn = 1:N_curves

    prof = getProf(indices(nn));

    % ------------------------------------------------------------------
    % Compute y-axis (altitude, raw or normalized)
    % ------------------------------------------------------------------
    if normalize_altitude

        alt_plot = (prof.altitude - min(prof.altitude)) ./ ...
            (max(prof.altitude) - min(prof.altitude));
        y_label_str = 'Normalized Altitude';

    else

        alt_plot = prof.altitude;
        y_label_str = 'Altitude ($m$)';

    end


    % ------------------------------------------------------------------
    % Panel 1: LWC
    % ------------------------------------------------------------------
    ax1 = subplot(1,3,1);
    plot(prof.lwc, alt_plot, 'Color', mySavedColors(clr_start + (nn-1), 'fixed'));
    hold on


    % ------------------------------------------------------------------
    % Panel 2: effective radius
    % Use the pre-computed re (all probes combined).
    % Fall back to re_CAS if re is all zeros (shouldn't happen in ORACLES).
    % ------------------------------------------------------------------
    ax2 = subplot(1,3,2);

    if any(prof.re > 0)
        plot(prof.re, alt_plot, 'Color', mySavedColors(clr_start + (nn-1), 'fixed'));
    else
        plot(prof.re_CAS, alt_plot, 'Color', mySavedColors(clr_start + (nn-1), 'fixed'));
    end
    hold on


    % ------------------------------------------------------------------
    % Panel 3: total droplet number concentration
    % ------------------------------------------------------------------
    ax3 = subplot(1,3,3);
    plot(prof.total_Nc, alt_plot, 'Color', mySavedColors(clr_start + (nn-1), 'fixed'));
    hold on


    legend_str{nn} = ['idx = ', num2str(indices(nn))];

end


% -------------------------------------------------------------------------
% Axis labels, titles, and formatting
% -------------------------------------------------------------------------

subplot(1,3,1)
grid on; grid minor
xlabel('LWC ($g/m^3$)',   'Interpreter', 'latex', 'FontSize', ax_fnt)
ylabel(y_label_str,        'Interpreter', 'latex', 'FontSize', ax_fnt)


subplot(1,3,2)
grid on; grid minor
xlabel('$r_e$ ($\mu m$)', 'Interpreter', 'latex', 'FontSize', ax_fnt)

% Title: show the thresholds used
prof_last = getProf(indices(N_curves));
if isfield(prof_last, 'LWC_threshold')
    title(['$LWC \geq$ ', num2str(prof_last.LWC_threshold), ' $g/m^{3}$', ...
        '   $N_c \geq$ ', num2str(prof_last.Nc_threshold), ' $cm^{-3}$'], ...
        'Interpreter', 'latex', 'FontSize', ttl_fnt)
end


subplot(1,3,3)
grid on; grid minor
xlabel('$N_c$ ($cm^{-3}$)', 'Interpreter', 'latex', 'FontSize', ax_fnt)
legend(legend_str, 'Interpreter', 'latex', 'Location', 'best', ...
    'Color', 'white', 'TextColor', 'k', 'FontSize', lgnd_fnt)


set(gcf, 'Position', [0 0 1200 625])

% Link y-axes so zooming one panel zooms all
linkaxes([ax1, ax2, ax3], 'y')


end
