%% Plot LWC, effective radius, and droplet number concentration vs optical depth
%  for ORACLES vertical profiles

% Optical depth (tau) is defined as 0 at cloud top and increasing toward
% cloud bottom. Ascending profiles are flipped so that all plots share the
% same orientation (tau increases downward on the y-axis).

% INPUTS:
% -------
%   vert_profiles          - struct array OR cell array of ORACLES vertical profiles
%                            (output of find_verticalProfiles_ORACLES, or ensemble_profiles)
%
%   indices                - vector of indices into vert_profiles to plot
%
%   normalize_opticalDepth - logical; if true, normalize each profile's tau
%                            to [0, 1] so that cloud top = 0 and cloud bottom = 1

% By Andrew John Buggee

%%

function [] = plot_LWC_and_re_and_Nc_vs_opticalDepth_ORACLES(vert_profiles, indices, normalize_opticalDepth)


if nargin ~= 3
    error([newline, 'Wrong number of inputs. Need 3: vertical profiles, indices to plot, ', ...
        'normalization flag', newline])
end


N_curves = length(indices);
legend_str = cell(1, N_curves);


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
    % Orient tau so that tau = 0 is always at cloud top.
    % Stored tau vectors always start at 0 (cloud top for descending,
    % cloud bottom for ascending). Flip ascending profiles.
    % ------------------------------------------------------------------
    is_ascending = (prof.altitude(end) - prof.altitude(1)) > 0;

    if is_ascending
        % Data starts at cloud bottom; flip so index 1 = cloud top
        if prof.tau(1)~=0
            % flip tau if it doesn't start at 0
            tau     = flipud(reshape(prof.tau,      [], 1));
        else
            tau = prof.tau;
        end
        lwc_plt = fliplr(prof.lwc);
        re_plt  = fliplr(prof.re);
        Nc_plt  = fliplr(prof.total_Nc);
    else
        % Data starts at cloud top; no flip needed
        tau     = reshape(prof.tau, [], 1);
        lwc_plt = prof.lwc;
        re_plt  = prof.re;
        Nc_plt  = prof.total_Nc;
    end

    % Fall back to re_CAS if re is all zeros
    if ~any(re_plt > 0)
        if is_ascending
            re_plt = fliplr(prof.re_CAS);
        else
            re_plt = prof.re_CAS;
        end
    end

    % ------------------------------------------------------------------
    % Compute y-axis (tau, raw or normalized)
    % ------------------------------------------------------------------
    if normalize_opticalDepth
        tau_plot = tau ./ tau(end);     % [0, 1]; 0 = cloud top, 1 = cloud bottom
        y_label_str = 'Normalized $\tau$';
    else
        tau_plot = tau;
        y_label_str = '$\tau$';
    end


    % ------------------------------------------------------------------
    % Determine line color: use a fixed palette for <= 6 curves,
    % cycle through all colors for larger sets
    % ------------------------------------------------------------------
    line_color = mySavedColors(nn, 'fixed');


    % ------------------------------------------------------------------
    % Panel 1: LWC
    % ------------------------------------------------------------------
    ax1 = subplot(1,3,1);
    plot(lwc_plt, tau_plot, 'Color', line_color);
    hold on


    % ------------------------------------------------------------------
    % Panel 2: effective radius
    % ------------------------------------------------------------------
    ax2 = subplot(1,3,2);
    plot(re_plt, tau_plot, 'Color', line_color);
    hold on


    % ------------------------------------------------------------------
    % Panel 3: total droplet number concentration
    % ------------------------------------------------------------------
    ax3 = subplot(1,3,3);
    plot(Nc_plt, tau_plot, 'Color', line_color);
    hold on


    legend_str{nn} = ['idx = ', num2str(indices(nn))];

end


% -------------------------------------------------------------------------
% Axis labels, titles, and formatting
% -------------------------------------------------------------------------

subplot(1,3,1)
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('LWC ($g/m^3$)',   'Interpreter', 'latex', 'FontSize', 20)
ylabel(y_label_str,        'Interpreter', 'latex', 'FontSize', 20)


subplot(1,3,2)
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('$r_e$ ($\mu m$)', 'Interpreter', 'latex', 'FontSize', 20)

% Title: show the thresholds used
prof_last = getProf(indices(N_curves));
if isfield(prof_last, 'LWC_threshold')
    title(['$LWC \geq$ ', num2str(prof_last.LWC_threshold), ' $g/m^{3}$', ...
        '   $N_c \geq$ ', num2str(prof_last.Nc_threshold), ' $cm^{-3}$'], ...
        'Interpreter', 'latex', 'FontSize', 16)
end

legend(legend_str, 'Interpreter', 'latex', 'FontSize', 18, 'Location', 'southeast')


subplot(1,3,3)
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('$N_c$ ($cm^{-3}$)', 'Interpreter', 'latex', 'FontSize', 20)


set(gcf, 'Position', [0 0 1200 650])

% Link y-axes so zooming one panel zooms all
linkaxes([ax1, ax2, ax3], 'y')


end
