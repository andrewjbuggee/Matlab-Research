%% Calculate ensemble statistics of vertical profiles from ORACLES in-situ measurements

% Loops through all ORACLES combined_microphysics netCDF files, finds vertical
% cloud profiles in each flight, and builds an ensemble. Computes statistics
% of effective radius, liquid water content, and number concentration as a
% function of normalized optical depth.
%
% Analogous to ensemble_vertical_statistics.m (VOCALS-REx version), adapted
% for the ORACLES P-3 dataset.

% By Andrew John Buggee

clear variables

%% -------------------------------------------------------------------------
%  -------------------- File locations ------------------------------------
%  -------------------------------------------------------------------------

which_computer = whatComputer;

if strcmp(which_computer, 'anbu8374')

    foldername_data = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/', ...
        'Hyperspectral_Cloud_Retrievals/ORACLES/oracles_data/combined_microphysics/'];

    foldername_2save = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/', ...
        'Hyperspectral_Cloud_Retrievals/ORACLES/oracles_data/'];

elseif strcmp(which_computer, 'andrewbuggee')

    foldername_data = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/', ...
        'Hyperspectral_Cloud_Retrievals/ORACLES/oracles_data/combined_microphysics/'];

end


% Find all combined_microphysics netCDF files
folder_contents = dir(foldername_data);

filename = cell(1, length(folder_contents));
backspace = 0;

for nn = 1:length(folder_contents)

    if length(folder_contents(nn).name) > 3 && ...
            strcmp(folder_contents(nn).name(1:14), 'Microphysics_P')

        filename{nn - backspace} = folder_contents(nn).name;

    else
        filename(nn - backspace) = [];
        backspace = backspace + 1;
    end

end


%% -------------------------------------------------------------------------
%  ------------ Define profile quality thresholds -------------------------
%  -------------------------------------------------------------------------

% LWC threshold: minimum liquid water content to define cloud boundaries
% The ORACLES dataset uses 0.05 g/m^3 as its cloud definition threshold.
inputs.LWC_threshold = 0.05;       % g/m^3

% Stop each profile at the maximum LWC value?
inputs.stop_at_max_LWC = false;

% Minimum total number concentration threshold
% The ORACLES dataset uses 10 cm^-3 as its cloud definition threshold.
inputs.Nc_threshold = 10;          % cm^-3

% Sort profiles by drizzle/precipitation?
inputs.sort_for_precip_driz = true;

% Precipitation is defined by the 2DS+HVPS-3 liquid water path (D > 50 µm)
% Keep profiles WITH precipitation/drizzle (true) or WITHOUT (false)?
inputs.keep_precip_drizzle_profiles = true;

% Rain water path threshold for precipitation classification
inputs.precip_driz_threshold = 5;  % g/m^2


%% -------------------------------------------------------------------------
%  ------------ Loop through all files, find vertical profiles -------------
%  -------------------------------------------------------------------------

ensemble_profiles = cell([]);

for nn = 1:length(filename)

    disp(['Processing file ', num2str(nn), ' of ', num2str(length(filename)), ...
        ': ', filename{nn}])

    % Read the ORACLES netCDF file
    oracles = readORACLES([foldername_data, filename{nn}]);

    % Find vertical cloud profiles in this flight
    vert_profs = find_verticalProfiles_ORACLES(oracles, inputs.LWC_threshold, ...
        inputs.stop_at_max_LWC, inputs.Nc_threshold, which_computer);

    if isempty(vert_profs)
        disp(['  No vertical profiles found in ', filename{nn}])
        continue
    end

    disp(['  Found ', num2str(length(vert_profs)), ' vertical profile(s)'])

    if inputs.sort_for_precip_driz

        % Sort profiles into precipitating and non-precipitating
        indexes_with_precip = sort_vert_profs_for_precipitation_ORACLES( ...
            vert_profs, inputs.precip_driz_threshold);

        if inputs.keep_precip_drizzle_profiles
            % Keep all profiles (precipitating and non-precipitating)
            indexes_2keep = 1:length(vert_profs);
        else
            % Keep only non-precipitating profiles (filter out drizzly ones)
            indexes_2keep = setxor((1:length(vert_profs)), indexes_with_precip);
        end

    else

        indexes_2keep = 1:length(vert_profs);

    end

    % Append selected profiles to ensemble
    l = length(ensemble_profiles);

    for mm = 1:length(indexes_2keep)

        if l == 0 && mm == 1
            ensemble_profiles{1} = vert_profs(indexes_2keep(mm));
        else
            ensemble_profiles{l + mm} = vert_profs(indexes_2keep(mm));
        end

    end


end

disp(['Total ensemble profiles collected: ', num2str(length(ensemble_profiles))])


%% -------------------------------------------------------------------------
%  ---------------------- Save ensemble profiles --------------------------
%  -------------------------------------------------------------------------

if inputs.sort_for_precip_driz

    if inputs.keep_precip_drizzle_profiles

        save([foldername_2save, 'ensemble_profiles_with_precip_from_', ...
            num2str(length(filename)), '_files_LWC-threshold_', ...
            num2str(inputs.LWC_threshold), '_Nc-threshold_', ...
            num2str(inputs.Nc_threshold), '_', char(datetime("today")), '.mat'], ...
            'ensemble_profiles', 'filename', 'inputs')

    else

        save([foldername_2save, 'ensemble_profiles_without_precip_from_', ...
            num2str(length(filename)), '_files_LWC-threshold_', ...
            num2str(inputs.LWC_threshold), '_Nc-threshold_', ...
            num2str(inputs.Nc_threshold), '_drizzleLWP-threshold_', ...
            num2str(inputs.precip_driz_threshold), '_', char(datetime("today")), '.mat'], ...
            'ensemble_profiles', 'filename', 'inputs')

    end

else

    save([foldername_2save, 'ensemble_profiles_from_', ...
        num2str(length(filename)), '_files_LWC-threshold_', ...
        num2str(inputs.LWC_threshold), '_Nc-threshold_', ...
        num2str(inputs.Nc_threshold), '_', char(datetime("today")), '.mat'], ...
        'ensemble_profiles', 'filename', 'inputs')

end


%% -------------------------------------------------------------------------
%  ----------- Plot histogram of optical depth and geometric thickness ----
%  -------------------------------------------------------------------------

cloud_optical_depth      = zeros(1, length(ensemble_profiles));
cloud_geometric_thickness = zeros(1, length(ensemble_profiles));

for nn = 1:length(ensemble_profiles)
    cloud_optical_depth(nn)       = ensemble_profiles{nn}.tau(end);
    cloud_geometric_thickness(nn) = abs(ensemble_profiles{nn}.altitude(1) - ...
        ensemble_profiles{nn}.altitude(end));
end

figure;

subplot(1,2,1)
histogram(cloud_geometric_thickness, 'NumBins', 20)
xlabel('Cloud Geometric Thickness (m)')
ylabel('Counts')
title([num2str(length(ensemble_profiles)), ' vertical profiles'])

subplot(1,2,2)
histogram(cloud_optical_depth, 'NumBins', 20)
xlabel('Cloud Optical Thickness')
ylabel('Counts')

set(gcf, 'Position', [0 0 1000 550])


%% -------------------------------------------------------------------------
%  ---- Bin all profiles into N vertical segments (normalized optical depth)
%  -------------------------------------------------------------------------

n_bins    = 30;
bin_edges = 0:1/n_bins:1;

% Cell array: rows = tau bins, cols = [re, lwc, Nc]
vertically_segmented_attributes = cell(n_bins, 3);

normalized_tau = cell(1, length(ensemble_profiles));

for nn = 1:length(ensemble_profiles)

    % Normalize optical depth to [0, 1]
    normalized_tau{nn} = ensemble_profiles{nn}.tau ./ ensemble_profiles{nn}.tau(end);

    % Orient data in optical-depth space (tau=0 at cloud top, tau=1 at bottom).
    % For ascending profiles, flip variables so index 1 = cloud top.
    dz_dt_prof = mean(diff(ensemble_profiles{nn}.altitude) ./ diff(ensemble_profiles{nn}.time));

    if dz_dt_prof > 0
        % Ascending: data starts at cloud bottom; flip to get cloud-top first
        re  = fliplr(ensemble_profiles{nn}.re);
        lwc = fliplr(ensemble_profiles{nn}.lwc);
        Nc  = fliplr(ensemble_profiles{nn}.total_Nc);
    else
        % Descending: data starts at cloud top; no flip needed
        re  = ensemble_profiles{nn}.re;
        lwc = ensemble_profiles{nn}.lwc;
        Nc  = ensemble_profiles{nn}.total_Nc;
    end

    % Bin variables by normalized optical depth
    for bb = 1:length(bin_edges) - 1

        if bb == 1
            index_segment = normalized_tau{nn} >= bin_edges(bb) & ...
                normalized_tau{nn} <= bin_edges(bb+1);
        else
            index_segment = normalized_tau{nn} > bin_edges(bb) & ...
                normalized_tau{nn} <= bin_edges(bb+1);
        end

        vertically_segmented_attributes{bb, 1} = ...
            [vertically_segmented_attributes{bb, 1}, re(index_segment)];

        vertically_segmented_attributes{bb, 2} = ...
            [vertically_segmented_attributes{bb, 2}, lwc(index_segment)];

        vertically_segmented_attributes{bb, 3} = ...
            [vertically_segmented_attributes{bb, 3}, Nc(index_segment)];

    end

end


%% -------------------------------------------------------------------------
%  ----------- Compute mean, median, std of re at each vertical bin -------
%  -------------------------------------------------------------------------

re_mean   = zeros(n_bins, 1);
re_median = zeros(n_bins, 1);
re_std    = zeros(n_bins, 1);
bin_center = zeros(n_bins, 1);

for bb = 1:n_bins
    re_mean(bb)   = mean(vertically_segmented_attributes{bb, 1});
    re_median(bb) = median(vertically_segmented_attributes{bb, 1});
    re_std(bb)    = std(vertically_segmented_attributes{bb, 1});
    bin_center(bb) = (bin_edges(bb+1) - bin_edges(bb))/2 + bin_edges(bb);
end


%% -------------------------------------------------------------------------
%  ----------- Compute mean, median, std of LWC at each vertical bin ------
%  -------------------------------------------------------------------------

lwc_mean   = zeros(n_bins, 1);
lwc_median = zeros(n_bins, 1);
lwc_std    = zeros(n_bins, 1);

for bb = 1:n_bins
    lwc_mean(bb)   = mean(vertically_segmented_attributes{bb, 2});
    lwc_median(bb) = median(vertically_segmented_attributes{bb, 2});
    lwc_std(bb)    = std(vertically_segmented_attributes{bb, 2});
end


%% -------------------------------------------------------------------------
%  ----------- Compute mean, median, std of Nc at each vertical bin -------
%  -------------------------------------------------------------------------

Nc_mean   = zeros(n_bins, 1);
Nc_median = zeros(n_bins, 1);
Nc_std    = zeros(n_bins, 1);

for bb = 1:n_bins
    Nc_mean(bb)   = mean(vertically_segmented_attributes{bb, 3});
    Nc_median(bb) = median(vertically_segmented_attributes{bb, 3});
    Nc_std(bb)    = std(vertically_segmented_attributes{bb, 3});
end


%% -------------------------------------------------------------------------
%  ---- Plot: median re, LWC, and Nc vs normalized optical depth (subplot) -
%  -------------------------------------------------------------------------

figure;

% Effective radius
subplot(1,3,1)
x = [re_median - re_std; flipud(re_median + re_std)];
y = [bin_center; flipud(bin_center)];
fill(x, y, mySavedColors(2,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
hold on
plot(re_median, bin_center, 'Color', mySavedColors(2,'fixed'))
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('$\langle r_e(z) \rangle$ $(\mu m)$', 'Interpreter', 'latex')
ylabel('Normalized Optical Depth')

% LWC
subplot(1,3,2)
x = [lwc_median - lwc_std; flipud(lwc_median + lwc_std)];
y = [bin_center; flipud(bin_center)];
fill(x, y, mySavedColors(2,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
hold on
plot(lwc_median, bin_center, 'Color', mySavedColors(2,'fixed'))
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('$\langle LWC(z) \rangle$ $(g/m^{3})$', 'Interpreter', 'latex')
ylabel('Normalized Optical Depth')
title(['Median profiles  |  LWC $\geq$ ', num2str(inputs.LWC_threshold), ' $g/m^3$  |  ', ...
    '$N_c \geq$ ', num2str(inputs.Nc_threshold), ' $cm^{-3}$'], 'Interpreter', 'latex')

% Number concentration
subplot(1,3,3)
x = [Nc_median - Nc_std; flipud(Nc_median + Nc_std)];
y = [bin_center; flipud(bin_center)];
fill(x, y, mySavedColors(2,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
hold on
plot(Nc_median, bin_center, 'Color', mySavedColors(2,'fixed'))
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('$\langle N_c(z) \rangle$ $(cm^{-3})$', 'Interpreter', 'latex')
ylabel('Normalized Optical Depth')

set(gcf, 'Position', [0 0 1200 625])


%% -------------------------------------------------------------------------
%  ---- Plot: mean and median re profiles with shaded standard deviation --
%  -------------------------------------------------------------------------

figure;

x = [re_mean - re_std; flipud(re_mean + re_std)];
y = [bin_center; flipud(bin_center)];
fill(x, y, mySavedColors(1,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
hold on
plot(re_mean, bin_center, 'Color', mySavedColors(1,'fixed'));

x = [re_median - re_std; flipud(re_median + re_std)];
y = [bin_center; flipud(bin_center)];
fill(x, y, mySavedColors(2,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
hold on
plot(re_median, bin_center, 'Color', mySavedColors(2,'fixed'))

legend('','Mean', '', 'Median', 'Location', 'best')
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('$\langle r_e(z) \rangle$', 'Interpreter', 'latex')
ylabel('Normalized Optical Depth')
set(gcf, 'Position', [0 0 650 600])


%% -------------------------------------------------------------------------
%  ---- Plot: mean re profile with adiabatic fit comparison ---------------
%  -------------------------------------------------------------------------

figure;

x = [re_mean - re_std; flipud(re_mean + re_std)];
y = [bin_center; flipud(bin_center)];
fill(x, y, mySavedColors(1,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
hold on
plot(re_mean, bin_center, 'Color', mySavedColors(1,'fixed'));

profile_type = 'subadiabatic_aloft';
re_fit = create_droplet_profile2([re_mean(1), re_mean(end)], bin_center, 'optical_depth', profile_type);
hold on
plot(re_fit, bin_center, 'Color', 'k', 'LineWidth', 1);

legend('','', 'subadiabatic', 'location', 'best')
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('$\langle r_e \rangle$', 'Interpreter', 'latex')
ylabel('Normalized Optical Depth')
title('Mean $r_e$ — ORACLES', 'Interpreter', 'latex')
set(gcf, 'Position', [0 0 650 600])


%% -------------------------------------------------------------------------
%  ---- Plot: mean and median LWC profiles ---------------------------------
%  -------------------------------------------------------------------------

figure;

x = [lwc_mean - lwc_std; flipud(lwc_mean + lwc_std)];
y = [bin_center; flipud(bin_center)];
fill(x, y, mySavedColors(6,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
hold on
plot(lwc_mean, bin_center, 'Color', mySavedColors(6,'fixed'));

x = [lwc_median - lwc_std; flipud(lwc_median + lwc_std)];
y = [bin_center; flipud(bin_center)];
fill(x, y, mySavedColors(4,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
hold on
plot(lwc_median, bin_center, 'Color', mySavedColors(4,'fixed'))

legend('','Mean', '', 'Median', 'Location', 'best')
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('$\langle LWC(z) \rangle$', 'Interpreter', 'latex')
ylabel('Normalized Optical Depth')
set(gcf, 'Position', [0 0 650 600])


%% -------------------------------------------------------------------------
%  ---- Plot: mean and median Nc profiles ----------------------------------
%  -------------------------------------------------------------------------

figure;

x = [Nc_mean - Nc_std; flipud(Nc_mean + Nc_std)];
y = [bin_center; flipud(bin_center)];
fill(x, y, mySavedColors(1,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
hold on
plot(Nc_mean, bin_center, 'Color', mySavedColors(1,'fixed'));

x = [Nc_median - Nc_std; flipud(Nc_median + Nc_std)];
y = [bin_center; flipud(bin_center)];
fill(x, y, mySavedColors(2,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
hold on
plot(Nc_median, bin_center, 'Color', mySavedColors(2,'fixed'))

legend('','Mean', '', 'Median', 'Location', 'best')
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('$\langle N_c(z) \rangle$ $(cm^{-3})$', 'Interpreter', 'latex')
ylabel('Normalized Optical Depth')
set(gcf, 'Position', [0 0 650 600])


%% -------------------------------------------------------------------------
%  ---- Segment by LWP regime ----------------------------------------------
%  -------------------------------------------------------------------------

LWP = zeros(1, length(ensemble_profiles));
for nn = 1:length(ensemble_profiles)
    dz_dt_prof = mean(diff(ensemble_profiles{nn}.altitude) ./ diff(ensemble_profiles{nn}.time));
    if dz_dt_prof > 0
        LWP(nn) = trapz(ensemble_profiles{nn}.altitude, ensemble_profiles{nn}.lwc);
    else
        LWP(nn) = -1 * trapz(ensemble_profiles{nn}.altitude, ensemble_profiles{nn}.lwc);
    end
end

figure; plot(LWP, '.'); grid on; grid minor
xlabel('Profile Index')
ylabel('LWP (g/m^2)')
title('Liquid Water Path — all ORACLES profiles')


% Define LWP bin edges for regime analysis
LWP_bin_edges = [0, 20, 40, 62, 100, 140];     % g/m^2

index_regime = false(length(LWP_bin_edges), length(LWP));
legend_str   = cell(1, length(LWP_bin_edges));

for nn = 1:length(LWP_bin_edges)

    if nn >= 1 && nn < length(LWP_bin_edges)
        index_regime(nn,:) = LWP >= LWP_bin_edges(nn) & LWP < LWP_bin_edges(nn+1);
        legend_str{nn} = [num2str(LWP_bin_edges(nn)), '$\leq LWP<$', num2str(LWP_bin_edges(nn+1))];
    else
        index_regime(nn,:) = LWP > LWP_bin_edges(nn);
        legend_str{nn} = ['$LWP>$', num2str(LWP_bin_edges(nn)), ' $g/m^{2}$'];
    end

end

changing_variable = 'LWP';
plot_multiple_median_ensemble_LWC_re_NC_vs_norm_optical_depth( ...
    ensemble_profiles, index_regime, changing_variable, legend_str)


%% -------------------------------------------------------------------------
%  ---- Segment by cloud optical depth regime ------------------------------
%  -------------------------------------------------------------------------

tau_total = zeros(1, length(ensemble_profiles));
for nn = 1:length(ensemble_profiles)
    tau_total(nn) = ensemble_profiles{nn}.tau(end);
end

tau_bin_edges = [0, 5, 9, 14, 21, 30];

index_regime = false(length(tau_bin_edges), length(tau_total));
legend_str   = cell(1, length(tau_bin_edges));

for nn = 1:length(tau_bin_edges)

    if nn >= 1 && nn < length(tau_bin_edges)
        index_regime(nn,:) = tau_total >= tau_bin_edges(nn) & tau_total < tau_bin_edges(nn+1);
        legend_str{nn} = [num2str(tau_bin_edges(nn)), '$\leq \tau <$', num2str(tau_bin_edges(nn+1))];
    else
        index_regime(nn,:) = tau_total > tau_bin_edges(nn);
        legend_str{nn} = ['$\tau >$', num2str(tau_bin_edges(nn))];
    end

end

changing_variable = '$\tau$';
plot_multiple_median_ensemble_LWC_re_NC_vs_norm_optical_depth( ...
    ensemble_profiles, index_regime, changing_variable, legend_str)
