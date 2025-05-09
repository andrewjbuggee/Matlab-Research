%% Retrieve Droplet profile using a look up table of precomputed radiances and interpolate on a finer grid


% By Anderw John Buggee

%%

clear variables

% Which computer are you working on?

which_computer = whatComputer;

if strcmp(which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------

    folderpath_2save = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/'];



elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    folderpath_2save = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'HySICS/Simulated_spectra/'];


elseif strcmp(which_computer,'curc')==true

    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------

    warning([newline, 'No folder to store things in!', newline])



end

%% LOAD FORWARD MODEL CALCULATIONS OF HYSICS TOA REFLECTANCE


% These are calculations similar to those that would be preformed by my
% retrieval algorithm. Droplet profiles are adiabatic.

forward_model_calcs = load([folderpath_2save, ...
    'forward_model_calcs_forRetrieval_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-01-May-2025_rev1.mat']);


%% 3D Interpolate the radiance calculations on a finer grid and compare with EMIT measurement


[Refl_model_fine, r_top_fine, r_bot_fine, tau_c_fine] = lookUp_table_3D_interpolate_HySICS(forward_model_calcs);



%% Load simulated measurements

simulated_measurements = load([folderpath_2save, ...
    'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-03-May-2025_rev1.mat']);


%% Find the l2-residual between the simulated measurements and the forward modeled reflectances
% If measurement uncertinaty is greater than 0, compute the relative l2-norm 


if simulated_measurements.inputs.measurement.uncert>0

        % Compute the rms difference between the measurements and the modeled
        % reflectances
        l2_residual = sqrt( sum( (repmat(reshape(simulated_measurements.Refl_model_with_noise, 1, 1, 1, []), length(r_top_fine), length(r_bot_fine), length(tau_c_fine))...
            - Refl_model_fine).^2, 4));


        % Compute the RMS of the synthetic measurement uncertainty
        rms_uncert = sqrt( mean( simulated_measurements.Refl_model_uncert.^2));

    else

        %Compute the l2 norm difference between the measurements and the modeled reflectances
        l2_residual = sqrt( sum( (repmat(simulated_measurements.Refl_model, 1, length(r_top_fine), length(r_bot_fine), length(tau_c_fine))...
            - Refl_model_fine).^2, 1));
        % reshape this array
        l2_residual = reshape(l2_residual, length(r_top_fine), length(r_bot_fine), []);


end




%% Find the states with the lowest l2-residual between the simulated measurements and the forward modeled reflectances

% find n smallest rms states
n_states = 50;

% create a mesh grid
[R_bot_fine, R_top_fine, Tau_c_fine] = meshgrid(r_bot_fine, r_top_fine, tau_c_fine);

% store the state values at each minimum
r_top_min = zeros(n_states, 1);
r_bot_min = zeros(n_states, 1);
tau_c_min = zeros(n_states, 1);

% store the rms value and the index
min_val = zeros(n_states, 1);
idx_min = zeros(n_states, 1);


% create a new array where the rms_residual can be used to determine the
% smallest values. We have to insert a nan each time
l2_residual_placeHolder = l2_residual;


for nn = 1:n_states

    % find the smallest rms residual value, omitting nans
    [min_val(nn), idx_min(nn)] = min(l2_residual_placeHolder, [], 'all', 'omitnan');

    r_top_min(nn) = R_top_fine(idx_min(nn));
    r_bot_min(nn) = R_bot_fine(idx_min(nn));
    tau_c_min(nn) = Tau_c_fine(idx_min(nn));

    

    % set the minimum value to nan and omit
    l2_residual_placeHolder(idx_min(nn)) = nan;


end

% for each retrieved state, compute the absolute difference between the retrieved
% state and the true state
% diff_abs_stateVector = abs([r_top_min, r_bot_min, tau_c_min] - ...
%  repmat([simulated_measurements.inputs.RT.r_top, simulated_measurements.inputs.RT.r_bot,simulated_measurements.inputs.RT.tau_c],...
%  n_states, 1));

l2_stateVector = sqrt( sum(([r_top_min, r_bot_min, tau_c_min] - ...
 repmat([simulated_measurements.inputs.RT.r_top, simulated_measurements.inputs.RT.r_bot,simulated_measurements.inputs.RT.tau_c],...
 n_states, 1)).^2, 2));


% get rid of the redundant l2_residual matrix and meshgrid matrices
clear l2_residual_placeHolder R_top_fine R_bot_fine Tau_c_fine

% save the reflectance estimates associated with the minimum rms value
% across all three variables
%min_Refl_model_fine = reshape(Refl_model_fine(r_top_fine==r_top_min(1), r_bot_fine==r_bot_min(1), tau_c_fine==tau_c_min(1),:), [],1);


disp([newline, 'State vector with minimum l2-residual between the associated forward model reflectances',...
    ' and the simulated measurements:', newline, '     Retrieved (LUT): [r_top, r_bot, tau_c] = ',...
    '[', num2str(r_top_min(1)), ',', num2str(r_bot_min(1)), ',', num2str(tau_c_min(1)), ']', newline,...
    '     Truth: [r_top, r_bot, tau_c] = [', num2str(simulated_measurements.inputs.RT.r_top), ',',...
    num2str(simulated_measurements.inputs.RT.r_bot), ',', num2str(simulated_measurements.inputs.RT.tau_c), ']'])


%% Create a 3D plot showing the retrievals with the lowest l2-residual
% Color each marker with the value of the l2-residual


figure; scatter3(r_bot_min, r_top_min, tau_c_min, 30, min_val, 'filled')

colorbar
% Create ylabel
ylabel('$r_{top}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create xlabel
xlabel('$r_{bot}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create xlabel
zlabel('$\tau_{c}$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% set the limits as the full range of the look-up table
ylim([r_top_fine(1), r_top_fine(end)])
xlim([r_bot_fine(1), r_bot_fine(end)])
zlim([tau_c_fine(1), tau_c_fine(end)])

% Create title
title('50 points with lowest rms', 'FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

grid on; grid minor

% set the figure size to be proportional to the length of the r_top and
% r_bot vectors
set(gcf, 'Position', [0 0 900 900])













%%

for nn = 1:1

    clear variables

    % define axes label font size
    axes_label_font_size = 40;

    % define axes tick label font size
    axes_tick_label_font_size = 25;

    % define colorbar font size
    cb_font_size = 40;

    % define contour label size
    contour_label_size = 25;





    % Create mesh grid
    [R_bot_fine, R_top_fine, Tau_c_fine] = meshgrid(r_bot_fine, r_top_fine, tau_c_fine);


    % Define synthetic model data

    % r_top_truth = 10.17;
    % r_bot_truth = 4.74;
    % tau_c_truth = 5.96;


    r_top_truth = 9.17;
    r_bot_truth = 5.74;
    tau_c_truth = 6.512;

    [~, idx_r_top] = min(abs(r_top_fine - r_top_truth));
    [~, idx_r_bot] = min(abs(r_bot_fine - r_bot_truth));
    [~, idx_tau_c] = min(abs(tau_c_fine - tau_c_truth));

    synthetic_measurement = reshape(Refl_model_fine(idx_r_top, idx_r_bot, idx_tau_c, :), [],1);


    % --- Create synthetic measurements with 5% uncertinaty ---
    % ---------------------------------------------------------
    % Add Gaussian Noise to the measurements

    % --- meausrement uncertainty ---
    % define this as a fraction of the measurement
    measurement_uncert_1 = 0.02;

    % Define a gaussian where the mean value is the true measurement, and twice
    % the standard deviation is the product of the measurement uncertainty and
    % the true measurements.
    % Remember: +/- 1*sigma = 68% of the area under the gaussian curve
    %           +/- 2*sigma = 95% of the area under the gaussian curve
    %           +/- 3*sigma = 99.7% of the area under the gaussian curve

    % Compute the new synethtic measurement with gaussian noise
    % *** Gaussian noise can be either positive or negative. Meaning, an
    % uncertainty of 5% implies the true value can lie anywhere between the
    % measured value +/- 5% of the measured value
    % define a
    synthetic_measurement_with_noise = synthetic_measurement + synthetic_measurement.*(measurement_uncert_1/3) .*...
        randn(length(inputs.bands2run), 1);

    % define the synthetic relfectance uncertainty
    synthetic_measurement_uncert = measurement_uncert_1 .* synthetic_measurement_with_noise;






    % Using an exact modeled estimate without noise
    use_l2_norm = true;

    l2_residual = [];
    rms_uncert = [];

    if use_l2_norm==false

        % Compute the rms difference between the measurements and the modeled
        % reflectances
        l2_residual = sqrt( mean( (repmat(reshape(synthetic_measurement_with_noise, 1, 1, 1, []), length(r_top_fine), length(r_bot_fine), length(tau_c_fine))...
            - Refl_model_fine).^2, 4));


        % Compute the RMS of the synthetic measurement uncertainty
        rms_uncert = sqrt( mean( synthetic_measurement_uncert.^2));

    else

        %Compute the l2 norm difference between the measurements and the modeled reflectances
        l2_residual = sqrt( sum( (repmat(reshape(synthetic_measurement_with_noise, 1, 1, 1, []), length(r_top_fine), length(r_bot_fine), length(tau_c_fine))...
            - Refl_model_fine).^2, 4));

        % Compute the l2 norm of the synthetic measurement uncertainty
        rms_uncert = sqrt( sum( synthetic_measurement_uncert.^2));

    end



    % Find the states with the lowest rms residul

    % find the smallest rms residual value, omitting nans
    [~, idx_min] = min(l2_residual, [], 'all', 'omitnan');

    r_top_min = R_top_fine(idx_min);




    % Create Contour plot of rms residual between true EMIT measurements and the libRadTran modeled measurements
    % --- (r_top - r_bot) versus tau  for the minimum r_top ----

    % define the optical depth slice you'd like to plot
    idx_rTop_5pct = r_top_fine == r_top_min(1);

    % Create figure
    figure;

    s1 = subplot(1,2,1);

    % set subplot position
    if strcmp(whatComputer, 'anbu8374')==true
        set(s1, 'Position', [0.0687 0.11 0.41 0.815])

    elseif strcmp(whatComputer, 'andrewbuggee')==true
        set(s1, 'Position', [0.0687 0.11 0.41 0.815])

    end


    % rms residual values to plot
    lvls = [0, 1:24];


    % Create filled contour
    % [c1,h1] = contourf(tau_c_fine, r_top_min(1)-r_bot_fine, reshape(rms_residual(idx_rTop_5pct,:, :)./rms_uncert, length(r_bot_fine),...
    %     length(tau_c_fine)),  lvls, 'LineWidth',4, 'EdgeColor', 'k');
    %clabel(c1,h1,'FontSize', contour_label_size,'FontWeight','bold');

    % Create contour
    [c1,h1] = contour(tau_c_fine, r_top_min(1)-r_bot_fine, reshape(l2_residual(idx_rTop_5pct,:, :)./rms_uncert, length(r_bot_fine),...
        length(tau_c_fine)),  lvls, 'LineWidth',4, 'EdgeColor', mySavedColors(20, 'fixed'));
    clabel(c1,h1,'FontSize', contour_label_size,'FontWeight','bold', 'Color', mySavedColors(20, 'fixed'));



    % Create ylabel
    ylabel('$r_{top}^{*} - r_{bot}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', axes_label_font_size);

    % Create xlabel
    xlabel('$\tau_c$','FontWeight','bold','Interpreter','latex', 'Fontsize', axes_label_font_size);

    % Create title
    if use_l2_norm==false

        title(['RMS Residual at global min $r_{top} = $', num2str(r_top_fine(idx_rTop_5pct)),...
            ' between Synthetic Measurements with ', num2str(100*measurement_uncert_1),...
            '\% uncertainty and LibRadTran'],'Interpreter','latex', 'FontSize', 16);

    else

        title(['RSS Residual at global min $r_{top} = $', num2str(r_top_fine(idx_rTop_5pct)),...
            ' between Synthetic Measurements with ', num2str(100*measurement_uncert_1),...
            '\% uncertainty and LibRadTran'],'Interpreter','latex', 'FontSize', 16);

    end

    idx_uncert = l2_residual./rms_uncert <= 1;
    percent_states_less_than_rms_uncert = sum(idx_uncert, 'all')/numel(l2_residual);
    disp([newline,'Percent of state space within the convergence region using 5% measurement uncertainty: ',...
        num2str(100*percent_states_less_than_rms_uncert),'%', newline])


    ylim([min(r_top_min(1)-r_bot_fine), max(r_top_min(1)-r_bot_fine)])

    grid on; grid minor




    % --- Create synthetic measurements with 1% uncertinaty ---
    % -------------------------------------
    % Add Gaussian Noise to the measurements

    % --- meausrement uncertainty ---
    % define this as a fraction of the measurement
    measurement_uncert_2 = 0.003;

    % Define a gaussian where the mean value is the true measurement, and twice
    % the standard deviation is the product of the measurement uncertainty and
    % the true measurements.
    % Remember: +/- 1*sigma = 68% of the area under the gaussian curve
    %           +/- 2*sigma = 95% of the area under the gaussian curve
    %           +/- 3*sigma = 99.7% of the area under the gaussian curve

    % Compute the new synethtic measurement with gaussian noise
    % *** Gaussian noise can be either positive or negative. Meaning, an
    % uncertainty of 5% implies the true value can lie anywhere between the
    % measured value +/- 5% of the measured value
    % define a
    synthetic_measurement_with_noise = synthetic_measurement + synthetic_measurement.*(measurement_uncert_2/3) .*...
        randn(length(inputs.bands2run), 1);

    % define the synthetic relfectance uncertainty
    synthetic_measurement_uncert = measurement_uncert_2 .* synthetic_measurement_with_noise;

    % Using an exact modeled estimate without noise
    use_l2_norm = true;

    l2_residual = [];
    rms_uncert = [];

    if use_l2_norm==false

        % Compute the rms difference between the measurements and the modeled
        % reflectances
        l2_residual = sqrt( mean( (repmat(reshape(synthetic_measurement_with_noise, 1, 1, 1, []), length(r_top_fine), length(r_bot_fine), length(tau_c_fine))...
            - Refl_model_fine).^2, 4));


        % Compute the RMS of the synthetic measurement uncertainty
        rms_uncert = sqrt( mean( synthetic_measurement_uncert.^2));

    else

        %Compute the l2 norm difference between the measurements and the modeled reflectances
        l2_residual = sqrt( sum( (repmat(reshape(synthetic_measurement_with_noise, 1, 1, 1, []), length(r_top_fine), length(r_bot_fine), length(tau_c_fine))...
            - Refl_model_fine).^2, 4));

        % Compute the l2 norm of the synthetic measurement uncertainty
        rms_uncert = sqrt( sum( synthetic_measurement_uncert.^2));

    end



    % Find the states with the lowest rms residul

    % find the smallest rms residual value, omitting nans
    [~, idx_min] = min(l2_residual, [], 'all', 'omitnan');

    r_top_min = R_top_fine(idx_min);


    % define the optical depth slice you'd like to plot
    idx_rTop_1pct = r_top_fine == r_top_min(1);



    s2 = subplot(1,2,2);

    % set subplot position
    if strcmp(whatComputer, 'anbu8374')==true
        set(s2, 'Position', [0.51 0.11 0.334659090909091 0.815])

    elseif strcmp(whatComputer, 'andrewbuggee')==true
        set(s2, 'Position', [0.56 0.11 0.41 0.815])

    end


    % rms residual values to plot
    lvls = [0, 1:4:32];


    % % Create filled contour
    % [c1,h1] = contourf(tau_c_fine, r_top_min(1)-r_bot_fine, reshape(rms_residual(idx_rTop_1pct,:, :)./rms_uncert, length(r_bot_fine),...
    %     length(tau_c_fine)),  lvls, 'LineWidth',4, 'EdgeColor', 'k');
    % clabel(c1,h1,'FontSize', contour_label_size,'FontWeight','bold');
    %
    % % Create colorbar
    % cb = colorbar();
    % % create colorbar label
    % ylabel(cb, '$\sqrt{ \Sigma{ \left(R(\vec{x}) - \vec{m} \right)^{2} }} / \sqrt{ \Sigma{ \left(\delta \vec{m} \right)^{2}}}$',...
    %     'FontSize', cb_font_size, 'Interpreter', 'latex')
    % clim([lvls(1), lvls(end)])



    % Create contour
    [c1,h1] = contour(tau_c_fine, r_top_min(1)-r_bot_fine, reshape(l2_residual(idx_rTop_1pct,:, :)./rms_uncert, length(r_bot_fine),...
        length(tau_c_fine)),  lvls, 'LineWidth',4, 'EdgeColor', mySavedColors(20, 'fixed'));
    clabel(c1,h1,'FontSize', contour_label_size,'FontWeight','bold', 'Color', mySavedColors(20, 'fixed'));




    % Create ylabel
    ylabel('$r_{top}^{*} - r_{bot}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', axes_label_font_size);

    % Create xlabel
    xlabel('$\tau_c$','FontWeight','bold','Interpreter','latex', 'Fontsize', axes_label_font_size);

    % Create title
    if use_l2_norm==false

        title(['RMS Residual at global min $r_{top} = $', num2str(r_top_fine(idx_rTop_1pct)),...
            ' between Synthetic Measurements with ', num2str(100*measurement_uncert_2),...
            '\% uncertainty and LibRadTran'],'Interpreter','latex', 'FontSize', 16);

    else

        title(['RSS Residual at global min $r_{top} = $', num2str(r_top_fine(idx_rTop_1pct)),...
            ' between Synthetic Measurements with ', num2str(100*measurement_uncert_2),...
            '\% uncertainty and LibRadTran'],'Interpreter','latex', 'FontSize', 16);

    end


    ylim([min(r_top_min(1)-r_bot_fine), max(r_top_min(1)-r_bot_fine)])

    grid on; grid minor



    % set the figure size to be proportional to the length of the r_top and
    % r_bot vectors

    %set the figure size
    if strcmp(whatComputer, 'anbu8374')==true
        set(gcf, 'Position', [0 0 2400 1200])


    elseif strcmp(whatComputer, 'andrewbuggee')==true
        set(gcf, 'Position', [0 0 1500 850])


    end


    idx_uncert = l2_residual./rms_uncert <= 1;
    percent_states_less_than_rms_uncert = sum(idx_uncert, 'all')/numel(l2_residual);
    disp([newline,'Percent of state space within the convergence region using 1% measurement uncertainty: ',...
        num2str(100*percent_states_less_than_rms_uncert),'%', newline])



    clear variables


end


% % ---------- Save figure --------------
% % save .fig file
% folderpath_figs = '/Users/andrewbuggee/Documents/CU-Boulder-ATOC/My Papers/Figures/';
% f = gcf;
% saveas(f,[folderpath_figs,'Fig 8 - relative l2-norm with wavelengths for synthetic data with 2% and 0.3% uncertainty.fig']);
%
%
% % save .png with 400 DPI resolution
% % remove title
% title('')
% folderpath_pngs = '/Users/andrewbuggee/Documents/CU-Boulder-ATOC/My Papers/Submission 1 Figures/';
% exportgraphics(f,[folderpath_pngs,'Fig 8 - relative l2-norm with wavelengths for synthetic data with 2% and 0.3% uncertainty.png'],'Resolution', 400);