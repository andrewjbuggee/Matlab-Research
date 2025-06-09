%% Create new versions of figures 7 and 8 for paper 1

% First, we need to interpolate on the forward model calculations to create
% a fine synthetic grid of reflectances

clear variables


% 3D Interpolate the radiance calculations on a finer grid and compare with EMIT measurement
which_computer = whatComputer();

% ------------ LOAD DATA SET -----------------------
if strcmp(which_computer, 'anbu8374')==true

    foldername = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/';

    % tauC = [5.5 : 0.5 : 7.5]
    filename = 'forward_model_calcs_forRetrieval_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-14-May-2025_rev1.mat';

    % tauC = [7.5 : 0.5 : 15]
    %     filename = 'forward_model_calcs_forRetrieval_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-19-May-2025_rev1.mat';
    % Load forward model cals over wide range of r_top, r_bot and tau
    fm = load([foldername,filename]);


elseif strcmp(which_computer, 'andrewbuggee')==true

    foldername = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/';

    % tauC = [5.5 : 0.5 : 7.5]
    %filename = 'forward_model_calcs_forRetrieval_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-14-May-2025_rev1.mat';

    % -----------------------
    % tauC = [7.5 : 0.5 : 15]
    % filename = 'forward_model_calcs_forRetrieval_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-19-May-2025_rev1.mat';

    % simulated fig 3.a on paper 1
    filename = 'forward_model_calcs_forRetrieval_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-05-Jun-2025_rev1.mat';
    % -----------------------

    % Load forward model cals over wide range of r_top, r_bot and tau
    fm = load([foldername,filename]);


end
% -----------------------------------------------------

% ---- Unpack the forward model simulated reflectances ---
num_rTop = length(fm.inputs.RT.r_top);
num_rBot = length(fm.inputs.RT.r_bot);
num_tauC = length(fm.inputs.RT.tau_c);
num_wl = length(fm.inputs.bands2run);
% Reshape the forward model into a 4-D array
% dim 1 = tau_c (rows)
% dim 2 = re_bot (cols)
% dim 3 = re_top (depth)
% dim 4 = wl (time?)
fm_refl = zeros(num_tauC, num_rBot, num_rTop, num_wl);
for nn = 1:num_wl

    fm_refl(:,:,:,nn) = reshape(fm.Refl_model(nn:num_wl:end), num_tauC, num_rBot, num_rTop);

end
% -----------------------------------------------------


% Meshgrid is defined on x,y,z space, not row, column, depth space
% In 3D space, z = depth, x = column, y = row
[R_bot, Tau_c, R_top] = meshgrid(fm.inputs.RT.r_bot, fm.inputs.RT.tau_c, fm.inputs.RT.r_top);

% Create the new fine grid to interpolate on
% define the discrete step length of each variable
d_r_top = 0.05;      % microns
d_r_bot = 0.05;      % microns
d_tau_c = 0.05;

r_top_fine = fm.inputs.RT.r_top(1):d_r_top:fm.inputs.RT.r_top(end);
r_bot_fine = fm.inputs.RT.r_bot(1):d_r_bot:fm.inputs.RT.r_bot(end);
tau_c_fine = fm.inputs.RT.tau_c(1):d_tau_c:fm.inputs.RT.tau_c(end);

% clear forward model structure to save memory
clear fm


[R_bot_fine, Tau_c_fine, R_top_fine] = meshgrid(r_bot_fine, tau_c_fine, r_top_fine);

Refl_model_fine = zeros(length(tau_c_fine), length(r_bot_fine), length(r_top_fine), num_wl);


tic
for wl = 1:size(fm_refl,4)

    Refl_model_fine(:,:,:,wl) = interp3(R_bot, Tau_c, R_top, fm_refl(:, :, :, wl),...
        R_bot_fine, Tau_c_fine, R_top_fine);


end
toc

% Save Refl_model_file and the rms_residual, because these calculations
% take a while!
save_filename = ['Interpolated_',filename];
save([foldername, save_filename],"Refl_model_fine", "r_top_fine", "r_bot_fine",...
    "tau_c_fine", "R_top_fine", "R_bot_fine", "Tau_c_fine", '-v7.3');


clear R_bot_fine R_top_fine Tau_c_fine Refl_model_fine fm_refl

%% Find the states with the lowest rms residul


% ------------ LOAD SIMULATED MEASUREMENT ------------------
if strcmp(which_computer, 'anbu8374')==true

    foldername = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/';

    %filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-14-May-2025_rev1.mat';
    filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-17-May-2025_rev1.mat';


    % Load forward model cals over wide range of r_top, r_bot and tau
    sim = load([foldername,filename]);


elseif strcmp(which_computer, 'andrewbuggee')==true

    foldername = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/';

    % filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-19-May-2025_rev1.mat';

    % simulated fig 3.a on paper 1
    filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-05-Jun-2025_rev1.mat';


    % Load forward model cals over wide range of r_top, r_bot and tau
    sim = load([foldername,filename]);

end
% -----------------------------------------------------


% find n smallest rms states
n_states = 20;

% store the state values at each minimum
r_top_min_1 = zeros(n_states, 1);
r_bot_min_1 = zeros(n_states, 1);
tau_c_min_1 = zeros(n_states, 1);

% store the rms value and the index
min_val = zeros(n_states, 1);
idx_min = zeros(n_states, 1);


rss_residual = sqrt( sum( (repmat(reshape(sim.Refl_model, 1, 1, 1, []),...
    length(tau_c_fine), length(r_bot_fine), length(r_top_fine))...
    - Refl_model_fine).^2, 4));


for nn = 1:n_states

    % find the smallest rms residual value, omitting nans
    [min_val(nn), idx_min(nn)] = min(rss_residual, [], 'all', 'omitnan');

    r_top_min_1(nn) = R_top_fine(idx_min(nn));
    r_bot_min_1(nn) = R_bot_fine(idx_min(nn));
    tau_c_min_1(nn) = Tau_c_fine(idx_min(nn));

    % set the minimum value to nan and omit
    rss_residual(idx_min(nn)) = nan;


end


% Save Refl_model_file and the rms_residual, because these calculations
% take a while!
save([foldername, save_filename],"r_top_min_1", "tau_c_min_1", "r_bot_min_1", '-append');







%% How does the LUT solution compare between a MODIS-like uncertainty system and a HySICS like uncertainty system?
% --- NO FIGURE! ---

clear variables


which_computer = whatComputer();

% ------------ LOAD DATA SET -----------------------
if strcmp(which_computer, 'anbu8374')==true

    foldername = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/';

    % tauC = [5.5 : 0.5 : 7.5]
    %filename = 'forward_model_calcs_forRetrieval_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-14-May-2025_rev1.mat';

    % tauC = [7.5 : 0.5 : 15]
    filename = 'forward_model_calcs_forRetrieval_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-19-May-2025_rev1.mat';
    % Load forward model cals over wide range of r_top, r_bot and tau
    fm = load([foldername,filename]);



    % Load a simulated measurement
    % r_top = 12, r_bot = 4, tau = 6
    sim_filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-14-May-2025_rev1.mat';

    % r_top = 8.5, r_bot = 6, tau_c = 5.9
    %sim_filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-17-May-2025_rev1.mat';

    % r_top = 8.5, r_bot = 6, tau_c = 5.9
    %sim_filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-18-May-2025_rev1.mat';

    % r_top = 8.5, r_bot = 6, tau_c = 5.9
    %sim_filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-18-May-2025_rev2.mat';

    % r_top = 12.5, r_bot = 4.3, tau_c = 6.9
    %sim_filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-18-May-2025_rev3.mat';


    sim = load([foldername, sim_filename]);

elseif strcmp(which_computer, 'andrewbuggee')==true


    % --- Load forward Model ---

    foldername = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/';

    % tauC = [5.5 : 0.5 : 7.5]
    %filename = 'forward_model_calcs_forRetrieval_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-14-May-2025_rev1.mat';
    filename = 'Interpolated_forward_model_calcs_forRetrieval_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-19-May-2025_rev1.mat';

    % tauC = [7.5 : 0.5 : 15]
    % filename = 'forward_model_calcs_forRetrieval_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-19-May-2025_rev1.mat';

    % Load forward model cals over wide range of r_top, r_bot and tau
    fm = load([foldername,filename]);


    % --- Load simulated measurement ---

    filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-19-May-2025_rev1.mat';


    % Load forward model cals over wide range of r_top, r_bot and tau
    sim = load([foldername,filename]);

end
% -----------------------------------------------------


% --- unpack the synthetic measurement ---
synthetic_measurement = sim.Refl_model;
r_top_true = sim.inputs.RT.r_top;
r_bot_true = sim.inputs.RT.r_bot;
tau_c_true = sim.inputs.RT.tau_c;

% delete the rest
clear sim


% ----- unpack the forward model array -----
fm_refl = fm.Refl_model_fine;
tau_c_fine = fm.tau_c_fine;
r_bot_fine = fm.r_bot_fine;
r_top_fine = fm.r_top_fine;
R_top_fine = fm.R_top_fine;
R_bot_fine = fm.R_bot_fine;
Tau_c_fine = fm.Tau_c_fine;

clear fm

% --- Create synthetic measurements with 3% uncertinaty ---
% ---------------------------------------------------------
% Add Gaussian Noise to the measurements

% --- Total uncertainty ---
% define this as a fraction
% define the total uncertainty (measurement + forward model) for a MODIS
% like instrument
measurement_uncert_1 = 0.03;

% define the smaller total uncertainty (measurement + foward model) to
% compare with
measurement_uncert_2 = 0.01;

% Define a gaussian where the mean value is the true measurement, and twice
% the standard deviation is the product of the measurement uncertainty and
% the true measurements.
% Remember: +/- 1*sigma = 68% of the area under the gaussian curve
%           +/- 2*sigma = 95% of the area under the gaussian curve
%           +/- 3*sigma = 99.7% of the area under the gaussian curve


% -----------------------------------------------------
% ----- Set up flags and empty arrays -----

% Using an exact modeled estimate without noise
num_iterations = 1000;
r_top_min_1 = zeros(num_iterations, 1);
r_bot_min_1 = zeros(num_iterations, 1);
tau_c_min_1 = zeros(num_iterations, 1);

r_top_min_2 = zeros(num_iterations, 1);
r_bot_min_2 = zeros(num_iterations, 1);
tau_c_min_2 = zeros(num_iterations, 1);

% -----------------------------------------------------



tic
parfor nn = 1:num_iterations


    rss_residual_1 = [];


    % Compute the new synethtic measurement with gaussian noise
    % *** Gaussian noise can be either positive or negative. Meaning, an
    % uncertainty of 5% implies the true value can lie anywhere between the
    % measured value +/- 5% of the measured value
    % define a
    synthetic_measurement_with_noise = synthetic_measurement + synthetic_measurement.*(measurement_uncert_1/3) .*...
        randn(length(synthetic_measurement), 1);




    %Compute the l2 norm difference between the measurements and the modeled reflectances
    rss_residual_1 = sqrt( sum( (repmat(reshape(synthetic_measurement_with_noise, 1, 1, 1, []),...
        length(tau_c_fine), length(r_bot_fine), length(r_top_fine))...
        - fm_refl).^2, 4));




    % Find the states with the lowest rms residul

    % find the smallest rms residual value, omitting nans
    [~, idx_min] = min(rss_residual_1, [], 'all', 'omitnan');

    r_top_min_1(nn) = R_top_fine(idx_min);
    r_bot_min_1(nn) = R_bot_fine(idx_min);
    tau_c_min_1(nn) = Tau_c_fine(idx_min);







    % ----------------------------------------------------------------
    % ***--- Now compute the synthetic measurement for a more accurate
    % system ---***
    % ----------------------------------------------------------------
    % Compute the new synethtic measurement with gaussian noise
    synthetic_measurement_with_noise_2 = synthetic_measurement + synthetic_measurement.*(measurement_uncert_2/3) .*...
        randn(length(synthetic_measurement), 1);





    rss_residual_2 = [];

    %Compute the l2 norm difference between the measurements and the modeled reflectances
    rss_residual_2 = sqrt( sum( (repmat(reshape(synthetic_measurement_with_noise_2, 1, 1, 1, []),...
        length(tau_c_fine), length(r_bot_fine), length(r_top_fine))...
        - fm_refl).^2, 4));




    % Find the states with the lowest rms residul

    % find the smallest rms residual value, omitting nans
    [~, idx_min] = min(rss_residual_2, [], 'all', 'omitnan');

    r_top_min_2(nn) = R_top_fine(idx_min);
    r_bot_min_2(nn) = R_bot_fine(idx_min);
    tau_c_min_2(nn) = Tau_c_fine(idx_min);




end
toc


% Compute the root-mean-square percent different between the two
% uncertainty scenarios and the true state vector

rms_state_vec_1 = 100.* sqrt( mean( ((repmat([r_top_true, r_bot_true, tau_c_true], num_iterations, 1) - ...
    [r_top_min_1, r_bot_min_1, tau_c_min_1])./...
    repmat([r_top_true, r_bot_true, tau_c_true], num_iterations, 1)).^2, 1));       % rms % difference for each variable

rms_state_vec_2 = 100.* sqrt( mean( ((repmat([r_top_true, r_bot_true, tau_c_true], num_iterations, 1) - ...
    [r_top_min_2, r_bot_min_2, tau_c_min_2])./...
    repmat([r_top_true, r_bot_true, tau_c_true], num_iterations, 1)).^2, 1));       % rms % difference for each variable





%% NEW FIGURE 8 (now figure 9) --- Subplot comparing l2 norm residual for synthetic HySICS data bewteen different measurement uncertainty scenarios

clear variables


which_computer = whatComputer();

% ------------ LOAD DATA SET -----------------------
if strcmp(which_computer, 'anbu8374')==true

    foldername = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/';

    filename = 'interpolated_forward_model_calcs_forRetrieval_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-14-May-2025_rev1.mat';

    % simulated
    % Load forward model cals over wide range of r_top, r_bot and tau
    fm = load([foldername,filename]);


    % ----------------------------
    % Load a simulated measurement
    % ----------------------------

    % r_top = 12, r_bot = 4, tau = 6
    %     sim_filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-14-May-2025_rev1.mat';

    % r_top = 8.5, r_bot = 6, tau_c = 5.9
    %sim_filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-17-May-2025_rev1.mat';

    % r_top = 9.5, r_bot = 4, tau_c = 6
    sim_filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-02-Jun-2025_rev6.mat';     % old rayliegh scattering model + adjusted CO2 column amount + surface albedo=0.04 + removed day of year + 10 layers instead of 250


    % r_top = 8.5, r_bot = 6, tau_c = 7.3
    % sim_filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-18-May-2025_rev1.mat';

    % r_top = 8.5, r_bot = 6, tau_c = 5.6
    % sim_filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-18-May-2025_rev2.mat';

    % r_top = 12.5, r_bot = 4.3, tau_c = 6.9
    %sim_filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-18-May-2025_rev4.mat';

    sim = load([foldername, sim_filename]);
    % ----------------------------
    % ----------------------------


elseif strcmp(which_computer, 'andrewbuggee')==true

    foldername = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/';


    % ----------------------------
    %    Load LUT calculations
    % ----------------------------
    % solar and viewing zenith angle of 0
    % filename = 'interpolated_forward_model_calcs_forRetrieval_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-14-May-2025_rev1.mat';


    % simulated calcs for MODIS obs on fig 3.a for paper 1
    filename = 'Interpolated_forward_model_calcs_forRetrieval_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-05-Jun-2025_rev1.mat';


    % Load forward model cals over wide range of r_top, r_bot and tau
    fm = load([foldername,filename]);


    % ----------------------------
    % Load a simulated measurement
    % ----------------------------

    % r_top = 9.5, r_bot = 4, tau_c = 6
    sim_filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-05-Jun-2025_rev1';     %


    sim = load([foldername, sim_filename]);



end
% -----------------------------------------------------


% --- unpack the synthetic measurement ---
synthetic_measurement = sim.Refl_model;
r_top_true = sim.inputs.RT.r_top;
r_bot_true = sim.inputs.RT.r_bot;
tau_c_true = sim.inputs.RT.tau_c;

% delete the rest
clear sim


% ----- unpack the forward model array -----
fm_refl = fm.Refl_model_fine;
tau_c_fine = fm.tau_c_fine;
r_bot_fine = fm.r_bot_fine;
r_top_fine = fm.r_top_fine;
R_top_fine = fm.R_top_fine;
R_bot_fine = fm.R_bot_fine;
Tau_c_fine = fm.Tau_c_fine;

clear fm

% --- Create synthetic measurements with 5% uncertinaty ---
% ---------------------------------------------------------
% Add Gaussian Noise to the measurements

% --- Total uncertainty ---
% define this as a fraction
% define the total uncertainty (measurement + forward model) for a MODIS
% like instrument
measurement_uncert_1 = 0.03;

% define the smaller total uncertainty (measurement + foward model) to
% compare with
measurement_uncert_2 = 0.01;

% Define a gaussian where the mean value is the true measurement, and twice
% the standard deviation is the product of the measurement uncertainty and
% the true measurements.
% Remember: +/- 1*sigma = 68% of the area under the gaussian curve
%           +/- 2*sigma = 95% of the area under the gaussian curve
%           +/- 3*sigma = 99.7% of the area under the gaussian curve





rss_residual = [];
rss_uncert = [];

% define axes label font size
axes_label_font_size = 40;

% define axes tick label font size
axes_tick_label_font_size = 25;

% define colorbar font size
cb_font_size = 40;

% define contour label size
contour_label_size = 25;






% define the synthetic relfectance uncertainty
synthetic_measurement_uncert = measurement_uncert_1 .* synthetic_measurement;





%Compute the l2 norm difference between the measurements and the modeled reflectances
rss_residual = sqrt( sum( (repmat(reshape(synthetic_measurement, 1, 1, 1, []),...
    length(tau_c_fine), length(r_bot_fine), length(r_top_fine))...
    - fm_refl).^2, 4));

% Compute the l2 norm of the synthetic measurement uncertainty
rss_uncert = sqrt( sum( synthetic_measurement_uncert.^2));





% Find the states with the lowest rms residul

% find the smallest rss residual value, omitting nans
[~, idx_min] = min(rss_residual, [], 'all', 'omitnan');

r_top_min_1 = R_top_fine(idx_min);
r_bot_min_1 = R_bot_fine(idx_min);
tau_c_min_1 = Tau_c_fine(idx_min);



% Create Contour plot of rms residual between simulated HySICS measurements and the libRadTran modeled measurements
% --- (r_top - r_bot) versus tau  for the minimum r_top ----

% define the optical depth slice you'd like to plot
idx_rTop_1 = r_top_fine == r_top_min_1(1);

% Create figure
figure;

s1 = subplot(1,2,1);

% set subplot position
if strcmp(which_computer, 'anbu8374')==true
    set(s1, 'Position', [0.0687 0.11 0.41 0.815])

elseif strcmp(which_computer, 'andrewbuggee')==true
    set(s1, 'Position', [0.0687 0.11 0.41 0.815])

end




% rms residual values to plot
lvls = [0, 1:24];


% Create filled contour
% [c1,h1] = contourf(tau_c_fine, r_top_min(1)-r_bot_fine, reshape(rms_residual(idx_rTop_5pct,:, :)./rms_uncert, length(r_bot_fine),...
%     length(tau_c_fine)),  lvls, 'LineWidth',4, 'EdgeColor', 'k');
%clabel(c1,h1,'FontSize', contour_label_size,'FontWeight','bold');

% Create contour
[c1,h1] = contour(tau_c_fine, r_top_min_1(1)-r_bot_fine, (rss_residual(:,:, idx_rTop_1)./rss_uncert)',...
    lvls, 'LineWidth',4, 'EdgeColor', mySavedColors(63, 'fixed'));
clabel(c1,h1,'FontSize', contour_label_size,'FontWeight','bold', 'Color', mySavedColors(63, 'fixed'));


% ---- Plot the true state vector ----
hold on
plot(tau_c_true, (r_top_true - r_bot_true), 'x', 'MarkerSize', 12, 'Color', ...
    mySavedColors(62, 'fixed'));

% ---- Plot the state vector associated with the Global Minimum RSS ----
hold on
plot(tau_c_min_1, (r_top_min_1 - r_bot_min_1), 'x', 'MarkerSize', 12, 'Color', ...
    mySavedColors(61, 'fixed'));



% Create ylabel
ylabel('$r_{top}^{*} - r_{bot}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', axes_label_font_size);

% Create xlabel
xlabel('$\tau_c$','FontWeight','bold','Interpreter','latex', 'Fontsize', axes_label_font_size);

% Create title
title(['RSS Residual at global min $r_{top} = $', num2str(r_top_fine(idx_rTop_1)),...
    ' between Synthetic Measurements with ', num2str(100*measurement_uncert_1),...
    '\% total uncertainty and LibRadTran'],'Interpreter','latex', 'FontSize', 16);



idx_uncert = rss_residual./rss_uncert <= 1;
percent_states_less_than_rms_uncert = sum(idx_uncert, 'all')/numel(rss_residual);
disp([newline,'Percent of state space within the convergence region using 5% measurement uncertainty: ',...
    num2str(100*percent_states_less_than_rms_uncert),'%', newline])


ylim([min(r_top_min_1(1)-r_bot_fine), max(r_top_min_1(1)-r_bot_fine)])

grid on; grid minor





% ----------------------------------------------------------------
% ***--- Now compute the synthetic measurement for a more accurate
% system ---***
% ----------------------------------------------------------------

% define the synthetic relfectance uncertainty
synthetic_measurement_uncert_2 = measurement_uncert_2 .* synthetic_measurement;



rss_residual = [];
rss_uncert = [];

%Compute the l2 norm difference between the measurements and the modeled reflectances
rss_residual = sqrt( sum( (repmat(reshape(synthetic_measurement, 1, 1, 1, []),...
    length(tau_c_fine), length(r_bot_fine), length(r_top_fine))...
    - fm_refl).^2, 4));

% Compute the l2 norm of the synthetic measurement uncertainty
rss_uncert = sqrt( sum( synthetic_measurement_uncert_2.^2));


% Find the states with the lowest rms residul

% find the smallest rms residual value, omitting nans
[~, idx_min] = min(rss_residual, [], 'all', 'omitnan');

r_top_min_2 = R_top_fine(idx_min);
r_bot_min_2 = R_bot_fine(idx_min);
tau_c_min_2 = Tau_c_fine(idx_min);


ylim([min(r_top_min_2(1)-r_bot_fine), max(r_top_min_2(1)-r_bot_fine)])


% define the optical depth slice you'd like to plot
idx_rTop_2 = r_top_fine == r_top_min_2(1);


% Get current y-axis limits
y_limits = ylim;
x_limits = xlim;

% Create shading for negative y-values
if y_limits(1) < 0
    area_x = [x_limits(1), x_limits(2), x_limits(2), x_limits(1)];
    area_y = [y_limits(1), y_limits(1), 0, 0];
    fill(area_x, area_y, mySavedColors(64, 'fixed'), 'EdgeColor', 'none', 'FaceAlpha', 0.15);
end





s2 = subplot(1,2,2);

% set subplot position
if strcmp(whatComputer, 'anbu8374')==true
    set(s2, 'Position', [0.53 0.11 0.41 0.815])

elseif strcmp(whatComputer, 'andrewbuggee')==true
    set(s2, 'Position', [0.56 0.11 0.41 0.815])

end


% rms residual values to plot
lvls = [0, 1:3:32];


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
[c1,h1] = contour(tau_c_fine, r_top_min_2(1)-r_bot_fine, (rss_residual(:,:, idx_rTop_2)./rss_uncert)',...
    lvls, 'LineWidth',4, 'EdgeColor', mySavedColors(63, 'fixed'));
clabel(c1,h1,'FontSize', contour_label_size,'FontWeight','bold', 'Color', mySavedColors(63, 'fixed'));



% ----- Plot the true state vector -----
hold on
plot(tau_c_true, (r_top_true - r_bot_true), 'x', 'MarkerSize', 12, 'Color', ...
    mySavedColors(62, 'fixed'));

% ---- Plot the state vector associated with the Global Minimum RSS ----
hold on
plot(tau_c_min_2, (r_top_min_2 - r_bot_min_2), 'x', 'MarkerSize', 12, 'Color', ...
    mySavedColors(61, 'fixed'));




% Create ylabel
ylabel('$r_{top}^{*} - r_{bot}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', axes_label_font_size);

% Create xlabel
xlabel('$\tau_c$','FontWeight','bold','Interpreter','latex', 'Fontsize', axes_label_font_size);

% Create title

title(['RSS Residual at global min $r_{top} = $', num2str(r_top_fine(idx_rTop_2)),...
    ' between Synthetic Measurements with ', num2str(100*measurement_uncert_2),...
    '\% total uncertainty and LibRadTran'],'Interpreter','latex', 'FontSize', 16);




ylim([min(r_top_min_2(1)-r_bot_fine), max(r_top_min_2(1)-r_bot_fine)])

grid on; grid minor


% Get current y-axis limits
y_limits = ylim;
x_limits = xlim;

% Create shading for negative y-values
if y_limits(1) < 0
    area_x = [x_limits(1), x_limits(2), x_limits(2), x_limits(1)];
    area_y = [y_limits(1), y_limits(1), 0, 0];
    fill(area_x, area_y, mySavedColors(64, 'fixed'), 'EdgeColor', 'none', 'FaceAlpha', 0.15);
end







% set the figure size to be proportional to the length of the r_top and
% r_bot vectors

%set the figure size
if strcmp(whatComputer, 'anbu8374')==true
    set(gcf, 'Position', [0 0 2400 1200])


elseif strcmp(whatComputer, 'andrewbuggee')==true
    set(gcf, 'Position', [0 0 1500 850])


end


idx_uncert = rss_residual./rss_uncert <= 1;
percent_states_less_than_rms_uncert = sum(idx_uncert, 'all')/numel(rss_residual);
disp([newline,'Percent of state space within the convergence region using 1% measurement uncertainty: ',...
    num2str(100*percent_states_less_than_rms_uncert),'%', newline])



clear variables



% ---------- Save figure --------------
% % save .fig file
% if strcmp(whatComputer, 'andrewbuggee')==true
%     folderpath_figs = '/Users/andrewbuggee/Documents/CU-Boulder-ATOC/My Papers/Paper 1/Figures/';
% elseif strcmp(whatComputer, 'anbu8374')==true
%     folderpath_figs = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/first_paper/figures_post_submission/';
% end
%
% f = gcf;
% saveas(f,[folderpath_figs,'Fig 9 - relative l2-norm with wavelengths for synthetic data with 3% and 1% uncertainty - ver 2.fig']);
%
%
% % save .png with 400 DPI resolution
% % remove title
% title('')
% if strcmp(whatComputer, 'andrewbuggee')==true
%     folderpath_pngs = '/Users/andrewbuggee/Documents/CU-Boulder-ATOC/My Papers/Paper 1/Submission 2 Figures/';
% elseif strcmp(whatComputer, 'anbu8374')==true
%     folderpath_pngs = '/Users/anbu8374/Documents/My Papers/Paper 1/Submission 2 Figures/';
% end
%
% exportgraphics(f,[folderpath_pngs,'Fig 9 - relative l2-norm with wavelengths for synthetic data with 3% and 1% uncertainty - ver 2.png'],'Resolution', 400);







%% How does the LUT solution compare between using 7 and 35 wavelengths
% --- NO FIGURE! ---

clear variables


which_computer = whatComputer();

% ------------ LOAD DATA SET -----------------------
if strcmp(which_computer, 'anbu8374')==true

    foldername = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/';

    filename = 'interpolated_forward_model_calcs_forRetrieval_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-14-May-2025_rev1.mat';

    % Load forward model cals over wide range of r_top, r_bot and tau
    fm = load([foldername,filename]);



    % Load a simulated measurement
    % r_top = 12, r_bot = 4, tau = 6
    sim_filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-14-May-2025_rev1.mat';

    % r_top = 8.5, r_bot = 6, tau_c = 5.9
    %sim_filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-17-May-2025_rev1.mat';

    % r_top = 8.5, r_bot = 6, tau_c = 5.9
    %sim_filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-18-May-2025_rev1.mat';

    % r_top = 8.5, r_bot = 6, tau_c = 5.9
    %sim_filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-18-May-2025_rev2.mat';

    % r_top = 12.5, r_bot = 4.3, tau_c = 6.9
    %sim_filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-18-May-2025_rev3.mat';


    sim = load([foldername, sim_filename]);

elseif strcmp(which_computer, 'andrewbuggee')==true

    error([newline, 'Where are these files?!', newline])

end
% -----------------------------------------------------

% Define 7 MODIS channels
% Let's now seperate out the interpolated relfectance at the seven MODIS
% wavelengths
wl_MODIS7_idx = [1, 3, 5, 6, 18, 22, 29];

% --- unpack the synthetic measurement ---
synthetic_measurement = sim.Refl_model;
synthetic_measurement_MODIS7 = synthetic_measurement(wl_MODIS7_idx);
r_top_true = sim.inputs.RT.r_top;
r_bot_true = sim.inputs.RT.r_bot;
tau_c_true = sim.inputs.RT.tau_c;

% delete the rest
clear sim


% ----- unpack the forward model array -----
fm_refl = fm.Refl_model_fine;
% Grab the modeled data at just the 7 MODIS bands
fm_refl_MODIS7 = fm_refl(:,:,:, wl_MODIS7_idx);
tau_c_fine = fm.tau_c_fine;
r_bot_fine = fm.r_bot_fine;
r_top_fine = fm.r_top_fine;
R_top_fine = fm.R_top_fine;
R_bot_fine = fm.R_bot_fine;
Tau_c_fine = fm.Tau_c_fine;

clear fm

% --- Create synthetic measurements using 7 MODIS Channels with 3% uncertinaty ---
% --------------------------------------------------------------------------------
% Add Gaussian Noise to the measurements

% --- meausrement uncertainty ---
% define this as a fraction
% This is the TOTAL uncertainty: radiometric + forward model
measurement_uncert_MODIS7 = 0.03;

% Define a gaussian where the mean value is the true measurement, and twice
% the standard deviation is the product of the measurement uncertainty and
% the true measurements.
% Remember: +/- 1*sigma = 68% of the area under the gaussian curve
%           +/- 2*sigma = 95% of the area under the gaussian curve
%           +/- 3*sigma = 99.7% of the area under the gaussian curve


% -----------------------------------------------------
% ----- Set up flags and empty arrays -----

% Using an exact modeled estimate without noise
num_iterations = 1000;
r_top_min_MODIS7 = zeros(num_iterations, 1);
r_bot_min_MODIS7 = zeros(num_iterations, 1);
tau_c_min_MODIS7 = zeros(num_iterations, 1);

r_top_min = zeros(num_iterations, 1);
r_bot_min = zeros(num_iterations, 1);
tau_c_min = zeros(num_iterations, 1);

% -----------------------------------------------------

tic
parfor nn = 1:num_iterations



    rss_residual_MODIS7 = [];



    % Compute the new synethtic measurement with gaussian noise
    % *** Gaussian noise can be either positive or negative. Meaning, an
    % uncertainty of 5% implies the true value can lie anywhere between the
    % measured value +/- 5% of the measured value
    % define a
    synthetic_measurement_with_noise_MODIS7 = synthetic_measurement_MODIS7 + synthetic_measurement_MODIS7.*(measurement_uncert_MODIS7/3) .*...
        randn(length(synthetic_measurement_MODIS7), 1);



    %Compute the l2 norm difference between the measurements and the modeled reflectances
    rss_residual_MODIS7 = sqrt( sum( (repmat(reshape(synthetic_measurement_with_noise_MODIS7, 1, 1, 1, []),...
        length(tau_c_fine), length(r_bot_fine), length(r_top_fine))...
        - fm_refl_MODIS7).^2, 4));



    % Find the states with the lowest rms residul

    % find the smallest rms residual value, omitting nans
    [~, idx_min] = min(rss_residual_MODIS7, [], 'all', 'omitnan');

    r_top_min_MODIS7(nn) = R_top_fine(idx_min);
    r_bot_min_MODIS7(nn) = R_bot_fine(idx_min);
    tau_c_min_MODIS7(nn) = Tau_c_fine(idx_min);






    % now define a new measurement using all 35 wavelengths
    synthetic_measurement_with_noise = synthetic_measurement + synthetic_measurement.*(measurement_uncert_MODIS7/3) .*...
        randn(length(synthetic_measurement), 1);



    rss_residual = [];

    %Compute the l2 norm difference between the measurements and the modeled reflectances
    rss_residual = sqrt( sum( (repmat(reshape(synthetic_measurement_with_noise, 1, 1, 1, []),...
        length(tau_c_fine), length(r_bot_fine), length(r_top_fine))...
        - fm_refl).^2, 4));




    % Find the states with the lowest rms residul

    % find the smallest rms residual value, omitting nans
    [~, idx_min] = min(rss_residual, [], 'all', 'omitnan');

    r_top_min(nn) = R_top_fine(idx_min);
    r_bot_min(nn) = R_bot_fine(idx_min);
    tau_c_min(nn) = Tau_c_fine(idx_min);





end
toc


% Compute the root-mean-square percent different between the two
% uncertainty scenarios and the true state vector

rms_state_vec_MODIS7 = 100.* sqrt( mean( ((repmat([r_top_true, r_bot_true, tau_c_true], num_iterations, 1) - ...
    [r_top_min_MODIS7, r_bot_min_MODIS7, tau_c_min_MODIS7])./...
    repmat([r_top_true, r_bot_true, tau_c_true], num_iterations, 1)).^2, 1));       % rms % difference for each variable

rms_state_vec_35 = 100.* sqrt( mean( ((repmat([r_top_true, r_bot_true, tau_c_true], num_iterations, 1) - ...
    [r_top_min, r_bot_min, tau_c_min])./...
    repmat([r_top_true, r_bot_true, tau_c_true], num_iterations, 1)).^2, 1));       % rms % difference for each variable





%% **NEW** FIGURE 7 (now Figure 8) --- Subplot comparing l2 norm residual for synthetic data bewteen 7 and 35 wavelengths employed in the retrieval

clear variables


which_computer = whatComputer();

% ------------ LOAD DATA SET -----------------------
if strcmp(which_computer, 'anbu8374')==true

    foldername = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/';

    filename = 'interpolated_forward_model_calcs_forRetrieval_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-14-May-2025_rev1.mat';

    % Load forward model cals over wide range of r_top, r_bot and tau
    fm = load([foldername,filename]);



    % Load a simulated measurement
    % r_top = 12, r_bot = 4, tau = 6
    %     sim_filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-14-May-2025_rev1.mat';

    % r_top = 8.5, r_bot = 6, tau_c = 5.9
    %     sim_filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-17-May-2025_rev1.mat';

    % r_top = 9.5, r_bot = 4, tau_c = 6
    %     sim_filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-02-Jun-2025_rev1.mat';
    %      sim_filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-02-Jun-2025_rev2.mat';   % old rayliegh scattering model
    %      sim_filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-02-Jun-2025_rev3.mat';   % old rayliegh scattering model + adjusted CO2 column amount
    %      sim_filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-02-Jun-2025_rev4.mat';   % old rayliegh scattering model + adjusted CO2 column amount + surface albedo=0.04
    %      sim_filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-02-Jun-2025_rev5.mat';     % old rayliegh scattering model + adjusted CO2 column amount + surface albedo=0.04 + removed day of year
    sim_filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-02-Jun-2025_rev6.mat';     % old rayliegh scattering model + adjusted CO2 column amount + surface albedo=0.04 + removed day of year + 10 layers instead of 250



    % r_top = 8.5, r_bot = 6, tau_c = 7.3
    %     sim_filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-18-May-2025_rev1.mat';

    % r_top = 8.5, r_bot = 6, tau_c = 5.6
    % sim_filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-18-May-2025_rev2.mat';

    % r_top = 12.5, r_bot = 4.3, tau_c = 6.9
    %sim_filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-18-May-2025_rev4.mat';

    sim = load([foldername, sim_filename]);

elseif strcmp(which_computer, 'andrewbuggee')==true


    foldername = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/';


    % ----------------------------
    %    Load LUT calculations
    % ----------------------------
    % solar and viewing zenith angle of 0
    % filename = 'interpolated_forward_model_calcs_forRetrieval_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-14-May-2025_rev1.mat';


    % simulated calcs for MODIS obs on fig 3.a for paper 1
    filename = 'Interpolated_forward_model_calcs_forRetrieval_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-05-Jun-2025_rev1.mat';


    % Load forward model cals over wide range of r_top, r_bot and tau
    fm = load([foldername,filename]);


    % ----------------------------
    % Load a simulated measurement
    % ----------------------------

    % r_top = 9.5, r_bot = 4, tau_c = 6
    % simulated calcs for MODIS obs on fig 3.a for paper 1
    sim_filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-05-Jun-2025_rev1';     %


    sim = load([foldername, sim_filename]);


end
% -----------------------------------------------------


% Define 7 MODIS channels
% Let's now seperate out the interpolated relfectance at the seven MODIS
% wavelengths
wl_MODIS7_idx = [1, 3, 5, 6, 18, 22, 29];

% --- unpack the synthetic measurement ---
synthetic_measurement = sim.Refl_model;
synthetic_measurement_MODIS7 = synthetic_measurement(wl_MODIS7_idx);
r_top_true = sim.inputs.RT.r_top;
r_bot_true = sim.inputs.RT.r_bot;
tau_c_true = sim.inputs.RT.tau_c;

% delete the rest
clear sim


% ----- unpack the forward model array -----
fm_refl = fm.Refl_model_fine;
% Grab the modeled data at just the 7 MODIS bands
fm_refl_MODIS7 = fm_refl(:,:,:, wl_MODIS7_idx);
tau_c_fine = fm.tau_c_fine;
r_bot_fine = fm.r_bot_fine;
r_top_fine = fm.r_top_fine;
R_top_fine = fm.R_top_fine;
R_bot_fine = fm.R_bot_fine;
Tau_c_fine = fm.Tau_c_fine;

clear fm

% --- Create synthetic measurements using 7 MODIS Channels with 5% uncertinaty ---
% --------------------------------------------------------------------------------
% Add Gaussian Noise to the measurements

% --- meausrement uncertainty ---
% define this as a fraction of the measurement
measurement_uncert_MODIS7 = 0.03;

% Define a gaussian where the mean value is the true measurement, and twice
% the standard deviation is the product of the measurement uncertainty and
% the true measurements.
% Remember: +/- 1*sigma = 68% of the area under the gaussian curve
%           +/- 2*sigma = 95% of the area under the gaussian curve
%           +/- 3*sigma = 99.7% of the area under the gaussian curve






% define axes label font size
axes_label_font_size = 40;

% define axes tick label font size
axes_tick_label_font_size = 25;

% define colorbar font size
cb_font_size = 40;

% define contour label size
contour_label_size = 25;



% define the synthetic relfectance uncertainty
synthetic_measurement_uncert_MODIS7 = measurement_uncert_MODIS7 .* synthetic_measurement_MODIS7;         % 1/sr




% Compute the l2 norm difference between the measurements and the modeled reflectances
rss_residual_MODIS7 = sqrt( sum( (repmat(reshape(synthetic_measurement_MODIS7, 1, 1, 1, []),...
    length(tau_c_fine), length(r_bot_fine), length(r_top_fine))...
    - fm_refl_MODIS7).^2, 4));

% Compute the l2 norm of the synthetic measurement uncertainty
rss_uncert_MODIS7 = sqrt( sum( synthetic_measurement_uncert_MODIS7.^2));


% Find the states with the lowest rms residul for the 7 MODIS wavelengths
% find the smallest rms residual value, omitting nans
[~, idx_min_MODIS7] = min(rss_residual_MODIS7, [], 'all', 'omitnan');

r_top_min_MODIS7 = R_top_fine(idx_min_MODIS7);
r_bot_min_MODIS7 = R_bot_fine(idx_min_MODIS7);
tau_c_min_MODIS7 = Tau_c_fine(idx_min_MODIS7);


% ------------- PLOT l2 norm using 7 MODIS wavelengths ------------
% Create Contour plot of rms residual between 7 MODIS wavelengths and
% Synthetic measurements
% --- (r_top - r_bot) versus tau  for the minimum r_top ----

% define the slice at some cloud top radii you'd like to plot
idx_rTop_MODIS7 = r_top_fine == r_top_min_MODIS7(1);

% Create figure
f = figure;

s1 = subplot(1,2,1);

% set subplot position
if strcmp(whatComputer, 'anbu8374')==true
    set(s1, 'Position', [0.0687 0.11 0.41 0.815])

elseif strcmp(whatComputer, 'andrewbuggee')==true
    set(s1, 'Position', [0.0687 0.11 0.41 0.815])

end


% rms residual values to plot
lvls = [0, 1:1:11];
%lvls = [0, 1:10];

% Create contour plot
[c1,h1] = contour(tau_c_fine, r_top_min_MODIS7(1)-r_bot_fine, (rss_residual_MODIS7(:,:, idx_rTop_MODIS7)...
    ./rss_uncert_MODIS7)',lvls, 'LineWidth',4, 'EdgeColor', mySavedColors(63, 'fixed'));
clabel(c1,h1,'FontSize',contour_label_size,'FontWeight','bold', 'Color', mySavedColors(63, 'fixed'));


% Create filled contour plot
% [c1,h1] = contourf(tau_c_fine, r_top_min_MODIS7(1)-r_bot_fine, reshape(rms_residual_MODIS7(idx_rTop_MODIS7,:, :)./rms_uncert_MODIS7, length(r_bot_fine),...
%     length(tau_c_fine)),  lvls, 'LineWidth', 3, 'EdgeColor', 'k');
%clabel(c1,h1,'FontSize',contour_label_size,'FontWeight','bold', 'Color', mySavedColors(9, 'fixed'));


% Set tick label font size
ax = gca(f);
ax.FontSize = axes_tick_label_font_size;

% Create ylabel
ylabel('$r_{top}^{*} - r_{bot}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', axes_label_font_size);

% Create xlabel
xlabel('$\tau_c$','FontWeight','bold','Interpreter','latex', 'Fontsize', axes_label_font_size);


title(['RSS Residual at $r_{top} = $', num2str(r_top_fine(idx_rTop_MODIS7)),...
    ' - 7 Synthetic Measurements with ', num2str(100*measurement_uncert_MODIS7),...
    '\% uncertainty'],'Interpreter','latex', 'FontSize', 16);

grid on; grid minor




% ----- Plot the true state vector -----
hold on
plot(tau_c_true, (r_top_true - r_bot_true), 'x', 'MarkerSize', 12, 'Color', ...
    mySavedColors(62, 'fixed'));

% ---- Plot the state vector associated with the Global Minimum RSS ----
hold on
plot(tau_c_min_MODIS7, (r_top_min_MODIS7 - r_bot_min_MODIS7), 'x', 'MarkerSize', 12, 'Color', ...
    mySavedColors(61, 'fixed'));




ylim([min(r_top_min_MODIS7(1)-r_bot_fine), max(r_top_min_MODIS7(1)-r_bot_fine)])


% Get current y-axis limits
y_limits = ylim;
x_limits = xlim;

% Create shading for negative y-values
if y_limits(1) < 0
    area_x = [x_limits(1), x_limits(2), x_limits(2), x_limits(1)];
    area_y = [y_limits(1), y_limits(1), 0, 0];
    fill(area_x, area_y, mySavedColors(64, 'fixed'), 'EdgeColor', 'none', 'FaceAlpha', 0.15);
end






% --- Create synthetic measurements using 35 MODIS Channels with 3% uncertinaty ---
% --------------------------------------------------------------------------------

% define the synthetic relfectance uncertainty
synthetic_measurement_uncert_35 = measurement_uncert_MODIS7 .* synthetic_measurement;

% --------------------------------------------------------------------------------

%Compute the l2 norm difference between the measurements and the modeled reflectances
rss_residual_35 = sqrt( sum( (repmat(reshape(synthetic_measurement, 1, 1, 1, []),...
    length(tau_c_fine), length(r_bot_fine), length(r_top_fine))...
    - fm_refl).^2, 4));


% Compute the l2 norm of the synthetic measurement uncertainty
rss_uncert_35 = sqrt( sum( synthetic_measurement_uncert_35.^2));


% Find the states with the lowest rms residul
% find the smallest rms residual value, omitting nans
[~, idx_min] = min(rss_residual_35, [], 'all', 'omitnan');

r_top_min_35 = R_top_fine(idx_min);
r_bot_min_35 = R_bot_fine(idx_min);
tau_c_min_35 = Tau_c_fine(idx_min);





% ------------- PLOT l2 norm using 35 MODIS wavelengths ------------
% Create Contour plot of rms residual between 35 MODIS wavelengths and
% Synthetic measurements
% --- (r_top - r_bot) versus tau  for the minimum r_top ----
% define the slice at some cloud top radii you'd like to plot
idx_rTop_35 = r_top_fine == r_top_min_35(1);

s2 = subplot(1,2,2);

% set subplot position
if strcmp(whatComputer, 'anbu8374')==true
    set(s2, 'Position', [0.53 0.11 0.41 0.815])

elseif strcmp(whatComputer, 'andrewbuggee')==true
    set(s2, 'Position', [0.56 0.11 0.41 0.815])

end


% rms residual values to plot
lvls = [0, 1:1:11];
%lvls = [0, 1:10];


% Create contour plot
[c1,h1] = contour(tau_c_fine, r_top_min_35(1)-r_bot_fine, (rss_residual_35(:,:, idx_rTop_35)...
    ./rss_uncert_35)',lvls, 'LineWidth',4, 'EdgeColor', mySavedColors(63, 'fixed'));
clabel(c1,h1,'FontSize', contour_label_size,'FontWeight','bold', 'Color', mySavedColors(63, 'fixed'));


% Create filled contour
% [c1,h1] = contourf(tau_c_fine, r_top_min(1)-r_bot_fine, reshape(rms_residual(idx_rTop,:, :)./rms_uncert, length(r_bot_fine),...
%     length(tau_c_fine)),  lvls, 'LineWidth', 3, 'EdgeColor', 'k');
%clabel(c1,h1,'FontSize', contour_label_size,'FontWeight','bold', 'Color', mySavedColors(9, 'fixed'));

% Create colorbar
% cb = colorbar();
% % create colorbar label
% ylabel(cb, '$\sqrt{ \Sigma{ \left(R(\vec{x}) - \vec{m} \right)^{2} }} / \sqrt{ \Sigma{ \left(\delta \vec{m} \right)^{2}}}$',...
%     'FontSize', cb_font_size, 'Interpreter', 'latex')
% clim([lvls(1), lvls(end)])




% Set tick label font size
ax = gca(f);
ax.FontSize = axes_tick_label_font_size;

% Create ylabel
ylabel('$r_{top}^{*} - r_{bot}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', axes_label_font_size);

% Create xlabel
xlabel('$\tau_c$','FontWeight','bold','Interpreter','latex', 'Fontsize', axes_label_font_size);





% ----- Plot the true state vector -----
hold on
plot(tau_c_true, (r_top_true - r_bot_true), 'x', 'MarkerSize', 12, 'Color', ...
    mySavedColors(62, 'fixed'));

% ---- Plot the state vector associated with the Global Minimum RSS ----
hold on
plot(tau_c_min_35, (r_top_min_35 - r_bot_min_35), 'x', 'MarkerSize', 12, 'Color', ...
    mySavedColors(61, 'fixed'));




title(['RSS Residual at $r_{top} = $', num2str(r_top_fine(idx_rTop_35)),...
    ' - 35 Synthetic Measurements with ', num2str(100*measurement_uncert_MODIS7),...
    '\% uncertainty'],'Interpreter','latex', 'FontSize', 16);


ylim([min(r_top_min_35(1)-r_bot_fine), max(r_top_min_35(1)-r_bot_fine)])

grid on; grid minor


% Get current y-axis limits
y_limits = ylim;
x_limits = xlim;

% Create shading for negative y-values
if y_limits(1) < 0
    area_x = [x_limits(1), x_limits(2), x_limits(2), x_limits(1)];
    area_y = [y_limits(1), y_limits(1), 0, 0];
    fill(area_x, area_y, mySavedColors(64, 'fixed'), 'EdgeColor', 'none', 'FaceAlpha', 0.15);
end






%set the figure size
% set subplot position
if strcmp(whatComputer, 'anbu8374')==true
    set(gcf, 'Position', [0 0 2400 1200])


elseif strcmp(whatComputer, 'andrewbuggee')==true
    set(gcf, 'Position', [0 0 1500 850])


end


idx_uncert_MODIS7 = rss_residual_MODIS7./rss_uncert_MODIS7 <= 1;
percent_states_less_than_rss_uncert_MODIS7 = sum(idx_uncert_MODIS7, 'all')/numel(rss_residual_MODIS7);
disp([newline,'Percent of state space within the convergence region using 7 MODIS channels: ',...
    num2str(100*percent_states_less_than_rss_uncert_MODIS7),'%', newline])



idx_uncert_35 = rss_residual_35./rss_uncert_35 <= 1;
percent_states_less_than_rss_uncert = sum(idx_uncert_35, 'all')/numel(rss_residual_35);
disp([newline,'Percent of state space within the convergence region using 35 EMIT channels: ',...
    num2str(100*percent_states_less_than_rss_uncert),'%', newline])



% clear variables




% % ---------- Save figure --------------
% % save .fig file
% if strcmp(whatComputer, 'andrewbuggee')==true
%     folderpath_figs = '/Users/andrewbuggee/Documents/CU-Boulder-ATOC/My Papers/Figures/';
% elseif strcmp(whatComputer, 'anbu8374')==true
%     folderpath_figs = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/first_paper/figures_post_submission/';
% end
%
% f = gcf;
% saveas(f,[folderpath_figs,'Fig 8 - relative l2-norm with wavelengths for synthetic data with ', ...
%     num2str(100*measurement_uncert_MODIS7), '% total uncertainty for 7 and 35 wavelengths - ver2.fig']);
%
%
% % save .png with 400 DPI resolution
% % remove title
% title('')
% if strcmp(whatComputer, 'andrewbuggee')==true
%     folderpath_pngs = '/Users/andrewbuggee/Documents/CU-Boulder-ATOC/My Papers/Submission 1 Figures/';
% elseif strcmp(whatComputer, 'anbu8374')==true
%     folderpath_pngs = '/Users/anbu8374/Documents/My Papers/Paper 1/Submission 2 Figures/';
% end
%
% exportgraphics(f,[folderpath_pngs,'Fig 8 - relative l2-norm with wavelengths for synthetic data with ', ...
%     num2str(100*measurement_uncert_MODIS7), '% total uncertainty for 7 and 35 wavelengths - ver2.png'],'Resolution', 400);




