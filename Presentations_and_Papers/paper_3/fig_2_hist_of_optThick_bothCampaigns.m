%% Figure 2 - Plot a histogram of all the optical thickness from the VOCALS-REx and ORACLES data


% ------------------------------------------------------------------------
% ------           Load  VOCALS-REx data ---
% -------------------------------------------------------------------------
clear variables

% *** Optical Depth Filter ***
tauC_min = 3;
tauC_max = 40;

% -------------------------------------------------------------------------
%  File locations
%  -------------------------------------------------------------------------

which_computer = whatComputer;

if strcmp(which_computer, 'anbu8374')



elseif strcmp(which_computer, 'andrewbuggee')

    foldername_data = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/'];

    foldername_save = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/', ...
        'Presentations_and_Papers/paper_3/'];

end


% -------------------------------------------------------------------------
%  Load ensemble profiles mat file
%  -------------------------------------------------------------------------

% Version from May 2026
% load([foldername_data,...
%     'ensemble_profiles_with_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25_04-Dec-2025.mat'])

% Version from November 2026 with optical depth correction
load([foldername_data,...
    'ensemble_profiles_with_precip_from_14_files_LWC-threshold_0.05_Nc-threshold_10_08-Nov-2026.mat'])

N_profiles = length(ensemble_profiles);
disp(['Loaded ', num2str(N_profiles), ' vertical profiles.'])

drizzle_LWP_threshold = 5;   % g/m^2
is_drizzle = false(1, N_profiles);
for nn = 1:N_profiles
    is_drizzle(nn) = ensemble_profiles{nn}.lwp_2DC >= drizzle_LWP_threshold;
end

tauC_nd  = [];
tauC_d   = [];
tauC_all = [];

num_removed = 0;
num_kept    = 0;

% step through all VOCALS-REx data

for nn = 1:N_profiles

    tau_prof = ensemble_profiles{nn}.tau;
    if max(tau_prof) > tauC_max || max(tau_prof) < tauC_min
        num_removed = num_removed + 1;
        continue
    else
        num_kept = num_kept + 1;
    end

    tauC_all = [tauC_all; max(tau_prof)];

    if is_drizzle(nn) == true

        tauC_d = [tauC_d; max(tau_prof)];

    else

        tauC_nd = [tauC_nd; max(tau_prof)];

    end


end

disp(['Mean cloud optical thickness for VOCALS-REx profiles: ', num2str(round(mean(tauC_all), 3))])
disp(['Mean cloud optical thickness for non-drizzling clouds from VOCALS-REx profiles: ', num2str(round(mean(tauC_nd), 3))])
disp(['Mean cloud optical thickness for drizzling clouds from VOCALS-REx profiles: ', num2str(round(mean(tauC_d), 3))])


% Create histogram of the drizzle and non-drizzle profiles

plt_clr_1   = mySavedColors(68,'fixed');  % bubble gum pink
plt_clr_2   = mySavedColors(67,'fixed');  % Broccoli
num_bins = 50;
alpha_val_1 = 0.3;
alpha_val_2 = 0.5;
ttl_fnt = 20;
ax_fnt = 20;

% Set the height and width as multiples of the IEEE paper requirements
w = 3.5; % inches
h = 3.5;    % inches

figure1 =     figure;
figure1.Units = 'inches';
figure1.Position = [1, 1, 2*w, 2*h]; % [left, bottom, width, height]
subplot(2,1,1)

histogram(tauC_nd, num_bins-15, 'FaceColor', plt_clr_1, 'FaceAlpha', alpha_val_1)
hold on
histogram(tauC_d, num_bins/2, 'FaceColor', plt_clr_2, 'FaceAlpha', alpha_val_2)
grid on; grid minor

legend('Non-drizzling', 'Drizzling', 'Interpreter', 'latex', 'Fontsize', ax_fnt,...
    'location', 'best')

xlabel('Cloud Optical Thickness', 'Interpreter', 'latex', 'FontSize', ax_fnt)
ylabel('Counts','Interpreter', 'latex', 'FontSize', ax_fnt)
% title(['VOCALS-REx - N = ', num2str(num_kept), ' - $\tau_{c} \geq$ ', num2str(tauC_min)],...
%     'Interpreter', 'latex','FontSize', ttl_fnt)
title(['VOCALS-REx - N = ', num2str(num_kept)],...
    'Interpreter', 'latex','FontSize', ttl_fnt)

% update x-axis tick labels
xticks([0, 10, 20, 30, 40, 50, 60, 70])
xticklabels({'0', '10', '20', '30', '40', '50', '60', '70'})
xlim([0, tauC_max])

% Create textbox
annotation(figure1,'textbox',[0.806725301587295 0.686952380952378 0.0781953333333357 0.0666],'String',{'(a)'},...
    'Interpreter','latex',...
    'FontSize',ax_fnt,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');



% -------------------------------------------------------------------------
% ------ ORACLES: load, separate drizzle, compute percentiles -------------
% -------------------------------------------------------------------------
clear ensemble_profiles



if strcmp(which_computer, 'anbu8374')

    foldername_data = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/', ...
        'Hyperspectral_Cloud_Retrievals/ORACLES/oracles_data/'];

elseif strcmp(which_computer, 'andrewbuggee')

    foldername_data = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/', ...
        'Hyperspectral_Cloud_Retrievals/ORACLES/oracles_data/'];

end

% % Version from May 2026
% load([foldername_data,...
%     'ensemble_profiles_with_precip_from_33_files_LWC-threshold_0.05_Nc-threshold_10_no_rEff_greaterThan50_microns_16-Mar-2026.mat'])

% Version from November 2026 with optical depth correction 
load([foldername_data,...
    'ensemble_profiles_from_33_files_LWC-threshold_0.05_Nc-threshold_10_08-Nov-2026.mat'])


% step through all ORACLES data
N_profiles = length(ensemble_profiles);
disp(['Loaded ', num2str(N_profiles), ' vertical profiles.'])

drizzle_LWP_threshold = 5;   % g/m^2
is_drizzle = false(1, N_profiles);

tauC_nd  = [];
tauC_d   = [];
tauC_all = [];


num_removed  = 0;
num_kept     = 0;
num_outliers = 0;

for nn = 1:N_profiles
    is_drizzle(nn) = ensemble_profiles{nn}.lwp_2DS_HVPS >= drizzle_LWP_threshold;
end

% --- Outliers excluded from the synthetic profile generator ---------------
% These are the same two ORACLES profiles listed in OUTLIERS in
% lasp-CU-paper-3/06_synthetic_profile_generator.py (lines 103-151), and
% excluded from Fig. 1 (fig_1_inSitu_dist_8_panel_vs_z.m). There,
% profile_idx is the 1-based position in ensemble_profiles of the
% 08-Nov-2026 ORACLES file loaded above (it is 'kk' in
% build_synthetic_profile_inputs_from_insitu_and_ERA5.m), so the indices
% carry over directly.
%
% Both are from 30 Aug 2017: median LWC of 1.9 and 2.3 g/m^3 in layers only
% 122 and 125 m deep, i.e. f_ad = LWP_obs/LWP_ad of 9.3 and 10.4 (Wood 2005
% definition). Every other profile in the generator's training set has
% f_ad <= 2.75. See the generator header for the full justification.
% Their tau_c (23.9 and 28.2) lies inside [tauC_min, tauC_max], so without
% this block they would be counted in the histogram.
% Empty these two arrays to restore them.
outlier_idx  = [111, 112];
outlier_date = {'2017-08-30', '2017-08-30'};

% Guard: make sure the indices still point at the intended profiles. If
% the ensemble file is rebuilt the indices can shift silently.
for oo = 1:length(outlier_idx)
    date_oo = char(ensemble_profiles{outlier_idx(oo)}.dateOfFlight, 'yyyy-MM-dd');
    if ~strcmp(date_oo, outlier_date{oo})
        error([newline, 'ORACLES outlier index ', num2str(outlier_idx(oo)), ...
            ' has flight date ', date_oo, ', expected ', outlier_date{oo}, '.', ...
            newline, 'The ensemble file may have changed; update outlier_idx.', newline])
    end
end

for nn = 1:N_profiles

    % skip the outliers first so they aren't also counted in num_removed
    if ismember(nn, outlier_idx)
        num_outliers = num_outliers + 1;
        continue
    end

    tau_prof = ensemble_profiles{nn}.tau;
    if max(tau_prof) > tauC_max || max(tau_prof) < tauC_min
        num_removed = num_removed + 1;
        continue
    else
        num_kept = num_kept + 1;
    end

    tauC_all = [tauC_all; max(tau_prof)];

    if is_drizzle(nn) == true

        tauC_d = [tauC_d; max(tau_prof)];

    else

        tauC_nd = [tauC_nd; max(tau_prof)];

    end


end

disp(['# ORACLES profiles kept   : ', num2str(num_kept)])
disp(['# ORACLES profiles removed: ', num2str(num_removed)])
disp(['# ORACLES outliers removed: ', num2str(num_outliers)])
disp(['Mean cloud optical thickness for ORALCES profiles: ', num2str(round(mean(tauC_all), 3))])
disp(['Mean cloud optical thickness for non-drizzling clouds from ORALCES profiles: ', num2str(round(mean(tauC_nd), 3))])
disp(['Mean cloud optical thickness for drizzling clouds from ORALCES profiles: ', num2str(round(mean(tauC_d), 3))])


subplot(2,1,2)

histogram(tauC_nd, num_bins, 'FaceColor', plt_clr_1, 'FaceAlpha', alpha_val_1)
hold on
histogram(tauC_d, num_bins/2, 'FaceColor', plt_clr_2, 'FaceAlpha', alpha_val_2)
grid on; grid minor


xlabel('Cloud Optical Thickness', 'Interpreter', 'latex', 'FontSize', ax_fnt)
ylabel('Counts','Interpreter', 'latex', 'FontSize', ax_fnt)
title(['ORACLES - N = ', num2str(num_kept), ' - $\tau_{c} \geq$ ', num2str(tauC_min)],...
    'Interpreter', 'latex','FontSize', ttl_fnt)
% update x-axis tick labels
xticks([0, 10, 20, 30, 40, 50, 60, 70])
xticklabels({'0', '10', '20', '30', '40', '50', '60', '70'})
xlim([0, tauC_max])

% Create textbox
annotation(figure1,'textbox',[0.806725301587295 0.313936507936503 0.0781953333333357 0.0665999999999995],'String',{'(b)'},...
    'Interpreter','latex',...
    'FontSize',ax_fnt,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');



% ** Paper Worthy **
% -------------------------------------
% ---------- Save figure --------------
% save .fig file
if strcmp(whatComputer,'anbu8374')==true
    error(['Where do I save the figure?'])
elseif strcmp(whatComputer,'andrewbuggee')==true
    folderpath_figs = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_3/saved_figures/';
end
saveas(figure1,[folderpath_figs,'Histogram of drizzling and non-drizzling cloud optical thickness',...
    ' for VOCALS-REx and ORACLES from in-situ data.fig']);


% save .png with 500 DPI resolution
% remove title
% title('');
exportgraphics(figure1,[folderpath_figs,'Histogram of drizzling and non-dizzling cloud optical thickness',...
    ' for VOCALS-REx and ORALCES from in-situ data.png'],...
    'Resolution', 500);
% -------------------------------------
% -------------------------------------





