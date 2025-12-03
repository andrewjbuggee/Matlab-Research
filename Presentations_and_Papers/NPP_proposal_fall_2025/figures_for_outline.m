%% Plots for NPP proposal

% By Andrew John Buggee

%% Create plot of droplet profiles from VOCALS-REx that show an increase in cloud droplet effective radius at cloud top

clear variables

% First, load ensemble data

% grab the filepath name according to which computer is being used

if strcmp(whatComputer, 'anbu8374')==true


    % Location of ensemble data
    folderpath = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/'];



    % --- non-precip profiles only, LWC>0.03, Nc>25 ----
    % load([folderpath,...
    % 'ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25_drizzleLWP-threshold_5_30-Oct-2025_rev1.mat'])

    load([folderpath,...
        'ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25_drizzleLWP-threshold_5_10-Nov-2025'])



elseif strcmp(whatComputer, 'andrewbuggee')==true

    % Location of ensemble data
    folderpath = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/VOCALS_REx/',...
        'vocals_rex_data/NCAR_C130/SPS_1/'];



    % --- non-precip profiles only, LWC>0.03, Nc>25 ----
    % load([folderpath,...
    % 'ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25_drizzleLWP-threshold_5_30-Oct-2025_rev1.mat'])

    load([folderpath,...
        'ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25_drizzleLWP-threshold_5_02-Nov-2025.mat'])

    % --- all profiles, LWC>0.03, Nc>25 ----
    % load([folderpath,...
    %     'ensemble_profiles_with_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25_02-Nov-2025.mat']);



end



% idx_2Plot = [2 ,6, 7, 9, 10, 11, 15, 16, 19, 20, 21, 28, 37, 43, 52, 57];
idx_2Plot = [2 ,6, 7, 9, 10, 11, ];



n_subPlots = 6;

ax = [];

for nn = 1:n_subPlots

     plot_LWC_and_re_and_Nc_vs_altitude(ensemble_profiles, idx_2Plot(nn), false)


end

%% OCI Weighting Functions


clear variables

% First, load ensemble data

% grab the filepath name according to which computer is being used

if strcmp(whatComputer, 'anbu8374')==true


    % Location of ensemble data
    folderpath = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/HySICS/weighting_functions/'];


    load([folderpath,...
        'disort_HySICS_reflectance_for_weightingFunctions_allBands_zOut_at_TOA_SZA-VZA_sim-ran-on-11-Jul-2025_rev1.mat'])

elseif strcmp(whatComputer, 'andrewbuggee')==true

    % Location of ensemble data
    folderpath = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/weighting_functions/'];


    load([folderpath,...
        'disort_HySICS_reflectance_for_weightingFunctions_allBands_zOut_at_TOA_SZA-VZA_sim-ran-on-11-Jul-2025_rev1.mat'])



end





% length of each independent variable
num_wl = length(inputs.bands2run);
num_tau_layers = length(tau_2run);

num_INP_files = num_wl*num_tau_layers;


% Let's fit a moving average to each weighting function and the renormalize

f = zeros(size(w));

N_mov_avg = 30;


tau_midPoint = tau_2run(1:end-1,:) + diff(tau_2run, 1, 1);

for ww = 1:num_wl

    % find the moving average
    % --- overlay a smoothed spline fit ---
    % Create smooth spline function
    %f=fit(diff(flipud(changing_variables(:,2)))/2 + flipud(tau), w, 'smoothingspline','SmoothingParam',0.95);
    f(:,ww) = movmean(w(:,ww), N_mov_avg);


    % renormalize!
    a = 1/trapz(tau_midPoint, f(:,ww));

    f(:,ww) = f(:,ww).*a;


end


% plot just the moving average weighting functions


% *** define which wavelengths to plot ***
% wl_2plot = inputs.RT.wavelengths2run(:,1);
% wl_2plot = [340:2.5:895, 940, 1038, 1250, 1378, 1615, 2130, 2260];
% wl_2plot = [340:50:895, 940, 1038, 1250, 1378, 1615, 2130, 2260];
wl_2plot = [400, 500, 600, 790, 850, 878, 940, 1038, 1250, 1378, 1615, 2130, 2260];
% wl_2plot = [900:10:1000];
% wl_2plot = [1300:10:1400];

lgnd_str = cell(numel(wl_2plot), 1);

figure;


if inputs.RT.monochromatic_calc==true

    for ww = 1:length(wl_2plot)

        [~,idx2plot] = min(abs(inputs.RT.wavelengths2run(:,1) - wl_2plot(ww)));

        plot(f(:,idx2plot), tau_midPoint, 'LineStyle', '-')

        hold on

        lgnd_str{ww} = ['$\lambda = $', num2str(round(inputs.RT.wavelengths2run(idx2plot, 1))), ' $nm$'];

    end



else

    tau = changing_variables(2:end,3);

    % plot(w, diff(flipud(changing_variables(:,2)))/2 + flipud(changing_variables(2:end,2)))
    plot(w, flipud(tau))

    % --- overlay a smoothed spline fit ---
    % Create smooth spline function
    % f=fit(diff(flipud(changing_variables(:,2)))/2 + flipud(tau), w, 'smoothingspline','SmoothingParam',0.95);
    f = movmean(w, N_mov_avg);

end




% Set up axes labels
set(gca, 'YDir','reverse')
grid on; grid minor
xlabel('$w(\tau)$','Interpreter','latex');
ylabel('$\tau$','Interpreter','latex')

% Create title
% title(['Weighting Function at ', num2str(changing_variables(1,1)), ' nm'],'Interpreter','latex')
title('Weighting Functions for OCI','Interpreter','latex')




set(gcf, 'Position',[0 0 1400 800])




% Create Legend
% legend(string(inputs.RT.wavelengths2run(:,1))','Interpreter','latex','Location','northwest','FontSize',22)
legend(lgnd_str,'Interpreter','latex','Location','northwest','FontSize',22,...
    'Color', 'white', 'TextColor', 'k')

% Create textbox with simulation properties

% Textbox
dim = [0.155714285714286 0.144548492431641 0.196462309701102 0.382951507568359];

if ischar(inputs.RT.sensor_altitude)==true
    texBox_str = {['$sza$ = ',num2str(inputs.RT.sza)],...
        ['$vza$ = ',num2str(inputs.RT.vza)],...
        ['$z_{out}$ = ', inputs.RT.sensor_altitude],...
        ['$Cloud\;top$ = ', num2str(inputs.RT.z_topBottom(1)), ' km'],...
        ['$Cloud\;base$ = ', num2str(inputs.RT.z_topBottom(2)), ' km'],...
        ['$r_{top}$ = ',num2str(round(inputs.RT.r_top)), ' $\mu m$'],...
        ['$r_{bot}$ = ',num2str(round(inputs.RT.r_bot)), ' $\mu m$'],...
        ['$\tau_0$ = ', num2str(inputs.RT.tau_c)],...
        ['$A_0$ = ', num2str(inputs.RT.surface_albedo)]};
else
    texBox_str = {['$sza$ = ',num2str(inputs.RT.sza)],...
        ['$vza$ = ',num2str(inputs.RT.vza)],...
        ['$z_{out}$ = ', num2str(inputs.RT.sensor_altitude), ' km'],...
        ['$Cloud\;top$ = ', num2str(inputs.RT.z_topBottom(1)), ' km'],...
        ['$Cloud\;base$ = ', num2str(inputs.RT.z_topBottom(2)), ' km'],...
        ['$r_{top}$ = ',num2str(round(inputs.RT.r_top)), ' $\mu m$'],...
        ['$r_{bot}$ = ',num2str(round(inputs.RT.r_bot)), ' $\mu m$'],...
        ['$\tau_0$ = ', num2str(inputs.RT.tau_c)],...
        ['$A_0$ = ', num2str(inputs.RT.surface_albedo)]};
end

t = annotation('textbox',dim,'string',texBox_str,'Interpreter','latex');
t.Color = 'black';
t.FontSize = 25;
t.FontWeight = 'bold';
t.EdgeColor = 'black';
t.FitBoxToText = 'on';






