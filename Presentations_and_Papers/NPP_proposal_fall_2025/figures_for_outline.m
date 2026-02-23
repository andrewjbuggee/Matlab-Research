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
% idx_2Plot = [2 ,6, 7, 9, 10, 11];
idx_2Plot = [2 ,6, 9, 10];
% idx_2Plot = [57];



n_subPlots = length(idx_2Plot);

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





%% Plot the HySICS retrieval along with the in-situ measurement

% ** only considering re_profile uncertainty **

clear variables


% Determine which computer you're using
which_computer = whatComputer();

% Load simulated measurements
if strcmp(which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------

    % define the folder where the Vocals-Rex in-situ derived measurements
    % are located
    folder_paths.retrieval = ['/Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
        'paper2_variableSweep/full_retrieval_logSpace_newCov_VR_meas_allBands_with_reProf_uncert_1/'];


elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------


    % define the folder where the Vocals-Rex in-situ derived measurements
    % are located
    % folder_paths.retrieval = ['/Users/andrewbuggee/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
    %     'paper2_variableSweep/full_retrieval_logSpace_newCov_VR_meas_allBands_3/'];

    % define the folder where retrievals are located
    folder_paths.retrieval = ['/Users/andrewbuggee/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
        'paper2_variableSweep/full_retrieval_logSpace_newCov_VR_meas_allBands_with_reProf_and_cloudTop_uncert_3/'];

    % mie folder location
    mie_folder_path = '/Users/andrewbuggee/Documents/libRadtran-2.0.6/Mie_Calculations/';



elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------


end


% First step through all files in the directory and remove invisble files
% or directories

% Grab filenames in drive
filenames_retrieval = dir(folder_paths.retrieval);
idx_2delete = [];
for nn = 1:length(filenames_retrieval)

    if strcmp(filenames_retrieval(nn).name(1), 'd')~=true

        idx_2delete = [idx_2delete, nn];

    end

end

% delete rows that don't have retrieval filenames
filenames_retrieval(idx_2delete) = [];


% ------------------------------------------------------------
% -- For full_retrieval_logSpace_newCov_VR_meas_allBands_3 ---
% ------------------------------------------------------------
% profile_indexes for paper = [3, 6, 7, 9, 18]
%plt_idx = 17;
% ------------------------------------------------------------


% -------------------------------------------------------------------------------
% -- For full_retrieval_logSpace_newCov_VR_meas_allBands_with_reProf_uncert_1 ---
% -------------------------------------------------------------------------------
% profile_indexes for paper = [3, 6, 7, 9, 18]
% plt_idx = 4;
% ------------------------------------------------------------


% ---------------------------------------------------------------------------------------------
% -- For full_retrieval_logSpace_newCov_VR_meas_allBands_with_reProf_and_cloudTop_uncert_3 ---
% ---------------------------------------------------------------------------------------------
% profile_indexes for paper = [1,2,4, 6, 7,]
plt_idx = 7;
% ------------------------------------------------------------












% ------------------------------------
% ---------*** MAKE PLOT ***----------
% ------------------------------------

% Load the data from the file
ds = load([folder_paths.retrieval, filenames_retrieval(plt_idx).name]);




ln_wdth = 1;
mkr_sz = 20;

C = mySavedColors(61:68, 'fixed');
C_idx_1 = 5;
C_idx_2 = 2;


% Plot the in-situ profile

fig1 = figure;

% Create axes
axes1 = axes('Parent',fig1,'Position',[0.13599062133646 0.11 0.719723664377824 0.815]);
hold(axes1,'on');

title('In-situ vs. Retrieved droplet profile', 'Interpreter','latex',...
    'FontSize', 26)

if isfield(ds.GN_inputs.measurement, 'tau')

    if ds.GN_inputs.measurement.tau(end) ~= ds.GN_inputs.measurement.tau_c

        % check the length
        if length(ds.GN_inputs.measurement.tau) ~= length(ds.GN_inputs.measurement.re_prof)

            warning([newline, 'Tau vector length doesnt match radius profile length', newline])

            if ds.GN_inputs.measurement.tau(length(ds.GN_inputs.measurement.re_prof)) == ds.GN_inputs.measurement.tau_c

                plot(ds.GN_inputs.measurement.re_prof, ds.GN_inputs.measurement.tau(1:length(ds.GN_inputs.measurement.re_prof)),...
                    'Marker','.','LineStyle','-', 'LineWidth',ln_wdth, 'MarkerSize', mkr_sz,...
                    'Color', 'k', 'MarkerFaceColor', 'k')


            end

        end

    else


        plot(ds.GN_inputs.measurement.re_prof, ds.GN_inputs.measurement.tau,...
            'Marker','.','LineStyle','-', 'LineWidth',ln_wdth, 'MarkerSize', mkr_sz,...
            'Color', 'k', 'MarkerFaceColor', 'k')

    end






elseif isfield(ds.GN_inputs.measurement, 'tau_prof')==true


    % check the length
    if length(ds.GN_inputs.measurement.tau_prof) == length(ds.GN_inputs.measurement.re_prof)


        % ** tau_prof starts at the bottom and re_prof starts at cloud top
        % **

        plot(ds.GN_inputs.measurement.re_prof, sort(ds.GN_inputs.measurement.tau_prof),...
            'Marker','.','LineStyle','-', 'LineWidth',ln_wdth, 'MarkerSize', mkr_sz,...
            'Color', 'k', 'MarkerFaceColor', 'k')

    else

        error([newline, 'Tau vector length doesnt match radius profile length', newline])

    end

else


    plot(ds.GN_inputs.measurement.re_prof, ds.GN_inputs.measurement.tau,...
        'Marker','.','LineStyle','-', 'LineWidth',ln_wdth, 'MarkerSize', mkr_sz,...
        'Color', 'k', 'MarkerFaceColor', 'k')

end






% flip y-axis and provide axes labels
set(gca,'YDir','reverse')

% -- update the axes labels --
ylabel('Optical Depth','interpreter','latex','FontSize',35);
xlabel('Droplet Effective Radius $$(\mu m)$$','Interpreter','latex', 'FontSize', 35)

grid on; grid minor; hold on;


% plot the retrieved profile
plot(ds.GN_outputs.re_profile, ds.GN_outputs.tau_vector,...
    'LineStyle','-', 'LineWidth',ln_wdth+3,'Color', C(C_idx_1,:))

% % Plot the retrieval uncertainty of the radius at cloud top
% errorbar(ds.GN_outputs.re_profile(1), ds.GN_outputs.tau_vector(1), sqrt(ds.GN_outputs.posterior_cov_log(1,1)),...
%     'horizontal', 'Color', C(C_idx,:), 'markersize', 20, 'Linewidth', 2)
% 
% % Plot the retrieval uncertainty of the radius at cloud bottom
% errorbar(ds.GN_outputs.re_profile(end), ds.GN_outputs.tau_vector(end), sqrt(ds.GN_outputs.posterior_cov_log(2,2)),...
%     'horizontal', 'Color', C(C_idx,:), 'markersize', 20, 'Linewidth', 2)
% 
% % Plot the retrieval uncertainty of the optical depth
% errorbar(ds.GN_outputs.re_profile(end), ds.GN_outputs.tau_vector(end), sqrt(ds.GN_outputs.posterior_cov_log(3,3)),...
%     'vertical', 'Color', C(C_idx,:), 'markersize', 20, 'Linewidth', 2)



% Label cloud top and cloud bottom
% Create textbox
annotation('textbox',[0.0109176753017712 0.862655122655124 0.0509999999999998 0.0777777777777779],...
    'String',{'Cloud','Top'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',22,...
    'FitBoxToText','off');

% Create textbox
annotation('textbox',[0.0109176753017712 0.0834920634920637 0.051 0.0777777777777779],...
    'String',{'Cloud','Bottom'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',22,...
    'FitBoxToText','off');



% Plot the TBLUT droplet estimate as a constant vertical line

xl0 = xline(ds.tblut_retrieval.minRe,'-',...
    [''], 'Fontsize',24,...
    'FontWeight', 'bold', 'Interpreter','latex','LineWidth',4,'Color', C(C_idx_2,:));
xl0.LabelVerticalAlignment = 'bottom';
xl0.LabelHorizontalAlignment = 'left';

% % Plot the optical depth TBLUT retrieval as a constant horizontal line
% yl0 = yline(ds.tblut_retrieval.minTau,':',...
%     ['TBLUT $\tau_{c} = $',num2str(round(ds.tblut_retrieval.minTau, 1))], 'Fontsize',24,...
%     'FontWeight', 'bold','Interpreter','latex','LineWidth',3,'Color', C(C_idx,:));
% yl0.LabelVerticalAlignment = 'bottom';
% yl0.LabelHorizontalAlignment = 'left';


% compute the LWP estimate using the TBLUT retrieval
con = physical_constants;
rho_h2o = con.density_h2o_liquid * 1e3;   % g/m^3

lwp_tblut = (2 * rho_h2o * (ds.tblut_retrieval.minRe/1e6) * ds.tblut_retrieval.minTau)/3; % g/m^2

% ** Compute the Wood-Hartmann LWP estimate asssuming Adiabatic **
lwp_tblut_WH = 5/9 * rho_h2o * ds.tblut_retrieval.minTau * (ds.tblut_retrieval.minRe/1e6);


% grab the hypersepctral retrieval estimate of LWP
% retrieved_LWP = ds.GN_outputs.LWP;        % g/m^2
% -------------------------------------------------------
% -------------------------------------------------------
% ** Compute new updated LWP calc ***
re_profile = create_droplet_profile2([ds.GN_outputs.retrieval(1,end), ds.GN_outputs.retrieval(2,end)],...
    ds.GN_inputs.RT.z, 'altitude', ds.GN_inputs.model.profile.type);


% define the z vector
z = linspace(ds.GN_inputs.RT.z_topBottom(2), ds.GN_inputs.RT.z_topBottom(1), length(re_profile)+1)';                 % km - altitude vector

% define the z midpoint at each layer and normalize it!
z_norm = z - z(1);
z_norm_mid = (diff(z_norm)/2 + z_norm(1:end-1));


% The radius input is defined as [r_start, r_end, r_step].
% where r_step is the interval between radii values (used only for
% vectors of radii). A 0 tells the code there is no step. Finally, the
% radius values have to be in increasing order.
ext_bulk_coeff_per_LWC = zeros(length(re_profile), 1);

for rr = 1:length(re_profile)

    mie_radius = [re_profile(rr), re_profile(rr), 0];    % microns

    size_distribution = {'gamma', ds.GN_inputs.RT.distribution_var(rr)};           % droplet distribution

    % Create a mie file
    [input_filename, output_filename] = write_mie_file('MIEV0', 'water',...
        mie_radius, 500, size_distribution, 'verbose', rr, round(re_profile(rr), 4), mie_folder_path);

    % run the mie file
    [~] = runMIE(mie_folder_path, input_filename,output_filename, which_computer);

    % Read the output of the mie file
    [mie,~,~] = readMIE(mie_folder_path, output_filename);

    ext_bulk_coeff_per_LWC(rr) = mie.Qext;       % km^-1 / (cm^3 / m^3)

end


% ** Assuming liquid water content increases linearly with depth **
z_kilometers_upper_boundary = z(2:end) - z(1);                     % kilometers - geometric depth at upper boundary of each cloud layer
dz_km = z(2) - z(1);           % kilometers

%slope = tau_c /(dz_km * sum(ext_bluk_coeff_per_LWC .* z_kilometers_midpoint ));     % g/m^3/m - slope of the lwc profile
slope = ds.GN_outputs.retrieval(3,end) /(dz_km * sum(ext_bulk_coeff_per_LWC .* z_kilometers_upper_boundary ));     % g/m^3/km - slope of the lwc profile

% solve for the linear liquid water content profile
%lwc = slope * z_kilometers_midpoint;                     % g/m^3 - grams of water per meter cubed of air
lwc = slope * z_kilometers_upper_boundary;                     % g/m^3 - grams of water per meter cubed of air


retrieved_LWP = trapz( 1e3 .* z_norm_mid, lwc);    % g/m^2


% -------------------------------------------------------
% -------------------------------------------------------



% What is the true LWP
LWP_true = ds.GN_inputs.measurement.lwp;   % g/m^2

% Print this information on the figure

dim = [0.443622045822601 0.330848760514647 0.37499026298523 0.122654526321976];
str = ['$LWP_{constant} = \,$',num2str(round(lwp_tblut,1)),' $g/m^{2}$', newline,...
        '$LWP_{adiabatic} = \,$',num2str(round(retrieved_LWP,1)),' $g/m^{2}$', newline...
        '$LWP_{in-situ} = \,$',num2str(round(LWP_true,1)),' $g/m^{2}$'];

annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',25,'FontWeight','bold');







% -- update the title --
title('Example 2',...
    'Fontsize', 25, 'Interpreter', 'latex');



% -- update the legend --
legend('In-Situ', 'Theoretical Adiabatic Droplet Profile', 'Vertically Constant Effective Radius', 'Interpreter','latex', 'Position',...
[0.486474657782653 0.774845679012346 0.460486869812011 0.0950617283950617], 'FontSize', 20,...
'Color', 'white', 'TextColor', 'k')


% set figure size
set(gcf,'Position',[0 0 800 810])

box(axes1,'on');
grid(axes1,'on');
axis(axes1,'ij');
hold(axes1,'off');

% set x axis limits
xlim([min(ds.GN_inputs.measurement.re_prof) - 0.5, max(ds.GN_inputs.measurement.re_prof) + 0.5])


% --- Add right y-axis showing altitude ---
% ds.GN_inputs.RT.z is the altitude vector (km), ordered cloud-top to cloud-base
% Map altitude to the right axis so it aligns with the optical depth left axis.
% Cloud top (tau=0, z=z_top) is at the top; cloud base (tau=tau_c, z=z_base) is at the bottom.

yyaxis(axes1, 'right')
if ds.GN_inputs.RT.z(1) > ds.GN_inputs.RT.z(end)

    z_top  = ds.GN_inputs.RT.z(1);    % km - cloud top altitude
    z_base = ds.GN_inputs.RT.z(end);  % km - cloud base altitude

else

    z_base  = ds.GN_inputs.RT.z(1);    % km - cloud top altitude
    z_top = ds.GN_inputs.RT.z(end);  % km - cloud base altitude

end

% Set limits so the right axis is oriented the same way as the left:
% cloud top at top (high altitude) and cloud base at bottom (low altitude).
ylim(axes1, [z_base, z_top])
set(axes1, 'YDir', 'normal', 'YColor', 'k')      % altitude increases upward, matching reversed tau axis

ylabel('Altitude (km)', 'Interpreter', 'latex', 'FontSize', 35, 'Color', 'k')

% Restore left axis as active so subsequent operations target optical depth
yyaxis(axes1, 'left')














% ** Paper Worthy **
% -------------------------------------
% ---------- Save figure --------------
% save .fig file
if strcmp(whatComputer,'anbu8374')==true
    error(['Where do I save the figure?'])
elseif strcmp(whatComputer,'andrewbuggee')==true
    folderpath_figs = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/NPP_proposal_fall_2025/saved_figures/';
end
saveas(fig1,[folderpath_figs,'HySICS retrieval with VR in-situ measurement and TBLUT effective radius - profile number ',...
    num2str(plt_idx), '.fig']);


% save .png with 500 DPI resolution
% remove title
% title('');
exportgraphics(fig1,[folderpath_figs,'HySICS retrieval with VR in-situ measurement and TBLUT effective radius - profile number ',...
    num2str(plt_idx), '.jpg'],'Resolution', 500);
% -------------------------------------
% -------------------------------------



