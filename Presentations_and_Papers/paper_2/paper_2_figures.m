%% Figures for Paper 2

% By Andrew John Buggee

%%  Plot droplet profile retrievals using simualted HySICS measurements with vocals-Rex in-situ measurements


clear variables


% Determine which computer you're using
which_computer = whatComputer();

% Load simulated measurements
if strcmp(which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------



elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------


    % define the folder where the Vocals-Rex in-situ derived measurements
    % are located
    folder_paths.retrieval = ['/Users/andrewbuggee/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
        'paper2_variableSweep/full_retrieval_logSpace_newCov_VR_meas_allBands_3/'];



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



plt_idx = 2;

plot_retrieved_prof_with_inSitu_paper2(folder_paths.retrieval, filenames_retrieval(plt_idx).name)




%% Compare full HySICS retrieval with exact knowledge of cloud top height with full retrieval assuming a 
% cloud top height of 1.5 km, and including a cloud top height unceratinty
% in the forward model covariance matrix


%%  Plot droplet profile retrievals using simualted HySICS measurements with vocals-Rex in-situ measurements



clear variables

% Determine which computer you're using
which_computer = whatComputer();

% Load simulated measurements
if strcmp(which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------



elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------


    % ** retrievals with exact knowledge of cloud top height **
    folder_paths.retrieval_exactCloudTop = ['/Users/andrewbuggee/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
        'paper2_variableSweep/full_retrieval_logSpace_newCov_VR_meas_allBands_3/'];


    % ** retrievals assuming cloud top height of 1.5 km and forward model
    % uncertainty **
    folder_paths.retrieval_assumedCloudTop = ['/Users/andrewbuggee/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
        'paper2_variableSweep/full_retrieval_logSpace_newCov_VR_meas_allBands_with_cloudTop_uncert_1/'];



elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------


end


% First step through all files in the directory and remove invisble files
% or directories

% Grab filenames in drive
filenames_retrieval_exactCloudTop = dir(folder_paths.retrieval_exactCloudTop);
idx_2delete = [];
for nn = 1:length(filenames_retrieval_exactCloudTop)

    if strcmp(filenames_retrieval_exactCloudTop(nn).name(1), 'd')~=true

        idx_2delete = [idx_2delete, nn];

    end

end

% delete rows that don't have retrieval filenames
filenames_retrieval_exactCloudTop(idx_2delete) = [];





% First step through all files in the directory and remove invisble files
% or directories

% Grab filenames in drive
filenames_retrieval_assumedCloudTop = dir(folder_paths.retrieval_assumedCloudTop);
idx_2delete = [];
for nn = 1:length(filenames_retrieval_assumedCloudTop)

    if strcmp(filenames_retrieval_assumedCloudTop(nn).name(1), 'd')~=true

        idx_2delete = [idx_2delete, nn];

    end

end

% delete rows that don't have retrieval filenames
filenames_retrieval_assumedCloudTop(idx_2delete) = [];











plt_idx = 70;

ln_wdth = 1;
mkr_sz = 20;

C = mySavedColors(61:68, 'fixed');
C_idx = 5;

% compute the LWP estimate using the TBLUT retrieval
rho_liquid_water = 10^6;        % g/m^3








% Plot the in-situ profile

figure;

subplot(1,2,1)


title('With exact cloud top', 'Interpreter','latex',...
    'FontSize', 26)


% Load the data for retrieval using exact cloud top height
ds = load([folder_paths.retrieval_exactCloudTop, filenames_retrieval_exactCloudTop(plt_idx).name]);


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

% flip y-axis and provide axes labels
set(gca,'YDir','reverse')
ylabel('$\tau$','interpreter','latex','FontSize',35);
xlabel('$r_{e}$ $$(\mu m)$$','Interpreter','latex', 'FontSize', 35)
grid on; grid minor; hold on;


% plot the retrieved profile
plot(ds.GN_outputs.re_profile, ds.GN_outputs.tau_vector,...
    'LineStyle','-', 'LineWidth',ln_wdth+2,'Color', C(C_idx,:))

% Plot the retrieval uncertainty of the radius at cloud top
errorbar(ds.GN_outputs.re_profile(1), ds.GN_outputs.tau_vector(1), sqrt(ds.GN_outputs.posterior_cov_log(1,1)),...
    'horizontal', 'Color', C(C_idx,:), 'markersize', 20, 'Linewidth', 2)

% Plot the retrieval uncertainty of the radius at cloud bottom
errorbar(ds.GN_outputs.re_profile(end), ds.GN_outputs.tau_vector(end), sqrt(ds.GN_outputs.posterior_cov_log(2,2)),...
    'horizontal', 'Color', C(C_idx,:), 'markersize', 20, 'Linewidth', 2)

% Plot the retrieval uncertainty of the optical depth
errorbar(ds.GN_outputs.re_profile(end), ds.GN_outputs.tau_vector(end), sqrt(ds.GN_outputs.posterior_cov_log(3,3)),...
    'vertical', 'Color', C(C_idx,:), 'markersize', 20, 'Linewidth', 2)



% Label cloud top and cloud bottom
% Create textbox
annotation('textbox',[0.02,0.865079365079366,0.051,0.077777777777778],...
    'String',{'Cloud','Top'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',22,...
    'FitBoxToText','off');

% Create textbox
annotation('textbox',[0.02,0.096825396825397,0.051,0.077777777777778],...
    'String',{'Cloud','Bottom'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',22,...
    'FitBoxToText','off');



% Plot the TBLUT droplet estimate as a constant vertical line

xl0 = xline(ds.tblut_retrieval.minRe,':',...
    ['TBLUT $r_{e} = $',num2str(round(ds.tblut_retrieval.minRe, 1)), '$\mu m$'], 'Fontsize',24,...
    'FontWeight', 'bold', 'Interpreter','latex','LineWidth',3,'Color', C(C_idx,:));
xl0.LabelVerticalAlignment = 'bottom';
xl0.LabelHorizontalAlignment = 'left';

% Plot the optical depth TBLUT retrieval as a constant horizontal line
yl0 = yline(ds.tblut_retrieval.minTau,':',...
    ['TBLUT $\tau_{c} = $',num2str(round(ds.tblut_retrieval.minTau, 1))], 'Fontsize',24,...
    'FontWeight', 'bold','Interpreter','latex','LineWidth',3,'Color', C(C_idx,:));
yl0.LabelVerticalAlignment = 'top';
yl0.LabelHorizontalAlignment = 'right';




lwp_tblut = (2*rho_liquid_water*(ds.tblut_retrieval.minRe/1e6) * ds.tblut_retrieval.minTau)/3; % g/m^2

% grab the hypersepctral retrieval estimate of LWP
retrieved_LWP = ds.GN_outputs.LWP;        % g/m^2

% What is the true LWP
LWP_true = ds.GN_inputs.measurement.lwp;   % g/m^2

% Print this information on the figure

dim = [0.141166666666667 0.71690449790349 0.293506310780843 0.157698676699684];
str = ['$LWP_{TBLUT} = \,$',num2str(round(lwp_tblut,1)),' $g/m^{2}$', newline,...
    '$LWP_{hyperspectral} = \,$',num2str(round(retrieved_LWP,1)),' $g/m^{2}$', newline...
    '$LWP_{true} = \,$',num2str(round(LWP_true,1)),' $g/m^{2}$'];

annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',25,'FontWeight','bold');




% plot the retrieved column water vapor if it was retireved
if size(ds.GN_outputs.retrieval, 1)>3

    acpw_true = ds.GN_inputs.measurement.actpw;             % kg/m^2 (mm)
    retrieved_CWV = ds.GN_outputs.retrieval(end, end);        % kg/m^2 (mm)

    % Print the simulated value and the retrieved value
    str = ['$acpw_{true} = \,$',num2str(round(acpw_true, 2)),' $mm$', newline,...
        '$acpw_{retrieved} = \,$',num2str(round(retrieved_CWV, 2)),' $mm$'];

else

    % plot the assumed column water vapor used in the forward model
    % plot the HySICS simulated above cloud column water vapor
    assumed_CWV = aboveCloud_CWV_simulated_hysics_spectra(ds.GN_inputs); % kg/m^2

    % print the simulated value and the foward model assumption
    str = ['$acpw_{forward \,model} = \,$',num2str(round(assumed_CWV, 2)),' $mm$', newline,...
        '$acpw_{MODIS} = \,$',num2str(modis_retrieved_aboveCloud_CWV),' $mm$'];

end


dim = [0.141166666666667 0.571252065991597 0.238002173105876 0.103351108611576];


annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',25,'FontWeight','bold');



legend('In-Situ', 'Hyperspectral Retrieval', 'Interpreter','latex', 'Location','best', 'FontSize', 20,...
            'Color', 'white', 'TextColor', 'k')



% ----------------------------------------------------------------
% ----------------------------------------------------------------
% ** Plot retrieval without exact knowledge of cloud top height **
% ----------------------------------------------------------------
% ----------------------------------------------------------------


subplot(1,2,2)


title('Assumed cloud top of 1.5 km', 'Interpreter','latex',...
    'FontSize', 26)


% Load the data for retrieval using exact cloud top height
ds = load([folder_paths.retrieval_assumedCloudTop, filenames_retrieval_assumedCloudTop(plt_idx).name]);


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

% flip y-axis and provide axes labels
set(gca,'YDir','reverse')
ylabel('$\tau$','interpreter','latex','FontSize',35);
xlabel('$r_{e}$ $$(\mu m)$$','Interpreter','latex', 'FontSize', 35)
grid on; grid minor; hold on;


% plot the retrieved profile
plot(ds.GN_outputs.re_profile, ds.GN_outputs.tau_vector,...
    'LineStyle','-', 'LineWidth',ln_wdth+2,'Color', C(C_idx,:))

% Plot the retrieval uncertainty of the radius at cloud top
errorbar(ds.GN_outputs.re_profile(1), ds.GN_outputs.tau_vector(1), sqrt(ds.GN_outputs.posterior_cov_log(1,1)),...
    'horizontal', 'Color', C(C_idx,:), 'markersize', 20, 'Linewidth', 2)

% Plot the retrieval uncertainty of the radius at cloud bottom
errorbar(ds.GN_outputs.re_profile(end), ds.GN_outputs.tau_vector(end), sqrt(ds.GN_outputs.posterior_cov_log(2,2)),...
    'horizontal', 'Color', C(C_idx,:), 'markersize', 20, 'Linewidth', 2)

% Plot the retrieval uncertainty of the optical depth
errorbar(ds.GN_outputs.re_profile(end), ds.GN_outputs.tau_vector(end), sqrt(ds.GN_outputs.posterior_cov_log(3,3)),...
    'vertical', 'Color', C(C_idx,:), 'markersize', 20, 'Linewidth', 2)



% Label cloud top and cloud bottom
% Create textbox
annotation('textbox',[0.02,0.865079365079366,0.051,0.077777777777778],...
    'String',{'Cloud','Top'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',22,...
    'FitBoxToText','off');

% Create textbox
annotation('textbox',[0.02,0.096825396825397,0.051,0.077777777777778],...
    'String',{'Cloud','Bottom'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',22,...
    'FitBoxToText','off');



% Plot the TBLUT droplet estimate as a constant vertical line

xl0 = xline(ds.tblut_retrieval.minRe,':',...
    ['TBLUT $r_{e} = $',num2str(round(ds.tblut_retrieval.minRe, 1)), '$\mu m$'], 'Fontsize',24,...
    'FontWeight', 'bold', 'Interpreter','latex','LineWidth',3,'Color', C(C_idx,:));
xl0.LabelVerticalAlignment = 'bottom';
xl0.LabelHorizontalAlignment = 'left';

% Plot the optical depth TBLUT retrieval as a constant horizontal line
yl0 = yline(ds.tblut_retrieval.minTau,':',...
    ['TBLUT $\tau_{c} = $',num2str(round(ds.tblut_retrieval.minTau, 1))], 'Fontsize',24,...
    'FontWeight', 'bold','Interpreter','latex','LineWidth',3,'Color', C(C_idx,:));
yl0.LabelVerticalAlignment = 'top';
yl0.LabelHorizontalAlignment = 'right';




lwp_tblut = (2*rho_liquid_water*(ds.tblut_retrieval.minRe/1e6) * ds.tblut_retrieval.minTau)/3; % g/m^2

% grab the hypersepctral retrieval estimate of LWP
retrieved_LWP = ds.GN_outputs.LWP;        % g/m^2

% What is the true LWP
LWP_true = ds.GN_inputs.measurement.lwp;   % g/m^2

% Print this information on the figure

dim = [0.141166666666667 0.71690449790349 0.293506310780843 0.157698676699684];
str = ['$LWP_{TBLUT} = \,$',num2str(round(lwp_tblut,1)),' $g/m^{2}$', newline,...
    '$LWP_{hyperspectral} = \,$',num2str(round(retrieved_LWP,1)),' $g/m^{2}$', newline...
    '$LWP_{true} = \,$',num2str(round(LWP_true,1)),' $g/m^{2}$'];

annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',25,'FontWeight','bold');




% plot the retrieved column water vapor if it was retireved
if size(ds.GN_outputs.retrieval, 1)>3

    acpw_true = ds.GN_inputs.measurement.actpw;             % kg/m^2 (mm)
    retrieved_CWV = ds.GN_outputs.retrieval(end, end);        % kg/m^2 (mm)

    % Print the simulated value and the retrieved value
    str = ['$acpw_{true} = \,$',num2str(round(acpw_true, 2)),' $mm$', newline,...
        '$acpw_{retrieved} = \,$',num2str(round(retrieved_CWV, 2)),' $mm$'];

else

    % plot the assumed column water vapor used in the forward model
    % plot the HySICS simulated above cloud column water vapor
    assumed_CWV = aboveCloud_CWV_simulated_hysics_spectra(ds.GN_inputs); % kg/m^2

    % print the simulated value and the foward model assumption
    str = ['$acpw_{forward \,model} = \,$',num2str(round(assumed_CWV, 2)),' $mm$', newline,...
        '$acpw_{MODIS} = \,$',num2str(modis_retrieved_aboveCloud_CWV),' $mm$'];

end


dim = [0.141166666666667 0.571252065991597 0.238002173105876 0.103351108611576];


annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',25,'FontWeight','bold');



legend('In-Situ', 'Hyperspectral Retrieval', 'Interpreter','latex', 'Location','best', 'FontSize', 20,...
            'Color', 'white', 'TextColor', 'k')







% set figure size
set(gcf,'Position',[0 0 1400 830])







%% Plot the new figure with cloud top height unceratinty included in the covariance matrix
 
clear variables


folder_path = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
    'HySICS/Droplet_profile_retrievals/paper2_variableSweep/test_logSpace_newCov_with_VR_inSitu_meas/'];

filename = ['dropletRetrieval_HySICS_636bands_0.3%_uncert_vocalsRex_recorded_',...
    '31-Oct-2008_9.0864UTC_vza_4_vaz_257_sza_31_saz_96_sim-ran-on-23-Dec-2025_1.mat'];

plot_retrieved_prof_with_inSitu_paper2(folder_path, filename)
