%% Make paneled figure with each optical-depth for a single total-column-water-vapor pair
% Compare the 4 retrievals that assume a total column water vapor amount with the retrieval for
% above cloud column water vapor


% *** DO THESE MATCH THE CALCULATIONS I MADE FOR GRC???? ******

% The files below were computed using the retrieval algorithm on 9-14-2025
% on the alpine cluster

clear variables



% load 5% uncertainty files
filenames_noACPW_startsWith = 'dropletRetrieval_noACPW_HySICS_35bands';
filenames_full_startsWith = 'dropletRetrieval_HySICS_66bands';

folder_path = ['/Users/andrewbuggee/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
    'paper2_variableSweep/rTop_10/test_folder_mimic_localMachine_results/'];


% load files
filenames = dir([folder_path, '*.mat']);

% Step through each file


% define the colors for each curve plotted
C = mySavedColors(61:(61+length(filenames)+1), 'fixed');

lgnd_str = cell(1, length(filenames) + 2);


figure;


for nn = 1:length(filenames)


    % Load a data set
    ds = load([folder_path, filenames(nn).name]);

    if nn==1


        % first, plot the simulated profile values as two lines
        xline(ds.GN_inputs.RT.r_bot, ':', ['Simulated $r_{bot}$'],...
            'Fontsize',20, 'Interpreter','latex','LineWidth',2,'Color', [147/255, 150/255, 151/255], 'LabelHorizontalAlignment','left',...
            'LabelVerticalAlignment','bottom');

        hold on

        yline(ds.GN_inputs.RT.r_top, ':', ['Simulated $r_{top}$'],...
            'Fontsize',20, 'Interpreter','latex','LineWidth',2,'Color', [147/255, 150/255, 151/255], 'LabelHorizontalAlignment','right',...
            'LabelVerticalAlignment','top');
        hold on


        % what was the assumed above cloud column water vapor path?
        simulated_CWV = ds.GN_inputs.measurement.actpw; % kg/m^2

        title(['Simulated profile - $acpw$ = ',num2str(round(simulated_CWV, 2)), ' $mm$'],...
            'Fontsize', 25, 'Interpreter', 'latex');

        % Skip the first two legend entries
        lgnd_str{1} = '';
        lgnd_str{2} = '';

    end



    % plot the retrieved droplet profile
    if nn<length(filenames)

        hold on


        % Plot the retrieval uncertainty of the radius at cloud top and
        % bottom
        e1 = errorbar(ds.GN_outputs.retrieval(2,end), ds.GN_outputs.retrieval(1,end), sqrt(ds.GN_outputs.posterior_cov(1,1))/2,...
            sqrt(ds.GN_outputs.posterior_cov(1,1))/2, sqrt(ds.GN_outputs.posterior_cov(2,2))/2,...
            sqrt(ds.GN_outputs.posterior_cov(2,2))/2, 'MarkerFaceColor', C(nn+1,:),...
            'MarkerEdgeColor', C(nn+1,:), 'Linewidth', 4, 'Marker', '.', 'MarkerSize', 40,...
            'Color', C(nn+1,:));

        e1.Bar.LineStyle = 'dotted';

        hold on



    else

        % give a different marker type for the retrieval using 66 bands
        % that also retrieved above cloud column water vapor


        errorbar(ds.GN_outputs.retrieval(2,end), ds.GN_outputs.retrieval(1,end), sqrt(ds.GN_outputs.posterior_cov(1,1))/2,...
            sqrt(ds.GN_outputs.posterior_cov(1,1))/2, sqrt(ds.GN_outputs.posterior_cov(2,2))/2,...
            sqrt(ds.GN_outputs.posterior_cov(2,2))/2, 'MarkerFaceColor', C(nn+1,:),...
            'MarkerEdgeColor', C(nn+1,:), 'Linewidth', 4, 'Marker', '.', 'MarkerSize', 40,...
            'Color', C(nn+1,:))

        hold on


    end



    % create the legend string
    if nn<length(filenames)

        % what was the assumed above cloud column water vapor path?
        assumed_CWV = aboveCloud_CWV_simulated_hysics_spectra(ds.GN_inputs); % kg/m^2

        lgnd_str{nn+2} = [num2str(numel(ds.GN_inputs.bands2run)), ' bands - assumed $acpw$ = ',...
            num2str(round(assumed_CWV,2)), ' $mm$'];





    else

        % create the string for the retrieval using CWV

        % what was the retrieved above cloud column water vapor path above
        % cloud?
        retrieved_CWV = ds.GN_outputs.retrieval(end, end);        % kg/m^2 (mm)

        lgnd_str{nn+2} = [num2str(numel(ds.GN_inputs.bands2run)), ' bands - retrieved $acpw$ = ',...
            num2str(round(retrieved_CWV, 2)), ' $mm$'];

    end

end


% Create a Legend with only the two black curves
legend(lgnd_str, 'Interpreter','latex', 'Location','northwest', 'FontSize', 20)

grid on; grid minor
ylabel('$r_{top}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',30)
xlabel('$r_{bot}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',30)

set(gcf,'Position',[0 0 950 750])





%% Make paneled figure with each optical-depth for a single total-column-water-vapor pair
% Compare the 4 retrievals that assume a total column water vapor amount with the retrieval for
% above cloud column water vapor


% *** DO THESE MATCH THE CALCULATIONS I MADE FOR GRC???? ******

% The files below were computed using the retrieval algorithm on 9-14-2025
% on my local machine



clear variables


% folder_path = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
%     'Hyperspectral_Cloud_Retrievals/HySICS/Droplet_profile_retrievals/'];

folder_path = ['/Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
    'paper2_variableSweep/test_folder_mimic_localMachine_results/'];

% define the mat files for each retreival to plot
filenames = {'dropletRetrieval_noACPW_HySICS_35bands_1%_uncert_rTop_9.2516_rBot_5.3192_tauC_6.1312_tcwv_20_tcwv-assumption_10_vza_4_vaz_257_sza_31_saz_96_sim-ran-on-13-Sep-2025_1.mat',...
    'dropletRetrieval_noACPW_HySICS_35bands_1%_uncert_rTop_9.2516_rBot_5.3192_tauC_6.1312_tcwv_20_tcwv-assumption_15_vza_4_vaz_257_sza_31_saz_96_sim-ran-on-13-Sep-2025_1.mat',...
    'dropletRetrieval_noACPW_HySICS_35bands_1%_uncert_rTop_9.2516_rBot_5.3192_tauC_6.1312_tcwv_20_tcwv-assumption_20_vza_4_vaz_257_sza_31_saz_96_sim-ran-on-13-Sep-2025_1.mat',...
    'dropletRetrieval_noACPW_HySICS_35bands_1%_uncert_rTop_9.2516_rBot_5.3192_tauC_6.1312_tcwv_20_tcwv-assumption_25_vza_4_vaz_257_sza_31_saz_96_sim-ran-on-13-Sep-2025_1.mat',...
    'dropletRetrieval_HySICS_66bands_1%_uncert_rTop_9.2516_rBot_5.3192_tauC_6.1312_tcwv_20_vza_4_vaz_257_sza_31_saz_96_sim-ran-on-15-Sep-2025_1'};




% Step through each file


% define the colors for each curve plotted
C = mySavedColors(61:(61+length(filenames)+1), 'fixed');

lgnd_str = cell(1, length(filenames) + 2);


figure;


for nn = 1:length(filenames)



    % Load a data set
    ds = load([folder_path, filenames{nn}]);

    if nn==1


        % first, plot the simulated profile values as two lines
        xline(ds.GN_inputs.measurement.r_bot, ':', ['Simulated $r_{bot}$'],...
            'Fontsize',20, 'Interpreter','latex','LineWidth',2,'Color', [147/255, 150/255, 151/255], 'LabelHorizontalAlignment','left',...
            'LabelVerticalAlignment','bottom');

        hold on

        yline(ds.GN_inputs.measurement.r_top, ':', ['Simulated $r_{top}$'],...
            'Fontsize',20, 'Interpreter','latex','LineWidth',2,'Color', [147/255, 150/255, 151/255], 'LabelHorizontalAlignment','right',...
            'LabelVerticalAlignment','top');
        hold on


        % what was the assumed above cloud column water vapor path?
        simulated_CWV = ds.GN_inputs.measurement.actpw; % kg/m^2

        title(['Simulated profile - $acpw$ = ',num2str(round(simulated_CWV, 2)), ' $mm$'],...
            'Fontsize', 25, 'Interpreter', 'latex');

        % Skip the first two legend entries
        lgnd_str{1} = '';
        lgnd_str{2} = '';

    end



    % plot the retrieved droplet profile
    if nn<length(filenames)

        hold on


        % Plot the retrieval uncertainty of the radius at cloud top and
        % bottom
        e1 = errorbar(ds.GN_outputs.retrieval(2,end), ds.GN_outputs.retrieval(1,end), sqrt(ds.GN_outputs.posterior_cov(1,1))/2,...
            sqrt(ds.GN_outputs.posterior_cov(1,1))/2, sqrt(ds.GN_outputs.posterior_cov(2,2))/2,...
            sqrt(ds.GN_outputs.posterior_cov(2,2))/2, 'MarkerFaceColor', C(nn+1,:),...
            'MarkerEdgeColor', C(nn+1,:), 'Linewidth', 4, 'Marker', '.', 'MarkerSize', 40,...
            'Color', C(nn+1,:));

        e1.Bar.LineStyle = 'dotted';

        hold on



    else

        % give a different marker type for the retrieval using 66 bands
        % that also retrieved above cloud column water vapor


        errorbar(ds.GN_outputs.retrieval(2,end), ds.GN_outputs.retrieval(1,end), sqrt(ds.GN_outputs.posterior_cov(1,1))/2,...
            sqrt(ds.GN_outputs.posterior_cov(1,1))/2, sqrt(ds.GN_outputs.posterior_cov(2,2))/2,...
            sqrt(ds.GN_outputs.posterior_cov(2,2))/2, 'MarkerFaceColor', C(nn+1,:),...
            'MarkerEdgeColor', C(nn+1,:), 'Linewidth', 4, 'Marker', '.', 'MarkerSize', 40,...
            'Color', C(nn+1,:))

        hold on


    end



    % create the legend string
    if nn<length(filenames)

        % what was the assumed above cloud column water vapor path?
        assumed_CWV = aboveCloud_CWV_simulated_hysics_spectra(ds.GN_inputs); % kg/m^2

        lgnd_str{nn+2} = [num2str(numel(ds.GN_inputs.bands2run)), ' bands - assumed $acpw$ = ',...
            num2str(round(assumed_CWV,2)), ' $mm$'];





    else

        % create the string for the retrieval using CWV

        % what was the retrieved above cloud column water vapor path above
        % cloud?
        retrieved_CWV = ds.GN_outputs.retrieval(end, end);        % kg/m^2 (mm)

        lgnd_str{nn+2} = [num2str(numel(ds.GN_inputs.bands2run)), ' bands - retrieved $acpw$ = ',...
            num2str(round(retrieved_CWV, 2)), ' $mm$'];

    end

end


% Create a Legend with only the two black curves
legend(lgnd_str, 'Interpreter','latex', 'Location','northwest', 'FontSize', 20)

grid on; grid minor
ylabel('$r_{top}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',30)
xlabel('$r_{bot}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',30)

set(gcf,'Position',[0 0 950 750])







%% Make paneled figure with each optical-depth for a single total-column-water-vapor pair
% Compare the 4 retrievals that assume a total column water vapor amount with the retrieval for
% above cloud column water vapor


% *** DO THESE MATCH THE CALCULATIONS I MADE FOR GRC???? ******

% The files below were computed using the retrieval algorithm on 9-14-2025
% on my local machine



clear variables


% folder_path = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
%     'Hyperspectral_Cloud_Retrievals/HySICS/Droplet_profile_retrievals/'];

folder_path = ['/Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
    'paper2_variableSweep/test_folder_mimic_localMachine_results/'];

% define the mat files for each retreival to plot
filenames = {'dropletRetrieval_noACPW_HySICS_35bands_0.3%_uncert_rTop_9.2516_rBot_5.3192_tauC_6.1312_tcwv_14_tcwv-assumption_10_vza_4_vaz_257_sza_31_saz_96_sim-ran-on-17-Sep-2025_1.mat',...
    'dropletRetrieval_noACPW_HySICS_35bands_0.3%_uncert_rTop_9.2516_rBot_5.3192_tauC_6.1312_tcwv_14_tcwv-assumption_15_vza_4_vaz_257_sza_31_saz_96_sim-ran-on-17-Sep-2025_1.mat',...
    'dropletRetrieval_noACPW_HySICS_35bands_0.3%_uncert_rTop_9.2516_rBot_5.3192_tauC_6.1312_tcwv_14_tcwv-assumption_20_vza_4_vaz_257_sza_31_saz_96_sim-ran-on-17-Sep-2025_1.mat',...
    'dropletRetrieval_noACPW_HySICS_35bands_0.3%_uncert_rTop_9.2516_rBot_5.3192_tauC_6.1312_tcwv_14_tcwv-assumption_25_vza_4_vaz_257_sza_31_saz_96_sim-ran-on-17-Sep-2025_1.mat',...
    'dropletRetrieval_HySICS_66bands_0.3%_uncert_rTop_9.2516_rBot_5.3192_tauC_6.1312_tcwv_14_vza_4_vaz_257_sza_31_saz_96_sim-ran-on-16-Sep-2025_1'};




% Step through each file


% define the colors for each curve plotted
C = mySavedColors(61:(61+length(filenames)+1), 'fixed');

lgnd_str = cell(1, length(filenames) + 2);


figure;


for nn = 1:length(filenames)



    % Load a data set
    ds = load([folder_path, filenames{nn}]);

    if nn==1


        % first, plot the simulated profile values as two lines
        xline(ds.GN_inputs.measurement.r_bot, ':', ['True $r_{bot}$'],...
            'Fontsize',20, 'Interpreter','latex','LineWidth',2,'Color', [87/255, 90/255, 91/255], 'LabelHorizontalAlignment','left',...
            'LabelVerticalAlignment','bottom');

        hold on

        yline(ds.GN_inputs.measurement.r_top, ':', ['True $r_{top}$'],...
            'Fontsize',20, 'Interpreter','latex','LineWidth',2,'Color', [87/255, 90/255, 91/255], 'LabelHorizontalAlignment','right',...
            'LabelVerticalAlignment','top');
        hold on


        % what was the assumed above cloud column water vapor path?
        simulated_CWV = ds.GN_inputs.measurement.actpw; % kg/m^2

        title(['True $acpw$ = ',num2str(round(simulated_CWV, 1)), ' $mm$'],...
            'Fontsize', 25, 'Interpreter', 'latex');

        % Skip the first two legend entries
        lgnd_str{1} = '';
        lgnd_str{2} = '';

    end



    % plot the retrieved droplet profile
    if nn<length(filenames)

        hold on


        % Plot the retrieval uncertainty of the radius at cloud top and
        % bottom
        e1 = errorbar(ds.GN_outputs.retrieval(2,end), ds.GN_outputs.retrieval(1,end), sqrt(ds.GN_outputs.posterior_cov(1,1))/2,...
            sqrt(ds.GN_outputs.posterior_cov(1,1))/2, sqrt(ds.GN_outputs.posterior_cov(2,2))/2,...
            sqrt(ds.GN_outputs.posterior_cov(2,2))/2, 'MarkerFaceColor', C(nn+1,:),...
            'MarkerEdgeColor', C(nn+1,:), 'Linewidth', 4, 'Marker', '.', 'MarkerSize', 40,...
            'Color', C(nn+1,:));

        e1.Bar.LineStyle = 'dotted';

        hold on



    else

        % give a different marker type for the retrieval using 66 bands
        % that also retrieved above cloud column water vapor


        errorbar(ds.GN_outputs.retrieval(2,end), ds.GN_outputs.retrieval(1,end), sqrt(ds.GN_outputs.posterior_cov(1,1))/2,...
            sqrt(ds.GN_outputs.posterior_cov(1,1))/2, sqrt(ds.GN_outputs.posterior_cov(2,2))/2,...
            sqrt(ds.GN_outputs.posterior_cov(2,2))/2, 'MarkerFaceColor', C(nn+1,:),...
            'MarkerEdgeColor', C(nn+1,:), 'Linewidth', 4, 'Marker', '.', 'MarkerSize', 40,...
            'Color', C(nn+1,:))

        hold on


    end



    % create the legend string
    if nn<length(filenames)

        % what was the assumed above cloud column water vapor path?
        assumed_CWV = aboveCloud_CWV_simulated_hysics_spectra(ds.GN_inputs); % kg/m^2

        lgnd_str{nn+2} = [num2str(numel(ds.GN_inputs.bands2run)), ' bands - assumed $acpw$ = ',...
            num2str(round(assumed_CWV,1)), ' $mm$'];





    else

        % create the string for the retrieval using CWV

        % what was the retrieved above cloud column water vapor path above
        % cloud?
        retrieved_CWV = ds.GN_outputs.retrieval(end, end);        % kg/m^2 (mm)

        lgnd_str{nn+2} = [num2str(numel(ds.GN_inputs.bands2run)), ' bands - retrieved $acpw$ = ',...
            num2str(round(retrieved_CWV, 1)), ' $mm$'];

    end

end


% Create a Legend with only the two black curves
legend(lgnd_str, 'Interpreter','latex', 'Location','best', 'FontSize', 20)

grid on; grid minor
ylabel('$r_{top}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',30)
xlabel('$r_{bot}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',30)

set(gcf,'Position',[0 0 950 750])









