%% Compare the 4 retrievals that assume a total column water vapor amount with the retrieval for
% above cloud column water vapor



% The files below were computed using the retrieval algorithm on 9-16-2025
% on my local machine



clear variables


% load 5% uncertainty files
filenames_noACPW_startsWith = 'dropletRetrieval_noACPW_HySICS_35bands_0.3%_uncert_';
filenames_full_startsWith = 'dropletRetrieval_HySICS_66bands_0.3%_uncert_';

if strcmp(whatComputer, 'anbu8374')==true

    folder_path = ['/Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
        'paper2_variableSweep/rTop_10/vza_4_vaz_257_sza_31_saz_96_subset_newRetrieval/'];

elseif strcmp(whatComputer, 'andrewbuggee')==true

    folder_path = ['/Users/andrewbuggee/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
        'paper2_variableSweep/rTop_10/vza_4_vaz_257_sza_31_saz_96_subset_newRetrieval2/'];

end



% what are the free parameters?
% r_top = 10;
% r_bot = 5;
% tau_c = [5,11,17,23];
% tau_c = 5;
% tcpw = [8, 14, 20];

r_top = 9.2516;
r_bot = 5.3192;
tau_c = [6.1312, 9.9317, 12.4356, 15.3172];
tcpw = [8.123, 10.6422, 13.6234, 15.8543, 18.9824];

tau_idx = 4;
tcpw_idx = 2;

% define the mat files for each retreival to plot
% Load all 5 retrievals
filenames = [dir([folder_path, filenames_noACPW_startsWith, 'rTop_', num2str(r_top),...
    '_rBot_', num2str(r_bot), '_tauC_', num2str(tau_c(tau_idx)),...
    '_tcwv_',  num2str(tcpw(tcpw_idx)), '*.mat']);...
    dir([folder_path, filenames_full_startsWith, 'rTop_', num2str(r_top),...
    '_rBot_', num2str(r_bot), '_tauC_', num2str(tau_c(tau_idx)),...
    '_tcwv_',  num2str(tcpw(tcpw_idx)), '*.mat']);];




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








%% Plot the retrieved r_top - r_bottom versus optical depth and color each marker according to the RSS residual

% The files below were computed using the retrieval algorithm on 9-16-2025
% on my local machine



clear variables


% load 5% uncertainty files
filenames_noACPW_startsWith = 'dropletRetrieval_noACPW_HySICS_35bands_0.3%_uncert_';
filenames_full_startsWith = 'dropletRetrieval_HySICS_66bands_0.3%_uncert_';

if strcmp(whatComputer, 'anbu8374')==true

    folder_path = ['/Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
        'paper2_variableSweep/rTop_10/vza_4_vaz_257_sza_31_saz_96_subset_newRetrieval/'];

elseif strcmp(whatComputer, 'andrewbuggee')==true

    folder_path = ['/Users/andrewbuggee/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
        'paper2_variableSweep/rTop_10/vza_4_vaz_257_sza_31_saz_96_subset_newRetrieval2/'];

end


% what are the free parameters?
% r_top = 10;
% r_bot = 5;
% % tau_c = [5,11,17,23];
% tau_c = 5;
% tcpw = [8, 14, 20];

r_top = 9.2516;
r_bot = 5.3192;
tau_c = [6.1312, 9.9317, 12.4356, 15.3172];
tcpw = [8.123, 10.6422, 13.6234, 15.8543, 18.9824];

tau_idx = 1;
tcpw_idx = 3;

% define the mat files for each retreival to plot
% Load all 5 retrievals
filenames = [dir([folder_path, filenames_noACPW_startsWith, 'rTop_', num2str(r_top),...
    '_rBot_', num2str(r_bot), '_tauC_', num2str(tau_c(tau_idx)),...
    '_tcwv_',  num2str(tcpw(tcpw_idx)), '*.mat']);...
    dir([folder_path, filenames_full_startsWith, 'rTop_', num2str(r_top),...
    '_rBot_', num2str(r_bot), '_tauC_', num2str(tau_c(tau_idx)),...
    '_tcwv_',  num2str(tcpw(tcpw_idx)), '*.mat']);];




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
        xline(ds.GN_inputs.measurement.tau_c, ':', ['True $\tau_c$'],...
            'Fontsize',20, 'Interpreter','latex','LineWidth',2,'Color', [87/255, 90/255, 91/255], 'LabelHorizontalAlignment','left',...
            'LabelVerticalAlignment','bottom');

        hold on

        yline((ds.GN_inputs.measurement.r_top - ds.GN_inputs.measurement.r_bot), ':', ['True $r_{top} - r_{bot}$'],...
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


        % plot the difference between r-top and r-bot versus cloud optical
        % depth
%         plot(ds.GN_outputs.retrieval(3,end), (ds.GN_outputs.retrieval(1,end)-ds.GN_outputs.retrieval(2,end)),...
%             'MarkerSize', 40, 'Marker', '.', 'MarkerFaceColor', C(nn+1,:),'MarkerEdgeColor', C(nn+1,:))
        e1 = errorbar(ds.GN_outputs.retrieval(3,end), (ds.GN_outputs.retrieval(1,end)-ds.GN_outputs.retrieval(2,end)),...
            1/2 * sqrt(ds.GN_outputs.posterior_cov(1,1) + ds.GN_outputs.posterior_cov(2,3)),...
            1/2 * sqrt(ds.GN_outputs.posterior_cov(1,1) + ds.GN_outputs.posterior_cov(2,3)),...
            sqrt(ds.GN_outputs.posterior_cov(3,3))/2,...
            sqrt(ds.GN_outputs.posterior_cov(3,3))/2, 'MarkerFaceColor', C(nn+1,:),...
            'MarkerEdgeColor', C(nn+1,:), 'Linewidth', 4, 'Marker', '.', 'MarkerSize', 40,...
            'Color', C(nn+1,:));

        e1.Bar.LineStyle = 'dotted';

        hold on



    else

        % give a different marker type for the retrieval using 66 bands
        % that also retrieved above cloud column water vapor

        %         plot(ds.GN_outputs.retrieval(3,end), (ds.GN_outputs.retrieval(1,end)-ds.GN_outputs.retrieval(2,end)),...
        %             'MarkerSize', 30, 'Marker', '*', 'MarkerFaceColor', C(nn+1,:),'MarkerEdgeColor', C(nn+1,:))


        e1 = errorbar(ds.GN_outputs.retrieval(3,end), (ds.GN_outputs.retrieval(1,end)-ds.GN_outputs.retrieval(2,end)),...
            1/2 * sqrt(ds.GN_outputs.posterior_cov(1,1) + ds.GN_outputs.posterior_cov(2,3)),...
            1/2 * sqrt(ds.GN_outputs.posterior_cov(1,1) + ds.GN_outputs.posterior_cov(2,3)),...
            sqrt(ds.GN_outputs.posterior_cov(3,3))/2,...
            sqrt(ds.GN_outputs.posterior_cov(3,3))/2, 'MarkerFaceColor', C(nn+1,:),...
            'MarkerEdgeColor', C(nn+1,:), 'Linewidth', 4, 'Marker', '.', 'MarkerSize', 40,...
            'Color', C(nn+1,:));

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
ylabel('$r_{top} - r_{bot}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',30)
xlabel('$\tau_c$ ', 'Interpreter','latex', 'FontSize',30)

set(gcf,'Position',[0 0 950 750])


