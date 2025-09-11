%% Retrieval analysis for subset of cases - for paper 2

% By Andrew John Buggee

%% Create panneled figure showing retrieval of droplet profile and optical depth for different ACPW

% *** 1% uncertainty ***

clear variables


% load 1% uncertainty files
filenames_startWith = 'dropletRetrieval_noACPW_HySICS_35bands_1%_uncert_rTop_10_rBot_5';

folder_path = ['/Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
    'paper2_variableSweep/rTop_10/vza_7_vaz_210_sza_10_saz_91_subset/'];

% % Access specific file or folder
% filenames_noACPW_1percent = dir(['//Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
%     'paper2_variableSweep/rTop_10/vza_7_vaz_210_sza_10_saz_91_subset/', filenames_startWith,...
%     '*.mat']);

% filenames_full_1percent = dir(['//Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
%     'paper2_variableSweep/rTop_10/vza_7_vaz_210_sza_10_saz_91_subset/',...
%     'dropletRetrieval_HySICS_66bands_1%*.mat']);

% what are the free parameters?
r_top = 10;
r_bot = 5;
tau_c = [5,11,17,23];
tcpw = [8, 14, 20];

acpw_true = [2.8932, 4.6291, 6.3650, 8.1010,...
    9.8369, 11.5728, 13.3088, 15.0447, 16.7806,...
    18.5165, 20.2525];

% length of each independent variable
num_tauC = length(tau_c);
num_tcpw = length(tcpw);




% Step through each file

% define the colors for each curve plotted
C = mySavedColors(61:65, 'fixed');

lgnd_str = cell(1, 7);


for pw = 1:length(tcpw)


    for tc = 1:length(tau_c)

        figure;


        % Load a data set
        fileNames = dir([folder_path, filenames_startWith, '_tauC_', num2str(tau_c(tc)),...
            '_tcwv_',  num2str(tcpw(pw)), '*.mat' ]);

        for nn = 1:length(fileNames)


            % Load a data set
            ds = load([fileNames(nn).folder, '/', fileNames(nn).name]);


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
                title(['Simulated profile - $acpw$ = ',num2str(round(ds.GN_inputs.measurement.actpw, 2)), ' $mm$'],...
                    'Fontsize', 25, 'Interpreter', 'latex');

                % Skip the first two legend entries
                lgnd_str{1} = '';
                lgnd_str{2} = '';

            end



            % plot the retrieved droplet profile
            if nn<length(fileNames)

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
            if nn<length(fileNames)

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




    end

end









%% Plot every optical-depth and total-column-water-vapor pair as a unique figure. Compare the 4 retrievals
% that assume a total column water vapor amount with the retrieval for
% above cloud column water vapor

% *** 5% uncertainty ***

clear variables


% load 5% uncertainty files
filenames_noACPW_startsWith = 'dropletRetrieval_noACPW_HySICS_35bands_5%_uncert_rTop_10_rBot_5';
filenames_full_startsWith = 'dropletRetrieval_HySICS_66bands_5%_uncert_rTop_10_rBot_5';

folder_path = ['/Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
    'paper2_variableSweep/rTop_10/vza_7_vaz_210_sza_10_saz_91_subset/'];

% % Access specific file or folder
% filenames_noACPW_1percent = dir(['//Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
%     'paper2_variableSweep/rTop_10/vza_7_vaz_210_sza_10_saz_91_subset/', filenames_startWith,...
%     '*.mat']);

% filenames_full_1percent = dir(['//Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
%     'paper2_variableSweep/rTop_10/vza_7_vaz_210_sza_10_saz_91_subset/',...
%     'dropletRetrieval_HySICS_66bands_1%*.mat']);

% what are the free parameters?
r_top = 10;
r_bot = 5;
tau_c = [5,11,17,23];
tcpw = [8, 14, 20];

acpw_true = [2.8932, 4.6291, 6.3650, 8.1010,...
    9.8369, 11.5728, 13.3088, 15.0447, 16.7806,...
    18.5165, 20.2525];

% length of each independent variable
num_tauC = length(tau_c);
num_tcpw = length(tcpw);




% Step through each file

% define the colors for each curve plotted
C = mySavedColors(61:66, 'fixed');

lgnd_str = cell(1, 7);


for pw = 1:length(tcpw)


    for tc = 1:length(tau_c)

        figure;


        % Load all 5 retrievals
        fileNames = [dir([folder_path, filenames_noACPW_startsWith, '_tauC_', num2str(tau_c(tc)),...
            '_tcwv_',  num2str(tcpw(pw)), '*.mat']);...
            dir([folder_path, filenames_full_startsWith, '_tauC_', num2str(tau_c(tc)),...
            '_tcwv_',  num2str(tcpw(pw)), '*.mat']);];

        for nn = 1:length(fileNames)


            % Load a data set
            ds = load([fileNames(nn).folder, '/', fileNames(nn).name]);


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
                title(['Simulated profile - $acpw$ = ',num2str(round(ds.GN_inputs.measurement.actpw, 2)), ' $mm$'],...
                    'Fontsize', 25, 'Interpreter', 'latex');

                % Skip the first two legend entries
                lgnd_str{1} = '';
                lgnd_str{2} = '';

            end



            % plot the retrieved droplet profile
            if nn<length(fileNames)

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
            if nn<length(fileNames)

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
        legend(lgnd_str, 'Interpreter','latex', 'Location','best', 'FontSize', 20)

        grid on; grid minor
        ylabel('$r_{top}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',30)
        xlabel('$r_{bot}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',30)

        set(gcf,'Position',[0 0 950 850])




    end

end













%% Make paneled figure with each optical-depth for a single total-column-water-vapor pair
% Compare the 4 retrievals that assume a total column water vapor amount with the retrieval for
% above cloud column water vapor

% *** 5% uncertainty ***

clear variables


% load 5% uncertainty files
filenames_noACPW_startsWith = 'dropletRetrieval_noACPW_HySICS_35bands_5%_uncert_rTop_10_rBot_5';
filenames_full_startsWith = 'dropletRetrieval_HySICS_66bands_5%_uncert_rTop_10_rBot_5';

folder_path = ['/Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
    'paper2_variableSweep/rTop_10/vza_7_vaz_210_sza_10_saz_91_subset2/'];

% % Access specific file or folder
% filenames_noACPW_1percent = dir(['//Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
%     'paper2_variableSweep/rTop_10/vza_7_vaz_210_sza_10_saz_91_subset/', filenames_startWith,...
%     '*.mat']);

% filenames_full_1percent = dir(['//Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
%     'paper2_variableSweep/rTop_10/vza_7_vaz_210_sza_10_saz_91_subset/',...
%     'dropletRetrieval_HySICS_66bands_1%*.mat']);

% what are the free parameters?
r_top = 10;
r_bot = 5;
tau_c = [5,11,17,23];
tcpw = [8, 14, 20];

acpw_true = [2.8932, 4.6291, 6.3650, 8.1010,...
    9.8369, 11.5728, 13.3088, 15.0447, 16.7806,...
    18.5165, 20.2525];

% length of each independent variable
num_tauC = length(tau_c);
num_tcpw = length(tcpw);




% Step through each file

% define the colors for each curve plotted
C = mySavedColors(61:66, 'fixed');




for pw = 1:length(tcpw)

    f = figure;


    for tc = 1:length(tau_c)

        subplot(1,4,tc)



        % Load all 5 retrievals
        fileNames = [dir([folder_path, filenames_noACPW_startsWith, '_tauC_', num2str(tau_c(tc)),...
            '_tcwv_',  num2str(tcpw(pw)), '*.mat']);...
            dir([folder_path, filenames_full_startsWith, '_tauC_', num2str(tau_c(tc)),...
            '_tcwv_',  num2str(tcpw(pw)), '*.mat']);];

        lgnd_str = cell(1, length(fileNames)+2);

        for nn = 1:length(fileNames)


            % Load a data set
            ds = load([fileNames(nn).folder, '/', fileNames(nn).name]);


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




                % Skip the first two legend entries
                lgnd_str{1} = '';
                lgnd_str{2} = '';

                % write titles and subtitles
                if tc==1

                    % what was the assumed above cloud column water vapor path?
                    title(['Simulated measurements with $5\%$ uncertainty - $acpw$ = ',num2str(round(ds.GN_inputs.measurement.actpw, 2)), ' $mm$'],...
                        'Fontsize', 25, 'Interpreter', 'latex');

                    % Create textbox
                    annotation(f,'textbox',...
                        [0.25 0.841666666666667 0.0369591836734694 0.0716666666666671],...
                        'String',{['$\tau_c = $ ', num2str(tau_c(tc))]},...
                        'Interpreter','latex',...
                        'FontSize',23,...
                        'FitBoxToText','off',...
                        'EdgeColor','none');

                elseif tc==2


                    % Create textbox
                    annotation(f,'textbox',...
                        [0.45 0.841666666666667 0.0369591836734694 0.0716666666666671],...
                        'String',{['$\tau_c = $ ', num2str(tau_c(tc))]},...
                        'Interpreter','latex',...
                        'FontSize',23,...
                        'FitBoxToText','off',...
                        'EdgeColor','none');

                elseif tc==3


                    % Create textbox
                    annotation(f,'textbox',...
                        [0.66 0.841666666666667 0.0369591836734694 0.0716666666666671],...
                        'String',{['$\tau_c = $ ', num2str(tau_c(tc))]},...
                        'Interpreter','latex',...
                        'FontSize',23,...
                        'FitBoxToText','off',...
                        'EdgeColor','none');

                elseif tc==4

                    % Create textbox
                    annotation(f,'textbox',...
                        [0.87 0.841666666666667 0.0369591836734694 0.0716666666666671],...
                        'String',{['$\tau_c = $ ', num2str(tau_c(tc))]},...
                        'Interpreter','latex',...
                        'FontSize',23,...
                        'FitBoxToText','off',...
                        'EdgeColor','none');



                end






            end



            % plot the retrieved droplet profile
            if nn<length(fileNames)

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
            if nn<length(fileNames)

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
        legend(lgnd_str, 'Interpreter','latex', 'FontSize', 18,...
            'Position',[0.166173414792696 0.16735419938909 0.142261290258291 0.153166666030884])

        grid on; grid minor
        ylabel('$r_{top}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',30)
        xlabel('$r_{bot}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',30)

        % set x and y lims
        % set figure limits
        dx = 2;
        dy = 1;

        xlim([ds.GN_inputs.measurement.r_bot - dx, ds.GN_inputs.measurement.r_bot + dx])
        ylim([ds.GN_inputs.measurement.r_top - dy, ds.GN_inputs.measurement.r_top + dy])



    end

    set(gcf,'Position',[0 0 2450 600])

end






%% Make paneled figure with each optical-depth for a single total-column-water-vapor pair
% Compare the 4 retrievals that assume a total column water vapor amount with the retrieval for
% above cloud column water vapor

% *** 1% uncertainty ***

clear variables


% load 5% uncertainty files
filenames_noACPW_startsWith = 'dropletRetrieval_noACPW_HySICS_35bands_1%_uncert_rTop_10_rBot_5';
filenames_full_startsWith = 'dropletRetrieval_HySICS_66bands_1%_uncert_rTop_10_rBot_5';

folder_path = ['/Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
    'paper2_variableSweep/rTop_10/vza_7_vaz_210_sza_10_saz_91_subset2/'];

% % Access specific file or folder
% filenames_noACPW_1percent = dir(['//Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
%     'paper2_variableSweep/rTop_10/vza_7_vaz_210_sza_10_saz_91_subset/', filenames_startWith,...
%     '*.mat']);

% filenames_full_1percent = dir(['//Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
%     'paper2_variableSweep/rTop_10/vza_7_vaz_210_sza_10_saz_91_subset/',...
%     'dropletRetrieval_HySICS_66bands_1%*.mat']);

% what are the free parameters?
r_top = 10;
r_bot = 5;
tau_c = [5,11,17,23];
tcpw = [8, 14, 20];

acpw_true = [2.8932, 4.6291, 6.3650, 8.1010,...
    9.8369, 11.5728, 13.3088, 15.0447, 16.7806,...
    18.5165, 20.2525];

% length of each independent variable
num_tauC = length(tau_c);
num_tcpw = length(tcpw);




% Step through each file

% define the colors for each curve plotted
C = mySavedColors(61:66, 'fixed');




for pw = 1:length(tcpw)

    f = figure;


    for tc = 1:length(tau_c)

        subplot(1,4,tc)



        % Load all 5 retrievals
        fileNames = [dir([folder_path, filenames_noACPW_startsWith, '_tauC_', num2str(tau_c(tc)),...
            '_tcwv_',  num2str(tcpw(pw)), '*.mat']);...
            dir([folder_path, filenames_full_startsWith, '_tauC_', num2str(tau_c(tc)),...
            '_tcwv_',  num2str(tcpw(pw)), '*.mat']);];

        lgnd_str = cell(1, length(fileNames)+2);

        for nn = 1:length(fileNames)


            % Load a data set
            ds = load([fileNames(nn).folder, '/', fileNames(nn).name]);


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




                % Skip the first two legend entries
                lgnd_str{1} = '';
                lgnd_str{2} = '';

                % write titles and subtitles
                if tc==1

                    % what was the assumed above cloud column water vapor path?
                    title(['Simulated measurements with $1\%$ uncertainty - $acpw$ = ',num2str(round(ds.GN_inputs.measurement.actpw, 2)), ' $mm$'],...
                        'Fontsize', 25, 'Interpreter', 'latex');

                    % Create textbox
                    annotation(f,'textbox',...
                        [0.25 0.841666666666667 0.0369591836734694 0.0716666666666671],...
                        'String',{['$\tau_c = $ ', num2str(tau_c(tc))]},...
                        'Interpreter','latex',...
                        'FontSize',23,...
                        'FitBoxToText','off',...
                        'EdgeColor','none');

                elseif tc==2


                    % Create textbox
                    annotation(f,'textbox',...
                        [0.45 0.841666666666667 0.0369591836734694 0.0716666666666671],...
                        'String',{['$\tau_c = $ ', num2str(tau_c(tc))]},...
                        'Interpreter','latex',...
                        'FontSize',23,...
                        'FitBoxToText','off',...
                        'EdgeColor','none');

                elseif tc==3


                    % Create textbox
                    annotation(f,'textbox',...
                        [0.66 0.841666666666667 0.0369591836734694 0.0716666666666671],...
                        'String',{['$\tau_c = $ ', num2str(tau_c(tc))]},...
                        'Interpreter','latex',...
                        'FontSize',23,...
                        'FitBoxToText','off',...
                        'EdgeColor','none');

                elseif tc==4

                    % Create textbox
                    annotation(f,'textbox',...
                        [0.87 0.841666666666667 0.0369591836734694 0.0716666666666671],...
                        'String',{['$\tau_c = $ ', num2str(tau_c(tc))]},...
                        'Interpreter','latex',...
                        'FontSize',23,...
                        'FitBoxToText','off',...
                        'EdgeColor','none');



                end






            end



            % plot the retrieved droplet profile
            if nn<length(fileNames)

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
            if nn<length(fileNames)

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
        legend(lgnd_str, 'Interpreter','latex', 'FontSize', 18,...
            'Position',[0.166173414792696 0.16735419938909 0.142261290258291 0.153166666030884])
        

        grid on; grid minor
        ylabel('$r_{top}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',30)
        xlabel('$r_{bot}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',30)

        % set x and y lims
        % set figure limits
        dx = 2;
        dy = 1;

        xlim([ds.GN_inputs.measurement.r_bot - dx, ds.GN_inputs.measurement.r_bot + dx])
        ylim([ds.GN_inputs.measurement.r_top - dy, ds.GN_inputs.measurement.r_top + dy])



    end

    
    if strcmp(whatComputer, 'anbu8374')==true

        set(gcf,'Position',[0 0 2450 600])

    elseif strcmp(whatComputer, 'andrewbuggee')==true

        set(gcf,'Position',[0 0 1200 600])

    end


end








%% Make paneled figure with each optical-depth for a single total-column-water-vapor pair
% Compare the 4 retrievals that assume a total column water vapor amount with the retrieval for
% above cloud column water vapor

% *** 0.1% uncertainty ***

clear variables


% load 5% uncertainty files
filenames_noACPW_startsWith = 'dropletRetrieval_noACPW_HySICS_35bands_0.1%_uncert_rTop_10_rBot_5';
filenames_full_startsWith = 'dropletRetrieval_HySICS_66bands_0.1%_uncert_rTop_10_rBot_5';

if strcmp(whatComputer, 'anbu8374')==true

folder_path = ['/Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
    'paper2_variableSweep/rTop_10/vza_7_vaz_210_sza_10_saz_91_subset2/'];

elseif strcmp(whatComputer, 'andrewbuggee')==true

    folder_path = ['/Users/andrewbuggee/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
    'paper2_variableSweep/rTop_10/vza_7_vaz_210_sza_10_saz_91_subset2/'];

end


% % Access specific file or folder
% filenames_noACPW_1percent = dir(['//Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
%     'paper2_variableSweep/rTop_10/vza_7_vaz_210_sza_10_saz_91_subset/', filenames_startWith,...
%     '*.mat']);

% filenames_full_1percent = dir(['//Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
%     'paper2_variableSweep/rTop_10/vza_7_vaz_210_sza_10_saz_91_subset/',...
%     'dropletRetrieval_HySICS_66bands_1%*.mat']);

% what are the free parameters?
r_top = 10;
r_bot = 5;
tau_c = [5,11,17,23];
tcpw = [8, 14, 20];

acpw_true = [2.8932, 4.6291, 6.3650, 8.1010,...
    9.8369, 11.5728, 13.3088, 15.0447, 16.7806,...
    18.5165, 20.2525];

% length of each independent variable
num_tauC = length(tau_c);
num_tcpw = length(tcpw);




% Step through each file

% define the colors for each curve plotted
C = mySavedColors(61:66, 'fixed');




for pw = 1:length(tcpw)

    f = figure;


    for tc = 1:length(tau_c)

        subplot(1,4,tc)



        % Load all 5 retrievals
        fileNames = [dir([folder_path, filenames_noACPW_startsWith, '_tauC_', num2str(tau_c(tc)),...
            '_tcwv_',  num2str(tcpw(pw)), '*.mat']);...
            dir([folder_path, filenames_full_startsWith, '_tauC_', num2str(tau_c(tc)),...
            '_tcwv_',  num2str(tcpw(pw)), '*.mat']);];

        lgnd_str = cell(1, length(fileNames)+2);

        for nn = 1:length(fileNames)


            % Load a data set
            ds = load([fileNames(nn).folder, '/', fileNames(nn).name]);


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




                % Skip the first two legend entries
                lgnd_str{1} = '';
                lgnd_str{2} = '';

                % write titles and subtitles
                if tc==1

                    % what was the assumed above cloud column water vapor path?
                    title(['Simulated measurements with $0.1\%$ uncertainty - $acpw$ = ',num2str(round(ds.GN_inputs.measurement.actpw, 2)), ' $mm$'],...
                        'Fontsize', 25, 'Interpreter', 'latex');

                    % Create textbox
                    annotation(f,'textbox',...
                        [0.25 0.841666666666667 0.0369591836734694 0.0716666666666671],...
                        'String',{['$\tau_c = $ ', num2str(tau_c(tc))]},...
                        'Interpreter','latex',...
                        'FontSize',23,...
                        'FitBoxToText','off',...
                        'EdgeColor','none');

                elseif tc==2


                    % Create textbox
                    annotation(f,'textbox',...
                        [0.45 0.841666666666667 0.0369591836734694 0.0716666666666671],...
                        'String',{['$\tau_c = $ ', num2str(tau_c(tc))]},...
                        'Interpreter','latex',...
                        'FontSize',23,...
                        'FitBoxToText','off',...
                        'EdgeColor','none');

                elseif tc==3


                    % Create textbox
                    annotation(f,'textbox',...
                        [0.66 0.841666666666667 0.0369591836734694 0.0716666666666671],...
                        'String',{['$\tau_c = $ ', num2str(tau_c(tc))]},...
                        'Interpreter','latex',...
                        'FontSize',23,...
                        'FitBoxToText','off',...
                        'EdgeColor','none');

                elseif tc==4

                    % Create textbox
                    annotation(f,'textbox',...
                        [0.87 0.841666666666667 0.0369591836734694 0.0716666666666671],...
                        'String',{['$\tau_c = $ ', num2str(tau_c(tc))]},...
                        'Interpreter','latex',...
                        'FontSize',23,...
                        'FitBoxToText','off',...
                        'EdgeColor','none');



                end






            end



            % plot the retrieved droplet profile
            if nn<length(fileNames)

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
            if nn<length(fileNames)

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
        legend(lgnd_str, 'Interpreter','latex', 'FontSize', 18,...
            'Position',[0.166173414792696 0.16735419938909 0.142261290258291 0.153166666030884])

        grid on; grid minor
        ylabel('$r_{top}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',30)
        xlabel('$r_{bot}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',30)

        % set x and y lims
        % set figure limits
        dx = 2;
        dy = 1;

        xlim([ds.GN_inputs.measurement.r_bot - dx, ds.GN_inputs.measurement.r_bot + dx])
        ylim([ds.GN_inputs.measurement.r_top - dy, ds.GN_inputs.measurement.r_top + dy])



    end


    if strcmp(whatComputer, 'anbu8374')==true

        set(gcf,'Position',[0 0 2450 600])

    elseif strcmp(whatComputer, 'andrewbuggee')==true

        set(gcf,'Position',[0 0 1200 600])

    end

end






%% Make paneled figure with each optical-depth for a single total-column-water-vapor pair
% Compare the 4 retrievals that assume a total column water vapor amount with the retrieval for
% above cloud column water vapor



clear variables

% define uncertainty
plt_uncert = 1;            % percent

% load 5% uncertainty files
filenames_noACPW_startsWith = ['dropletRetrieval_noACPW_HySICS_35bands_',...
    num2str(plt_uncert), '%_uncert_rTop_10_rBot_5'];
filenames_full_startsWith = ['dropletRetrieval_HySICS_66bands_',...
    num2str(plt_uncert), '%_uncert_rTop_10_rBot_5'];

folder_path = ['/Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
    'paper2_variableSweep/rTop_10/vza_4_vaz_257_sza_31_saz_96_subset_lowUncert/'];

% % Access specific file or folder
% filenames_noACPW_1percent = dir(['//Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
%     'paper2_variableSweep/rTop_10/vza_7_vaz_210_sza_10_saz_91_subset/', filenames_startWith,...
%     '*.mat']);

% filenames_full_1percent = dir(['//Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
%     'paper2_variableSweep/rTop_10/vza_7_vaz_210_sza_10_saz_91_subset/',...
%     'dropletRetrieval_HySICS_66bands_1%*.mat']);

% what are the free parameters?
r_top = 10;
r_bot = 5;
tau_c = [5,11,17,23];
tcpw = [8, 14, 20];

acpw_true = [2.8932, 4.6291, 6.3650, 8.1010,...
    9.8369, 11.5728, 13.3088, 15.0447, 16.7806,...
    18.5165, 20.2525];

% length of each independent variable
num_tauC = length(tau_c);
num_tcpw = length(tcpw);




% Step through each file

% define the colors for each curve plotted
C = mySavedColors(61:66, 'fixed');




for pw = 1:length(tcpw)

    f = figure;


    for tc = 1:length(tau_c)

        subplot(1,4,tc)



        % Load all 5 retrievals
        fileNames = [dir([folder_path, filenames_noACPW_startsWith, '_tauC_', num2str(tau_c(tc)),...
            '_tcwv_',  num2str(tcpw(pw)), '*.mat']);...
            dir([folder_path, filenames_full_startsWith, '_tauC_', num2str(tau_c(tc)),...
            '_tcwv_',  num2str(tcpw(pw)), '*.mat']);];

        lgnd_str = cell(1, length(fileNames)+2);

        for nn = 1:length(fileNames)


            % Load a data set
            ds = load([fileNames(nn).folder, '/', fileNames(nn).name]);


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




                % Skip the first two legend entries
                lgnd_str{1} = '';
                lgnd_str{2} = '';

                % write titles and subtitles
                if tc==1

                    % what was the assumed above cloud column water vapor path?
                    title(['Simulated measurements with $', num2str(plt_uncert), '\%$ uncertainty - $acpw$ = ',...
                        num2str(round(ds.GN_inputs.measurement.actpw, 2)), ' $mm$'],...
                        'Fontsize', 25, 'Interpreter', 'latex');

                    % Create textbox
                    annotation(f,'textbox',...
                        [0.25 0.841666666666667 0.0369591836734694 0.0716666666666671],...
                        'String',{['$\tau_c = $ ', num2str(tau_c(tc))]},...
                        'Interpreter','latex',...
                        'FontSize',23,...
                        'FitBoxToText','off',...
                        'EdgeColor','none');

                elseif tc==2


                    % Create textbox
                    annotation(f,'textbox',...
                        [0.45 0.841666666666667 0.0369591836734694 0.0716666666666671],...
                        'String',{['$\tau_c = $ ', num2str(tau_c(tc))]},...
                        'Interpreter','latex',...
                        'FontSize',23,...
                        'FitBoxToText','off',...
                        'EdgeColor','none');

                elseif tc==3


                    % Create textbox
                    annotation(f,'textbox',...
                        [0.66 0.841666666666667 0.0369591836734694 0.0716666666666671],...
                        'String',{['$\tau_c = $ ', num2str(tau_c(tc))]},...
                        'Interpreter','latex',...
                        'FontSize',23,...
                        'FitBoxToText','off',...
                        'EdgeColor','none');

                elseif tc==4

                    % Create textbox
                    annotation(f,'textbox',...
                        [0.87 0.841666666666667 0.0369591836734694 0.0716666666666671],...
                        'String',{['$\tau_c = $ ', num2str(tau_c(tc))]},...
                        'Interpreter','latex',...
                        'FontSize',23,...
                        'FitBoxToText','off',...
                        'EdgeColor','none');



                end






            end



            % plot the retrieved droplet profile
            if strcmp(fileNames(nn).name(1:30), 'dropletRetrieval_noACPW_HySICS')==true

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



            elseif strcmp(fileNames(nn).name(1:23), 'dropletRetrieval_HySICS')==true

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
            if nn<length(fileNames)

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
        legend(lgnd_str, 'Interpreter','latex', 'FontSize', 18,...
            'Position',[0.166173414792696 0.16735419938909 0.142261290258291 0.153166666030884])

        grid on; grid minor
        ylabel('$r_{top}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',30)
        xlabel('$r_{bot}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',30)

        % set x and y lims
        % set figure limits
        dx = 2;
        dy = 1;

        xlim([ds.GN_inputs.measurement.r_bot - dx, ds.GN_inputs.measurement.r_bot + dx])
        ylim([ds.GN_inputs.measurement.r_top - dy, ds.GN_inputs.measurement.r_top + dy])



    end

    set(gcf,'Position',[0 0 2450 600])

end



