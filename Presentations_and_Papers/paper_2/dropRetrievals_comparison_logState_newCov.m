%% There are 840 files with different r_top, r_bot_ tau_c and acpw combinations
% Make plots with the same r_top and r_bot but where the colors represent
% different tau_c and acpw combinations. Solid markers represent full
% retrieval and dashed markers represent the assumed ACPW retrieval



clear variables


% Determine which computer you're using
which_computer = whatComputer();

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------

    % ---- Define where the retrievals are stored ---
    folder_paths.drive.fullRetrieval = ['/Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
        'paper2_variableSweep/newRetrieval_logTransform_newCov/'];

    folder_paths.drive.noACPW_10 = ['/Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
        'paper2_variableSweep/rTop_10/vza_4_vaz_257_sza_31_saz_96_subset_logState_1/'];


elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

     % ---- Define where the retrievals are stored ---
    folder_paths.drive.fullRetrieval = ['/Users/andrewbuggee/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
        'paper2_variableSweep/newRetrieval_logTransform_newCov/'];

    folder_paths.drive.noACPW_10 = ['/Users/andrewbuggee/MATLAB-Drive/HySICS/Droplet_profile_retrievals/paper2_variableSweep/rTop_10/vza_4_vaz_257_sza_31_saz_96_subset_logState_1/'];



elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------

end


% First step through all files in the directory and remove invisble files
% or directories

% Grab filenames in drive
filenames_fullRetrieval = dir(folder_paths.drive.fullRetrieval);
idx_2delete = [];
for nn = 1:length(filenames_fullRetrieval)

    if strcmp(filenames_fullRetrieval(nn).name(1), 'd')~=true

        idx_2delete = [idx_2delete, nn];

    end

end

% delete rows that don't have retrieval filenames
filenames_fullRetrieval(idx_2delete) = [];


% define the 4 variables and the values that were simulated
r_top = [7.3418, 8.4296, 9.2516, 10.8153, 11.6297];
r_bot = [4.6147, 5.3192, 6.8372, 7.0135];
tau_c = [6.1312, 9.9317, 12.4356, 15.3172, 18.2468, 21.3466, 24.1697];
tcpw = [5.3741, 8.123, 10.6422, 13.6234, 15.8543, 18.9824];
acpw = [3.1097, 4.7003, 6.1580, 7.8831, 9.1740, 10.9840];

% Find all files within filenames_fullRetrieval with the same r_top and
% r_bot
idx_profiles = cell(length(r_top), length(r_top));

% Define the linewidth
ln_wdth = 2;

% define marker width for the circle
circ_size = 35;
% define size for square and diamond
sq_dmnd_size = 12;

for rt = 1:length(r_top)

    for rb = 1:length(r_bot)


        % Find files matching current r_top and r_bot
        filePattern = sprintf('dropletRetrieval_HySICS_66bands_0.3%%_uncert_rTop_%.4f_rBot_%.4f_*.mat', r_top(rt), r_bot(rb));
        matchingFiles = dir(fullfile(folder_paths.drive.fullRetrieval, filePattern));



        % define the colors for each curve plotted
        C = mySavedColors(1:(1+length(tau_c) +1), 'fixed');

        lgnd_str = cell(1, length(tau_c) + length(tcpw) + 2);

        figure;



        % Process each matching file
        for mf = 1:length(matchingFiles)


            % Load the data from the file
            ds = load(fullfile(matchingFiles(mf).folder, matchingFiles(mf).name));


            if mf==1

                % first, plot the simulated profile values as two lines
                xline(ds.GN_inputs.measurement.r_bot, ':', ['True $r_{bot}$'],...
                    'Fontsize',20, 'Interpreter','latex','LineWidth',2,'Color', [147/255, 150/255, 151/255], 'LabelHorizontalAlignment','left',...
                    'LabelVerticalAlignment','bottom');

                hold on

                yline(ds.GN_inputs.measurement.r_top, ':', ['True $r_{top}$'],...
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



            % determine the optical deph
            opticalDepth = ds.GN_inputs.measurement.tau_c; % Extract optical depth from loaded data
            % extract the acpw
            aboveCloud_pw = ds.GN_inputs.measurement.actpw;

            [~, idx_tauC] = min(abs(tau_c - opticalDepth));
            [~,idx_acpw] = min(abs(acpw - aboveCloud_pw));

            if idx_acpw==1


                if isfield(ds.GN_outputs, 'posterior_cov')==true

                    e1 = errorbar(ds.GN_outputs.retrieval(2,end), ds.GN_outputs.retrieval(1,end), sqrt(ds.GN_outputs.posterior_cov(1,1))/2,...
                        sqrt(ds.GN_outputs.posterior_cov(1,1))/2, sqrt(ds.GN_outputs.posterior_cov(2,2))/2,...
                        sqrt(ds.GN_outputs.posterior_cov(2,2))/2, 'Marker','.', 'MarkerSize', circ_size,...
                        'MarkerFaceColor', C(idx_tauC,:), 'MarkerEdgeColor', C(idx_tauC,:), 'Linewidth', ln_wdth,...
                        'Color', C(idx_tauC,:));

                elseif isfield(ds.GN_outputs, 'posterior_cov_lin')==true

                    e1 = errorbar(ds.GN_outputs.retrieval(2,end), ds.GN_outputs.retrieval(1,end), sqrt(ds.GN_outputs.posterior_cov_lin(1,1))/2,...
                        sqrt(ds.GN_outputs.posterior_cov_lin(1,1))/2, sqrt(ds.GN_outputs.posterior_cov_lin(2,2))/2,...
                        sqrt(ds.GN_outputs.posterior_cov_lin(2,2))/2,'Marker','.', 'MarkerSize', circ_size,...
                        'MarkerFaceColor', C(idx_tauC,:), 'MarkerEdgeColor', C(idx_tauC,:), 'Linewidth', ln_wdth,...
                        'Color', C(idx_tauC,:));

                end

            elseif idx_acpw==4


                if isfield(ds.GN_outputs, 'posterior_cov')==true

                    e1 = errorbar(ds.GN_outputs.retrieval(2,end), ds.GN_outputs.retrieval(1,end), sqrt(ds.GN_outputs.posterior_cov(1,1))/2,...
                        sqrt(ds.GN_outputs.posterior_cov(1,1))/2, sqrt(ds.GN_outputs.posterior_cov(2,2))/2,...
                        sqrt(ds.GN_outputs.posterior_cov(2,2))/2, 'Marker','square', 'MarkerSize', sq_dmnd_size,...
                        'MarkerFaceColor', C(idx_tauC,:), 'MarkerEdgeColor', C(idx_tauC,:), 'Linewidth', ln_wdth,...
                        'Color', C(idx_tauC,:));

                elseif isfield(ds.GN_outputs, 'posterior_cov_lin')==true

                    e1 = errorbar(ds.GN_outputs.retrieval(2,end), ds.GN_outputs.retrieval(1,end), sqrt(ds.GN_outputs.posterior_cov_lin(1,1))/2,...
                        sqrt(ds.GN_outputs.posterior_cov_lin(1,1))/2, sqrt(ds.GN_outputs.posterior_cov_lin(2,2))/2,...
                        sqrt(ds.GN_outputs.posterior_cov_lin(2,2))/2,'Marker','square', 'MarkerSize', sq_dmnd_size,...
                        'MarkerFaceColor', C(idx_tauC,:), 'MarkerEdgeColor', C(idx_tauC,:), 'Linewidth', ln_wdth,...
                        'Color', C(idx_tauC,:));

                end


            elseif idx_acpw==6


                if isfield(ds.GN_outputs, 'posterior_cov')==true

                    e1 = errorbar(ds.GN_outputs.retrieval(2,end), ds.GN_outputs.retrieval(1,end), sqrt(ds.GN_outputs.posterior_cov(1,1))/2,...
                        sqrt(ds.GN_outputs.posterior_cov(1,1))/2, sqrt(ds.GN_outputs.posterior_cov(2,2))/2,...
                        sqrt(ds.GN_outputs.posterior_cov(2,2))/2, 'Marker','diamond', 'MarkerSize', sq_dmnd_size,...
                        'MarkerFaceColor', C(idx_tauC,:), 'MarkerEdgeColor', C(idx_tauC,:), 'Linewidth', ln_wdth,...
                        'Color', C(idx_tauC,:));

                elseif isfield(ds.GN_outputs, 'posterior_cov_lin')==true

                    e1 = errorbar(ds.GN_outputs.retrieval(2,end), ds.GN_outputs.retrieval(1,end), sqrt(ds.GN_outputs.posterior_cov_lin(1,1))/2,...
                        sqrt(ds.GN_outputs.posterior_cov_lin(1,1))/2, sqrt(ds.GN_outputs.posterior_cov_lin(2,2))/2,...
                        sqrt(ds.GN_outputs.posterior_cov_lin(2,2))/2,'Marker','diamond', 'MarkerSize', sq_dmnd_size,...
                        'MarkerFaceColor', C(idx_tauC,:), 'MarkerEdgeColor', C(idx_tauC,:), 'Linewidth', ln_wdth,...
                        'Color', C(idx_tauC,:));

                end


            end





        end


        %         % Create a Legend with only the two black curves
        % legend(lgnd_str, 'Interpreter','latex', 'Location','northwest', 'FontSize', 20,...
        %     'Color', 'white', 'TextColor', 'k')

        grid on; grid minor
        ylabel('$r_{top}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',30)
        xlabel('$r_{bot}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',30)

        set(gcf,'Position',[0 0 1250 950])

        % xlim([3.5, 6])
        % ylim([8, 9.8])




    end



end




%% Show 6 examples of the full drolet profile compared with the assumed ACPW retrieval





clear variables


% Determine which computer you're using
which_computer = whatComputer();

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------

    % ---- Define where the retrievals are stored ---
    folder_paths.drive.fullRetrieval = ['/Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
        'paper2_variableSweep/newRetrieval_logTransform_newCov/'];

    folder_paths.drive.noACPW_10 = ['/Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
        'paper2_variableSweep/newRetrieval_logTransform_newCov_noACPW_10/'];

    folder_paths.drive.noACPW_15 = ['/Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
        'paper2_variableSweep/newRetrieval_logTransform_newCov_noACPW_15/'];

    folder_paths.drive.noACPW_20 = ['/Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
        'paper2_variableSweep/newRetrieval_logTransform_newCov_noACPW_20/'];

    folder_paths.drive.noACPW_25 = ['/Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
        'paper2_variableSweep/newRetrieval_logTransform_newCov_noACPW_25/'];


elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    % ---- Define where the retrievals are stored ---
    folder_paths.drive.fullRetrieval = ['/Users/andrewbuggee/MATLAB-Drive/HySICS/Droplet_profile_retrievals/paper2_variableSweep/newRetrieval_logTransform_newCov/'];

    folder_paths.drive.noACPW_10 = ['/Users/andrewbuggee/MATLAB-Drive/HySICS/Droplet_profile_retrievals/paper2_variableSweep/newRetrieval_logTransform_newCov_noACPW_10'];

    folder_paths.drive.noACPW_15 = ['/Users/andrewbuggee/MATLAB-Drive/HySICS/Droplet_profile_retrievals/paper2_variableSweep/newRetrieval_logTransform_newCov_noACPW_15/'];

    folder_paths.drive.noACPW_20 = ['/Users/andrewbuggee/MATLAB-Drive/HySICS/Droplet_profile_retrievals/paper2_variableSweep/newRetrieval_logTransform_newCov_noACPW_20/'];

    folder_paths.drive.noACPW_25 = ['/Users/andrewbuggee/MATLAB-Drive/HySICS/Droplet_profile_retrievals/paper2_variableSweep/newRetrieval_logTransform_newCov_noACPW_25/'];



elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------

end


% First step through all files in the directory and remove invisble files
% or directories

% Grab filenames in drive
filenames_fullRetrieval = dir(folder_paths.drive.fullRetrieval);
idx_2delete = [];
for nn = 1:length(filenames_fullRetrieval)

    if strcmp(filenames_fullRetrieval(nn).name(1), 'd')~=true

        idx_2delete = [idx_2delete, nn];

    end

end

% delete rows that don't have retrieval filenames
filenames_fullRetrieval(idx_2delete) = [];



% % Grab filenames in drive
% filenames_noACPW_10 = dir(folder_paths.drive.noACPW_10);
% idx_2delete = [];
% for nn = 1:length(filenames_noACPW_10)
%
%     if strcmp(filenames_noACPW_10(nn).name(1), 'd')~=true
%
%         idx_2delete = [idx_2delete, nn];
%
%     end
%
% end
%
% % delete rows that don't have retrieval filenames
% filenames_noACPW_10(idx_2delete) = [];




% define the 4 variables and the values that were simulated
r_top = [7.3418, 8.4296, 9.2516, 10.8153, 11.6297];
r_bot = [4.6147, 5.3192, 6.8372, 7.0135];
tau_c = [6.1312, 9.9317, 12.4356, 15.3172, 18.2468, 21.3466, 24.1697];
tcpw = [5.3741, 8.123, 10.6422, 13.6234, 15.8543, 18.9824];
acpw = [3.1097, 4.7003, 6.1580, 7.8831, 9.1740, 10.9840];

% Find all files within filenames_fullRetrieval with the same r_top and
% r_bot
% idx_profiles_2_plot = [200, 463, 744, 598, 81, 315];
% idx_profiles_2_plot = [666   806   551    30   714   785];
idx_profiles_2_plot = [666   806   551    30   785];
% idx_profiles_2_plot = randi([1, length(filenames_fullRetrieval)], 1, 6);

% Define the linewidth
ln_wdth = 2;

% define marker width for the circle
circ_size = 30;
% define size for square and diamond
sq_dmnd_size = 12;

% define the colors for each curve plotted
C = mySavedColors(61:(61+5 +1), 'fixed');

% lgnd_str = cell(1, length(tau_c) + length(tcpw) + 2);
lgnd_str = cell(1, 6);

figure;

for nn = 1:length(idx_profiles_2_plot)


    subplot(3,2,nn)



    % Load the data from the file
    ds = load(fullfile(filenames_fullRetrieval(idx_profiles_2_plot(nn)).folder,...
        filenames_fullRetrieval(idx_profiles_2_plot(nn)).name));



    % first, plot the simulated profile values as two lines
    xline(ds.GN_inputs.measurement.r_bot, ':', ['True $r_{bot}$'],...
        'Fontsize',20, 'Interpreter','latex','LineWidth',2,'Color', [147/255, 150/255, 151/255], 'LabelHorizontalAlignment','left',...
        'LabelVerticalAlignment','bottom');

    hold on

    yline(ds.GN_inputs.measurement.r_top, ':', ['True $r_{top}$'],...
        'Fontsize',20, 'Interpreter','latex','LineWidth',2,'Color', [147/255, 150/255, 151/255], 'LabelHorizontalAlignment','right',...
        'LabelVerticalAlignment','top');
    hold on


    % what was the assumed above cloud column water vapor path?
    simulated_ACPW = ds.GN_inputs.measurement.actpw; % kg/m^2
    simulated_tauC = ds.GN_inputs.measurement.tau_c; % 
    simulated_rTop = ds.GN_inputs.measurement.r_top; % 
    simulated_rBot = ds.GN_inputs.measurement.r_bot; % 

    % title(['Simulated profile - $acpw$ = ',num2str(round(simulated_ACPW, 2)), ' $mm$',...
    %     ', $\tau_{c} = $', num2str(simulated_tauC)],...
    %     'Fontsize', 15, 'Interpreter', 'latex');

    title(['$r_{top} = $', num2str(round(simulated_rTop,1)),', $r_{bot} = $', num2str(round(simulated_rBot,1)),...
        ', $\tau_{c} = $', num2str(round(simulated_tauC,1)), ', $acpw$ = ',num2str(round(simulated_ACPW, 1)), ' $mm$'],...
        'Fontsize', 15, 'Interpreter', 'latex');

    % Skip the first two legend entries
    lgnd_str{1} = '';
    lgnd_str{2} = '';





    % ----------------------------------------------------------------
    % --------------- First, plot the whole retrieval ----------------
    % ----------------------------------------------------------------
    hold on


    % if isfield(ds, 'GN_outputs')==true && isfield(ds.GN_outputs, 'posterior_cov')==true
    % 
    %     e1 = errorbar(ds.GN_outputs.retrieval(2,end), ds.GN_outputs.retrieval(1,end), sqrt(ds.GN_outputs.posterior_cov(1,1))/2,...
    %         sqrt(ds.GN_outputs.posterior_cov(1,1))/2, sqrt(ds.GN_outputs.posterior_cov(2,2))/2,...
    %         sqrt(ds.GN_outputs.posterior_cov(2,2))/2, 'MarkerFaceColor', C(1,:),...
    %         'MarkerEdgeColor', C(1,:), 'Linewidth', ln_wdth, 'Marker', '.', 'MarkerSize', circ_size,...
    %         'Color', C(1,:));
    % 
    % elseif isfield(ds, 'GN_outputs')==true && isfield(ds.GN_outputs, 'posterior_cov_lin')==true
    % 
    %     e1 = errorbar(ds.GN_outputs.retrieval(2,end), ds.GN_outputs.retrieval(1,end), sqrt(ds.GN_outputs.posterior_cov_lin(1,1))/2,...
    %         sqrt(ds.GN_outputs.posterior_cov_lin(1,1))/2, sqrt(ds.GN_outputs.posterior_cov_lin(2,2))/2,...
    %         sqrt(ds.GN_outputs.posterior_cov_lin(2,2))/2, 'MarkerFaceColor', C(1,:),...
    %         'MarkerEdgeColor', C(1,:), 'Linewidth', ln_wdth, 'Marker', '.', 'MarkerSize', circ_size,...
    %         'Color', C(1,:));
    % 
    % end

    hold on


    % determine which TCPW is being plotted
    [~,idx_acpw] = min(abs(acpw - ds.GN_inputs.measurement.actpw));
    tcpw_file = tcpw(idx_acpw);
    % Define the matching file name
    if idx_acpw == 2

        filePattern = sprintf(['dropletRetrieval_noACPW_HySICS_35bands_0.3%%_uncert_rTop_%.4f_rBot_%.4f',...
            '_tauC_%.4f_tcwv_%.3f_*.mat'], ds.GN_inputs.measurement.r_top, ds.GN_inputs.measurement.r_bot,...
            ds.GN_inputs.measurement.tau_c, tcpw_file);

    else

        filePattern = sprintf(['dropletRetrieval_noACPW_HySICS_35bands_0.3%%_uncert_rTop_%.4f_rBot_%.4f',...
            '_tauC_%.4f_tcwv_%.4f_*.mat'], ds.GN_inputs.measurement.r_top, ds.GN_inputs.measurement.r_bot,...
            ds.GN_inputs.measurement.tau_c, tcpw_file);

    end



    % ----------------------------------------------------------------
    % Next, find and plot the correspond file for assumed ACPW of 10mm
    % ----------------------------------------------------------------

    matchingFiles_10 = dir(fullfile(folder_paths.drive.noACPW_10, filePattern));

    if isempty(matchingFiles_10)==false

        % Load the data from the file
        ds_10 = load(fullfile(matchingFiles_10(1).folder, matchingFiles_10(1).name));

        hold on


        if isfield(ds_10, 'GN_outputs')==true && isfield(ds_10.GN_outputs, 'posterior_cov')==true

            e2 = errorbar(ds_10.GN_outputs.retrieval(2,end), ds_10.GN_outputs.retrieval(1,end), sqrt(ds_10.GN_outputs.posterior_cov(1,1))/2,...
                sqrt(ds_10.GN_outputs.posterior_cov(1,1))/2, sqrt(ds_10.GN_outputs.posterior_cov(2,2))/2,...
                sqrt(ds_10.GN_outputs.posterior_cov(2,2))/2, 'MarkerFaceColor', C(2,:),...
                'MarkerEdgeColor', C(2,:), 'Linewidth', ln_wdth, 'Marker', '.', 'MarkerSize', circ_size,...
                'Color', C(2,:));

        elseif isfield(ds_10, 'GN_outputs')==true && isfield(ds_10.GN_outputs, 'posterior_cov_lin')==true

            e2 = errorbar(ds_10.GN_outputs.retrieval(2,end), ds_10.GN_outputs.retrieval(1,end), sqrt(ds_10.GN_outputs.posterior_cov_lin(1,1))/2,...
                sqrt(ds_10.GN_outputs.posterior_cov_lin(1,1))/2, sqrt(ds_10.GN_outputs.posterior_cov_lin(2,2))/2,...
                sqrt(ds_10.GN_outputs.posterior_cov_lin(2,2))/2, 'MarkerFaceColor', C(2,:),...
                'MarkerEdgeColor', C(2,:), 'Linewidth', ln_wdth, 'Marker', '.', 'MarkerSize', circ_size,...
                'Color', C(2,:));

        end

        e2.Bar.LineStyle = 'dotted';


    end

    lgnd_str{3} = 'assumed $acpw = 10 mm$';



    % ----------------------------------------------------------------
    % Next, find and plot the correspond file for assumed ACPW of 15mm
    % ----------------------------------------------------------------

    matchingFiles_15 = dir(fullfile(folder_paths.drive.noACPW_15, filePattern));

    if isempty(matchingFiles_15)==false

        % Load the data from the file
        ds_15 = load(fullfile(matchingFiles_15(1).folder, matchingFiles_15(1).name));

        hold on


        if isfield(ds_15, 'GN_outputs')==true && isfield(ds_15.GN_outputs, 'posterior_cov')==true

            e3 = errorbar(ds_15.GN_outputs.retrieval(2,end), ds_15.GN_outputs.retrieval(1,end), sqrt(ds_15.GN_outputs.posterior_cov(1,1))/2,...
                sqrt(ds_15.GN_outputs.posterior_cov(1,1))/2, sqrt(ds_15.GN_outputs.posterior_cov(2,2))/2,...
                sqrt(ds_15.GN_outputs.posterior_cov(2,2))/2, 'MarkerFaceColor', C(3,:),...
                'MarkerEdgeColor', C(3,:), 'Linewidth', ln_wdth, 'Marker', '.', 'MarkerSize', circ_size,...
                'Color', C(3,:));

        elseif isfield(ds_15, 'GN_outputs')==true && isfield(ds_15.GN_outputs, 'posterior_cov_lin')==true

            e3 = errorbar(ds_15.GN_outputs.retrieval(2,end), ds_15.GN_outputs.retrieval(1,end), sqrt(ds_15.GN_outputs.posterior_cov_lin(1,1))/2,...
                sqrt(ds_15.GN_outputs.posterior_cov_lin(1,1))/2, sqrt(ds_15.GN_outputs.posterior_cov_lin(2,2))/2,...
                sqrt(ds_15.GN_outputs.posterior_cov_lin(2,2))/2, 'MarkerFaceColor', C(3,:),...
                'MarkerEdgeColor', C(3,:), 'Linewidth', ln_wdth, 'Marker', '.', 'MarkerSize', circ_size,...
                'Color', C(3,:));

        end

        e3.Bar.LineStyle = 'dotted';


    end

    lgnd_str{4} = 'assumed $acpw = 15 mm$';



    % ----------------------------------------------------------------
    % Next, find and plot the correspond file for assumed ACPW of 20mm
    % ----------------------------------------------------------------

    matchingFiles_20 = dir(fullfile(folder_paths.drive.noACPW_20, filePattern));

    if isempty(matchingFiles_20)==false

        % Load the data from the file
        ds_20 = load(fullfile(matchingFiles_20(1).folder, matchingFiles_20(1).name));

        hold on


        if isfield(ds_20, 'GN_outputs')==true && isfield(ds_20.GN_outputs, 'posterior_cov')==true

            e4 = errorbar(ds_20.GN_outputs.retrieval(2,end), ds_20.GN_outputs.retrieval(1,end), sqrt(ds_20.GN_outputs.posterior_cov(1,1))/2,...
                sqrt(ds_20.GN_outputs.posterior_cov(1,1))/2, sqrt(ds_20.GN_outputs.posterior_cov(2,2))/2,...
                sqrt(ds_20.GN_outputs.posterior_cov(2,2))/2, 'MarkerFaceColor', C(4,:),...
                'MarkerEdgeColor', C(4,:), 'Linewidth', ln_wdth, 'Marker', '.', 'MarkerSize', circ_size,...
                'Color', C(4,:));

        elseif isfield(ds_20, 'GN_outputs')==true && isfield(ds_20.GN_outputs, 'posterior_cov_lin')==true

            e4 = errorbar(ds_20.GN_outputs.retrieval(2,end), ds_20.GN_outputs.retrieval(1,end), sqrt(ds_20.GN_outputs.posterior_cov_lin(1,1))/2,...
                sqrt(ds_20.GN_outputs.posterior_cov_lin(1,1))/2, sqrt(ds_20.GN_outputs.posterior_cov_lin(2,2))/2,...
                sqrt(ds_20.GN_outputs.posterior_cov_lin(2,2))/2, 'MarkerFaceColor', C(4,:),...
                'MarkerEdgeColor', C(4,:), 'Linewidth', ln_wdth, 'Marker', '.', 'MarkerSize', circ_size,...
                'Color', C(4,:));

        end

        e4.Bar.LineStyle = 'dotted';


    end

    lgnd_str{5} = 'assumed $acpw = 20 mm$';



    % -------------------------------------------------------------------
    % Lastly, find and plot the correspond file for assumed ACPW of 25mm
    % -------------------------------------------------------------------

    matchingFiles_25 = dir(fullfile(folder_paths.drive.noACPW_25, filePattern));

    if isempty(matchingFiles_25)==false

        % Load the data from the file
        ds_25 = load(fullfile(matchingFiles_25(1).folder, matchingFiles_25(1).name));

        hold on


        if isfield(ds_25, 'GN_outputs')==true && isfield(ds_25.GN_outputs, 'posterior_cov')==true

            e5 = errorbar(ds_25.GN_outputs.retrieval(2,end), ds_25.GN_outputs.retrieval(1,end), sqrt(ds_25.GN_outputs.posterior_cov(1,1))/2,...
                sqrt(ds_25.GN_outputs.posterior_cov(1,1))/2, sqrt(ds_25.GN_outputs.posterior_cov(2,2))/2,...
                sqrt(ds_25.GN_outputs.posterior_cov(2,2))/2, 'MarkerFaceColor', C(5,:),...
                'MarkerEdgeColor', C(5,:), 'Linewidth', ln_wdth, 'Marker', '.', 'MarkerSize', circ_size,...
                'Color', C(5,:));

        elseif isfield(ds_25, 'GN_outputs')==true && isfield(ds_25.GN_outputs, 'posterior_cov_lin')==true 

            e5 = errorbar(ds_25.GN_outputs.retrieval(2,end), ds_25.GN_outputs.retrieval(1,end), sqrt(ds_25.GN_outputs.posterior_cov_lin(1,1))/2,...
                sqrt(ds_25.GN_outputs.posterior_cov_lin(1,1))/2, sqrt(ds_25.GN_outputs.posterior_cov_lin(2,2))/2,...
                sqrt(ds_25.GN_outputs.posterior_cov_lin(2,2))/2, 'MarkerFaceColor', C(5,:),...
                'MarkerEdgeColor', C(5,:), 'Linewidth', ln_wdth, 'Marker', '.', 'MarkerSize', circ_size,...
                'Color', C(5,:));

        end

        e5.Bar.LineStyle = 'dotted';


    end

    lgnd_str{6} = 'assumed $acpw = 25 mm$';



    %         % Create a Legend with only the two black curves
    % legend(lgnd_str, 'Interpreter','latex', 'Location','northwest', 'FontSize', 20,...
    %     'Color', 'white', 'TextColor', 'k')

    grid on; grid minor
    ylabel('$r_{top}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',15)
    xlabel('$r_{bot}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',15)

    if nn==1

        xlim([ds.GN_inputs.measurement.r_bot - 2.75, ds.GN_inputs.measurement.r_bot + 2.75])
        ylim([ds.GN_inputs.measurement.r_top - 0.75, ds.GN_inputs.measurement.r_top + 0.75])

    else

        xlim([ds.GN_inputs.measurement.r_bot - 1.8, ds.GN_inputs.measurement.r_bot + 1.8])
        ylim([ds.GN_inputs.measurement.r_top - 0.75, ds.GN_inputs.measurement.r_top + 0.75])

    end


end



    set(gcf,'Position',[0 0 1250 950])



% plot legend
legend(lgnd_str, 'Interpreter','latex', 'Location','best', 'FontSize', 20,...
    'Color', 'white', 'TextColor', 'k')





% create a figure with subplots by stepping through optical depth
% and plotting the full droplet profile retrieval along with the
% assumed ACPW retreivals



%% Plot droplet profile retrieval results using the new algorithm that works in log space

% Using the new measurements and:
% (1) Transformed the algorithm to operate in log space
% (2) the new Jacobian calcualtion that should be better at ensuring
% the change or some small perturbation of each state variable is linear
% (3) Using a smaller maximum scalar value to compute the change in the
% state vector, New value is 3 instead of 20

% Plot full retrieval only

% *** Using HySICS data ***

clear variables


% Determine which computer you're using
which_computer = whatComputer();

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------

    % ---- Define where the retrievals are stored ---
    folder_paths.drive = ['/Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
        'paper2_variableSweep/rTop_10/vza_4_vaz_257_sza_31_saz_96_subset_logState_1/'];


elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    % ---- Define where the retrievals are stored ---
    % folder_paths.drive = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
    %     'HySICS/Droplet_profile_retrievals/paper2_variableSweep/test_logSpace_newCov/'];

    folder_paths.drive = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'HySICS/Droplet_profile_retrievals/paper2_variableSweep/test_logSpace_newCov_2/'];


elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------

end



% Grab filenames in drive
filenames = dir(folder_paths.drive);
idx_2delete = [];
for nn = 1:length(filenames)

    if strcmp(filenames(nn).name(1), 'd')~=true

        idx_2delete = [idx_2delete, nn];

    end

end

% delete rows that don't have retrieval filenames
filenames(idx_2delete) = [];




% Step through each file

% define the colors for each curve plotted
C = mySavedColors(1:(1+length(filenames) +1), 'fixed');

lgnd_str = cell(1, length(filenames) + 2);


figure;


for nn = 1:length(filenames)


    % Load a data set
    ds = load([filenames(nn).folder, '/', filenames(nn).name]);

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




    if strcmp(filenames(nn).name(1:23), 'dropletRetrieval_noACPW')==true

        % plot the retrieved droplet profile using 35 bands and assuming an
        % ACPW

        hold on


        % Plot the retrieval uncertainty of the radius at cloud top and
        % bottom
        if isfield(ds.GN_outputs, 'posterior_cov')==true

            e1 = errorbar(ds.GN_outputs.retrieval(2,end), ds.GN_outputs.retrieval(1,end), sqrt(ds.GN_outputs.posterior_cov(1,1))/2,...
                sqrt(ds.GN_outputs.posterior_cov(1,1))/2, sqrt(ds.GN_outputs.posterior_cov(2,2))/2,...
                sqrt(ds.GN_outputs.posterior_cov(2,2))/2, 'MarkerFaceColor', C(idx_tauC,:),...
                'MarkerEdgeColor', C(idx_tauC,:), 'Linewidth', 4, 'Marker', '.', 'MarkerSize', 40,...
                'Color', C(idx_tauC,:));

        elseif isfield(ds.GN_outputs, 'posterior_cov_lin')==true

            e1 = errorbar(ds.GN_outputs.retrieval(2,end), ds.GN_outputs.retrieval(1,end), sqrt(ds.GN_outputs.posterior_cov_lin(1,1))/2,...
                sqrt(ds.GN_outputs.posterior_cov_lin(1,1))/2, sqrt(ds.GN_outputs.posterior_cov_lin(2,2))/2,...
                sqrt(ds.GN_outputs.posterior_cov_lin(2,2))/2, 'MarkerFaceColor', C(idx_tauC,:),...
                'MarkerEdgeColor', C(idx_tauC,:), 'Linewidth', 4, 'Marker', '.', 'MarkerSize', 40,...
                'Color', C(idx_tauC,:));

        end


        e1.Bar.LineStyle = 'dotted';

        hold on


        % what was the assumed above cloud column water vapor path?
        assumed_CWV = aboveCloud_CWV_simulated_hysics_spectra(ds.GN_inputs); % kg/m^2

        lgnd_str{nn+2} = [num2str(numel(ds.GN_inputs.bands2run)), ' bands - assumed $acpw$ = ',...
            num2str(round(assumed_CWV,2)), ' $mm$ - ', num2str(ds.GN_inputs.measurement.uncertainty(1)*100),...
            '\% meas uncert'];



    elseif strcmp(filenames(nn).name(1:23), 'dropletRetrieval_HySICS')==true

        % give a different marker type for the retrieval using 66 bands
        % that also retrieved above cloud column water vapor

        hold on


        if isfield(ds.GN_outputs, 'posterior_cov')==true

            e1 = errorbar(ds.GN_outputs.retrieval(2,end), ds.GN_outputs.retrieval(1,end), sqrt(ds.GN_outputs.posterior_cov(1,1))/2,...
                sqrt(ds.GN_outputs.posterior_cov(1,1))/2, sqrt(ds.GN_outputs.posterior_cov(2,2))/2,...
                sqrt(ds.GN_outputs.posterior_cov(2,2))/2, 'MarkerFaceColor', C(idx_tauC,:),...
                'MarkerEdgeColor', C(idx_tauC,:), 'Linewidth', 4, 'Marker', '.', 'MarkerSize', 40,...
                'Color', C(idx_tauC,:));

        elseif isfield(ds.GN_outputs, 'posterior_cov_lin')==true

            e1 = errorbar(ds.GN_outputs.retrieval(2,end), ds.GN_outputs.retrieval(1,end), sqrt(ds.GN_outputs.posterior_cov_lin(1,1))/2,...
                sqrt(ds.GN_outputs.posterior_cov_lin(1,1))/2, sqrt(ds.GN_outputs.posterior_cov_lin(2,2))/2,...
                sqrt(ds.GN_outputs.posterior_cov_lin(2,2))/2, 'MarkerFaceColor', C(idx_tauC,:),...
                'MarkerEdgeColor', C(idx_tauC,:), 'Linewidth', 4, 'Marker', '.', 'MarkerSize', 40,...
                'Color', C(idx_tauC,:));

        end

        hold on


        % create the string for the retrieval using CWV

        % what was the retrieved above cloud column water vapor path above
        % cloud?
        retrieved_CWV = ds.GN_outputs.retrieval(end, end);        % kg/m^2 (mm)

        lgnd_str{nn+2} = [num2str(numel(ds.GN_inputs.bands2run)), ' bands - retrieved $acpw$ = ',...
            num2str(round(retrieved_CWV, 2)), ' $mm$ - ', num2str(ds.GN_inputs.measurement.uncertainty(1)*100),...
            '\% meas uncert'];


    else

        error([newline, 'There are no files to read.', newline])


    end





end




% Create a Legend with only the two black curves
legend(lgnd_str, 'Interpreter','latex', 'Location','northwest', 'FontSize', 20,...
    'Color', 'white', 'TextColor', 'k')

grid on; grid minor
ylabel('$r_{top}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',30)
xlabel('$r_{bot}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',30)

set(gcf,'Position',[0 0 950 750])

xlim([3.5, 6])
ylim([8, 9.8])
