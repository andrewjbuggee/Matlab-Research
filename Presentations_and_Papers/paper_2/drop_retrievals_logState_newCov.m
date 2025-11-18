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
    error('Where is the folder?')


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
        e1 = errorbar(ds.GN_outputs.retrieval(2,end), ds.GN_outputs.retrieval(1,end), sqrt(ds.GN_outputs.posterior_cov(1,1))/2,...
            sqrt(ds.GN_outputs.posterior_cov(1,1))/2, sqrt(ds.GN_outputs.posterior_cov(2,2))/2,...
            sqrt(ds.GN_outputs.posterior_cov(2,2))/2, 'MarkerFaceColor', C(nn+1,:),...
            'MarkerEdgeColor', C(nn+1,:), 'Linewidth', 4, 'Marker', '.', 'MarkerSize', 40,...
            'Color', C(nn+1,:));

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


        errorbar(ds.GN_outputs.retrieval(2,end), ds.GN_outputs.retrieval(1,end), sqrt(ds.GN_outputs.posterior_cov(1,1))/2,...
            sqrt(ds.GN_outputs.posterior_cov(1,1))/2, sqrt(ds.GN_outputs.posterior_cov(2,2))/2,...
            sqrt(ds.GN_outputs.posterior_cov(2,2))/2, 'MarkerFaceColor', C(nn+1,:),...
            'MarkerEdgeColor', C(nn+1,:), 'Linewidth', 4, 'Marker', '.', 'MarkerSize', 40,...
            'Color', C(nn+1,:))

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
legend(lgnd_str, 'Interpreter','latex', 'Location','northwest', 'FontSize', 10,...
    'Color', 'white', 'TextColor', 'k')

grid on; grid minor
ylabel('$r_{top}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',30)
xlabel('$r_{bot}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',30)

set(gcf,'Position',[0 0 950 750])

xlim([3.5, 6])
ylim([8, 9.8])





%% Plot droplet profile retrieval results using the new algorithm that works in log space
% And show the retrieval results using the new covariance matrix

% Using the new measurements and:
% (1) Transformed the algorithm to operate in log space
% (2) Used in-situ CDP and radiosonde measurements to compute new
% prior covariance matrix
% (3) the new Jacobian calcualtion that should be better at ensuring
% the change or some small perturbation of each state variable is linear
% (4) Using a smaller maximum scalar value to compute the change in the
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

    % % ---- Define where the log state retrievals are stored ---
    % folder_paths.drive_1 = ['/Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
    %     'paper2_variableSweep/rTop_10/vza_4_vaz_257_sza_31_saz_96_subset_logState_1/'];
    % 
    % % ---- Define where the log state retrievals with new prior covariance are stored ---
    % folder_paths.drive_2 = ['/Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
    %     'paper2_variableSweep/rTop_10/vza_4_vaz_257_sza_31_saz_96_subset_logState_newCov_1/'];


        % ---- Define where the log state retrievals are stored ---
    folder_paths.drive_1 = ['/Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
        'paper2_variableSweep/rTop_10/vza_4_vaz_257_sza_31_saz_96_subset_logState_2/'];

    % ---- Define where the log state retrievals with new prior covariance are stored ---
    folder_paths.drive_2 = ['/Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
        'paper2_variableSweep/rTop_10/vza_4_vaz_257_sza_31_saz_96_subset_logState_newCov_2/'];


elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    % ---- Define where the retrievals are stored ---
    error('Where is the folder?')


elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------

end



% Grab filenames for the log-state retrieval
filenames_log = dir(folder_paths.drive_1);
idx_2delete = [];
for nn = 1:length(filenames_log)

    if strcmp(filenames_log(nn).name(1), 'd')~=true

        idx_2delete = [idx_2delete, nn];

    end

end

% delete rows that don't have retrieval filenames
filenames_log(idx_2delete) = [];



% Grab filenames for the log-state and new covariance retreival
filenames_log_newCov = dir(folder_paths.drive_2);
idx_2delete = [];
for nn = 1:length(filenames_log_newCov)

    if strcmp(filenames_log_newCov(nn).name(1), 'd')~=true

        idx_2delete = [idx_2delete, nn];

    end

end

% delete rows that don't have retrieval filenames
filenames_log_newCov(idx_2delete) = [];





% Step through each file

% define the colors for each curve plotted
C = mySavedColors(1:(1+length(filenames_log) +1), 'fixed');

lgnd_str = cell(1, length(filenames_log) + length(filenames_log_newCov) + 2);


figure;


for nn = 1:length(filenames_log)


    % Load a data set
    ds = load([filenames_log(nn).folder, '/', filenames_log(nn).name]);

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


        % % what was the assumed above cloud column water vapor path?
        % simulated_CWV = ds.GN_inputs.measurement.actpw; % kg/m^2

        title('Comparison between log-state and log-state + new prior cov',...
            'Fontsize', 25, 'Interpreter', 'latex');

        % Skip the first two legend entries
        lgnd_str{1} = '';
        lgnd_str{2} = '';

    end



    hold on


    errorbar(ds.GN_outputs.retrieval(2,end), ds.GN_outputs.retrieval(1,end), sqrt(ds.GN_outputs.posterior_cov(1,1))/2,...
        sqrt(ds.GN_outputs.posterior_cov(1,1))/2, sqrt(ds.GN_outputs.posterior_cov(2,2))/2,...
        sqrt(ds.GN_outputs.posterior_cov(2,2))/2, 'MarkerFaceColor', C(nn+1,:),...
        'MarkerEdgeColor', C(nn+1,:), 'Linewidth', 4, 'Marker', '.', 'MarkerSize', 40,...
        'Color', C(nn+1,:))

    hold on


    % create the string for the retrieval using CWV

    % what was the retrieved above cloud column water vapor path above
    % cloud?
    retrieved_CWV = ds.GN_outputs.retrieval(end, end);        % kg/m^2 (mm)

    lgnd_str{nn+2} = [num2str(numel(ds.GN_inputs.bands2run)), ' bands - log-state - retrieved $acpw$ = ',...
        num2str(round(retrieved_CWV, 2)), ' $mm$ - ', num2str(ds.GN_inputs.measurement.uncertainty(1)*100),...
        '\% meas uncert'];








end









% In case the two folders don't have the same number of files...
% define the colors for each curve plotted
C = mySavedColors(1:(1+length(filenames_log_newCov) +1), 'fixed');


% Step through each file


for nn = 1:length(filenames_log_newCov)


    % Load a data set
    ds = load([filenames_log_newCov(nn).folder, '/', filenames_log_newCov(nn).name]);




    % plot the retrieved droplet profile using 35 bands and assuming an
    % ACPW

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


    % what was the assumed above cloud column water vapor path?
    assumed_CWV = aboveCloud_CWV_simulated_hysics_spectra(ds.GN_inputs); % kg/m^2

    lgnd_str{nn+2+length(filenames_log)} = [num2str(numel(ds.GN_inputs.bands2run)), ' bands - log-state w/ new Cov ',...
        '- retrieved $acpw$ = ',...
        num2str(round(retrieved_CWV, 2)), ' $mm$ - ', num2str(ds.GN_inputs.measurement.uncertainty(1)*100),...
        '\% meas uncert'];


end




% Create a Legend with only the two black curves
legend(lgnd_str, 'Interpreter','latex', 'Location','northwest', 'FontSize', 10,...
    'Color', 'white', 'TextColor', 'k')

grid on; grid minor
ylabel('$r_{top}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',30)
xlabel('$r_{bot}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',30)

set(gcf,'Position',[0 0 950 750])

xlim([3.5, 6])
ylim([8, 9.8])






%% Plot the retrieved above cloud precipitable water versus the true value for all measurements


clear variables




% Determine which computer you're using
which_computer = whatComputer();

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------


    % ---- Define where the log state retrievals are stored ---
    folder_paths.drive_1 = ['/Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
        'paper2_variableSweep/rTop_10/vza_4_vaz_257_sza_31_saz_96_subset_logState_2/'];

    % ---- Define where the log state retrievals with new prior covariance are stored ---
    folder_paths.drive_2 = ['/Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
        'paper2_variableSweep/rTop_10/vza_4_vaz_257_sza_31_saz_96_subset_logState_newCov_2/'];


elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    % ---- Define where the retrievals are stored ---
    error('Where is the folder?')


elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------

end




% Grab filenames for the log-state retrieval
filenames = dir(folder_paths.drive_2);
idx_2delete = [];
for nn = 1:length(filenames)

    if strcmp(filenames(nn).name(1), 'd')~=true

        idx_2delete = [idx_2delete, nn];

    end

end
% delete rows that don't have retrieval filenames
filenames(idx_2delete) = [];



% what are the free parameters?
r_top = 9.2516;
r_bot = 5.3192;
tau_c = [6.1312, 9.9317, 12.4356, 15.3172];
tcpw = [8.123, 10.6422, 13.6234, 15.8543, 18.9824];


% length of each independent variable
num_rTop = length(r_top);
num_rBot = length(r_bot);
num_tauC = length(tau_c);
num_tcpw = length(tcpw);


% Step through each file

% For some simulated above cloud total column precipitable water amount,
% let's show the retrieval of r-top and r-bot for different optical depths.
% But somehow I want to show a comparison between the simulated values and
% the retrieved values for each independent variable.

% But that isn't how my data is arranged. I have a retrieval set where
% every profile has a radius at clod top of 10 microns. So I will have to
% step through the other variables for now.

% How do I look at the difference between the retrieved values and the true
% value for 4 different variables? Could I have panels for different values
% of the above cloud precipitable water?

% The y axis is the percent difference between the true r-top value and the
% retrieved

% The x axis is the percent difference between the true r-bot value and the
% retrieved

% Each curve represents a different optical depth

% define the colors for each curve plotted
C = mySavedColors(61:(61+num_tauC-1), 'fixed');

% lgnd_str = cell(1, length(filenames) + 2);


% one-to-one line
x = linspace(3,13,100);


figure;


for nn = 1:length(filenames)

    if nn==1

        % make a 1-to-1 line
        plot(x, x, 'k-', 'LineWidth',1)

        grid on; grid minor
        hold on

    end


    % Load a data set
    ds = load([filenames(nn).folder, '/', filenames(nn).name]);

    if isfield(ds, 'GN_outputs')

        % acpw_true = [acpw_true, ds.GN_inputs.measurement.actpw];

        [~, idx_tauC] = min(abs(ds.GN_inputs.measurement.tau_c - tau_c));

        plot(ds.GN_inputs.measurement.actpw, ds.GN_outputs.retrieval(4,end), '.', 'MarkerSize', 35,...
            'Color', C(idx_tauC, :))

        hold on

    end


end

% legend(["one to one", "$\tau_c = $"+string((round(tau_c, 2)))], 'Interpreter','latex', 'Location','northwest', 'FontSize', 20)
legend("one to one", 'Interpreter','latex', 'Location','northwest', 'FontSize', 20)

title('Simulated HySICS measurements', 'Interpreter','latex', 'FontSize',30)
ylabel('$acpw$ retrieved (mm)', 'Interpreter','latex', 'FontSize',30)
xlabel('$acpw$ true (mm)', 'Interpreter','latex', 'FontSize',30)

set(gcf,'Position',[0 0 750 750])

xlim([3,13])
ylim([3,13])







%% Plot the retrieved optical depth versus the true optical depth for all measurements



clear variables




% Determine which computer you're using
which_computer = whatComputer();

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------


    % ---- Define where the log state retrievals are stored ---
    folder_paths.drive_1 = ['/Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
        'paper2_variableSweep/rTop_10/vza_4_vaz_257_sza_31_saz_96_subset_logState_2/'];

    % ---- Define where the log state retrievals with new prior covariance are stored ---
    folder_paths.drive_2 = ['/Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
        'paper2_variableSweep/rTop_10/vza_4_vaz_257_sza_31_saz_96_subset_logState_newCov_2/'];


elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    % ---- Define where the retrievals are stored ---
    error('Where is the folder?')


elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------

end




% Grab filenames for the log-state retrieval
filenames = dir(folder_paths.drive_2);
idx_2delete = [];
for nn = 1:length(filenames)

    if strcmp(filenames(nn).name(1), 'd')~=true

        idx_2delete = [idx_2delete, nn];

    end

end
% delete rows that don't have retrieval filenames
filenames(idx_2delete) = [];


% what are the free parameters?
r_top = 9.2516;
r_bot = 5.3192;
tau_c = [6.1312, 9.9317, 12.4356, 15.3172];
tcpw = [8.123, 10.6422, 13.6234, 15.8543, 18.9824];

acpw_true = [4.7003, 6.158, 7.8831, 9.1740, 10.9840];


% length of each independent variable
num_rTop = length(r_top);
num_rBot = length(r_bot);
num_tauC = length(tau_c);
num_tcpw = length(tcpw);

num_files = num_rTop * num_rBot * num_tauC * num_tcpw;

% Step through each file

% For some simulated above cloud total column precipitable water amount,
% let's show the retrieval of r-top and r-bot for different optical depths.
% But somehow I want to show a comparison between the simulated values and
% the retrieved values for each independent variable.

% But that isn't how my data is arranged. I have a retrieval set where
% every profile has a radius at clod top of 10 microns. So I will have to
% step through the other variables for now.

% How do I look at the difference between the retrieved values and the true
% value for 4 different variables? Could I have panels for different values
% of the above cloud precipitable water?

% The y axis is the percent difference between the true r-top value and the
% retrieved

% The x axis is the percent difference between the true r-bot value and the
% retrieved

% Each curve represents a different optical depth

% define the colors for each curve plotted
C = mySavedColors(61:(61+num_tcpw-1), 'fixed');

% lgnd_str = cell(1, length(filenames) + 2);



figure;


for nn = 1:length(filenames)

    if nn==1

        % make a 1-to-1 line
        plot(tau_c, tau_c, 'k-', 'LineWidth',1)

        grid on; grid minor
        hold on

    end


    % Load a data set
    ds = load([filenames(nn).folder, '/', filenames(nn).name]);

    if isfield(ds, 'GN_outputs')

        % acpw_true = [acpw_true, ds.GN_inputs.measurement.actpw];

        [~, idx_actpw] = min(abs(ds.GN_inputs.measurement.actpw - acpw_true));

        plot(ds.GN_inputs.measurement.tau_c, ds.GN_outputs.retrieval(3,end), '.', 'MarkerSize', 20,...
            'Color', C(idx_actpw, :))

        hold on

    end


end

legend(["one to one", "$acpw = $"+string((round(acpw_true, 2)))], 'Interpreter','latex', 'Location','northwest', 'FontSize', 20)

title('$r_{top} = 10$ simulated HySICS measurements', 'Interpreter','latex', 'FontSize',30)
ylabel('$\tau_c$ retrieved', 'Interpreter','latex', 'FontSize',30)
xlabel('$\tau_c$ true', 'Interpreter','latex', 'FontSize',30)

set(gcf,'Position',[0 0 950 750])


