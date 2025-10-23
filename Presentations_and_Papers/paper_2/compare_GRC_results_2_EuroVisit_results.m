%% Plot the GRC results and overlay the three retrievals with measurement uncertainty of 0%, 0.3% and 1% using the EuroVisit algorithm

% the one file with retrieved column water vapor

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
    folder_paths.HySICS_retrievals = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'HySICS/Droplet_profile_retrievals/'];

    folder_paths.drive = '/Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/testGRC_results/';


elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    % ---- Define where the retrievals are stored ---
    folder_paths.HySICS_retrievals = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/HySICS/Droplet_profile_retrievals/'];


    folder_paths.drive = '/Users/andrewbuggee/MATLAB-Drive/HySICS/Droplet_profile_retrievals/testGRC_results/';


elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------

end


% define the mat files for each retreival to plot
filenames = {'dropletRetrieval_HySICS_35bands_withNoise_10mm-totalCWV_sim-ran-on-12-Jul-2025_rev1.mat',...
    'dropletRetrieval_HySICS_35bands_withNoise_15mm-totalCWV_sim-ran-on-12-Jul-2025_rev1.mat',...
    'dropletRetrieval_HySICS_35bands_withNoise_20mm-totalCWV_sim-ran-on-13-Jul-2025_rev1.mat',...
    'dropletRetrieval_HySICS_35bands_withNoise_25mm-totalCWV_sim-ran-on-12-Jul-2025_rev1.mat',...
    'dropletRetrieval_HySICS_66bands_withNoise_cwvRetrieval_sim-ran-on-12-Jul-2025_rev1.mat'};

% define the mat files for each retreival to plot
% filenames = {'dropletRetrieval_HySICS_35bands_withNoise_10mm-totalCWV_sim-ran-on-12-Jul-2025_rev1.mat',...
%     'dropletRetrieval_HySICS_35bands_withNoise_15mm-totalCWV_sim-ran-on-12-Jul-2025_rev1.mat',...
%     'dropletRetrieval_HySICS_35bands_withNoise_20mm-totalCWV_sim-ran-on-13-Jul-2025_rev1.mat',...
%     'dropletRetrieval_HySICS_66bands_withNoise_cwvRetrieval_sim-ran-on-12-Jul-2025_rev1.mat'};

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
C = mySavedColors(61:(61+length(filenames) +1), 'fixed');

lgnd_str = cell(1, length(filenames) + 2);


figure;


for nn = 1:length(filenames)


    % Load a data set
    ds = load([folder_paths.HySICS_retrievals, filenames{nn}]);

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
        simulated_CWV = aboveCloud_CWV_simulated_hysics_spectra(ds.simulated_measurements.inputs); % kg/m^2

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


% --- *** Overlay new retrieval results *** ---
% define a new set of colors
C2 = mySavedColors(61+length(filenames) +1 : 61+length(filenames) +1 + length(filenames)+1, 'fixed');
lgnd_str2 = cell(1, length(filenames));

% Load a data set
for nn = 1:length(filenames)

    ds = load([filenames(nn).folder, '/', filenames(nn).name]);

    errorbar(ds.GN_outputs.retrieval(2,end), ds.GN_outputs.retrieval(1,end), sqrt(ds.GN_outputs.posterior_cov(1,1))/2,...
        sqrt(ds.GN_outputs.posterior_cov(1,1))/2, sqrt(ds.GN_outputs.posterior_cov(2,2))/2,...
        sqrt(ds.GN_outputs.posterior_cov(2,2))/2, 'MarkerFaceColor', C2(nn+1,:),...
        'MarkerEdgeColor', C2(nn+1,:), 'Linewidth', 4, 'Marker', '.', 'MarkerSize', 40,...
        'Color', C2(nn+1,:))

    hold on


    % what was the retrieved above cloud column water vapor path above
    % cloud?
    retrieved_CWV = ds.GN_outputs.retrieval(end, end);        % kg/m^2 (mm)

    lgnd_str2{nn} = ['*new retrieval* - ', num2str(numel(ds.GN_inputs.bands2run)), ' bands - ',...
        num2str(round(ds.GN_inputs.measurement.uncertainty(1)*100,1)), '$\%$ uncert - retrieved $acpw$ = ',...
        num2str(round(retrieved_CWV, 2)), ' $mm$'];

end



% ----------------------------------------------


% Create a Legend with only the two black curves
legend([lgnd_str, lgnd_str2], 'Interpreter','latex', 'Location','northwest', 'FontSize', 20)

grid on; grid minor
ylabel('$r_{top}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',30)
xlabel('$r_{bot}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',30)

set(gcf,'Position',[0 0 950 750])







%% Plot the GRC results and overlay the three retrievals with measurement uncertainty of 0%, 0.3% and 1% using the EuroVisit algorithm

% Plot full retrieval with ACPW only

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
    folder_paths.drive = '/Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/testGRC_results/';


elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    % ---- Define where the retrievals are stored ---
    folder_paths.drive = '/Users/andrewbuggee/MATLAB-Drive/HySICS/Droplet_profile_retrievals/testGRC_results/';


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
C = mySavedColors(61:(61+length(filenames) +1), 'fixed');

lgnd_str = cell(1, length(filenames) + 2);


figure;


for nn = 1:length(filenames)


    % Load a data set
    ds = load([filenames(nn).folder, '/', filenames(nn).name]);

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


    % give a different marker type for the retrieval using 66 bands
    % that also retrieved above cloud column water vapor


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

    lgnd_str{nn+2} = ['*new retrieval* - ', num2str(numel(ds.GN_inputs.bands2run)), ' bands - ',...
        num2str(round(ds.GN_inputs.measurement.uncertainty(1)*100,1)), '$\%$ uncert - retrieved $acpw$ = ',...
        num2str(round(retrieved_CWV, 2)), ' $mm$'];



end








% ----------------------------------------------


% Create a Legend with only the two black curves
legend(lgnd_str, 'Interpreter','latex', 'Location','northwest', 'FontSize', 20)

grid on; grid minor
ylabel('$r_{top}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',30)
xlabel('$r_{bot}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',30)

set(gcf,'Position',[0 0 950 750])









%% Plot the GRC results and overlay the three retrievals with measurement uncertainty of 0%, 0.3% and 1% using the EuroVisit algorithm

% Using the new measurements and deriving all 15 retrievals 

% Plot full retrieval and no-ACPW retrieval

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
    folder_paths.drive = '/Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/testGRC_results2/';


elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    % ---- Define where the retrievals are stored ---
    folder_paths.drive = '/Users/andrewbuggee/MATLAB-Drive/HySICS/Droplet_profile_retrievals/testGRC_results2/';


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
C = mySavedColors(61:(61+length(filenames) +1), 'fixed');

lgnd_str = cell(1, length(filenames) + 2);


figure;


for nn = 1:length(filenames)


    % Load a data set
    ds = load([filenames(nn).folder, '/', filenames(nn).name]);

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








%% Plot the GRC results and overlay the three retrievals with measurement uncertainty of 0%, 0.3% and 1% using the EuroVisit algorithm

% Using the new measurements and deriving all 15 retrievals 

% Plot full retrieval and no-ACPW retrieval

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
    folder_paths.drive = '/Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/testGRC_results3/';


elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    % ---- Define where the retrievals are stored ---
    folder_paths.drive = '/Users/andrewbuggee/MATLAB-Drive/HySICS/Droplet_profile_retrievals/testGRC_results3/';


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
C = mySavedColors(55:(55+length(filenames) +1), 'fixed');

lgnd_str = cell(1, length(filenames) + 2);


figure;


for nn = 1:length(filenames)


    % Load a data set
    ds = load([filenames(nn).folder, '/', filenames(nn).name]);

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








%% Plot the GRC results and overlay the three retrievals with measurement uncertainty of 0%, 0.3% and 1% using the EuroVisit algorithm

% Using the new measurements and the new a priori values and a priori
% covariance matrix. 

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
    folder_paths.drive = '/Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/testGRC_results4/';


elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    % ---- Define where the retrievals are stored ---
    folder_paths.drive = '/Users/andrewbuggee/MATLAB-Drive/HySICS/Droplet_profile_retrievals/testGRC_results4/';


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
C = mySavedColors(61:(61+length(filenames) +1), 'fixed');

lgnd_str = cell(1, length(filenames) + 2);


figure;


for nn = 1:length(filenames)


    % Load a data set
    ds = load([filenames(nn).folder, '/', filenames(nn).name]);

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






%% Plot the 3 spectrum to view differences - should be different since there is added noise



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
    folder_paths.simulated_measurements = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/testGRC_results/'];


elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    % ---- Define where the retrievals are stored ---
    folder_paths.simulated_measurements = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/testGRC_results'];


elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------

end


% Grab filenames in drive
filenames = dir(folder_paths.simulated_measurements);
idx_2delete = [];
for nn = 1:length(filenames)

    if strcmp(filenames(nn).name(1), 's')~=true

        idx_2delete = [idx_2delete, nn];

    end

end

% delete rows that don't have retrieval filenames
filenames(idx_2delete) = [];



% Step through each file

% define the colors for each curve plotted
C = mySavedColors(61:(61+length(filenames) -1), 'fixed');

lgnd_str = cell(1, length(filenames));


figure;


for nn = 1:length(filenames)


    % Load a data set
    ds = load([filenames(nn).folder, '/', filenames(nn).name]);

    if isfield(ds, "Refl_model_with_noise")==true

        plot(mean(ds.inputs.RT.wavelengths2run,2), ds.Refl_model_with_noise, 'Color', C(nn,:))

    else

        plot(mean(ds.inputs.RT.wavelengths2run,2), ds.Refl_model, 'Color', C(nn,:))

    end

    hold on

    lgnd_str{nn} = [num2str(numel(ds.inputs.bands2run)), ' bands - ',...
        num2str(round(ds.inputs.measurement.uncert*100,1)), '$\%$ uncert'];



end


% Create a Legend with only the two black curves
legend(lgnd_str, 'Interpreter','latex', 'Location','northwest', 'FontSize', 20)

grid on; grid minor
ylabel('$r_{top}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',30)
xlabel('$r_{bot}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',30)

set(gcf,'Position',[0 0 950 750])