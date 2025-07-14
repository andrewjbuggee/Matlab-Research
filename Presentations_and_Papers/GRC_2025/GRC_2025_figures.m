%% Figures for GRC Poster!!

% By Andrew John Buggee

%% Plot droplet profile retrieval using 35 wavelengths with different forward model assumptions
% along with the retrieval using 66 wavelengths and the retrieval of CWV

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


elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    % ---- Define where the retrievals are stored ---
    folder_paths.HySICS_retrievals = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/HySICS/Droplet_profile_retrievals/'];




elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------

end


% define the mat files for each retreival to plot
% filenames = {'dropletRetrieval_HySICS_35bands_withNoise_10mm-totalCWV_sim-ran-on-12-Jul-2025_rev1.mat',...
%     'dropletRetrieval_HySICS_35bands_withNoise_15mm-totalCWV_sim-ran-on-12-Jul-2025_rev1.mat',...
%     'dropletRetrieval_HySICS_35bands_withNoise_20mm-totalCWV_sim-ran-on-13-Jul-2025_rev1.mat',...
%     'dropletRetrieval_HySICS_35bands_withNoise_25mm-totalCWV_sim-ran-on-12-Jul-2025_rev1.mat',...
%     'dropletRetrieval_HySICS_66bands_withNoise_cwvRetrieval_sim-ran-on-12-Jul-2025_rev1.mat'};


% define the mat files for each retreival to plot
filenames = {'dropletRetrieval_HySICS_35bands_withNoise_10mm-totalCWV_sim-ran-on-12-Jul-2025_rev1.mat',...
    'dropletRetrieval_HySICS_35bands_withNoise_15mm-totalCWV_sim-ran-on-12-Jul-2025_rev1.mat',...
    'dropletRetrieval_HySICS_35bands_withNoise_20mm-totalCWV_sim-ran-on-13-Jul-2025_rev1.mat',...
    'dropletRetrieval_HySICS_66bands_withNoise_cwvRetrieval_sim-ran-on-12-Jul-2025_rev1.mat'};




% Step through each file

% define the colors for each curve plotted
C = mySavedColors(61:(61+length(filenames)+1), 'fixed');

lgnd_str = cell(1, 4*length(filenames) + 1);


figure;


for nn = 1:length(filenames)


    % Load a data set
    ds = load([folder_paths.HySICS_retrievals, filenames{nn}]);

    if nn==1

        % create a droplet profile from simulated measurement inputs
        re_sim = create_droplet_profile2([ds.GN_inputs.RT.r_top, ds.GN_inputs.RT.r_bot],...
            ds.GN_outputs.tau_vector, 'optical_depth', ds.GN_inputs.model.profile.type);

        tau_sim = linspace(0, ds.GN_inputs.RT.tau_c, 100);

        % first, plot the simulated profile
        plot(re_sim, tau_sim, 'Color', 'k','LineStyle','-', 'LineWidth',5)
        hold on


        % what was the assumed above cloud column water vapor path?
        simulated_CWV = aboveCloud_CWV_simulated_hysics_spectra(ds.simulated_measurements.inputs); % kg/m^2

        lgnd_str{nn} = ['Simulated profile - cwv = ',...
            num2str(round(simulated_CWV, 2)), ' $kg/m^2$'];

    end



    % plot the retrieved droplet profile
    if nn<length(filenames)

        plot(ds.GN_outputs.re_profile, ds.GN_outputs.tau_vector', 'Color', C(nn+1,:),'LineStyle',':', 'LineWidth',3)

        hold on

        % plot two markers at cloud bottom and cloud top
        plot([ds.GN_outputs.re_profile(1), ds.GN_outputs.re_profile(end)], [ds.GN_outputs.tau_vector(1), ds.GN_outputs.tau_vector(end)],...
            '.', 'MarkerSize', 23, 'Color', C(nn+1,:))


        % Plot the retrieval uncertainty of the radius at cloud top
        errorbar(ds.GN_outputs.re_profile(1), ds.GN_outputs.tau_vector(1), sqrt(ds.GN_outputs.posterior_cov(1,1)),...
            'horizontal', 'Color', C(nn+1,:), 'markersize', 20, 'Linewidth', 2)

        % Plot the retrieval uncertainty of the radius at cloud bottom
        errorbar(ds.GN_outputs.re_profile(end), ds.GN_outputs.tau_vector(end), sqrt(ds.GN_outputs.posterior_cov(2,2)),...
            'horizontal', 'Color', C(nn+1,:), 'markersize', 20, 'Linewidth', 2)

        % Plot the retrieval uncertainty of the optical depth
        errorbar(ds.GN_outputs.re_profile(end), ds.GN_outputs.tau_vector(end), sqrt(ds.GN_outputs.posterior_cov(3,3)),...
            'vertical', 'Color', C(nn+1,:), 'markersize', 20, 'Linewidth', 2)

    else

        % give a different marker type for the retrieval using 66 bands
        % that also retrieved above cloud column water vapor


        plot(ds.GN_outputs.re_profile, ds.GN_outputs.tau_vector', 'Color', C(nn+1,:),'LineStyle','--', 'LineWidth',3)

        hold on

        % plot two markers at cloud bottom and cloud top
        plot([ds.GN_outputs.re_profile(1), ds.GN_outputs.re_profile(end)], [ds.GN_outputs.tau_vector(1), ds.GN_outputs.tau_vector(end)],...
            '.', 'MarkerSize', 23, 'Color', C(nn+1,:))


        % Plot the retrieval uncertainty of the radius at cloud top
        errorbar(ds.GN_outputs.re_profile(1), ds.GN_outputs.tau_vector(1), sqrt(ds.GN_outputs.posterior_cov(1,1)),...
            'horizontal', 'Color', C(nn+1,:), 'markersize', 20, 'Linewidth', 2)

        % Plot the retrieval uncertainty of the radius at cloud bottom
        errorbar(ds.GN_outputs.re_profile(end), ds.GN_outputs.tau_vector(end), sqrt(ds.GN_outputs.posterior_cov(2,2)),...
            'horizontal', 'Color', C(nn+1,:), 'markersize', 20, 'Linewidth', 2)

        % Plot the retrieval uncertainty of the optical depth
        errorbar(ds.GN_outputs.re_profile(end), ds.GN_outputs.tau_vector(end), sqrt(ds.GN_outputs.posterior_cov(3,3)),...
            'vertical', 'Color', C(nn+1,:), 'markersize', 20, 'Linewidth', 2)


    end



    % create the legend string
    if nn<length(filenames)

        % what was the assumed above cloud column water vapor path?
        assumed_CWV = aboveCloud_CWV_simulated_hysics_spectra(ds.GN_inputs); % kg/m^2

        lgnd_str{nn+1 + 4*(nn-1)} = [num2str(numel(ds.GN_inputs.bands2run)), ' bands - assumed cwv = ',...
            num2str(round(assumed_CWV,2)), ' $kg/m^2$'];

        % Leave blanks inbetween
        lgnd_str{nn+1 + 4*(nn-1) +1} = '';
        lgnd_str{nn+1 + 4*(nn-1) +2} = '';
        lgnd_str{nn+1 + 4*(nn-1) +3} = '';
        lgnd_str{nn+1 + 4*(nn-1) +4} = '';


    else

        % create the string for the retrieval using CWV

        % what was the retrieved above cloud column water vapor path above
        % cloud?
        retrieved_CWV = ds.GN_outputs.retrieval(end, end);        % kg/m^2 (mm)

        lgnd_str{nn+1 + 4*(nn-1)} = [num2str(numel(ds.GN_inputs.bands2run)), ' bands - retrieved cwv = ',...
            num2str(round(retrieved_CWV, 2)), ' $kg/m^2$'];

    end

end


% flip the axes
set(gca, 'YDir','reverse')



% Create a Legend with only the two black curves
legend(lgnd_str, 'Interpreter','latex', 'Location','northwest', 'FontSize', 25)



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

% Plot the z-space in meters on the right axis
yyaxis right
ylim([0, abs(ds.GN_inputs.RT.z(end) - ds.GN_inputs.RT.z(1))])
set(gca,'YColor','black')
ylabel('Altitude within cloud $(m)$', 'Interpreter','latex','FontSize',30);
yyaxis left

ylabel('Optical Depth', 'Interpreter','latex', 'FontSize',30)
xlabel('Effective Radius ($\mu m$)', 'Interpreter','latex', 'FontSize',30)


grid on; grid minor

ylim([-0.1, 1.2*ds.GN_inputs.RT.tau_c])


set(gcf,'Position',[0 0 1200 630])



%% Plot the retrieved droplet size at cloud top and bottom for the 4 profiles with different amounts of assumed CWV and
% the one file with retrieved column water vapor



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


elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    % ---- Define where the retrievals are stored ---
    folder_paths.HySICS_retrievals = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/HySICS/Droplet_profile_retrievals/'];




elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------

end


% define the mat files for each retreival to plot
% filenames = {'dropletRetrieval_HySICS_35bands_withNoise_10mm-totalCWV_sim-ran-on-12-Jul-2025_rev1.mat',...
%     'dropletRetrieval_HySICS_35bands_withNoise_15mm-totalCWV_sim-ran-on-12-Jul-2025_rev1.mat',...
%     'dropletRetrieval_HySICS_35bands_withNoise_20mm-totalCWV_sim-ran-on-13-Jul-2025_rev1.mat',...
%     'dropletRetrieval_HySICS_35bands_withNoise_25mm-totalCWV_sim-ran-on-12-Jul-2025_rev1.mat',...
%     'dropletRetrieval_HySICS_66bands_withNoise_cwvRetrieval_sim-ran-on-12-Jul-2025_rev1.mat'};

% define the mat files for each retreival to plot
filenames = {'dropletRetrieval_HySICS_35bands_withNoise_10mm-totalCWV_sim-ran-on-12-Jul-2025_rev1.mat',...
    'dropletRetrieval_HySICS_35bands_withNoise_15mm-totalCWV_sim-ran-on-12-Jul-2025_rev1.mat',...
    'dropletRetrieval_HySICS_35bands_withNoise_20mm-totalCWV_sim-ran-on-13-Jul-2025_rev1.mat',...
    'dropletRetrieval_HySICS_66bands_withNoise_cwvRetrieval_sim-ran-on-12-Jul-2025_rev1.mat'};


% Step through each file

% define the colors for each curve plotted
C = mySavedColors(61:(61+length(filenames)+1), 'fixed');

lgnd_str = cell(1, length(filenames) + 2);


figure;


for nn = 1:length(filenames)


    % Load a data set
    ds = load([folder_paths.HySICS_retrievals, filenames{nn}]);

    if nn==1


        % first, plot the simulated profile values as two lines
        xline(ds.GN_inputs.RT.r_bot, '-', ['Simulated $r_{bot}$'],...
            'Fontsize',20, 'Interpreter','latex','LineWidth',2,'Color','black', 'LabelHorizontalAlignment','left',...
            'LabelVerticalAlignment','bottom');

        hold on

        yline(ds.GN_inputs.RT.r_top, '-', ['Simulated $r_{top}$'],...
            'Fontsize',20, 'Interpreter','latex','LineWidth',2,'Color','black', 'LabelHorizontalAlignment','right',...
            'LabelVerticalAlignment','top');
        hold on


        % what was the assumed above cloud column water vapor path?
        simulated_CWV = aboveCloud_CWV_simulated_hysics_spectra(ds.simulated_measurements.inputs); % kg/m^2

        title(['Simulated profile - cwv = ',num2str(round(simulated_CWV, 2)), ' $kg/m^2$'],...
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
        errorbar(ds.GN_outputs.retrieval(2,end), ds.GN_outputs.retrieval(1,end), sqrt(ds.GN_outputs.posterior_cov(1,1))/2,...
            sqrt(ds.GN_outputs.posterior_cov(1,1))/2, sqrt(ds.GN_outputs.posterior_cov(2,2))/2,...
            sqrt(ds.GN_outputs.posterior_cov(2,2))/2, 'MarkerFaceColor', C(nn+1,:),...
            'MarkerEdgeColor', C(nn+1,:),'markersize', 15, 'Linewidth', 2, 'Marker', '.', 'MarkerSize', 35,...
            'Color', C(nn+1,:))

        hold on



    else

        % give a different marker type for the retrieval using 66 bands
        % that also retrieved above cloud column water vapor


        errorbar(ds.GN_outputs.retrieval(2,end), ds.GN_outputs.retrieval(1,end), sqrt(ds.GN_outputs.posterior_cov(1,1))/2,...
            sqrt(ds.GN_outputs.posterior_cov(1,1))/2, sqrt(ds.GN_outputs.posterior_cov(2,2))/2,...
            sqrt(ds.GN_outputs.posterior_cov(2,2))/2, 'MarkerFaceColor', C(nn+1,:),...
            'MarkerEdgeColor', C(nn+1,:),'markersize', 15, 'Linewidth', 2, 'Marker', 'square', 'MarkerSize', 12,...
            'Color', C(nn+1,:))

        hold on


    end



    % create the legend string
    if nn<length(filenames)

        % what was the assumed above cloud column water vapor path?
        assumed_CWV = aboveCloud_CWV_simulated_hysics_spectra(ds.GN_inputs); % kg/m^2

        lgnd_str{nn+2} = [num2str(numel(ds.GN_inputs.bands2run)), ' bands - assumed cwv = ',...
            num2str(round(assumed_CWV,2)), ' $kg/m^2$'];





    else

        % create the string for the retrieval using CWV

        % what was the retrieved above cloud column water vapor path above
        % cloud?
        retrieved_CWV = ds.GN_outputs.retrieval(end, end);        % kg/m^2 (mm)

        lgnd_str{nn+2} = [num2str(numel(ds.GN_inputs.bands2run)), ' bands - retrieved cwv = ',...
            num2str(round(retrieved_CWV, 2)), ' $kg/m^2$'];

    end

end


% Create a Legend with only the two black curves
legend(lgnd_str, 'Interpreter','latex', 'Location','northwest', 'FontSize', 20)

grid on; grid minor
ylabel('$r_{top}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',30)
xlabel('$r_{bot}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',30)

set(gcf,'Position',[0 0 950 750])



%% Plot 35 wavelengths over a cloud Hysics Spectrum
% Reflectance spectra example over cloudy pixel

clear variables




% Determine which computer you're using

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(whatComputer,'anbu8374')==true

    % ------ Folders on my Mac Desktop --------

    % Define the Simulated HySICS data folder path

    folderpath = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/';

    % filename
    filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-02-Jun-2025_ALL_BANDS.mat';



elseif strcmp(whatComputer,'andrewbuggee')==true

    % ------ Folders on my Macbook --------

    % Define the Simulated HySICS data folder path

    folderpath = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/';

    % filename
    filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_allBands_20mm-totalColumnWaterVapor_sim-ran-on-08-Jul-2025_rev1.mat';


end

ds = load([folderpath, filename]);

% --- Create plot ---

figure;
plot(mean(ds.spec_response.wavelength, 2), ds.Refl_model, 'Color', 'k')
xlabel('Wavelength ($nm$)', Interpreter='latex', FontSize=30)
ylabel('Reflectance ($1/sr$)', Interpreter='latex', FontSize=30)
grid on; grid minor




% --- shows the spectral bands used in the hyperspectral retireval ---
bands2run = [49, 57, 69, 86, 103, 166, 169, 171, 174, 217, 220,...
    222, 224, 227, 237, 288, 290, 293, 388, 390, 393,...
    426, 434, 436, 570, 574, 577, 579, 582, 613, 616,...
    618, 620, 623, 625]';


% define the color of the filled patch

for bb = 1:length(bands2run)

    % hold on
    % % plot the bands used as transparent area
    % x = [ds.spec_response.wavelength(bands2run(bb), :),...
    %     fliplr(emit.spec_response.wavelength(inputs.bands2run(bb), wl_range))];
    % %y = [1,1, 0,0];
    % y = [1e5,1e5, 1e-15,1e-15];
    % fill(x,y, C, 'EdgeAlpha', 0, 'FaceAlpha', 1)

    hold on
    xline(mean(ds.spec_response.wavelength(bands2run(bb), :)), 'Color', mySavedColors(61, 'fixed'),...
        'linewidth', 2)

end

% Define x limits
xlim([300, 2400])

% set figure size
set(gcf, 'Position', [0 0 1250 700])

% create title
title('Wavelengths used to retrieve droplet profile and above cloud column water vapor', ...
    'FontSize', 25, 'Interpreter','latex')

% Create legend
legend('HySICS Reflectance', 'Wavelengths used in retrieval', 'Location', 'best',...
    'Interpreter', 'latex', 'FontSize', 25, 'Position',[0.690049780273438 0.7386 0.25 0.129])





%% Plot 66 wavelengths over a cloudy spectrum

% Reflectance spectra example over cloudy pixel

clear variables



% Determine which computer you're using

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(whatComputer,'anbu8374')==true

    % ------ Folders on my Mac Desktop --------

    % Define the Simulated HySICS data folder path

    folderpath = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/';

    % filename
    filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-02-Jun-2025_ALL_BANDS.mat';



elseif strcmp(whatComputer,'andrewbuggee')==true

    % ------ Folders on my Macbook --------

    % Define the Simulated HySICS data folder path

    folderpath = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/';

    % filename
    filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_allBands_20mm-totalColumnWaterVapor_sim-ran-on-08-Jul-2025_rev1.mat';



end

ds = load([folderpath, filename]);

% --- Create plot ---

figure;
plot(mean(ds.spec_response.wavelength, 2), ds.Refl_model, 'Color', 'k')
xlabel('Wavelength ($nm$)', Interpreter='latex', FontSize=30)
ylabel('Reflectance ($1/sr$)', Interpreter='latex', FontSize=30)
grid on; grid minor





% The same 35 spectral channels above that avoid water vapor and other
% gaseous absorbers, AND 31 bands in the wings of water vapor absorption
% features for a total of 66 bands
bands2run = [49, 57, 69, 86, 103, 166, 169, 171, 174, 180, 188,...
    198, 217, 220, 245, 249, 254, 264, 222, 224, 227, 237, 288, 290, 293,...
    346, 351, 354, 360, 365, 367, 372, 379, 388, 390, 393, 426, 434, 436,...
    462, 468, 469, 520, 524, 525, 526, 527, 530, 531, 533, 535, 537, 539,...
    543, 547, 570, 574, 577, 579, 582, 613, 616,618, 620, 623, 625]';

% define the color of the filled patch

for bb = 1:length(bands2run)

    % hold on
    % % plot the bands used as transparent area
    % x = [ds.spec_response.wavelength(bands2run(bb), :),...
    %     fliplr(emit.spec_response.wavelength(inputs.bands2run(bb), wl_range))];
    % %y = [1,1, 0,0];
    % y = [1e5,1e5, 1e-15,1e-15];
    % fill(x,y, C, 'EdgeAlpha', 0, 'FaceAlpha', 1)

    hold on
    xline(mean(ds.spec_response.wavelength(bands2run(bb), :)), 'Color', mySavedColors(61, 'fixed'),...
        'linewidth', 2)

end

% Define x limits
xlim([300, 2400])


% set figure size
set(gcf, 'Position', [0 0 1250 700])

% create title
title('Wavelengths used to retrieve droplet profile and above cloud column water vapor', ...
    'FontSize', 25, 'Interpreter','latex')

% Create legend
legend('HySICS Reflectance', 'Wavelengths used in retrieval', 'Location', 'best',...
    'Interpreter', 'latex', 'FontSize', 25, 'Position',[0.690049780273438 0.7386 0.25 0.129])



