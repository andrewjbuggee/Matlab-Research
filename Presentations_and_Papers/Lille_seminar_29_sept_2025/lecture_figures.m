%% Reflectance spectra example over cloudy pixel over the HySICS range, with MODIS bands overlaid

clear variables

% --- Load EMIT Data ---

% -------------------------------------
% ------- PICK EMIT DATA SET  --------
% -------------------------------------

emitDataFolder = '17_Jan_2024_coast/';

% -------------------------------------


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
    filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-02-Jun-2025_rev1.mat';



end

ds = load([folderpath, filename]);

% --- Create plot ---

figure;
plot(mean(ds.spec_response.wavelength, 2), ds.Refl_model, 'Color', 'k')
xlabel('Wavelength ($nm$)', Interpreter='latex', FontSize=30)
ylabel('Reflectance ($1/sr$)', Interpreter='latex', FontSize=30)
grid on; grid minor


% set figure size
set(gcf, 'Position', [0 0 1250 500])




% --- shows the spectral bands used in the hyperspectral retireval ---
% bands2run = [49, 57, 69, 86, 103, 166, 169, 171, 174, 217, 220,...
%     222, 224, 227, 237, 288, 290, 293, 388, 390, 393,...
%     426, 434, 436, 570, 574, 577, 579, 582, 613, 616,...
%     618, 620, 623, 625]';

% The same 35 spectral channels above that avoid water vapor and other
% gaseous absorbers, AND 31 bands in the wings of water vapor absorption
% features for a total of 66 bands
% bands2run = [49, 57, 69, 86, 103, 166, 169, 171, 174, 180, 188,...
%     198, 217, 220, 245, 249, 254, 264, 222, 224, 227, 237, 288, 290, 293,...
%     346, 351, 354, 360, 365, 367, 372, 379, 388, 390, 393, 426, 434, 436,...
%     462, 468, 469, 520, 524, 525, 526, 527, 530, 531, 533, 535, 537, 539,...
%     543, 547, 570, 574, 577, 579, 582, 613, 616,618, 620, 623, 625]';

% plot the MODIS bands within this range
bands2run = modisBands(1:19);

% define the color of the filled patch

for bb = 1:length(bands2run)

    hold on
    % plot the bands used as transparent area
    x = [bands2run(bb,1), bands2run(bb,2), bands2run(bb,2), bands2run(bb,1)];
    y = [1,1, 0,0];
%     y = [1e5,1e5, 1e-15,1e-15];
    fill(x,y, mySavedColors(65, 'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 1)

%     hold on
%     xline(mean(ds.spec_response.wavelength(bands2run(bb), :)), 'Color', mySavedColors(61, 'fixed'),...
%         'linewidth', 2)

end

% Define x limits
xlim([300, 2400])

% Create legend
legend('HySICS Reflectance', 'MODIS Channels', 'Location', 'best',...
    'Interpreter', 'latex', 'FontSize', 25, 'Position',[0.690049780273438 0.8386 0.25 0.129])




% ---------- Save figure --------------
% % save .fig file
% folderpath_figs = '/Users/andrewbuggee/Documents/CU-Boulder-ATOC/My Papers/Paper 1/Figures/';
% f = gcf;
% saveas(f,[folderpath_figs,'Fig 7 - reflectance for cloudy scene with 35 wavelengths used in LUT analysis.fig']);
% 
% 
% % save .png with 400 DPI resolution
% % remove title
% title('')
% folderpath_pngs = '/Users/andrewbuggee/Documents/CU-Boulder-ATOC/My Papers/Paper 1/Resubmission Figures/';
% exportgraphics(f,[folderpath_pngs,'Fig 7 - reflectance for cloudy scene with 35 wavelengths used in LUT analysis.png'],'Resolution', 400);
% 




%% Plot the retrieved droplet size at cloud top and bottom for the 4 profiles with different amounts of assumed CWV and
% the one file with retrieved column water vapor

% *** Using EMIT data ***

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
    folder_paths.HySICS_retrievals = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'EMIT/EMIT_data/27_Jan_2024/Droplet_profile_retrievals/'];


elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    % ---- Define where the retrievals are stored ---
    folder_paths.HySICS_retrievals = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'EMIT/EMIT_data/27_Jan_2024/Droplet_profile_retrievals/'];




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
filenames = {'35bands_assumed-15mm_totalCWV_ran-on-17-Jul-2025_rev1.mat',...
    '35bands_assumed-20mm_totalCWV_ran-on-17-Jul-2025_rev1.mat',...
    '35bands_assumed-30mm_totalCWV_ran-on-17-Jul-2025_rev1.mat',...
        '64bands_ran-on-16-Jul-2025_rev6_useForGRC.mat'};


% Step through each file

% define the colors for each curve plotted
C = mySavedColors(61:(61+length(filenames)+2), 'fixed');

lgnd_str = cell(1, length(filenames) + 1);

load('35bands_assumed-26mm_totalCWV_ran-on-16-Jul-2025_rev5.mat', 'tblut_retrieval')


figure;


for nn = 1:length(filenames)


    % Load a data set
    ds = load([folder_paths.HySICS_retrievals, filenames{nn}]);

    if nn==1


        xline(tblut_retrieval.minTau, ':', ['MODIS $\tau_{c}$'],...
            'Fontsize',20, 'Interpreter','latex','LineWidth',2,'Color', [87/255, 90/255, 91/255], 'LabelHorizontalAlignment','right',...
            'LabelVerticalAlignment','top');

        hold on

        yline(tblut_retrieval.minRe, ':', ['MODIS $r_{e}$'],...
            'Fontsize',20, 'Interpreter','latex','LineWidth',2,'Color', [87/255, 90/255, 91/255], 'LabelHorizontalAlignment','right',...
            'LabelVerticalAlignment','top');
        
        hold on


        % Show the MODIS retrieved above cloud column water vapor

        title(['Coincident MODIS retrieval: $acpw$ = 26 $ mm$'],...
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
        e1 = errorbar(ds.GN_outputs.retrieval(3,end), ds.GN_outputs.retrieval(1,end),...
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


        e2 = errorbar(ds.GN_outputs.retrieval(3,end), ds.GN_outputs.retrieval(1,end),...
            1/2 * sqrt(ds.GN_outputs.posterior_cov(1,1) + ds.GN_outputs.posterior_cov(2,3)),...
            1/2 * sqrt(ds.GN_outputs.posterior_cov(1,1) + ds.GN_outputs.posterior_cov(2,3)),...
            sqrt(ds.GN_outputs.posterior_cov(3,3))/2,...
            sqrt(ds.GN_outputs.posterior_cov(3,3))/2, 'MarkerFaceColor', C(nn+2,:),...
            'MarkerEdgeColor', C(nn+2,:), 'Linewidth', 4, 'Marker', '.', 'MarkerSize', 40,...
            'Color', C(nn+2,:));

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
legend(lgnd_str, 'Interpreter','latex', 'Location','best', 'FontSize', 20)


grid on; grid minor
ylabel('$r_{top}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',30)
xlabel('$\tau_c$ ', 'Interpreter','latex', 'FontSize',30)

set(gcf,'Position',[0 0 950 750])


