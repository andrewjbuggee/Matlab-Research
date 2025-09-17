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

