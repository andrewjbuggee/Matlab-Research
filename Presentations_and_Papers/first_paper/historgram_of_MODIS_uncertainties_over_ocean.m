%% For Zhibo Zhang and other reviewers - Histogram of various pixel-level properties for MODIS observations 
% over the ocean observing liquid water clouds with an optical depth of at
% least 3


clear variables


% Determine which computer you're using

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(whatComputer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------

    % ***** Define the MODIS Folder *****

    modisFolder = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'MODIS_Cloud_Retrieval/MODIS_data/'];


elseif strcmp(whatComputer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    % ----- Define the MODIS folder name -----

    modisFolder = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'MODIS_Cloud_Retrieval/MODIS_data/'];


elseif strcmp(whatComputer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------

    % Define the MODIS folder name

    modisFolder = '/projects/anbu8374/MODIS_data/';


end


% Loop through three MODIS scenes used in this analysis and store retrieved
% effective radius and optical thickness

% ----- November 9th at decimal time 0.611 (14:40) -----
% ----- November 11th at decimal time 0.604 (14:30) -----
% ----- November 11th at decimal time 0.784 (18:50) -----   
modisData = {'2008_11_09/', '2008_11_11_1430/', '2008_11_11_1850/'};

% variables to keep
sig_Refl = [];

effRad_uncert_17_percent = [];
effRad_uncert_17 = [];

optThickness_uncert_17_precent = [];
optThickness_uncert_17 = [];

% limit to only pixels over ocean with an optical depth of at least 3
tau_min = 3;


for nn = 1:length(modisData)

    [modis,L1B_fileName] = retrieveMODIS_data([modisFolder, modisData{nn}]);

    % this stacks columns one on top of the another
    modis_lat = double(modis.geo.lat(:));
    modis_long = double(modis.geo.long(:));
    
    isOcean = land_or_ocean(modis_lat, modis_long, 15, false);

    % rearrance isOcean back to the same size as the MODIS swath
    % reshape will restack the columns
    isOcean = reshape(isOcean, size(modis.cloud.effRadius17,1), size(modis.cloud.effRad_uncert_17, 2));
    
    % Show pixels over ocean
    figure; 
    geoscatter(modis_lat(isOcean), modis_long(isOcean), 100, modis.cloud.effRadius17(isOcean), '.');
    colorbar


    % grab indexes for all liquid water clouds
    idx_liquidWater = modis.cloud.phase==1;
    
    % ------------- effective radius -------------------
    % grab effective radius uncertainty over ocean and convert from percent uncertainty to microns
    temp_effRad_uncert_17 = modis.cloud.effRadius17(isOcean) .* modis.cloud.effRad_uncert_17(isOcean) .* 0.01;          % microns

    % find values whose uncertainty is less than 0 or nan (applied to non-cloudy
    % pixels) and with an optical depth of less than 3
    pix_lessThan0_or_Nan_or_tauLessThan = temp_effRad_uncert_17<0 | isnan(temp_effRad_uncert_17) | modis.cloud.optThickness17(isOcean) < tau_min ...
                                          | modis.cloud.phase(isOcean)~=1;

    % remove values found above
    temp_effRad_uncert_17(pix_lessThan0_or_Nan_or_tauLessThan) = [];

    % store the effective radius index for warm cloud pixels with an
    % optical depth of at least 3, and a positive non-NAN effective radius
    % uncertainty
    effRad_uncert_17 = [effRad_uncert_17; temp_effRad_uncert_17];


    % grab effective radius uncertainty over ocean 
    temp_effRad_uncert_17_percent = modis.cloud.effRad_uncert_17(isOcean);          % percent of retrieved effective radius
    
    % find values whose uncertainty is less than 0 or nan (applied to non-cloudy
    % pixels) and with an optical depth of less than 3
    pix_lessThan0_or_Nan_or_tauLessThan = temp_effRad_uncert_17_percent<0 | isnan(temp_effRad_uncert_17_percent) | modis.cloud.optThickness17(isOcean) < tau_min;

    % remove values found above
    temp_effRad_uncert_17_percent(pix_lessThan0_or_Nan_or_tauLessThan) = [];

    % store remaining values in global variable
    effRad_uncert_17_percent = [effRad_uncert_17_percent; temp_effRad_uncert_17_percent];



    
    % ------------- optical depth -------------------
    % grab optical depth uncertainty over ocean and convert from percent
    % uncertainty to opical depth
    temp_optThickness_uncert_17 = modis.cloud.optThickness17(isOcean) .* modis.cloud.optThickness_uncert_17(isOcean) .* 0.01;          % optical depth

    % find values whose uncertainty is less than 0 or nan (applied to non-cloudy
    % pixels) and with an optical depth of less than 3
    pix_lessThan0_or_Nan_or_tauLessThan = temp_optThickness_uncert_17<0 | isnan(temp_optThickness_uncert_17) | modis.cloud.optThickness17(isOcean) < tau_min;

    % remove values found above
    temp_optThickness_uncert_17(pix_lessThan0_or_Nan_or_tauLessThan) = [];

    % store remaining values in a global variable
    optThickness_uncert_17 = [optThickness_uncert_17; temp_optThickness_uncert_17];        % microns


    % grab optical depth uncertainty over ocean 
    temp_optThickness_uncert_17_percent = modis.cloud.optThickness_uncert_17(isOcean);           % percent of retrieved optical depth
    
    % find values whose uncertainty is less than 0 or nan (applied to non-cloudy
    % pixels) and with an optical depth of less than 3
    pix_lessThan0_or_Nan_or_tauLessThan = temp_optThickness_uncert_17_percent<0 | isnan(temp_optThickness_uncert_17_percent) | modis.cloud.optThickness17(isOcean) < tau_min;

    % remove values found above
    temp_optThickness_uncert_17_percent(pix_lessThan0_or_Nan_or_tauLessThan) = [];

    % store remaining values in global variable
    optThickness_uncert_17_precent = [optThickness_uncert_17_precent; temp_optThickness_uncert_17_percent];
    

    
  

end


% create histogram of COD and CER uncertainty for pixels over ocean
figure; 

subplot(1,2,1)
histogram(effRad_uncert_17)
grid on; grid minor
ylabel('Counts', 'interpreter', 'latex')
xlabel('$\delta r_e$ $(\mu m)$', 'interpreter', 'latex')
set(gca, 'YScale', 'log')

subplot(1,2,2)
histogram(optThickness_uncert_17)
grid on; grid minor
ylabel('Counts', 'interpreter', 'latex')
xlabel('$\delta \tau_c$', 'interpreter', 'latex')
set(gca, 'YScale', 'log')

set(gcf, 'Position', [0,0, 1300, 750])


figure; 

subplot(1,2,1)
histogram(effRad_uncert_17_percent)
grid on; grid minor
ylabel('Counts', 'interpreter', 'latex')
xlabel('$\delta r_e$ $(\%)$', 'interpreter', 'latex')
set(gca, 'YScale', 'log')

subplot(1,2,2)
histogram(optThickness_uncert_17_precent)
grid on; grid minor
ylabel('Counts', 'interpreter', 'latex')
xlabel('$\delta \tau_c$  (\%)', 'interpreter', 'latex')
set(gca, 'YScale', 'log')

set(gcf, 'Position', [0,0, 1300, 750])

% print the mean values of uncertainty as a percent for each variable
disp([newline, 'Average CER uncertainty for pixels over ocean: ', ...
    num2str(mean(effRad_uncert_17_percent)), '%, with a standard deviation of: ',...
    num2str(std(effRad_uncert_17_percent)), '%', newline])
disp([newline, 'Average COD uncertainty for pixels over ocean: ', ...
    num2str(mean(optThickness_uncert_17_precent)), '%, with a standard deviation of: ',...
    num2str(std(optThickness_uncert_17_precent)), '%', newline])

% print the mean values of uncertainty for each variable
disp([newline, 'Average CER uncertainty for pixels over ocean: ', ...
    num2str(mean(effRad_uncert_17)), ' microns, with a standard deviation of: ',...
    num2str(std(effRad_uncert_17)), ' microns', newline])
disp([newline, 'Average COD uncertainty for pixels over ocean: ', ...
    num2str(mean(optThickness_uncert_17)), ', with a standard deviation of: ',...
    num2str(std(optThickness_uncert_17)), newline])

% Create plot of retireved effective radius versus retrieval uncertainaty
% and include a one-to-one line



clear variables