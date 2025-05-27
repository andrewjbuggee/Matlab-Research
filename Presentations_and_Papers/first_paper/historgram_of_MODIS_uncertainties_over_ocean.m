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

refl_uncert_precent = [];
refl_uncert = [];

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

    % Show pixels over ocean colored by the effective radius
    % figure;
    % geoscatter(modis_lat(isOcean), modis_long(isOcean), 100, modis.cloud.effRadius17(isOcean), '.');
    % colorbar



    % grab indexes for all liquid water clouds
    % 0 -- cloud mask undetermined
    % 1 -- clear sky
    % 2 -- liquid water cloud
    % 3 -- ice cloud
    % 4 -- undetermined phase cloud (but retrieval is attempted as  liquid water)
    idx_liquidWater = modis.cloud.phase==2;

    % only keep pixels with an optical depth of at least 3
    idx_tau = modis.cloud.optThickness17 >= tau_min;

    % grab all pixels over the ocean with liquid water clouds and an
    % optical depth of at least 3
    idx_master = logical(isOcean .* idx_liquidWater .* idx_tau);

    % --------------------------------------------------
    % ------------- effective radius -------------------
    % --------------------------------------------------
    % grab effective radius uncertainty over ocean and convert from percent uncertainty to microns
    temp_effRad_uncert_17 = modis.cloud.effRadius17(idx_master) .* modis.cloud.effRad_uncert_17(idx_master) .* 0.01;          % microns

    % find values whose uncertainty is less than 0 or nan (applied to non-cloudy
    % pixels) 
    pix_lessThan0_or_Nan = temp_effRad_uncert_17<0 | isnan(temp_effRad_uncert_17);

    % remove values found above
    temp_effRad_uncert_17(pix_lessThan0_or_Nan) = [];

    % store the uncertainty of effective radius in microns
    effRad_uncert_17 = [effRad_uncert_17; temp_effRad_uncert_17];     % microns

    
    % grab effective radius uncertainty over ocean
    temp_effRad_uncert_17_percent = modis.cloud.effRad_uncert_17(idx_master);          % percent of retrieved effective radius

    % find values whose uncertainty is less than 0 or nan (applied to non-cloudy
    % pixels) and with an optical depth of less than 3
    pix_lessThan0_or_Nan = temp_effRad_uncert_17_percent<0 | isnan(temp_effRad_uncert_17_percent);

    % remove values found above
    temp_effRad_uncert_17_percent(pix_lessThan0_or_Nan) = [];

    % store the uncertainty of effective radius as a percent
    effRad_uncert_17_percent = [effRad_uncert_17_percent; temp_effRad_uncert_17_percent];   % percent



    % -----------------------------------------------
    % ------------- optical depth -------------------
    % -----------------------------------------------
    % grab optical depth uncertainty over ocean and convert from percent
    % uncertainty to opical depth
    temp_optThickness_uncert_17 = modis.cloud.optThickness17(idx_master) .* modis.cloud.optThickness_uncert_17(idx_master) .* 0.01;          % optical depth

    % find values whose uncertainty is less than 0 or nan (applied to non-cloudy
    % pixels) and with an optical depth of less than 3
    pix_lessThan0_or_Nan = temp_optThickness_uncert_17<0 | isnan(temp_optThickness_uncert_17);

    % remove values found above
    temp_optThickness_uncert_17(pix_lessThan0_or_Nan) = [];

    % store remaining values in a global variable
    optThickness_uncert_17 = [optThickness_uncert_17; temp_optThickness_uncert_17];        % optical thickness uncertainty


    % grab optical depth uncertainty over ocean
    temp_optThickness_uncert_17_percent = modis.cloud.optThickness_uncert_17(idx_master);           % percent of retrieved optical depth

    % find values whose uncertainty is less than 0 or nan (applied to non-cloudy
    % pixels) and with an optical depth of less than 3
    pix_lessThan0_or_Nan = temp_optThickness_uncert_17_percent<0 | isnan(temp_optThickness_uncert_17_percent);

    % remove values found above
    temp_optThickness_uncert_17_percent(pix_lessThan0_or_Nan) = [];

    % store remaining values in global variable as a percent
    optThickness_uncert_17_precent = [optThickness_uncert_17_precent; temp_optThickness_uncert_17_percent];      % optical thickness unceratinty as a percent


    % --------------------------------------------------
    % ---------------- Reflectance ---------------------
    % --------------------------------------------------

    temp_Refl_uncert = [];
    temp_Refl_uncert_percent = [];
    % grab reflectance uncertainty from the first seven 
    % MODIS channels and convert it from a percentage to a
    % reflectance
    for ww = 1:size(modis.EV1km.reflectance,3)
        
        temp_R = modis.EV1km.reflectance(:,:,ww);
        temp_U = modis.EV1km.reflectanceUncert(:,:,ww);

        temp_Refl_uncert(:,ww) = temp_R(idx_master) .* temp_U(idx_master).* 0.01;          % 1/sr - reflectance

        temp_Refl_uncert_percent(:,ww) = temp_U(idx_master);          % 1/sr - reflectance

    end

    % find values whose uncertainty is less than 0 
    pix_lessThan0 = temp_Refl_uncert<0;


    % remove values found above
    temp_Refl_uncert(pix_lessThan0) = NaN;

    % store remaining values in a global variable
    refl_uncert = [refl_uncert; temp_Refl_uncert];        % 1/sr - reflectance uncertainty


    % grab reflectance uncertainty from the first seven 
    % MODIS channels and convert it from a percentage to a
    % reflectance

    % find values whose uncertainty is less than 0 or nan (applied to non-cloudy
    % pixels) and with an optical depth of less than 3
    pix_lessThan0 = temp_Refl_uncert_percent<0;

    % remove values found above
    temp_Refl_uncert_percent(pix_lessThan0) = NaN;

    % store remaining values in a global variable
    refl_uncert_precent = [refl_uncert_precent; temp_Refl_uncert_percent];        % 1/sr - reflectance uncertainty



end

clear modis

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


% create histogram of COD and CER uncertainty as a percent for pixels over ocean

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




% print the mean values of reflectance uncertainty reflectance uncertainty
% percent
disp([newline, 'Average Reflectance uncertainty for pixels over ocean per channel: ', newline])
for ww = 1:size(refl_uncert,2)

    % grab MODIS wavelenghts
    sr = modis_terra_specResponse_func_2(ww,true);
    modis_wl_center(ww) = mean(sr.wavelength);

    if ww == 6
        % *** One data set is faulty, with uncertanties in the 6th
        % channel (1.6 microns) at 30%! Ignore these values.
        pix_lessThan20 = refl_uncert_precent(:,ww)<20;
        % remove these values
        temp_Refl_percent = refl_uncert_precent(pix_lessThan20, ww);
        temp_Refl = refl_uncert(pix_lessThan20, ww);

        disp(['     Band ', num2str(ww), ' - center wl ',num2str(modis_wl_center(ww)),...
            ' (nm) - ', num2str(mean(temp_Refl)), ' 1/sr or, ',...
            num2str(mean(temp_Refl_percent)), '%'])

    else

        disp(['     Band ', num2str(ww), ' - center wl ',num2str(modis_wl_center(ww)),...
            ' (nm) - ', num2str(mean(refl_uncert(:,ww))), ' 1/sr or, ',...
            num2str(mean(refl_uncert_precent(:,ww))), '%'])

    end

end

disp([newline, 'Standard Deviation of Reflectance uncertainty for pixels over ocean per channel: ', newline])
for ww = 1:size(refl_uncert,2)

    % grab MODIS wavelenghts
    sr = modis_terra_specResponse_func_2(ww,true);
    modis_wl_center(ww) = mean(sr.wavelength);

    if ww == 6
        % *** One data set is faulty, with uncertanties in the 6th
        % channel (1.6 microns) at 30%! Ignore these values.
        pix_lessThan20 = refl_uncert_precent(:,ww)<20;
        % remove these values
        temp_Refl_percent = refl_uncert_precent(pix_lessThan20, ww);
        temp_Refl = refl_uncert(pix_lessThan20, ww);

        disp(['     Band ', num2str(ww), ' - center wl ',num2str(modis_wl_center(ww)),...
            ' (nm) - ', num2str(std(temp_Refl)), ' 1/sr or, ',...
            num2str(std(temp_Refl_percent)), '%'])

    else

        disp(['     Band ', num2str(ww), ' - center wl ',num2str(modis_wl_center(ww)),...
            ' (nm) - ', num2str(std(refl_uncert(:,ww))), ' 1/sr or, ',...
            num2str(std(refl_uncert_precent(:,ww))), '%'])

    end
end




% Create plot of retireved effective radius versus retrieval uncertainaty
% and include a one-to-one line



clear variables