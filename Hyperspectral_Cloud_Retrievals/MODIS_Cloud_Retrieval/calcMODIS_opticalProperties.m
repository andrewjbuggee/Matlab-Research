%% ----- Estimate MODIS Cloud Optical Properties -----


clear variables;
% By Andrew J. Buggee

%% ----- Extract MODIS data of Interest -----


% Define data folders and files that you'd like to read

% add the libradtran directory to the path

% this one is for my LASP computer
% modisINP_folderName = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/MODIS_08_25_2021/';

% this one is for my personal laptop
% modisINP_folderName = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/',...
%     'Hyperspectral-Cloud-Droplet-Retrieval-Research/LibRadTran/libRadtran-2.0.4/MODIS_08_25_2021/'];

% define the files names
if strcmp(whatComputer, 'anbu8374')==true

    folderName = '/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/MODIS_Cloud_Retrieval/MODIS_data/2023_04_13/';

elseif strcmp(whatComputer, 'andrewbuggee') == true

    folderName = '/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/MODIS_Cloud_Retrieval/MODIS_data/2023_04_13/';
end


[modis,L1B_1km_fileName] = retrieveMODIS_data(folderName);

%% ----- Create a structure defining inputs of the problem -----

% this is a built-in function that is defined at the bottom of this script
inputs = create_modis_inputs(folderName, L1B_1km_fileName);

disp('Check inputs to make sure they are what you want!!')


%% ----- Find suitable Pixels! -----

% find pixels within the modis data set that fit our needs
% if we've alraedy done this long calculation, we just load the pixel set
% from a saved .mat file

% lets check to see if there is a suitable pixels .mat file in our folder

pixels_file_flag = isfile([inputs.savedCalculations_folderName,'suitablePixels.mat']);

if inputs.flags.findSuitablePixels == true && pixels_file_flag == false
    
    pixels = findSuitablePixel(modis,inputs);
    
    % save these pixels in a .mat file, along with the inputs
    save([inputs.savedCalculations_folderName,'suitablePixels.mat'],'pixels','inputs');
    
    % now we want only a random subset of pixels to write INP files. inputs
    % defines how many pixels we should grab from the suitable pixels file.
    % This function will also grab the geometry for each pixel and store
    % this information in pixels2use
    
    pixels2use = subset_suitablePixels(inputs,modis);
    
elseif pixels_file_flag == true && inputs.flags.loadPixelSet == false
    
    % we don't need to load the entire pixels file into out workspace. But
    % if we do load a subset of pixels, we need a way to trace back to what
    % pixels we used. Save the pixels used to the data folder listed in the
    % inputs
    
    pixels2use = subset_suitablePixels(inputs,modis);

elseif pixels_file_flag == true && inputs.flags.loadPixelSet == true
    
    error([newline, 'I dont know what you want me to do', newline]);  
end

%% Plot MODIS measured relfectance at 650nm with selected pixels

figure; geoscatter(modis.geo.lat(:), modis.geo.long(:), 10, reshape(modis.EV1km.reflectance(:,:,1),[],1),'.');
cb = colorbar;
set(get(cb, 'label'), 'string', 'Reflectance $(1/sr)$','Interpreter','latex', 'Fontsize',22)
set(gca, 'FontSize',25)
set(gca, 'FontWeight', 'bold')
set(gcf, 'Position', [0 0 800 800])
title('Reflectance $650 \; nm$ band','Interpreter','latex', 'FontSize', 40)
hold on

% convert row column into index
pixels2use.res1km.index = sub2ind(pixels2use.res1km.size, pixels2use.res1km.row, pixels2use.res1km.col);
geoscatter(modis.geo.lat(pixels2use.res1km.index),modis.geo.long(pixels2use.res1km.index),...
    500,"red", '.')



%% ----- Create .INP files for MODIS Geometry -----


if inputs.flags.writeINPfiles == true
    % which pixels on the MODIS array are we using for the gemoetry of the
    % problem? The pixels that we found to be suitable!
    
    %names.inp = write_INP_4_MODIS_hdf(inputs,pixels2use,modis);
    [names.inp, inputs] = write_INP_file_4MODIS_homogenous(inputs, pixels2use, modis);
    
    % now lets write the output names
    
    names.out = writeOutputNames(names.inp);
else
    
    % if the files already exist, just grab the names!
    [names.inp, inputs] = getMODIS_INPnames_withClouds(modis.solar,inputs,pixels2use);
    names.out = writeOutputNames(names.inp);
end

%% ----- Run uvspec and calculate Reflectance Function for Model -----

% geometry stays the same, but we calculate the radiative transfer equation
% for different cloud values (tau,re)

if inputs.flags.runUVSPEC == true
    
    % 1st output - R is the reflectance integrated over a bandwidth
    % 2nd output - Rl is the reflectance at each spectral bin
    tic
    [R,~] = runReflectanceFunction(inputs,names, inputs.spec_response);
    toc
    
elseif inputs.flags.runUVSPEC == false
    
    load([inputs.savedCalculations_folderName,inputs.saveCalculations_fileName] ,'inputs','R');
    
end
%% ----- Compute the Reflectance Function for the MODIS Observations -----

% We don't have to calculate the reflectance function of MODIS if we don't
% want to. They proivde it as an output in their data

% We want to grab modis reflectances at 1km!!

modisR = grab_modis_reflectance(modis,inputs, pixels2use);


%% ----- Compare Reflectance Fucntion of MODIS with Theoretical Calculations (Grid Search) -----

% first grid search is on a coarse grid
% we want to minimize two the reflectance for two wavelengths

% if interpGridScalFactor is 10, then 9 rows will be interpolated to be 90
% rows, and 10 columns will be interpolated to be 100 columns

minVals = leastSquaresGridSearch(modisR, R, inputs);

[truth_estimate_table] = gatherTruthEstimateVals(modis, minVals, inputs, pixels2use); % containts truth ad estimates and difference


%% ----- Make Plots -----


% plot relfectance curves with lines of constant radius
% if there are more than 3 pixels, this function will plot three random
% pixels from the set
plotReflectanceCurves_singleBand(R,inputs,pixels2use, modis);

% plot reflectance contours where x and y are tau and r

plotReflectanceContours(R,inputs,pixels2use)

% Plot MODIS re values against my calculated re values

plot_effRadius_modis_estimates(truth_estimate_table, inputs)


plot_tau_modis_estimates(truth_estimate_table, inputs)

% Plot both re and tau on two panels within the same figure
plot_re_tau_modis_vs_myEstimates(truth_estimate_table)

% plot the ratio of my retireval to the MODIS retrieval of droplet size and
% optical thickness
plot_myRetrieval_vs_modisRetrieval_hist(truth_estimate_table)


% plot the reflectance for bands 1 and 7 and show the MODIS measurement
plot2ReflectanceFuncBands(modis,R, inputs, pixels2use, 'king')



% Plot and compare the MODIS measured reflectance and my estiamte of
% reflectance using LibRadTran
compare_my_reflectance_with_MODIS(modisR,R,inputs, pixels2use)
