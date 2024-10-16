%% Detecting Thermodynamic Phase using EMIT measurements

% By Andrew John Buggee

clear variables
%% Define the EMIT data file to use by defining the data folder

% -------------------------------------
% ------- PICK EMIT DATA SET  --------
% -------------------------------------

%emitDataFolder = '17_Jan_2024_coast/';

% 27 january has overlap with MODIS observations
emitDataFolder = '27_Jan_2024/';
% -------------------------------------


%% Load EMIT data and define folders 

[emitDataPath, folder2save] = define_EMIT_dataPath_and_saveFolders();

[emit,L1B_fileName] = retrieveEMIT_data([emitDataPath, emitDataFolder]);

%% Plot the EMIT data on a geo scatter plot

figure; geoscatter(emit.radiance.geo.lat(:), emit.radiance.geo.long(:), 10, reshape(emit.radiance.measurements(:,:,17),[],1),'.');
cb = colorbar;
set(get(cb, 'label'), 'string', 'Radiance $(\mu W/nm/cm^{2}/sr)$','Interpreter','latex', 'Fontsize',22)
set(gca, 'FontSize',25)
set(gca, 'FontWeight', 'bold')
set(gcf, 'Position', [0 0 800 800])
title('Radiance $500 \; nm$ band','Interpreter','latex', 'FontSize', 40)
hold on





%% Define the pixels to use


% find all pixels over ocean that haven't been masked out
% first find masked pixels
nan_mask = isnan(emit.radiance.measurements(:,:,1));
% next find pixels over ocean
coastal_res = 40;    % 1 is low resolution, 10 is decently high resolution
make_plot = 1;  %0 = no plot, 1 = plot
isOcean = land_or_ocean(double(emit.radiance.geo.lat(:)), double(emit.radiance.geo.long(:)),...
    coastal_res,make_plot);
% reshape ocean so it's the same size as combined_mask
isOcean = reshape(isOcean, size(emit.radiance.geo.lat,1), size(emit.radiance.geo.long,2));

% grab the pixels that are NOT nan and that are over ocean
combined_mask = logical(~nan_mask .* isOcean);

% find the linear indexes for each pixel that meets the above criteria
lin_idx = find(combined_mask);

% convert this to the row and column value

for nn = 1:numel(lin_idx)

    [r,c] = ind2sub(size(combined_mask), lin_idx(nn));

    pixels2use(nn).row = r;
    pixels2use(nn).col = c;

end

% Grab the pixel indices
pixels2use = grab_pixel_indices(pixels2use, [size(emit.radiance.measurements,1),...
    size(emit.radiance.measurements, 2)]);

%% Remove data that is not needed

emit = remove_unwanted_emit_data(emit, pixels2use);


%% Create an input structure that helps write the INP files

% this is a built-in function that is defined at the bottom of this script
inputs = create_emit_inputs_hyperspectral_top_bottom(emitDataFolder, folder2save, L1B_fileName, emit);
%inputs = create_emit_inputs_hyperspectral_top_middle(emitDataFolder, folder2save, L1B_fileName, emit);

% *** Check Inputs ***

%% Define the spectral response function of EMIT for the desired Bands

% create the spectral response functions
emit.spec_response = create_EMIT_specResponse(emit, inputs);


%% Define the solar source file name and read in the solar source data

% ********* IMPORTANT *************
% The source flux is integrated with the EMIT spectral response function

% define the source file using the input resolution
inputs = define_source_for_EMIT(inputs, emit);


%% Convert radiance measurements to TOA reflectance for the desired pixels

emit = convert_EMIT_radiance_2_reflectance(emit, inputs);


%% Compute the radiance measurement uncertainty 

emit.radiance.uncertainty = compute_EMIT_radiance_uncertainty(emit);


%% Compute the reflectance uncertainty

emit.reflectance.uncertainty = compute_EMIT_reflectance_uncertainty(emit, inputs);

%% Compute the spectral shape parameter (Knap et al. 2002)

wl_1640 = emit.radiance.wavelength>1638 & emit.radiance.wavelength<1642;
wl_1670 = emit.radiance.wavelength>1668 & emit.radiance.wavelength<1672;

S = 100* (emit.reflectance.value(wl_1670,:) - emit.reflectance.value(wl_1640,:))./...
            emit.reflectance.value(wl_1640, :);


%% Plot the spectral shape parameter on a geoscatter plot

figure; geoscatter(emit.radiance.geo.lat, emit.radiance.geo.long,...
    60, reshape(S,[],1),'.');
cb = colorbar;
set(get(cb, 'label'), 'string', 'Spectral Shape Parameter','Interpreter','latex', 'Fontsize',22)
set(gca, 'FontSize',25)
set(gca, 'FontWeight', 'bold')
set(gcf, 'Position', [0 0 1000 700])
title('EMIT Spectral Shape Parameter','Interpreter','latex', 'FontSize', 40)



%%








%% Load MODIS data over the same region

% Load modis data and create input structure


% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(whatComputer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------

    % Define the MODIS folder name

    % ----- November 9th at decimal time 0.611 (14:40) -----
    modisFolder = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/MODIS_Cloud_Retrieval/MODIS_data/2024_01_27/';


elseif strcmp(whatComputer,'andrewbuggee')==true


    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    % Define the MODIS folder name

end


[modis,L1B_fileName] = retrieveMODIS_data(modisFolder);

%% Plot MODIS data on a geoscatter plot

lim_buffer = 2;      % degrees

figure; geoscatter(modis.geo.lat(:), modis.geo.long(:), 10, reshape(modis.EV1km.reflectance(:,:,3),[],1),'.');
cb = colorbar;
set(get(cb, 'label'), 'string', 'Reflectance $(1/sr)$','Interpreter','latex', 'Fontsize',22)
set(gca, 'FontSize',25)
set(gca, 'FontWeight', 'bold')
set(gcf, 'Position', [0 0 1000 700])
title('Reflectance $469 \; nm$ band','Interpreter','latex', 'FontSize', 40)
hold on
geolimits([min(emit.radiance.geo.lat, [], 'all')-lim_buffer max(emit.radiance.geo.lat, [], 'all')+lim_buffer],...
    [min(emit.radiance.geo.long, [], 'all')-lim_buffer max(emit.radiance.geo.long, [], 'all')+lim_buffer])


%% find all MODIS pixels over ocean with a cloud

% find all pixels with an optical depth greater than 0.5
tau_mask = modis.cloud.optThickness17>0.5;

% next find pixels over ocean
coastal_res = 40;    % 1 is low resolution, 10 is decently high resolution
make_plot = 1;  %0 = no plot, 1 = plot
isOcean = land_or_ocean(double(modis.geo.lat(:)), double(modis.geo.long(:)),...
    coastal_res,make_plot);
% reshape ocean so it's the same size as combined_mask
isOcean = reshape(isOcean, size(modis.geo.lat,1), size(modis.geo.long,2));

% Combine these two masks
combined_mask_modis = logical(tau_mask .* isOcean);

%% Plot the MODIS phase retrieval for pixels with a cloud over ocean

% ----- MODIS Cloud Phase Index Values -----
% 0 -- cloud mask undetermined                                                       
% 1 -- clear sky                                                                     
% 2 -- liquid water cloud                                                            
% 3 -- ice cloud                                                                     
% 4 -- undetermined phase cloud (but retrieval is attempted as  liquid water)   

lim_buffer = 0;      % degrees

figure; geoscatter(modis.geo.lat(combined_mask_modis), modis.geo.long(combined_mask_modis),...
    100, reshape(modis.cloud.phase(combined_mask_modis),[],1),'.');
cb = colorbar;
set(get(cb, 'label'), 'string', 'Phase Index','Interpreter','latex', 'Fontsize',22)
set(gca, 'FontSize',25)
set(gca, 'FontWeight', 'bold')
set(gcf, 'Position', [0 0 1000 700])
title('MODIS Cloud Phase Index','Interpreter','latex', 'FontSize', 40)
hold on
geolimits([min(emit.radiance.geo.lat, [], 'all')-lim_buffer max(emit.radiance.geo.lat, [], 'all')+lim_buffer],...
    [min(emit.radiance.geo.long, [], 'all')-lim_buffer max(emit.radiance.geo.long, [], 'all')+lim_buffer])
