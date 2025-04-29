%% ----- read the spectral response function for the HySICS instrument -----

% Response functions are stored in the text file labeled
% 'ILS Map.txt'

% The center wavelengths are stored in the text file labeled
% 'Center Wavelength.txt'

% These files contain the spectral response functions from 407 to 1813 nm

% INPUTS
%   (1) band_number: the MODIS band number tells the code which relative
%   spectral response to read in.

%   (2) sourceFile - name of uvspec solar flux file (string) - there are
%   three options to choose from:
%       (a) 'kurudz_1.0nm.dat' - These data range from 250 - 10000 nm and
%       are spaced by 1 nm. These data were taken from LBLRTM 5.21
%       (http://www.meto.umd.edu/~bobe/LBLRTM/). The original Kurudz [1992]
%       data were converted to mW / (m2 nm) and averaged over 1nm intervals
%       centered around the given wavelength.

%       (b) 'kurudz_0.1nm.dat' - These data range from 250 - 10000 nm and
%       are spaced by 0.1 nm. These data were taken from LBLRTM 5.21
%       (http://www.meto.umd.edu/~bobe/LBLRTM/). The original Kurudz [1992]
%       data were converted to mW / (m2 nm) and averaged over 0.1nm intervals
%       centered around the given wavelength.

%       (c) 'atlas_plus_modtran' - Use this if you wish to have solar
%       flux values from 200 - 250 nm. These data range from 200 to 800 nm
%       with 0.05 nm resolution. You do not need an extension if using this
%       file!

%       (d)
%       'hybrid_reference_spectrum_p025nm_resolution_c2022-11-30_with_unc.dat'
%       - this is a hybrind reference spectrum downloaded from LASP's
%       LISIRD tool (https://lasp.colorado.edu/lisird/data/tsis1_hsrs_p1nm)
%       These data range from 202 to 2730 nm
%       These data have a smapling resolution of 0.005 nm.  The original
%       data were converted to mW / (m2 nm)

%       (e)
%       'hybrid_reference_spectrum_p1nm_resolution_c2022-11-30_with_unc.dat'
%       - this is a hybrind reference spectrum downloaded from LASP's
%       LISIRD tool (https://lasp.colorado.edu/lisird/data/tsis1_hsrs_p1nm)
%       These data range from 202 to 2730 nm
%       These data have a smapling resolution of 0.025 nm.  The original
%       data were converted to mW / (m2 nm)


%       (f)
%       'hybrid_reference_spectrum_1nm_resolution_c2022-11-30_with_unc.dat'
%       - this is a hybrind reference spectrum downloaded from LASP's
%       LISIRD tool (https://lasp.colorado.edu/lisird/data/tsis1_hsrs_p1nm)
%       These data range from 202 to 2730 nm
%       These data have a smapling resolution of 0.1 nm.  The original
%       data were converted to mW / (m2 nm)



% By Andrew John Buggee

%%

function spec_response = create_HySICS_specResponse(band_number, sourceFile)

if nargin ~=2

    error([newline, 'Need two inputs: the band numbers, and the source file.', newline])

end

% Check inputs

if isnumeric(band_number)

else
    error('input must be a numeric entry')

end

if length(band_number)>636
    error('HySICS has 636 spectral bands. You requested more than 460.')
end

if sum(band_number>636)>0
    error('You requested a band number that is higher than 460, which doesnt exist')
end



% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(whatComputer,'anbu8374')==true

    % ------ Folders on my Mac Desktop --------

    folderPath = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/';



elseif strcmp(whatComputer,'andrewbuggee')==true

    % ------ Folders on my Macbook --------

    folderPath = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'HySICS/'];



elseif strcmp(whatComputer,'curc')==true

    % ------ Folders on the CU Supercomputer /projects folder --------

    error([newline, 'Not on the CU supercomputer yet'])


end




% define the filenames
filename_specResponse = 'ILS Map.txt';
filename_centerWavelength = 'Center Wavelength.txt';


% ------------------------------------------------------
% -------- Reading .txt file using textscan ------------
% ------------------------------------------------------

% ---- Open the center wavelength File ----
file_id = fopen([folderPath,filename_centerWavelength], 'r');   % 'r' tells the function to open the file for reading
% Use the textscan() function instead, which allows us to define comments to ignore
format_spec = '%f';                                  % two floating point numbers
wl_center = textscan(file_id, format_spec, 'Delimiter',',',...
    'MultipleDelimsAsOne',1, 'CommentStyle',';');
wl_center = wl_center{1};





% ---- Open the spectral resposne function File ----
file_id = fopen([folderPath,filename_specResponse], 'r');   % 'r' tells the function to open the file for reading

% determine the file formatting properties
file_prop = detectImportOptions([folderPath,filename_specResponse]);

% define the columns according to the band number and define the wavelength
% range for each spectral response function
native_resolution = 0.1;     % nm

wavelength = cell(1, length(wl_center));

for nn = 1:length(wl_center)

    var_names(nn) = string({['Band ', num2str(nn)]});

    % The displacement steps are in 0.1 nm, so span a range from -20 to +20 nm around the center wavelength
    wavelength{nn} = (-20:0.1:20) + wl_center(nn);    % nm

end

% read in table
data = readtable(filename_specResponse, file_prop);

fclose(file_id);

% reset the variable names
data.Properties.VariableNames = var_names;






% ------------------------------------------
% ----------- IMPORTANT ADDITION -----------
% ------------------------------------------

% The above files only include center wavelengths and the associated
% spectral response fucntions from 407 - 1813 nm, but HySICS has a spectral
% range from 350 to 2300 nm. Using the file above, fit a spline the FWHM's
% in the range provided to extend this to the full spectral range. The
% spectral response functions are then defined using the FWHM and the
% center wavelengths

% Use the mean difference to find the center wavelengths for the remaining
% spectral range
mean_wl_diff = mean(diff(wl_center));        % nm
full_center_wl = [fliplr(wl_center(1)-mean_wl_diff:-mean_wl_diff:350)'; wl_center;...
    (wl_center(end)+mean_wl_diff:mean_wl_diff:2300)'];


% --------------------------------------------------
% ------------- Read in the source file ------------
% --------------------------------------------------

% read in the desired source file wavelength grid
[~, source_wavelength] = read_solar_flux_file([full_center_wl(1)-25, full_center_wl(end)+25], sourceFile);  % nm

sourceFile_wavelength_resolution = source_wavelength(2) - source_wavelength(1);              % nm
% define the wavelength grid of the spectral response on the same
% wavelength grid as the source file
source_wavelength_grid = -20:sourceFile_wavelength_resolution:20;     % nm


% --------------------------------------------------
% ---- Compute the FWHM from the file provided -----
% --------------------------------------------------

% For a normal distribution, FWHM = 2*sqrt(2*log(2))*std
fwhm = zeros(length(wl_center), 1);

parfor nn = 1:length(wl_center)

    f = fit(wavelength{nn}', data.(string(['Band ', num2str(nn)])), 'gauss1');

    std = f.c1/sqrt(2);

    fwhm(nn) = 2*sqrt(2*log(2))*std;


end

% Now let's extrapolate the full spectral range of HySICS
fwhm_extrapolated = interp1(wl_center, fwhm, full_center_wl, "pchip", "extrap");

% compute the spectral response on a fine grid
specResponse_wavelength_grid = -21:(sourceFile_wavelength_resolution/20):21;  % nm




spec_response.value = zeros(length(band_number), length(source_wavelength_grid));
spec_response.wavelength = zeros(length(band_number), length(source_wavelength_grid));

% step through each spectral band for which a spectral response function is
% to be computed
for nn = 1:length(band_number)

    % In libRadtran, everything is computed on the source wavelength grid.
    % So if the source is defined at integer values of wavelength (100,
    % 101, ...), then the computations must be done at those wavelengths

    % So, let's define a continuous spectral response function and define
    % the output for values at the wavelength values of the source grid. 

    % compute the standard deviation from the FWHM
    sigma = fwhm_extrapolated(band_number(nn))/(2*sqrt(2*log(2)));      % std
        
    % define the high-resolution grid to compute the spectral response
    % function
    wl_grid_hiRes = specResponse_wavelength_grid + full_center_wl(band_number(nn));

    % Find source wavelength closest to each center wavelength of the
    % HySICS spectrometer
    [~, center_wavelength_idx] = min(abs(source_wavelength - full_center_wl(band_number(nn))));

    % define the wavelength grid for each spectral channel
    % This should be defined as part of the source wavelength grid
    % The displacement should span a range from -20 to +20 nm around the center wavelength
    % the center wavelength of the output grid is NOT aligned with the peak
    % of the spectral response function. It's the value closest to the
    % center wavelength on the source wavelength grid
    spec_response.wavelength(nn, :) = source_wavelength_grid + source_wavelength(center_wavelength_idx);

    % compute the gaussian spectral response function
    response_hiRes = pdf('Normal', wl_grid_hiRes, full_center_wl(band_number(nn)), sigma);

    % interpolate on the lower resolution source grid to get the values of
    % the spectral response function
    spec_response.value(nn, :) = interp1(wl_grid_hiRes, response_hiRes, spec_response.wavelength(nn, :));


end



%% Read in the correct response function

% ****** For when I have the response functions across to the full
% spectrum ******



% % define an empty cell aray
% spec_response = cell(1, length(band_number));
% 
% % define a place holder cell aray
% spec_response_temporary = cell(1, length(band_number));
% 
% for nn = 1:length(band_number)
% 
% 
%     % Select this spectral response function
%     data2keep = data{:, band_number(nn)};
% 
%     % only keep the non-zero values
%     index_nonZero = find(data2keep);
% 
%     spec_response_temporary{nn}(:,1) = wavelength{band_number(nn)}(index_nonZero);
%     spec_response_temporary{nn}(:,2) = data2keep(index_nonZero);
% 
% 
%     % ------ Check the wavelength resolution desired -------
% 
%     if native_resolution~=sourceFile_wavelength_resolution
%         % then we linear interpolate!
% 
%         % round the wavelength grid to the source file resolution
%         if sourceFile_wavelength_resolution==1
% 
%             specResponse_wavelength_grid = ceil(spec_response_temporary{nn}(1,1)):sourceFile_wavelength_resolution:...
%                 floor(spec_response_temporary{nn}(end,1));
% 
%         elseif sourceFile_wavelength_resolution==0.1
% 
%             specResponse_wavelength_grid = round(spec_response_temporary{nn}(2,1),1):sourceFile_wavelength_resolution:...
%                 spec_response_temporary{nn}(end-1,1);
% 
%         else
% 
%             error([newline, 'I dont know how to work with the wavelength resolution', newline])
%         end
% 
% 
%         new_spec_response = interp1(spec_response_temporary{nn}(:,1), spec_response_temporary{nn}(:,2), specResponse_wavelength_grid);
% 
% 
%         spec_response{nn}(:,1) = specResponse_wavelength_grid;
%         spec_response{nn}(:,2) = new_spec_response;
% 
%     else
% 
%         % if the resolution of the spectral response functions are
%         % equal to the resolution of the source file, keep the
%         % temporary spectral response function
% 
%         spec_response{nn}(:,1) = spec_response_temporary{nn}(:,1);
%         spec_response{nn}(:,2) = spec_response_temporary{nn}(:,2);
%     end
% 
% 
% 
% 
% 
% end



end