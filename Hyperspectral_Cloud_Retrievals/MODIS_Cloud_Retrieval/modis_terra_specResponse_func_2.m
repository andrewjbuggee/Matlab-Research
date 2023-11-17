%% ----- read the spectral response functions from the MODIS Characterization Support Team -----

% Response functions are stored in the rich text file (.rtf) labeled
% 'modis_spectral_response_func.rtf'

% This file contains the spectral response functions below 3 microns


% INPUTS
%   (1) band_number: the MODIS band number tells the code which relative
%   spectral response to read in. Only 1-7 are valid at the moment

%   (2) linearInterp: a flag that is either true or false. If true, the
%   relative spectral response function will use linear interpolation to
%   provide a set of weightings on a grid of monotonically increasing
%   wavelengths


% By Andrew John Buggee

%%

function spec_response = modis_terra_specResponse_func_2(band_number, linearInterp)

% Check inputs 

    if isnumeric(band_number)
        
    else
        error('input must be a numeric entry')
        
    end
    
    if band_number<1 || band_number>7
        error([newline,'I can only read first 7 MODIS spectral channels', newline])
    end

    

% Define the Folder Name

if strcmp(whatComputer, 'andrewbuggee')==true
    foldername = ['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/',...
        'MODIS_Cloud_Retrieval/MODIS_data/relative_spectral_response_terra/'];
elseif strcmp(whatComputer, 'anbu8374')==true
    error([newline, 'Define the folder for this compute!', newline])
end


% define the filename
filename = ['rsr.',num2str(band_number),'.inb.final'];



% open the file for reading
file_id = fopen([foldername,filename], 'r');   % 'r' tells the function to open the file for reading

format_spec = '%f %f %f %f';                                  % two floating point numbers
data = textscan(file_id, format_spec, 'Delimiter',' ',...
'MultipleDelimsAsOne',1, 'CommentStyle','#');

%% Read in the correct response function

% Store the wavelength data in a vector. Values associated with non-zero
% responses will be thrown out

spec_response.wavelength = data{3};       % nm

% Store the spectral response

spec_response.value = data{4};   % unitless

% There may be some spectral response values less than 0. Throw these away.

idx = spec_response.value<0;

spec_response.wavelength(idx) = [];
spec_response.value(idx) = [];



if linearInterp==true
    % The data is usually very long, and the wavelengths don't always increase
    % monotonically. Let's use linear interpolation to provide a smaller,
    % monotonic set of relative spectral responses


    % Extract the unique wavelength grid from the multiples
    [unique_wavelength_vector, idx_original, idx_new] = unique(spec_response.wavelength, 'first');
    % grab the spectral response values corresponding to the unique
    % wavelength vector above
    unique_spectral_response = spec_response.value(idx_original);

    % create a wavelength vector to interpolat with
    wavelength_vector_interp = ceil(unique_wavelength_vector(1)):1:floor(unique_wavelength_vector(end));
    
    % interpolate using the new wavelength vector
    interp_method = 'linear';
    new_spec_response = interp1(unique_wavelength_vector, unique_spectral_response, wavelength_vector_interp, interp_method);

    % rewrite the structure
    spec_response.wavelength = wavelength_vector_interp;
    spec_response.value = new_spec_response;

end





end


