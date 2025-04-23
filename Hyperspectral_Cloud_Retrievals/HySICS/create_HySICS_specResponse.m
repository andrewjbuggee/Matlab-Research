%% ----- read the spectral response function for the HySICS instrument -----

% Response functions are stored in the text file labeled
% 'ILS Map.txt'

% The center wavelengths are stored in the text file labeled
% 'Center Wavelength.txt'

% These files contain the spectral response functions from 407 to 1813 nm

% INPUTS
%   (1) band_number: the MODIS band number tells the code which relative
%   spectral response to read in.

%   (2) wavelength_resolution: this is the resolution of the wavelength
%   vector (nm). If it is different from the native resolution of the
%   files, the function will perform linear interpolation

% By Andrew John Buggee

%%

function spec_response = create_HySICS_specResponse(band_number, sourceFile_wavelength_resolution)


% Check inputs

if isnumeric(band_number)

else
    error('input must be a numeric entry')

end

if length(band_number)>460
    error('HySICS has 460 spectral bands. You requested more than 460.')
end

if sum(band_number>460)>0
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





%% Read in the correct response function


% define an empty cell aray
spec_response = cell(1, length(band_number));

% define a place holder cell aray
spec_response_temporary = cell(1, length(band_number));

for nn = 1:length(band_number)


    % Select this spectral response function
    data2keep = data{:, band_number(nn)};

    % only keep the non-zero values
    index_nonZero = find(data2keep);

    spec_response_temporary{nn}(:,1) = wavelength{band_number(nn)}(index_nonZero);
    spec_response_temporary{nn}(:,2) = data2keep(index_nonZero);


    % ------ Check the wavelength resolution desired -------

    if native_resolution~=sourceFile_wavelength_resolution
        % then we linear interpolate!

        % round the wavelength grid to the source file resolution
        if sourceFile_wavelength_resolution==1

            new_wavelength_grid = ceil(spec_response_temporary{nn}(1,1)):sourceFile_wavelength_resolution:...
                floor(spec_response_temporary{nn}(end,1));

        elseif sourceFile_wavelength_resolution==0.1

            new_wavelength_grid = round(spec_response_temporary{nn}(2,1),1):sourceFile_wavelength_resolution:...
                spec_response_temporary{nn}(end-1,1);

        else

            error([newline, 'I dont know how to work with the wavelength resolution', newline])
        end


        new_spec_response = interp1(spec_response_temporary{nn}(:,1), spec_response_temporary{nn}(:,2), new_wavelength_grid);


        spec_response{nn}(:,1) = new_wavelength_grid;
        spec_response{nn}(:,2) = new_spec_response;

    else

        % if the resolution of the spectral response functions are
        % equal to the resolution of the source file, keep the
        % temporary spectral response function

        spec_response{nn}(:,1) = spec_response_temporary{nn}(:,1);
        spec_response{nn}(:,2) = spec_response_temporary{nn}(:,2);
    end





end



end