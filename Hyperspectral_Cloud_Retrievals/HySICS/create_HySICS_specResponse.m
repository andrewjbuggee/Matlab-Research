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

function spec_response = create_HySICS_specResponse(band_number, wavelength_resolution)


% Check inputs

if isnumeric(band_number)

else
    error('input must be a numeric entry')

end

if length(band_number)>36
    error('MODIS has 36 spectral bands. You requested more than 36. You may pull any of of them')
end

if sum(band_number>36)>0
    error('You requested a band number that is higher than 36, which doesnt exist')
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




% The displacement steps are in 0.1 nm, so span a range from -20 to +20 nm around the center wavelength
sr_length = length(-20:0.1:20);

% ---- Open the center wavelength File ----
file_id = fopen([folderPath,filename_specResponse], 'r');   % 'r' tells the function to open the file for reading
% Use the textscan() function instead, which allows us to define comments to ignore
format_spec = '%f';                                  % two floating point numbers
specResponse = textscan(file_id, format_spec, 'Delimiter',{',','\n'},...
    'Whitespace', '\b\t', 'LineEnding', {'\n'  '\r'  '\r\n'}, 'CommentStyle',';');

fclose(file_id)

% determine the file formatting properties
file_prop = detectImportOptions([folderPath,filename_specResponse]);

% define the columns according to the MODIS band number
var_names = string({'Wavelength (nm)', 'Band 8', 'Band 9', 'Band 3',...
    'Band 10', 'Band 11', 'Band 12', 'Band 4', 'Band 1',...
    'Band 13', 'Band 14', 'Band 15', 'Band 2', 'Band 16', ...
    'Band 5', 'Band 6', 'Band 7'}); %#ok<STRCLQT>

% read in table
data = readtable(filename_specResponse, file_prop);

% reset the variable names
data.Properties.VariableNames = var_names;


% The last column needs to be fixed. It's read in as a cell vector of
% string characters



% % open the file for reading
% file_id = fopen(filename, 'r');   % 'r' tells the function to open the file for reading
%
% format_spec = '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f';                                  % two floating point numbers
% data = textscan(file_id, format_spec, 'Delimiter',' ',...
% 'MultipleDelimsAsOne',1, 'CommentStyle','/', 'EndOfLine','\');

%% Read in the correct response function

% Store the wavelength data in a vector. Values associated with non-zero
% responses will be thrown out

wavelength = data.("Wavelength (nm)");

% define an empty cell aray
spec_response = cell(1, length(band_number));

% define a place holder cell aray
spec_response_temporary = cell(1, length(band_number));

for nn = 1:length(band_number)

    for bb = 1:(length(data.Properties.VariableNames)-1)

        if strcmp(['Band ',num2str(band_number(nn))], data.Properties.VariableNames{bb+1})==true

            % Select this spectral response function
            data2keep = data{:,bb+1};

            % only keep the non-zero values
            index_nonZero = find(data2keep);

            spec_response_temporary{nn}(:,1) = wavelength(index_nonZero);
            spec_response_temporary{nn}(:,2) = data2keep(index_nonZero);


            % ------ Check the wavelength resolution desired -------
            native_resolution = spec_response_temporary{nn}(2,1) - spec_response_temporary{nn}(1,1);

            if native_resolution~=wavelength_resolution
                % then we linear interpolate!
                new_wavelength = spec_response_temporary{nn}(1,1):wavelength_resolution:spec_response_temporary{nn}(end,1);
                new_spec_response = interp1(spec_response_temporary{nn}(:,1), spec_response_temporary{nn}(:,2), new_wavelength);


                spec_response{nn}(:,1) = new_wavelength;
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

end



end