%% ----- Aggregate all MODIS data -----


% by Andrew J. Buggee

%%


function [emit,L1B_fileNames] = retrieveEMIT_data(folderName)


% --- Add foldername to the matlab path ---
addpath(folderName)

% -----------------------------------------------------
% ----- Check to see if the folder name is valid ------
% -----------------------------------------------------


files = dir([folderName,'*.nc']); % find all files that end in .hdf

% check to see if we found any files!
if isempty(files)==true
    error([newline,'There are not files in the folder provided!', newline])
end

L1B_fileNames = cell(1,length(files));

for ii = 1:length(files)

    file_ii = files(ii).name;



    if strcmp(file_ii(10:12), 'RAD')==true

        % read in the radiance data
        emit = readEMIT_L1B_radiance_data(file_ii);

    elseif strcmp(file_ii(10:12), 'OBS')==true

        % read in the observation data
        warning([newline, 'I dont know how to read observation data', newline])


    end


end

% Record the time the data was taken
% Time recorded in [hours, minutes]
emit.time = [str2double(file_ii(27:28)), str2double(file_ii(29:30))];        % UTC time of data recording
emit.time_decimal = emit.time(1) + emit.time(2)/60;                        % UTC time in decimal format



end
