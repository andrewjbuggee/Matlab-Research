%% ----- Aggregate all MODIS data -----


% by Andrew J. Buggee

%%


function [emit,L1B_fileNames] = retrieveEMIT_data(folderName)


% --- Add foldername to the matlab path ---
addpath(folderName)

% -----------------------------------------------------
% ----- Check to see if the folder name is valid ------
% -----------------------------------------------------


files = dir([folderName, '*.nc']);       % find all files that end in .hdf

% check to see if we found any files!
if isempty(files)==true
    error([newline,'There are no files in the folder provided!', newline])
end

L1B_fileNames = cell(1,length(files));

for ii = 1:length(files)

    L1B_fileNames{ii} = files(ii).name;



    if strcmp(L1B_fileNames{ii}(10:12), 'RAD')==true
    

        % read in the radiance data
        emit.radiance = readEMIT_L1B_radiance_data(L1B_fileNames{ii});

    elseif strcmp(L1B_fileNames{ii}(10:12), 'OBS')==true

        % read in the observation data
        emit.obs = readEMIT_L1B_observation_data(L1B_fileNames{ii});


    end


end


%%

% Record the time the data was taken
% Time recorded in [hours, minutes]
emit.time = [str2double(L1B_fileNames{ii}(27:28)), str2double(L1B_fileNames{ii}(29:30))];        % UTC time of data recording
emit.time_decimal = emit.time(1) + emit.time(2)/60;                        % UTC time in decimal format

% Convert the date to day of year
d = datetime(str2double(L1B_fileNames{1}(18:21)), str2double(L1B_fileNames{1}(22:23)),...
    str2double(L1B_fileNames{1}(24:25)));
emit.day_of_year = day(d, 'dayofyear');


end
