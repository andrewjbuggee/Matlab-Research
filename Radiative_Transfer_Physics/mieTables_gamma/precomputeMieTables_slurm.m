%% Function to run all the mie calculations required to construct new precompute mie tables

% By Andrew John Buggee

%%

function precomputeMieTables_slurm(folder_path)

clear inputs


which_computer = whatComputer;

%% Start parallel pool

start_parallel_pool(which_computer)


%% Find all INP files and run them in a parallel for loop!


% Find all .INP files within the input defined folder path
inpFiles = dir(fullfile(folder_path, '*.INP'));

% step through each .INP file in a parallel for loop and run the .INP file
% using runMie
outputFileName = cell(size(inpFiles));

tic
parfor nn = 1:length(inpFiles)

    % Can i set the output name? I don't think it's used when output_user
    % netcdf is enabled
    outputFileName{nn} = ['OUTPUT_', inpFiles(nn).name(1:end-4)];


    % run the mie file
    runMIE(folder_path, inpFiles(nn).name, outputFileName{nn}, which_computer);

end

disp([newline, 'Total time to run ', num2str(length(inpFiles)), ' files was: ', toc, ' seconds.', newline])

end