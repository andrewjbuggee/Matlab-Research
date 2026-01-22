%% Function to run all the mie calculations required to construct new precompute mie tables

% By Andrew John Buggee

%%

function precomputeMieTables_slurm(folder_path)

clear inputs


which_computer = whatComputer;

%% Start parallel pool

% *** Start parallel pool ***
% Is parpool running?
p = gcp('nocreate');
if isempty(p)==true

    % first read the local number of workers avilabile.
    p = parcluster('local');
    % start the cluster with the number of workers available
    if p.NumWorkers>64
        % Likely the amilan128c partition with 2.1 GB per core
        % Leave some cores for overhead
        parpool(p.NumWorkers - 8);

    elseif p.NumWorkers<=64 && p.NumWorkers>10

        % Leave a core for overhead
        parpool(p.NumWorkers -1);

    elseif p.NumWorkers<=10

        % Leave a core for overhead
        parpool(p.NumWorkers -1);


    end

end

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