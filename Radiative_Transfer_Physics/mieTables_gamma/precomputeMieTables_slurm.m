%% Function to run all the mie calculations required to construct new precompute mie tables

% By Andrew John Buggee

%%

function precomputeMieTables_slurm(folder_path)


% Find all .INP files within the input defined folder path
inpFiles = dir(fullfile(folder_path, '*.INP'));

% step through each .INP file in a parallel for loop and run the .INP file
% using runMie

parfor nn = 1:length(inpFiles)

    % run the mie file
    runMIE(mie_folder_path,input_filename,output_filename, computer_name);

end


end