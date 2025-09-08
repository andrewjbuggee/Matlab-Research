% Clear old parallel pool job files

% By Andrew John Buggee

%%

function [] = clear_parPool_jobFiles(which_computer)


% Is parpool running?
p = gcp('nocreate');

% Clean up any existing parallel pool and its associated job files
if ~isempty(p)
    delete(p);
end

% Clean up any stray job files
localJobStorageLocation = fullfile(getenv('HOME'), '.matlab', 'local_cluster_jobs');
if exist(localJobStorageLocation, 'dir')==7
    rmdir(localJobStorageLocation, 's');
end



    % Clean up any stale job files before starting (especially for CURC)
    if strcmp(which_computer,'curc')==true
        try
            % Clear parallel preferences to remove stale job references
            parallel.Settings.clearAll();

            % Optional: Clean up job files directory
            job_dir = fullfile(prefdir, 'local_cluster_jobs');
            if exist(job_dir, 'dir')
                fprintf('Cleaning up old job files...\n');
                rmdir(job_dir, 's');
            end
        catch ME
            % If cleanup fails, just warn and continue
            warning([newline, 'Could not clean up old job files: %s', ME.message, newline]);
        end
    end


    



end