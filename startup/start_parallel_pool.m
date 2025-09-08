%% Start parallel pool

% By Andrew John Buggee

%%

function [] = start_parallel_pool(which_computer)


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

% *** Start parallel pool ***

if isempty(p)==true


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




    % first read the local number of workers avilabile.
    p = parcluster('local');


    % Find the folder where the mie calculations are stored
    % find the folder where the water cloud files are stored.
    if strcmp(which_computer,'anbu8374')==true

        % -----------------------------------------
        % -------------- Mac Desktop --------------
        % -----------------------------------------


        % *** Start parallel pool ***


        % start the cluster with the number of workers available
        if p.NumWorkers>64
            % Likely the amilan128c partition with 2.1 GB per core
            % Leave some cores for overhead
            parpool(p.NumWorkers - 8);

        elseif p.NumWorkers<=64 && p.NumWorkers>10

            parpool(p.NumWorkers);

        elseif p.NumWorkers<=10

            parpool(p.NumWorkers);

        end





    elseif strcmp(which_computer,'andrewbuggee')==true



        % --------------------------------------
        % --------------- Macbook --------------
        % --------------------------------------


        % *** Start parallel pool ***


        % start the cluster with the number of workers available
        if p.NumWorkers>64
            % Likely the amilan128c partition with 2.1 GB per core
            % Leave some cores for overhead
            parpool(p.NumWorkers - 8);

        elseif p.NumWorkers<=64 && p.NumWorkers>10

            parpool(p.NumWorkers);

        elseif p.NumWorkers<=10

            parpool(p.NumWorkers);

        end







    elseif strcmp(which_computer,'curc')==true


        % ----------------------------------------------
        % --------------  CU Super Computer ------------
        % ----------------------------------------------



        % Add a small delay to let file system sync (cluster-specific)
        pause(1);


        % *** Start parallel pool ***


        try
            % start the cluster with the number of workers available
            if p.NumWorkers>64
                % Likely the amilan128c partition with 2.1 GB per core
                % Leave some cores for overhead
                fprintf('Starting parallel pool with %d workers...\n', p.NumWorkers - 8);
                parpool(p.NumWorkers - 8);

            elseif p.NumWorkers<=64 && p.NumWorkers>10
                fprintf('Starting parallel pool with %d workers...\n', p.NumWorkers);
                parpool(p.NumWorkers);

            elseif p.NumWorkers<=10
                fprintf('Starting parallel pool with %d workers...\n', p.NumWorkers);
                parpool(p.NumWorkers);

            end

        catch ME
            % If parpool fails, try once more after cleanup
            warning([newline,'First parpool attempt failed: %s', ME.message, newline]);
            fprintf('Attempting cleanup and retry...\n');

            % Force cleanup - use alternative methods
            delete(gcp('nocreate'));  % Delete any existing pool
            
            % Clear job files directory
            job_dir = fullfile(prefdir, 'local_cluster_jobs');
            if exist(job_dir, 'dir')
                rmdir(job_dir, 's');
            end
            pause(2);
            
            % Retry parpool
            if p.NumWorkers>64
                parpool(p.NumWorkers - 8);
            elseif p.NumWorkers<=64 && p.NumWorkers>10
                parpool(p.NumWorkers);
            elseif p.NumWorkers<=10
                parpool(p.NumWorkers);
            end
        end





    end


end






