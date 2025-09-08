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


% *** Start parallel pool ***

if isempty(p)==true


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







    end


end






