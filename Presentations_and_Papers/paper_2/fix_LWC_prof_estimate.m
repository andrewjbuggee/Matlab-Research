%% Fix the LWC profile estimate for each GN_inputs structure

clear variables

% Determine which computer you're using
which_computer = whatComputer();

% Load simulated measurements
if strcmp(which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------



elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------


    % define the folder where the Vocals-Rex in-situ derived measurements
    % are located
    folder_paths_1.retrieval = ['/Users/andrewbuggee/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
        'paper2_variableSweep/full_retrieval_logSpace_newCov_VR_meas_allBands_3/'];

    folder_paths_1.ensemble_VR_profs = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/',...
        'ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25_drizzleLWP-threshold_5_10-Nov-2025.mat'];


    % *** updated ensemble profiles file with corrected LWP ** (values were
    % negative but correct)
    % folder_paths_1.ensemble_VR_profs = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
    %     'Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/',...
    %     'ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25_drizzleLWP-threshold_5_04-Dec-2025.mat'];



elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------


end


% First step through all files in the directory and remove invisble files
% or directories

% Grab filenames in drive
filenames_retrieval = dir(folder_paths_1.retrieval);
idx_2delete = [];
for nn = 1:length(filenames_retrieval)

    if strcmp(filenames_retrieval(nn).name(1), 'd')~=true

        idx_2delete = [idx_2delete, nn];

    end

end


% delete rows that don't have retrieval filenames
filenames_retrieval(idx_2delete) = [];


% load the ensemble profiles
load(folder_paths_1.ensemble_VR_profs)

for nn = 1:length(filenames_retrieval)

    % Load the data from the file
    ds = load([filenames_retrieval(nn).folder, '/',filenames_retrieval(nn).name]);

    % first, plot the simulated profile values as two lines
    if isfield(ds, 'GN_inputs')==true

        % find the ensemble profile that was used to generate the
        % measurement
        VR_date = extractBetween(filenames_retrieval(nn).name, 'recorded_', '-2008');
        % append 2008
        VR_date{1} = [VR_date{1}, '-2008'];
        % find the time the measurement was recorded
        VR_time = extractBetween(filenames_retrieval(nn).name, '2008_', 'UTC');

        % Find the ensemble profile with the same date
        for pp = 1:length(ensemble_profiles)

            if strcmp(VR_date{1}, string(ensemble_profiles{pp}.dateOfFlight))==true

                % check it the times match
                ensem_times = ensemble_profiles{pp}.time_utc;

                if str2double(VR_time{1}) >= min(ensem_times) && str2double(VR_time{1}) <= max(ensem_times)

                    % We found a match!
                    % add the tau vector to the GN_inputs and resave the
                    % file
                    % determine direction of flight
                    if (ensemble_profiles{pp}.altitude(2) - ensemble_profiles{pp}.altitude(1)) < 0
                        % the plane is descending
                        ds.GN_inputs.measurement.lwp = (-1) * trapz(ensemble_profiles{pp}.altitude, ensemble_profiles{pp}.lwc);

                    else

                        % the plane is ascending
                        ds.GN_inputs.measurement.lwp = trapz(ensemble_profiles{pp}.altitude, ensemble_profiles{pp}.lwc);

                    end
                    
                    % open all variables and save them
                    inputs_tblut = ds.inputs_tblut;
                    Refl_model_tblut = ds.Refl_model_tblut;
                    tblut_retrieval = ds.tblut_retrieval;
                    acpw_retrieval = ds.acpw_retrieval;
                    inputs_acpw = ds.inputs_acpw;
                    Refl_model_acpw = ds.Refl_model_acpw;
                    folder_paths = ds.folder_paths;
                    GN_inputs = ds.GN_inputs;
                    GN_outputs = ds.GN_outputs;

                    save([folder_paths_1.retrieval, filenames_retrieval(nn).name], "inputs_tblut", "Refl_model_tblut", ...
                        "tblut_retrieval", "acpw_retrieval", "inputs_acpw", "Refl_model_acpw", "folder_paths",...
                        "GN_inputs", "GN_outputs");

                    % clear variables
                    clear("inputs_tblut", "Refl_model_tblut", ...
                        "tblut_retrieval", "acpw_retrieval", "inputs_acpw", "Refl_model_acpw", "folder_paths",...
                        "GN_inputs", "GN_outputs");


                end





            end




        end



    end


end

