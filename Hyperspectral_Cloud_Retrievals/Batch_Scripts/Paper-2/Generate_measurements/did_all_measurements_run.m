%% Check to see which files haven't been run

%% For the full variable space

clear variables

meas_path = ['/Users/andrewbuggee/MATLAB-Drive/HySICS/Simulated_spectra/',...
    'paper2_variableSweep/rTop_10/vza_7_vaz_210_sza_10_saz_91/'];

retrieval_path = ['/Users/andrewbuggee/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
    'paper2_variableSweep/rTop_10/vza_7_vaz_210_sza_10_saz_91/'];

meas_path_remaining = ['/Users/andrewbuggee/MATLAB-Drive/HySICS/Simulated_spectra/',...
    'paper2_variableSweep/rTop_10/vza_7_vaz_210_sza_10_saz_91_remaining/'];

% meas_files = dir(['/Users/andrewbuggee/MATLAB-Drive/HySICS/Simulated_spectra/',...
%     'paper2_variableSweep/rTop_10/vza_7_vaz_210_sza_10_saz_91/*.mat']);

% step through each file name and check to see if there is a corresponding
% droplet retrieval file. If not, store the filename

r_top = 10;
r_bot = 3:10;
tau_c = 5:3:29;
tcpw = 5:3:35;

for rb = 1:numel(r_bot)
    for tt = 1:numel(tau_c)
        for wv = 1:numel(tcpw)


            searchPattern = ['*', 'rTop_10_rBot_', num2str(r_bot(rb)), '_tauC_',...
                num2str(tau_c(tt)),'_tcwv_',num2str(tcpw(wv)), '*'];

            % Use dir to find files matching the pattern in the current directory
            matchingFiles = dir([retrieval_path, searchPattern]);

            % if empty, copy file from original folder to new folder
            if isempty(matchingFiles)

                source_file = dir([meas_path, searchPattern]);
                copyfile([source_file.folder, '/', source_file.name], meas_path_remaining);

            end



        end


    end

end



%% For the subset variable space

clear variables


% Access specific file or folder
meas_path = ['/Users/anbu8374/MATLAB-Drive/HySICS/Simulated_spectra/paper2_variableSweep/',...
    'rTop_10/vza_7_vaz_210_sza_10_saz_91_subset/'];

retrieval_path = ['/Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
    'paper2_variableSweep/rTop_10/vza_7_vaz_210_sza_10_saz_91_subset/'];





% step through each file name and check to see if there is a corresponding
% droplet retrieval file for four difference uncertainties: 0.001%, 0.1%,
% 1%, and 5%. If not, store the file name and percentage.

r_top = 10;
r_bot = 5;
tau_c = [5,11,17,23];
tcpw = [8, 14, 20];

meas_uncert = [0.001, 0.1, 1, 5];
tcpw_assumed = [10, 15, 20, 25];


no_retrieval = {};
missing_files = {};

for tt = 1:numel(tau_c)
    for wv = 1:numel(tcpw)
        for uu = 1:numel(meas_uncert)


            searchPattern = ['*', num2str(meas_uncert(uu)), '%_uncert_rTop_10_rBot_5_tauC_',...
                num2str(tau_c(tt)),'_tcwv_',num2str(tcpw(wv)), '*'];

            % Use dir to find files matching the pattern in the current directory
            matchingFiles = dir([retrieval_path, searchPattern]);


            % there should be 5 files, 1 for the full retrieval, and 4 for
            % the retrieval without column water vapor
            if meas_uncert(uu)==0.1 || meas_uncert(uu)==1 || meas_uncert(uu)==5

                if length(matchingFiles)==4

                    % Check to see if each one has a retrieval
                    for nn = 1:length(matchingFiles)

                        ds = load([matchingFiles(nn).folder, '/', matchingFiles(nn).name]);

                        if isfield(ds, 'GN_outputs')~=1

                            no_retrieval = [no_retrieval; matchingFiles(nn).name];

                        elseif isfield(ds.GN_outputs, 'retrieval')~=1

                            no_retrieval = [no_retrieval; matchingFiles(nn).name];


                        end

                    end

                elseif length(matchingFiles)<=4

                    % check to see which file is missing
                    missing_files = [missing_files, matchingFiles(:).name];

                end


                % There should be 5 files here
            elseif meas_uncert(uu)==0.001

                if length(matchingFiles)==5

                    % Check to see if each one has a retrieval
                    for nn = 1:length(matchingFiles)

                        ds = load([matchingFiles(nn).folder, '/', matchingFiles(nn).name]);

                        if isfield(ds, 'GN_outputs')~=1

                            no_retrieval = [no_retrieval; matchingFiles(nn).name];

                        elseif isfield(ds.GN_outputs, 'retrieval')~=1

                            no_retrieval = [no_retrieval; matchingFiles(nn).name];


                        end

                    end

                elseif length(matchingFiles)<4

                    % check to see which file is missing
                    missing_files = [missing_files, matchingFiles(:).name];

                end

            end







        end


    end


end






