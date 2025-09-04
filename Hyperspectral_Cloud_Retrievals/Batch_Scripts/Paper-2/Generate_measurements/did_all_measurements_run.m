%% Check to see which files haven't been run

%%

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




