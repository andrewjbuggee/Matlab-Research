%% Make slit function files for libRadTran

% -t type of slitfunction
%       1. triangular (default)
%       2. rectangular
%       3. Gaussian

% -f full width at half maximum, in nm

% -r spectral resolution, in nm

% -n number of fwhm (in nm) spanned by the slit function. Only applicable with Gaussian
% (type 3) slit function. Default value is 4.

% -h Print help message.

%%

function [filename] = run_makeSlitFunction(slitFunc_type, fwhm, spectral_res, n_fwhm, computer_name)


% --- Point to locaiton of make_slitfunction program ---

% To run uvspec in the command line we have to point to its full location.
% To do this we will check to see what computer we are using



if strcmp('anbu8374',computer_name)

    make_slitfunction_folderName = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/bin/';

elseif strcmp('andrewbuggee',computer_name)

    make_slitfunction_folderName = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/'...
        'Hyperspectral-Cloud-Droplet-Retrieval/',...
        'LibRadTran/libRadtran-2.0.4/bin/'];

elseif strcmp('curc', computer_name)
    % location of the mie program
    make_slitfunction_folderName = '/projects/anbu8374/software/libRadtran-2.0.5/bin/';

    % if running on the CURC supercomputer, you need to load the modules
    % everytime you run a job. That's because each time you run the
    % function 'system', a new unique terminal window is open. Each time a
    % new terminal window is open, the modules need to be loaded.
    cmnd_modules = ['ml purge', newline, 'ml gcc/11.2.0', newline,...
        'ml netcdf/4.8.1', newline, 'ml perl/5.36.0', newline, 'ml texlive/2021',...
        newline, 'export PATH=/projects/$USER/software/libRadtran-2.0.5/:$PATH',...         % add libRadtran to the path
        newline, 'export PATH=/projects/$USER/software/libRadtran-2.0.5/data/:$PATH',...    % add libRadtran data files to the path
        newline, 'export PATH=/projects/$USER/software/libRadtran-2.0.5/bin/:$PATH',...     % add uvspec location to the path
        newline, 'export PATH=', folderName_INP_OUT,':$PATH',...                            % add inp/out file locations to the path
        newline, 'export GSL_BIN=/projects/$USER/software/gsl-2.6/bin',...                  % define the binary file location for the GSL pacakge for libRadtran to find
        newline, 'export GSL_LIB=/projects/$USER/software/gsl-2.6/lib',...                  % define the library folder for the GSL package for libRadtran to fine
        newline, 'export GSL_INC=/projects/$USER/software/gsl-2.6/include',...              % define pointer for libRadtran to find the GSL pacakge
        newline, 'export LD_LIBRARY_PATH=$GSL_LIB:$LD_LIBRARY_PATH',...                     % define the install library path for libRadtran
        newline, 'export INSTALL_DIR=/projects/$USER/software/libRadtran-2.0.5',...         % define the install directory
        newline, 'export PATH=$GSL_BIN:$PATH'];                                             % add the GSL binary files location to the path

end


% using the function 'system' runs commans in the terminal window
cmnd1 = ['cd ', make_slitfunction_folderName];






cmnd2 = ['(',make_slitfunction_folderName,'make_slitfunction',...
    ' -f ',num2str(fwhm), ' -r ', num2str(spectral_res),...
    ' -n ', num2str(n_fwhm), ' -t ', num2str(slitFunc_type), ')'];
% a successful command will return a status of 0
% an unsuccessful command will return a status of 1



% run all commands in the terminal window
if strcmp('curc', computer_name)

    [status] = system([cmnd_modules, ' ; ', cmnd1, ' ; ', cmnd2]);
    %[status] = system([cmnd_modules, ' ; ', cmnd2]);

else
    [status] = system([cmnd1, ' ; ', cmnd2]);
end

if status ~= 0
    error(['Status returned value of ',num2str(status)])
end






end
