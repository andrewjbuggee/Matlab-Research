%% --- Run Intput Files with UVSPEC ---
% ======================================
% The purpose of this script is to run uvspec with the input files specfied
% by the user. The user will be required to provide to folder location of
% the input file, and the input file name. The user will also need to
% provide the output file name. This output file will be saved in the
% same folder as the input file. The input file will be fed into the
% command line in order to run uvspec.



% --- By Andrew J. Buggee ---
%% Creating the .INP file

function runUVSPEC_ver2(folderName_INP_OUT,inputName,outputName, computer_name)
%% ---- A Few Checks are Needed ----

if iscell(inputName)==true && iscell(outputName)==false
    error('inputName is a cell array, while outputName is not')

elseif iscell(inputName)==false && iscell(outputName)==true
    error('outputName is a cell array, while inputName is not')

elseif iscell(inputName)==true && iscell(outputName)==true
    if length(inputName)~=length(outputName)
        error('The number of input files doesns equal the number of output files')
    end
end

% -------------------------------------------




%% Running the .INP file from the command line


% First we need to determine how many files we need to run
if iscell(inputName)==true
    numFiles2Run = length(inputName);
elseif ischar(inputName)==true
    numFiles2Run = 1;
else
    error('I dont understand the input file')
end



% to run a .INP file from the command line, we run uvspec by pointing the
% command line to its full location. AND the command line directory must be
% in the folder where the file you are running sits. For example, if the
% file you wish to run lives in the folder: directory/folder/runMe.file
% and the program runnning the file is in the folder:
% directory/program/fancyProgram then the code in the command line to run
% this file is:
%
% cd directory/folder/
% directory/program/fancyProgram < runMe.file > output.OUT





% --- Now we Can Run the Files ----


% --- Point to locaiton of uvspec program ---

% To run uvspec in the command line we have to point to its full location.
% To do this we will check to see what computer we are using



if strcmp('anbu8374',computer_name)

    uvspec_folderName = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/bin/';

elseif strcmp('andrewbuggee',computer_name)

    uvspec_folderName = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/'...
        'Hyperspectral-Cloud-Droplet-Retrieval/',...
        'LibRadTran/libRadtran-2.0.4/bin/'];

elseif strcmp('curc', computer_name)
    % location of the mie program
    uvspec_folderName = '/projects/anbu8374/software/libRadtran-2.0.5/bin/';

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
cmnd1 = ['cd ', uvspec_folderName];



if numFiles2Run==1

    if ischar(inputName)==true
        % cmnd2 = [uvspec_folderName,'uvspec ',...
        %            '< ',inputName,' > ', outputName];

        cmnd2 = ['(',uvspec_folderName,'uvspec ',...
            '< ',folderName_INP_OUT,inputName,' > ', folderName_INP_OUT, outputName,'.OUT',...
            ')>& ', folderName_INP_OUT,'errMsg.txt'];
        % a successful command will return a status of 0
        % an unsuccessful command will return a status of 1

    elseif iscell(inputName)==true

        cmnd2 = ['(',uvspec_folderName,'uvspec ',...
            '< ', folderName_INP_OUT, inputName{1},' > ', folderName_INP_OUT, outputName{1},'.OUT',...
            ')>&', folderName_INP_OUT, 'errMsg.txt'];

    else

        error([newline, 'I dont understand the INP filename structure', newline])

    end


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



elseif numFiles2Run>1


    for ii = 1:numFiles2Run


        % cmnd2 = [uvspec_folderName,'uvspec ',...
        %            '< ',inputName,' > ', outputName];

        cmnd2 = ['(',uvspec_folderName,'uvspec ',...
            '< ', folderName_INP_OUT, inputName{ii},' > ', folderName_INP_OUT,...
            outputName{ii},'.OUT',')>& errMsg.txt'];
        % a successful command will return a status of 0
        % an unsuccessful command will return a status of 1

        if strcmp('curc', computer_name)

            [status] = system([cmnd_modules, ' ; ', cmnd1, ' ; ', cmnd2]);

        else
            [status] = system([cmnd1, ' ; ', cmnd2]);
        end


        if status ~= 0
            error(['Status returned value of ',num2str(status)])
        end
    end

else

    error('I Dont understand the file names you are trying to run')

end



end






