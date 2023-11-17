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

function [inputSettings] = runUVSPEC(folderName,inputName,outputName)
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
% -----

% Determine which computer you're using
computer_name = whatComputer;

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(computer_name,'anbu8374')==true
    

    folderSolar = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/solar_flux/';
    
elseif strcmp(computer_name,'andrewbuggee')==true
    
    folderSolar = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
        'LibRadTran/libRadtran-2.0.4/data/solar_flux/'];
end






%% ----- Lets Read the input file -----

% Lets determine the input settings

% First we need to determine how many files we need to run

if iscell(inputName)==true
    numFiles2Run = length(inputName);
elseif ischar(inputName)==true
    numFiles2Run = 1;
else
    error('I dont understand the input file')
end

if numFiles2Run==1
    
    textFile = fileread([folderName,inputName]);
    
    expr1 = '[^\n]*rte_solver [^\n]*';
    expr2 = '[^\n]*umu [^\n]*';
    expr3 = '[^\n]*phi [^\n]*';
    expr4 = '[^\n]*sza [^\n]*';
    expr5 = '[^\n]*phi0 [^\n]*';
    expr6 = '[^\n]*zout [^\n]*';
    expr7 = '[^\n]*source [^\n]*';
    expr8 = '[^\n]*wavelength [^\n]*';
    
    match1 = regexp(textFile,expr1,'match'); % find rte_solver typ
    match2 = regexp(textFile,expr2,'match'); % find consine of viewing angle vector
    match3 = regexp(textFile,expr3,'match'); % find azimuth viewing angle vector
    match4 = regexp(textFile,expr4,'match'); % find the solar zenith angle
    match5 = regexp(textFile,expr5,'match'); % find the solar azimuth angle
    match6 = regexp(textFile,expr6,'match'); % find the sensor altitude
    match7 = regexp(textFile,expr7,'match'); % find the source file
    match8 = regexp(textFile,expr8,'match'); % find the wavelength range in order to trim the source file
    
    index1_space1 = regexp(match1{1},'\s[a-z]+'); % find the spaces
    index1_space2 = regexp(match1{1},'[a-zA-Z_0-9]\s+'); % the second index should be the last letter in the solver type
    
    index2_space1 = regexp(match2{1},'\s[0123456789]+'); % find a space that is followed by a number
    index2_space2 = regexp(match2{1},'[0123456789]\s+'); % find a space that comes after a number
    
    index3_space1 = regexp(match3{1},'\s[0123456789]+'); % find a space that is followed by a number
    index3_space2 = regexp(match3{1},'[0123456789]\s+'); % find a space that comes after a number
    
    index4_space1 = regexp(match4{1},'\s[0123456789]+'); % There is only 1 value for the solar zenith angle
    index4_space2 = regexp(match4{1},'[0123456789]\s+'); % find a space that comes after a number
    
    index5_space1 = regexp(match5{1},'\s[0123456789]+'); % There is only 1 value for the solar zenith angle
    index5_space2 = regexp(match5{1},'[0123456789]\s+'); % find a space that comes after a number - but since the varaible name, phi0, ends in a 0, we have to manually choose the second space
    
    index6_space1 = regexp(match6{1},'\s[0123456789]+'); % There is only 1 value for the solar zenith angle
    index6_space2 = regexp(match6{1},'[0123456789]\s+'); % find a space that comes after a number
    
    % don't let the lack of a source stop you!
    if isempty(match7)==true

    else
        index7_space1 = regexp(match7{1},'\s[a-z]'); % find the spaces
        index7_space2 = regexp(match7{1},'[a-z]\s'); % Brackets treat the symbol literally. number of decimals tells us how many values there are in the vector
        index7_file1 = regexp(match7{1},'flux[/][a-z]'); % find the locaition a letter follows two dots and a forward slash
        index7_file2 = regexp(match7{1},'[.]dat');
    end

    % don't let a lack of wavelengths stop you!
    if isempty(match8)==false
        index8_space1 = regexp(match8{1},'\s[0123456789]+'); % There is only 1 value for the solar zenith angle
        index8_space2 = regexp(match8{1},'[0123456789]\s+'); % find a space that comes after a number
    end

    
    % determine the rte_solver type
    rte_solver = match1{1}(index1_space1(1)+1:index1_space2(2));
    
    % determine the umu vector
    umuStr = cell(1,length(index2_space1));
    
    for ii = 1:length(index2_space1)
        umuStr{ii} = match2{1}(index2_space1(ii)+1:index2_space2(ii));
    end
    
    umuVec = str2double(umuStr);
    
    % determine the phi vector
    phiStr = cell(1,length(index3_space1));
    
    for ii = 1:length(index3_space1)
        phiStr{ii} = match3{1}(index3_space1(ii)+1:index3_space2(ii));
    end
    
    phiVec = str2double(phiStr);
    
    % find the solar zenith angle
    sza = match4{1}(index4_space1+1:index4_space2);
    sza = str2double(sza);
    
    % find the solar azimuth angle
    saz = match5{1}(index5_space1+1:index5_space2(2));
    saz = str2double(saz);
    
    % find the sensor altitude
    
    
    if isempty(index6_space1)==false
        zout = match6{1}(index6_space1(1)+1:index6_space2(2)-1); % this would be for a numeric value
        zout = str2double(zout);
    elseif isempty(index6_space1)==true
        indexString1 = regexp(match6{1},'\s[a-z][a-z][a-z]'); % this would be for a string value like toa or boa
        indexString2 = regexp(match6{1},'[a-z]\s');
        zout_str = match6{1}(indexString1(1)+1:indexString2(2));
        
        if strcmp(zout_str,'toa')==true
            zout = 100;
        elseif strcmp(zout_str,'sur')==true
            zout = 0;
        end
    end
    
    % find the wavelength range of the output file
    if isempty(match8)==false
        wavelength_str= cell(1,length(index8_space1));

        for ii = 1:length(index8_space1)
            wavelength_str{ii} = match8{1}(index8_space1(ii)+1:index8_space2(ii));
        end
        wavelength = str2double(wavelength_str);

    else
        wavelength = [];
    end

    % ---------------------------------------------------------
    % ------ determine the solar or thermal source file -------
    % we want to store the source flux as a vector
    % all solar source files will be located in the folder: /Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval-Research/LibRadTran/libRadtran-2.0.4/data/solar_flux
    % all thermal source files will be located in the foler:

    if isempty(match7)==false
        if strcmp('solar',match7{1}(index7_space1(1)+1:index7_space2(2)))

            if length(match7{1})<=54
                % This happens when the input is simple 'source solar' with no
                % specified file
                fileSolar = 'internal';
                source = [];

            else
                fileSolar = match7{1}(index7_file1(1)+5:index7_file2(1)+3);

                % open the file for reading
                file_id = fopen([folderSolar,fileSolar], 'r');   % 'r' tells the function to open the file for reading

                format_spec = '%f %f';                                  % two floating point numbers
                source_data = textscan(file_id, format_spec, 'Delimiter',' ',...
                    'MultipleDelimsAsOne',1, 'CommentStyle','#');

                % now we clip source to match the length of our wavelength vector
                % if we run a monochromatic calculation, we do the following first.
                % Then, for multispectral calculations
                if length(wavelength)==1
                    indexSource = source_data{1}==round(wavelength); % can only have integer values for wavelength
                    source(:,1) = source_data{1}(indexSource);
                    source(:,2) = source_data{2}(indexSource);

                elseif length(wavelength)>1

                    indexSource = source_data{1}>=wavelength(1) & source_data{1}<=wavelength(2);
                    source(:,1) = source_data{1}(indexSource);
                    source(:,2) = source_data{2}(indexSource);

                elseif isempty(wavelength)==true
                    source(:,1) = source_data{1};
                    source(:,2) = source_data{2};
                end
            end
        end

    else

        source = [];

    end

    % Pull all input settings into a cell array
    % first lets give them headers and labels:
    
    inputSettings{1,1} = 'Solver Type';
    inputSettings{1,2} = 'Cos(zva)';
    inputSettings{1,3} = 'Azimuthal Angle';
    inputSettings{1,4} = 'Solar Zenith Angle';
    inputSettings{1,5} = 'Solar Azimuthal Angle';
    inputSettings{1,6} = 'Sensor Altitude (km)';
    inputSettings{1,7} = 'Source Wavelength (nm) and Irradiance';
    
    inputSettings{2,1} = rte_solver;
    inputSettings{2,2} = umuVec;
    inputSettings{2,3} = phiVec;
    inputSettings{2,4} = sza;
    inputSettings{2,5} = saz;
    inputSettings{2,6} = zout;
    inputSettings{2,7} = source;
    
elseif numFiles2Run>1
    
    inputSettings = cell(numFiles2Run+1,7);
    
    inputSettings{1,1} = 'Solver Type';
    inputSettings{1,2} = 'Cos(zva)';
    inputSettings{1,3} = 'Azimuthal Angle';
    inputSettings{1,4} = 'Solar Zenith Angle';
    inputSettings{1,5} = 'Solar Azimuthal Angle';
    inputSettings{1,6} = 'Sensor Altitude (km)';
    inputSettings{1,7} = 'Source Wavelength (nm) and Irradiance';
    
    for jj=1:numFiles2Run
        
        textFile = fileread([folderName,inputName{jj}]);
        
        expr1 = '[^\n]*rte_solver [^\n]*';
        expr2 = '[^\n]*umu [^\n]*';
        expr3 = '[^\n]*phi [^\n]*';
        expr4 = '[^\n]*sza [^\n]*';
        expr5 = '[^\n]*phi0 [^\n]*';
        expr6 = '[^\n]*zout [^\n]*';
        expr7 = '[^\n]*source [^\n]*';
        expr8 = '[^\n]*wavelength [^\n]*';
        
        match1 = regexp(textFile,expr1,'match'); % find rte_solver typ
        match2 = regexp(textFile,expr2,'match'); % find consine of viewing angle vector
        match3 = regexp(textFile,expr3,'match'); % find azimuth viewing angle vector
        match4 = regexp(textFile,expr4,'match'); % find the solar zenith angle
        match5 = regexp(textFile,expr5,'match'); % find the solar azimuth angle
        match6 = regexp(textFile,expr6,'match'); % find the sensor altitude
        match7 = regexp(textFile,expr7,'match'); % find the source file
        match8 = regexp(textFile,expr8,'match'); % find the wavelength range in order to trim the source file
        
        index1_space1 = regexp(match1{1},'\s[a-z]+'); % find the spaces
        index1_space2 = regexp(match1{1},'[a-z]\s+'); % the second index should be the last letter in the solver type
        
        index2_space1 = regexp(match2{1},'\s[0123456789]+'); % find a space that is followed by a number
        index2_space2 = regexp(match2{1},'[0123456789]\s+'); % find a space that comes after a number
        
        index3_space1 = regexp(match3{1},'\s[0123456789]+'); % find a space that is followed by a number
        index3_space2 = regexp(match3{1},'[0123456789]\s+'); % find a space that comes after a number
        
        index4_space1 = regexp(match4{1},'\s[0123456789]+'); % There is only 1 value for the solar zenith angle
        index4_space2 = regexp(match4{1},'[0123456789]\s+'); % find a space that comes after a number
        
        index5_space1 = regexp(match5{1},'\s[0123456789]+'); % There is only 1 value for the solar zenith angle
        index5_space2 = regexp(match5{1},'[0123456789]\s+'); % find a space that comes after a number
        
        index6_space1 = regexp(match6{1},'\s[0123456789]+'); % There is only 1 value for the solar zenith angle
        index6_space2 = regexp(match6{1},'[0123456789]\s+'); % find a space that comes after a number
        
        index7_space1 = regexp(match7{1},'\s[a-z]'); % find the spaces
        index7_space2 = regexp(match7{1},'[a-z]\s'); % Brackets treat the symbol literally. number of decimals tells us how many values there are in the vector
        index7_file1 = regexp(match7{1},'flux[/][a-z]'); % find the locaition a letter follows two dots and a forward slash
        index7_file2 = regexp(match7{1},'[.]dat');
        
        index8_space1 = regexp(match8{1},'\s[0123456789]+'); % There is only 1 value for the solar zenith angle
        index8_space2 = regexp(match8{1},'[0123456789]\s+'); % find a space that comes after a number
        
        
        % determine the rte_solver type
        rte_solver = match1{1}(index1_space1(1)+1:index1_space2(2));
        
        % determine the umu vector
        umuStr = cell(1,length(index2_space1));
        
        for ii = 1:length(index2_space1)
            umuStr{ii} = match2{1}(index2_space1(ii)+1:index2_space2(ii));
        end
        
        umuVec = str2double(umuStr);
        
        % determine the phi vector
        phiStr = cell(1,length(index3_space1));
        
        for ii = 1:length(index3_space1)
            phiStr{ii} = match3{1}(index3_space1(ii)+1:index3_space2(ii));
        end
        
        phiVec = str2double(phiStr);
        
        % find the solar zenith angle
        sza = match4{1}(index4_space1+1:index4_space2);
        sza = str2double(sza);
        
        % find the solar azimuth angle
        saz = match5{1}(index5_space1+1:index5_space2(2));
        saz = str2double(saz);
        
        % find the sensor altitude
        
        
        if isempty(index6_space1)==false
            zout = match6{1}(index6_space1(1)+1:index6_space2(2)-1); % this would be for a numeric value
            zout = str2double(zout);
        elseif isempty(index6_space1)==true
            indexString1 = regexp(match6{1},'\s[a-z][a-z][a-z]'); % this would be for a string value like toa or boa
            indexString2 = regexp(match6{1},'[a-z]\s');
            zout_str = match6{1}(indexString1(1)+1:indexString2(2));
            
            if strcmp(zout_str,'toa')==true
                zout = 100;
            elseif strcmp(zout_str,'boa')==true
                zout = 0;
            end
        end
        
        % find the wavelength range of the output file
        wavelength_str= cell(1,length(index8_space1));
        
        for ii = 1:length(index8_space1)
            wavelength_str{ii} = match8{1}(index8_space1(ii)+1:index8_space2(ii));
        end
        wavelength = str2double(wavelength_str);
        
        % determine the solar or thermal source file
        % we want to store the source flux as a vector
        % all solar source files will be located in the folder: /Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval-Research/LibRadTran/libRadtran-2.0.4/data/solar_flux
        % all thermal source files will be located in the foler:
        if strcmp('solar',match7{1}(index7_space1(1)+1:index7_space2(2)))
            % find which computer you're running on
            if strcmp(folderName(1:15),'/Users/anbu8374')
                folderSolar = ['/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/solar_flux/'];
                
            elseif strcmp(folderName(1:19),'/Users/andrewbuggee')
                folderSolar = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
                    'LibRadTran/libRadtran-2.0.4/data/solar_flux/'];
            else
                error('Dont reconginze the libradtran solar flux folder')
            end
            fileSolar = match7{1}(index7_file1(1)+5:index7_file2(1)+3);
            sourceFile = fileread([folderSolar,fileSolar]);
            
            exprSource = '[^\n]*[\d][\d] [^\n]*'; % look for the new lines with atleast two digits in a row
            matchSource = regexp(sourceFile,exprSource,'match'); % find rte_solver typ
            
            source = zeros(length(matchSource),2);
            for ii = 1:length(matchSource)
                source(ii,:) = str2num(matchSource{ii});
            end
            
            % now we clip source to match the length of our wavelength vector
            % if we run a monochromatic calculation, we do the following first.
            % Then, for multispectral calculations
            if length(wavelength)==1
                indexSource = source(:,1)==round(wavelength); % can only have integer values for wavelength
                source = source(repmat(indexSource,1,2));
                source = reshape(source,size(source,1)/2,[]);
                
            elseif length(wavelength)>1
                
                indexSource = source(:,1)>=wavelength(1) & source(:,1)<=wavelength(2);
                source = source(repmat(indexSource,1,2));
                source = reshape(source,size(source,1)/2,[]);
            end
        end
        
        % Pull all input settings into a cell array
        inputSettings{jj+1,1} = rte_solver;
        inputSettings{jj+1,2} = umuVec;
        inputSettings{jj+1,3} = phiVec;
        inputSettings{jj+1,4} = sza;
        inputSettings{jj+1,5} = saz;
        inputSettings{jj+1,6} = zout;
        inputSettings{jj+1,7} = source;
        
    end
end



%% Running the .INP file from the command line

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

usrName = whatComputer;

if strcmp('anbu8374',usrName)
    
    uvspec_folderName = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/bin/';
    
elseif strcmp('andrewbuggee',usrName)
    
    uvspec_folderName = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/'...
        'Hyperspectral-Cloud-Droplet-Retrieval/',...
        'LibRadTran/libRadtran-2.0.4/bin/'];
    
else
    
    error('I dont know where the uvspec function is! Tell me what folder its in, please')
    
end


% using the function 'system' runs commans in the terminal window
cmnd1 = ['cd ', folderName];



if numFiles2Run==1
    
    
    % cmnd2 = [uvspec_folderName,'uvspec ',...
    %            '< ',inputName,' > ', outputName];
    
    cmnd2 = ['(',uvspec_folderName,'uvspec ',...
        '< ',inputName,' > ', outputName,'.OUT',')>& errMsg.txt'];
    % a successful command will return a status of 0
    % an unsuccessful command will return a status of 1
    
    
    
    [status] = system([cmnd1, ' ; ', cmnd2]);
    if status ~= 0
        error(['Status returned value of ',num2str(status)])
    end
    
elseif numFiles2Run>1
    
    
    for ii = 1:numFiles2Run
        
        
        % cmnd2 = [uvspec_folderName,'uvspec ',...
        %            '< ',inputName,' > ', outputName];
        
        cmnd2 = ['(',uvspec_folderName,'uvspec ',...
            '< ',inputName{ii},' > ', outputName{ii},'.OUT',')>& errMsg.txt'];
        % a successful command will return a status of 0
        % an unsuccessful command will return a status of 1
        
        [status] = system([cmnd1, ' ; ', cmnd2]);
        if status ~= 0
            error(['Status returned value of ',num2str(status)])
        end
    end
    
else
    
    error('I Dont understand the file names you are trying to run')
    
end



end






