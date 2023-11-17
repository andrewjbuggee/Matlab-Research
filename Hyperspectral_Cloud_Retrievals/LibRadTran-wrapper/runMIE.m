%% --- Run Intput Files with MIE ---
% ======================================
% The purpose of this script is to run MIE with the input files specfied
% by the user. The user will be required to provide to folder location of
% the input file, and the input file name. The user will also need to
% provide the output file name. This output file will be saved in the
% same folder as the input file. The input file will be fed into the
% command line in order to run uvspec.



% --- By Andrew J. Buggee ---
%%

function [droplet_settings] = runMIE(folderName,inputName,outputName)
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

% We only care about finding the effective radius and the droplet size
% distribution

if numFiles2Run==1
    
    textFile = fileread([folderName,inputName]);
    
    expr1 = '[^\n]*r_eff [^\n]*';
    expr2 = '[^\n]*distribution [^\n]*';
    
    
    match1 = regexp(textFile,expr1,'match'); % find effective radius
    match2 = regexp(textFile,expr2,'match'); % find distribution type
    
    
    index1_space1 = regexp(match1{1},'\s[0123456789]+'); % find a space that is followed by a number
    index1_space2 = regexp(match1{1},'[0123456789]\s+'); % find a space that comes after a number
    
    
    % determine the effective particle radius
    r_eff_str = cell(1,length(index1_space1));
    
    for ii = 1:length(index1_space1)
        
        r_eff_str{ii} = match1{1}(index1_space1(ii)+1:index1_space2(ii));
        
    end
    
    
    r_eff = str2double(r_eff_str); % convert to type double
    
    if length(r_eff) == 3
        r_eff = r_eff(1):r_eff(3):r_eff(2); % if the length is 3, then the entries are: fist,last,step
    elseif length(r_eff)==2 || length(r_eff)>3
        error('r_eff vector set up incorrectly')
    end
    
    
    droplet_settings{1} = r_eff; % store effective radius in the droplet settings
    
    
    % Find the distribution type
    
    if isempty(match2)==false
        matchGamma = regexp(match2{1},'gamma','match'); % find the space that is follwed by a g (for gamma) or an l (for log-normal)
        matchLog = regexp(match2{1},'lognormal','match');
        
        if isempty(matchGamma)==false
            distType = 'gamma';
            
            indexAlpha = regexp(match2{1},'\s[0123456789]+'); % find the space before the alpha parameter
            indexAlpha2 = regexp(match2{1},'[0123456789]\s+'); % find the last number before a space
            alpha = match2{1}(indexAlpha+1:indexAlpha2);
            
            droplet_settings{2} = distType;
            droplet_settings{3} = alpha;
            
        elseif isempty(matchLog)==false
            distType = 'log normal';
            
            indexSigma = regexp(match2{1},'\s[0123456789]+'); % find the space before the alpha parameter
            indexSigma2 = regexp(match2{1},'[0123456789]\s+'); % find the last number before a space
            sigma = match2{1}(indexSigma+1:indexSigma2);
            
            droplet_settings{2} = distType;
            droplet_settings{3} = sigma;
            
        else
            error('Something is wrong with the distribution type')
            
            
        end
    else
        distType = 'Homogenous';
        
        droplet_settings{2} = distType;
    end
    
    
    
    
    
    
    
elseif numFiles2Run>1
    
    droplet_settings = cell(numFiles2Run,3);
    
    for jj=1:numFiles2Run
        
        textFile = fileread([folderName,inputName{jj}]);
        
        expr1 = '[^\n]*r_eff [^\n]*';
        expr2 = '[^\n]distribution [^\n]*';
        
        
        match1 = regexp(textFile,expr1,'match'); % find effective radius
        match2 = regexp(textFile,expr2,'match'); % find distribution type
        
        
        index1_space1 = regexp(match1{1},'\s[0123456789]+'); % find a space that is followed by a number
        index1_space2 = regexp(match1{1},'[0123456789]\s+'); % find a space that comes after a number
        
        
        % determine the effective particle radius
        r_eff_str = cell(1,length(index1_space1));
        
        for ii = 1:length(index1_space1)
            
            r_eff_str{ii} = match1{1}(index1_space1(ii)+1:index1_space2(ii));
            
        end
        
        
        r_eff = str2double(r_eff_str); % convert to type double
        
        droplet_settings{jj,1} = r_eff; % store effective radius in the droplet settings
        
        
        % Find the dsitribution type
        
        if isempty(match2)==false
            matchGamma = regexp(match2{1},'gamma','match'); % find the space that is follwed by a g (for gamma) or an l (for log-normal)
            matchLog = regexp(match2{1},'lognormal','match');
            
            if isempty(matchGamma)==false
                distType = 'gamma';
                
                indexAlpha = regexp(match2{1},'\s[0123456789]+'); % find the space before the alpha parameter
                indexAlpha2 = regexp(match2{1},'[0123456789]\s+'); % find the last number before a space
                alpha = match2{1}(indexAlpha+1:indexAlpha2);
                
                droplet_settings{jj,2} = distType;
                droplet_settings{jj,3} = alpha;
                
            elseif isempty(matchLog)==false
                distType = 'log normal';
                
                indexSigma = regexp(match2{1},'\s[0123456789]+'); % find the space before the alpha parameter
                indexSigma2 = regexp(match2{1},'[0123456789]\s+'); % find the last number before a space
                sigma = match2{1}(indexSigma+1:indexSigma2);
                
                droplet_settings{jj,2} = distType;
                droplet_settings{jj,3} = sigma;
                
            else
                error('Something is wrong with the distribution type')
                
                
            end
        else
            distType = 'Homogenous';
            
            droplet_settings{jj,2} = distType;
            droplet_settings{jj,3} = [];
        end
        
        
        
        
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
% To run mie in the command line we have to point to its full location.
% To do this we will check to see what computer we are using

usrName = whatComputer;

if strcmp('anbu8374',usrName)
    
        mie_folderName = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/bin/';
    
elseif strcmp('andrewbuggee',usrName)
    
    mie_folderName = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/'...
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
    
    cmnd2 = ['(',mie_folderName,'mie ',...
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
        
        cmnd2 = ['(',mie_folderName,'mie ',...
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






