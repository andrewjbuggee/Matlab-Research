%% ----- WHATS CHANGING?! -----

% this fun little code checks to see what is changing in between different
% .INP files for LibRadTran! The differences should all be located in the
% inputSettings file. So this code will look through for changes


% By Andrew J. Buggee
%%

function [changes] = whatsChanging_INP(inputSettings)

% for there to even be multiple settings, inputSettings must be a cell
% array. If its not, we should reject the call to this function

if size(inputSettings,1)>1 && size(inputSettings,2)>1
    
    % if this is true, then we have a cell array as an input
    
else
    
    error('This is not an array! Nothing can change without an array.')
end


% first lets unpack the settings file

rte_solver = inputSettings(1,:); % solver type across files
umuVec = [inputSettings{2,:}]; % cos(viewing angle) across files
phiVec = [inputSettings{3,:}]; % observer azimuth angle across files
szaVec = [inputSettings{4,:}]; % solar zenith angle across files
sazVec = [inputSettings{5,:}]; % solar azimuth across files
zoutVec = [inputSettings{6,:}]; % solar azimuth across files
source = [inputSettings{7,:}]; % wavelength band across different files - (wavelength_file1,irradiance_file1, wavelength_file2,irradiance_file2,...)


% rearrange the source inputs into wavelength and irradiance components

index_w = 1:2:(size(source,2)-1);
index_i = 2:2:size(source,2);
    
source_w = source(:,index_w);
source_i = source(:,index_i);

%% ----- Check each input Settings -----
% we now have to check each setting to see if anything changed between
% files

% lets check the solver type first

numU = size(unique(rte_solver),1)*size(unique(rte_solver),2);


if numU~=1
    changes{1,1} = 'rte_solver';
    changes{1,2} = numU;
    
end


% ----------------------------------
% lets check the cos(observer angle)

numU = size(unique(umuVec),1)*size(unique(umuVec),2);


if numU~=1
    % if changes has already been made, simply add onto it. Otherwise,
    % make it
    if exist('changes','var') == true
        numRow = size(changes,1);
        numCol = size(changes,2);
        
        changes{numRow+1,1} = 'umu';
        changes{numCol+1,2} = numU;
    else
        changes{1,1} = 'umu';
        changes{1,2} = numU;
    end
 
end


% ---------------------------------------
% lets check the observer azimuthal angle

numU = size(unique(phiVec),1)*size(unique(phiVec),2);


if numU~=1
    % if changes has already been made, simply add onto it. Otherwise,
    % make it
    if exist('changes','var') == true
        numRow = size(changes,1);
        numCol = size(changes,2);
        
        changes{numRow+1,1} = 'phi';
        changes{numCol+1,2} = numU;
    else
        changes{1,1} = 'phi';
        changes{1,2} = numU;
    end
 
end


% ----------------------------------
% lets check the solar zenith angle

numU = size(unique(szaVec),1)*size(unique(szaVec),2);


if numU~=1
    % if changes has already been made, simply add onto it. Otherwise,
    % make it
    if exist('changes','var') == true
        numRow = size(changes,1);
        numCol = size(changes,2);
        
        changes{numRow+1,1} = 'sza';
        changes{numCol+1,2} = numU;
    else
        changes{1,1} = 'sza';
        changes{1,2} = numU;
    end
 
end


% ------------------------------------
% lets check the solar azimuthal angle

numU = size(unique(sazVec),1)*size(unique(sazVec),2);


if numU~=1
    % if changes has already been made, simply add onto it. Otherwise,
    % make it
    if exist('changes','var') == true
        numRow = size(changes,1);
        numCol = size(changes,2);
        
        changes{numRow+1,1} = 'saz';
        changes{numCol+1,2} = numU;
    else
        changes{1,1} = 'saz';
        changes{1,2} = numU;
    end
 
end

% ----------------------------------
% lets check the observer altitude

numU = size(unique(zoutVec),1)*size(unique(zoutVec),2);


if numU~=1
    % if changes has already been made, simply add onto it. Otherwise,
    % make it
    if exist('changes','var') == true
        numRow = size(changes,1);
        numCol = size(changes,2);
        
        changes{numRow+1,1} = 'zout';
        changes{numCol+1,2} = numU;
    else
        changes{1,1} = 'zout';
        changes{1,2} = numU;
    end
 
end


% ----------------------------------
% lets check the source wavelength

% in this case we deal with vectors, so we check to see that we have only 1
% vector. Otherwise, if we have multiple, then we must have a changing
% source. So we only need to check the number of rows 



numU = size(unique(source_w','rows'),1);


if numU~=1
    % if changes has already been made, simply add onto it. Otherwise,
    % make it
    if exist('changes','var') == true
        numRow = size(changes,1);
        numCol = size(changes,2);
        
        changes{numRow+1,1} = 'source_wavelength';
        changes{numCol+1,2} = numU;
    else
        changes{1,1} = 'source_wavelength';
        changes{1,2} = numU;
    end
 
end


% ----------------------------------
% lets check the source irradiance

% in this case we deal with vectors, so we check to see that we have only 1
% vector. Otherwise, if we have multiple, then we must have a changing
% source. So we only need to check the number of rows 



numU = size(unique(source_i','rows'),1);


if numU~=1
    % if changes has already been made, simply add onto it. Otherwise,
    % make it
    if exist('changes','var') == true
        numRow = size(changes,1);
        numCol = size(changes,2);
        
        changes{numRow+1,1} = 'source_irrad';
        changes{numCol+1,2} = numU;
    else
        changes{1,1} = 'source_irrad';
        changes{1,2} = numU;
    end
 
end


%% Now we need to collect all of the changing variables

if exist('changes','var') ~= true
    
    changes = 0; % we need someway of telling matlab that there are no chagnes between files
    
end
    
    


end
