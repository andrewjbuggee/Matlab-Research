% read in data from LibRadTran atmospheric data files


% INPUTS
%   (1) filename:
%       (a) afglus.dat - AFGL atmospheric constituent profile. U.S. standard atmosphere 1976. ( AFGL-TR-86-0110) 
%       (b) afglms.dat - AFGL atmospheric constituent profile. midlatitude summer. ( AFGL-TR-86-0110)    
%       (c) afglmw.dat - AFGL atmospheric constituent profile. midlatitude winter. ( AFGL-TR-86-0110)   
%       (d) afglt.dat - AFGL atmospheric constituent profile. tropical. ( AFGL-TR-86-0110)


% By Andrew John Buggee

%%


function atm = read_atmos_prof_data(filename)

% Point to the folder in the computer that is being used


if strcmp('anbu8374',whatComputer)
    
        foldername = ['/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/atmmod/'];
    
elseif strcmp('andrewbuggee',whatComputer)
    
    foldername = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
             'LibRadTran/libRadtran-2.0.4/data/atmmod/'];
        
    
end



%% Read in data


delimeter = ' ';
headerLine = 2; % 0 if no header

data = importdata([foldername,filename],delimeter,headerLine);

if isstruct(data)==false
    % If this is false, then we may have read in the data incorrectly.
    % Try a tab delimeted file
    clear data
    delimeter = '\t';

    data = importdata([foldername,filename],delimeter,headerLine);

    % Delete the first column if it consists only of nans
    if all(data.data(:,1))==true
        data.data(:,1) = [];
    end

end

% Load the columns into structures
atm.altitude = data.data(:,1);      % km - altitude above surface
atm.pressure = data.data(:,2);      % mb - pressure above surface
atm.temp = data.data(:,3);          % K -   Temperature profile
atm.air_molec = data.data(:,4);     % #/cm^3 - number density of air molecules
atm.O3 = data.data(:,5);            % #/cm^3 - number density of Ozone
atm.O2 = data.data(:,6);            % #/cm^3 - number density of molecular oxygen
atm.H2O = data.data(:,7);           % #/cm^3 - number density of water vapor molecules
atm.CO2 = data.data(:,8);           % #/cm^3 - number density of carbon dioxide
atm.NO2 = data.data(:,9);           % #/cm^3 - number density of nitrogen dioxide






end