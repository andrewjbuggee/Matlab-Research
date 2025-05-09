%% ----- Create new .INP files from MODIS hdf data -----

% Inputs:
%   - EV: modis data. you must input an array, thus you need to pick a
%   specific band. A data cube will not work
%   - solar: the solar strucutre that contains the solar locaiton
%   - sensor: the sensor structure that contains locaiton of the sensor


% By Andrew J. Buggee

%%

function [inpNames] = write_INP_4_MODIS_singlePixel_withProfile(wc_profile_fileName,INP_folderName,modis,pixel_row,pixel_col,retrieval)


% what computer are we using?

% a template file has been set up to be edited and saved as a new file
% determine which computer is being used
userName = whatComputer;

if strcmp(userName,'anbu8374')
    
    libRadTran_path = ['/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4'];
    templateFolder = ['/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/INP_template_files/'];
elseif strcmp(userName,'andrewbuggee')
    
    libRadTran_path = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval-Research/',...
        'LibRadTran/libRadtran-2.0.4'];
    templateFolder = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/',...
        'Hyperspectral-Cloud-Droplet-Retrieval-Research/LibRadTran/libRadtran-2.0.4/INP_template_files/'];
else
    error('I dont recognize this computer user name')
end



addpath(libRadTran_path);


newFolder = [libRadTran_path,'/',INP_folderName]; % where the newly created .inp files will be saved




% we always edit the template file
oldFile = 'band_sza_saz__oceanSurface_adiabatic_profile_template.INP';


% define the expressions that you wish to edit in the template file

oldExpr = {'wc_file 1D ../data/wc/WC_profile_adiabatic_rbot00_rtop00_T00.DAT', 'wavelength 0.00000 0.000000',...
    'sza 0000.0','phi0 0000.0','umu 0000.0','phi 0000.0'};


% --- compute the forward model at our current estimate ---
r_top = retrieval(1);
r_bottom = retrieval(2);
tau_c = retrieval(3);



% lets determine the geometry of the pixel in question

sza = modis.solar.zenith(pixel_row,pixel_col);
saz = modis.solar.azimuth(pixel_row,pixel_col);

% we need the cosine of the zenith viewing angle
umu = round(cosd(double(modis.sensor.zenith(pixel_row,pixel_col))),3); % values are in degrees
phi = modis.sensor.azimuth(pixel_row,pixel_col);


% ----- Define bands 2 run uvspec ------------
bands2run = 1:7;



% step through each band, each effective raidus and each optical depth
inpNames = cell(1,length(bands2run));


% ----- lets edit the newExpression string -----


% some new expressions change in the for loop, and others are fixed like
% the geometry of the chosen pixel
% the new expressions have to be the same length. For geometry that means
% we always need 4 numerals and a decimal point

% lets fix the solar zenith angle first
str = ['sza ',num2str(sza),'.0'];

if length(str) < length(oldExpr{3})
    
    lengthDiff = length(oldExpr{3}) - length(str);
    for ii = 1:(lengthDiff) % the minus one accounts for the decimal point
        str = [str,'0'];
    end
    
elseif length(str) > length(oldExpr{3})
    
    error('new expression is greater than the old in expression in length')
    
end

newExpr{3} = str;

% now lets write the new solar azimuth angle. Modis defines the azimuth
% angle as [0,180] and [-180,0, whereas libradtran defines the azimuth
% angle as [0,360]. So we need to make this adjustment
if saz<0
    saz = saz+360;
end
str = ['phi0 ',num2str(saz),'.0'];

if length(str) < length(oldExpr{4})
    
    lengthDiff = length(oldExpr{4}) - length(str);
    for ii = 1:(lengthDiff) % the minus one accounts for the decimal point
        str = [str,'0'];
    end
    
elseif length(str) > length(oldExpr{4})
    
    error('new expression is greater than the old in expression in length')
    
end

newExpr{4} = str;


% now lets write the cosine of the zentih viewing angle. LibRadTran defines
% this as looking down on the earth, measuring upwelling radiation if umu>0
% and being on the earth looking up, measuring downwelling radition, if
% umu<0. If umu is 0, it implies a device looking horizontally. MODIS
% defines the sensor zenith angle to be between 0 and 180, where 0 is at
% zenith. This is equivelant to the libradtran method since cos(0)=1
% implies a device looking down. A sensor zenith of 180 implies looking up,
% and cos(180) = -1, which is the same definition libradtran uses.

if umu==1 || umu==-1 || umu==0
    str = ['umu ',num2str(umu),'.0'];
else
    str = ['umu ',num2str(umu)];
end

if length(str) < length(oldExpr{5})
    
    lengthDiff = length(oldExpr{5}) - length(str);
    for ii = 1:(lengthDiff) % the minus one accounts for the decimal point
        str = [str,'0'];
    end
    
elseif length(str) > length(oldExpr{5})
    
    % if the new string is greater than the old, we will assume there are
    % zeros at the end that can be delted, or that the extra precision is
    % not neccessary
    lengthDiff = length(str) - length(oldExpr{5});
    
    str = str(1:(end-lengthDiff));
    
    
end

newExpr{5} = str;

% now lets write the azimuth viewing angle. Modis defines the azimuth
% angle as [0,180] and [-180,0], whereas libradtran defines the azimuth
% angle as [0,360]. So we need to make this adjustment

if phi<0
    phi = phi+360;
end

str = ['phi ',num2str(phi),'.0'];

if length(str) < length(oldExpr{6})
    
    lengthDiff = length(oldExpr{6}) - length(str);
    for ii = 1:(lengthDiff) % the minus one accounts for the decimal point
        str = [str,'0'];
    end
    
elseif length(str) > length(oldExpr{6})
    
    error('new expression is greater than the old in expression in length')
    
end

newExpr{6} = str;




% create the begining of the file name string
fileBegin = ['wc_profile_pixel_',num2str(pixel_row),'r_',num2str(pixel_col),'c_sza_',num2str(sza),'_saz_',num2str(saz),'_band_'];


for bb = 1: length(bands2run)
    
    modis_band_num = bands2run(bb);
    
    % create the new expression for the wavelength band of interest
    if modis_band_num<=2
        str = ['wavelength ',num2str(modis.EV.m250.bands.lowerBound(modis_band_num)),'.0 ',...
            num2str(modis.EV.m250.bands.upperBound(modis_band_num)),'.0'];
    elseif modis_band_num>2
        str = ['wavelength ',num2str(modis.EV.m500.bands.lowerBound(modis_band_num-2)),'.0 ',...
            num2str(modis.EV.m500.bands.upperBound(modis_band_num-2)),'.0'];
    end
    
    if length(str) < length(oldExpr{2})
        
        lengthDiff = length(oldExpr{2}) - length(str);
        for ii = 1:(lengthDiff) % the minus one accounts for the decimal point
            str = [str,'0'];
        end
        
    elseif length(str) > length(oldExpr{2})
        
        error('new expression is greater than the old in expression in length')
        
    end
    
    newExpr{2} = str;
    
    
    % redefine the old file each time
    inpNames{bb} = [fileBegin,num2str(modis_band_num),'_rbot_',num2str(round(r_bottom)),'_rtop_',num2str(round(r_top)),'_T_',num2str(round(tau_c)),'.INP'];
    
    
    % lets define the new expression that defines the water cloud file
    newExpr{1} = ['wc_file 1D ../data/wc/',wc_profile_fileName];
    
    
    
    
    
    
    % ------ THIS WONT RUN IN PARALELL -----
    % because we are passing temporary variables out of the
    % loop and into this function below. This is not allowed in
    % parfor
    
    % for checking the strings are all the same size
%     for cc = 1:length(oldExpr)
%         disp([length(oldExpr{cc}),length(newExpr{cc})])
%     end
    
    edit_INP_DAT_files(templateFolder,newFolder,oldFile,inpNames{bb},oldExpr,newExpr);
    
end



end












