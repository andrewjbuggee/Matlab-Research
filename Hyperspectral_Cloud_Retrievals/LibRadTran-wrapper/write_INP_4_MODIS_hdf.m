%% ----- Create new .INP files from MODIS hdf data -----

% Inputs:
%   - EV: modis data. you must input an array, thus you need to pick a
%   specific band. A data cube will not work
%   - solar: the solar strucutre that contains the solar locaiton
%   - sensor: the sensor structure that contains locaiton of the sensor


% By Andrew J. Buggee

%%

function [inpNames,pixels2use] = write_INP_4_MODIS_hdf(inputs,pixels2use,modis)


% what computer are we using?

% a template file has been set up to be edited and saved as a new file
% determine which computer is being used
userName = whatComputer;

if strcmp(userName,'anbu8374')
    
    libRadTran_path = ['/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4'];
    templateFolder = ['/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/INP_template_files/'];
elseif strcmp(userName,'andrewbuggee')
    
    libRadTran_path = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
        'LibRadTran/libRadtran-2.0.4'];
    templateFolder = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/',...
        'Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/INP_template_files/'];
else
    error('I dont recognize this computer user name')
end



addpath(libRadTran_path);

% for each spectral bin, we have an image on the ground composed of 2030 *
% 1354 pixels. The swath on the ground is large enough that the solar
% zenith and solar azimuth change significantly. Ideally at each spectral bin, and
% each pixel, we could calculate the reflectance function and how it varies
% with changing tau and change re.

re = inputs.re;
tau_c = inputs.tau_c;
bands2run = inputs.bands2run;
pixel_row = pixels2use.res1km.row; % for everything I need in this code, we use the 1km resolution pixel locations
pixel_col = pixels2use.res1km.col; %
newFolder = [libRadTran_path,'/',inputs.INP_folderName]; % where the newly created .inp files will be saved




% we always edit the template file
oldFile = 'band_sza_saz__oceanSurface_template.INP';


% define the expressions that you wish to edit in the template file

oldExpr = {'wc_file 1D ../data/wc/WC_r04_T01_dist_.DAT', 'wavelength 0.00000 0.000000',...
    'sza 0000.0','phi0 0000.0','umu 0000.0','phi 0000.0'};



% we have 4 for loops to step throug

%   1) pixel
%   2) modis band
%   3) re
%   4) tau_c


% step through each band, each effective raidus and each optical depth
inpNames = cell(length(pixel_row), length(re),length(tau_c),length(bands2run));

for pp = 1:length(pixel_row)
    
    % lets determine the geometry of the pixel in question
    
    sza = modis.solar.zenith(pixel_row(pp),pixel_col(pp));
    saz = modis.solar.azimuth(pixel_row(pp),pixel_col(pp));
    
    % we need the cosine of the zenith viewing angle
    umu = round(cosd(double(modis.sensor.zenith(pixel_row(pp),pixel_col(pp)))),3); % values are in degrees
    phi = modis.sensor.azimuth(pixel_row(pp),pixel_col(pp));
    
    
    
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
    fileBegin = ['pixel_',num2str(pixel_row(pp)),'r_',num2str(pixel_col(pp)),'c_sza_',num2str(sza),'_saz_',num2str(saz),'_band_'];
    
    
    for bb = 1: length(bands2run)
        
        modis_band_num = inputs.bands2run(bb);
        modis_band = modisBands(modis_band_num);        % nm - grab the wavelenghts

        % create the new expression for the wavelength band of interest
        str = ['wavelength ',num2str(modis_band(2)),'.0 ',num2str(modis_band(3)),'.0'];
        
        if length(str) < length(oldExpr{2})
            
            lengthDiff = length(oldExpr{2}) - length(str);
            for ii = 1:(lengthDiff) % the minus one accounts for the decimal point
                str = [str,'0'];
            end
            
        elseif length(str) > length(oldExpr{2})
            
            error('new expression is greater than the old in expression in length')
            
        end
        
        newExpr{2} = str;
        
        % now lets step through radius values
        
        for rr = 1:length(re)
            
            % now lets step through tau values
            
            for tt = 1:length(tau_c)
                
                % redefine the old file each time
                inpNames{pp,rr,tt,bb} = [fileBegin,num2str(modis_band_num),'_r_',num2str(re(rr)),'_T_',num2str(tau_c(tt)),'.INP'];
                
                
                % lets define the new expressions to substitute the old ones
                
                
                % ----- Write a wc file for each re and tau -----
                
                % for now, lets code every cloud to be at the same
                % altitude, with the same thickness.
                
                % ------ FUTURE WORK -------
                % Retireve the cloud top pressure and infer a more accurate
                % cloud top height and thickness is set to be the same values 
                % described in: "Overview of the MODIS Collection 6 Cloud Optical 
                % Property (MOD06) Retrieval Look-up Tables" Amarasinghe
                % et. al 2017
                z_topBottom = [9, 8];                     % km - cloud top and bottom altitude above ground
                lambda = modis_band(1);                     % nm - center wavelength
                
                wc_filename = write_wc_file(re(rr), tau_c(tt), z_topBottom, lambda,inputs.clouds.distribution_type,...
                    inputs.clouds.distribution_variance, inputs.clouds.vert_homogeneity,...
                    inputs.clouds.wc_parameterization);
                % open up the cell structure
                wc_filename = wc_filename{1};
                % lets start with inserting the proper water cloud file
                
                if strcmp(inputs.clouds.distribution_type,'mono')==true
                    
                    if re(rr)<10 && tau_c(tt)<10
                        % newExpr{1} = ['wc_file 1D ../data/wc/WC_r0',num2str(re(rr)),'_T0',num2str(tau_c(tt)),'.DAT'];
                        newExpr{1} = ['wc_file 1D ../data/wc/',wc_filename(1:4),'0',wc_filename(5:7), '00',wc_filename(8:end)];
                    elseif re(rr)>=10 && tau_c(tt)<10
                        %newExpr{1} = ['wc_file 1D ../data/wc/WC_r',num2str(re(rr)),'_T0',num2str(tau_c(tt)),'.DAT'];
                        newExpr{1} = ['wc_file 1D ../data/wc/',wc_filename(1:8),'00',wc_filename(9:end)];
                    elseif re(rr)>=10 && tau_c(tt)>=10
                        %newExpr{1} = ['wc_file 1D ../data/wc/WC_r',num2str(re(rr)),'_T',num2str(tau_c(tt)),'.DAT'];
                        newExpr{1} = ['wc_file 1D ../data/wc/',wc_filename(1:10),'_',wc_filename(11:end)];
                    elseif re(rr)<10 && tau_c(tt)>=10
                        %newExpr{1} = ['wc_file 1D ../data/wc/WC_r0',num2str(re(rr)),'_T',num2str(tau_c(tt)),'.DAT'];
                        newExpr{1} = ['wc_file 1D ../data/wc/',wc_filename(1:4),'00',wc_filename(5:end)];
                    end
                    
                elseif strcmp(inputs.clouds.distribution_type, 'gamma')==true
                    
                    if re(rr)<10 && tau_c(tt)<10
                        % newExpr{1} = ['wc_file 1D ../data/wc/WC_r0',num2str(re(rr)),'_T0',num2str(tau_c(tt)),'.DAT'];
                        newExpr{1} = ['wc_file 1D ../data/wc/',wc_filename(1:4),'0',wc_filename(5:7), '0',wc_filename(8:end)];
                    elseif re(rr)>=10 && tau_c(tt)<10
                        %newExpr{1} = ['wc_file 1D ../data/wc/WC_r',num2str(re(rr)),'_T0',num2str(tau_c(tt)),'.DAT'];
                        newExpr{1} = ['wc_file 1D ../data/wc/',wc_filename(1:8),'0',wc_filename(9:end)];
                    elseif re(rr)>=10 && tau_c(tt)>=10
                        %newExpr{1} = ['wc_file 1D ../data/wc/WC_r',num2str(re(rr)),'_T',num2str(tau_c(tt)),'.DAT'];
                        newExpr{1} = ['wc_file 1D ../data/wc/',wc_filename];
                    elseif re(rr)<10 && tau_c(tt)>=10
                        %newExpr{1} = ['wc_file 1D ../data/wc/WC_r0',num2str(re(rr)),'_T',num2str(tau_c(tt)),'.DAT'];
                        newExpr{1} = ['wc_file 1D ../data/wc/',wc_filename(1:4),'0',wc_filename(5:end)];
                    end
                    
                else
                    
                    error([newline, 'I dont recognize the distribution string.', newline])
                end
                
                
                
                
                
                
                
                % ------ THIS WONT RUN IN PARALELL -----
                % because we are passing temporary variables out of the
                % loop and into this function below. This is not allowed in
                % parfor
                
                edit_INP_DAT_files(templateFolder,newFolder,oldFile,inpNames{pp,rr,tt,bb},oldExpr,newExpr);
                
            end
        end
        
    end
    
    
    
    
end





