%% ----- Create wc.dat files with droplet profile -----

% this function writes .dat files for clouds with a specific droplet
% profile

% By Andrew J. Buggee
%%

function [newFile] = write_waterCloud_files_with_profile(r_top,r_bottom,tau_c,profile_type, wavelength_tau_c)

% determine what computer is being used in order to find the water cloud
% folder
computer_name = whatComputer;
if strcmp(computer_name,'andrewbuggee')==true
    
    oldFolder = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval-Research/LibRadTran/libRadtran-2.0.4/data/wc/'];
    
elseif strcmp(computer_name,'anbu8374')==true
    
    oldFolder = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/wc/';
    
end

% we save all water cloud files in the same locaiton as the template
newFolder = oldFolder;

% define the template that we will be editing
oldFile = 'file2Edit_profile';

% define the old expressions that will be edited
oldExpr = {'      2.000   1.0000  2.5000', '      2.100   1.0000  2.5000',...
    '      2.200   1.0000  2.5000', '      2.300   1.0000  2.5000',...
    '      2.400   1.0000  2.5000', '      2.500   1.0000  2.5000',...
    '      2.600   1.0000  2.5000', '      2.700   1.0000  2.5000',...
    '      2.800   1.0000  2.5000', '      2.900   1.0000  2.5000'};

% create a cell array that is the same length as oldExpr, since we create
% as many new expressions as there are old ones to replace
newExpr = cell(1,length(oldExpr));


% no matter the profile, we are create a droplet profile that extends 1 km,
% from the base of the cloud at 2km above the ground to the top of the
% cloud at 3 km above the ground
% we only define the values at the base of each layer
h = 1000; % meters - cloud thickness
dh = 100;
z0 = 2000; % meters - cloud base location

% Set up the z-vector in the same manner as the watercloud files set up the
% z coordinate. That is, an entry for a LWC of 1 at a height of 2000 meters
% implies that LWC is 1 from 200 meters to the next height entry in the
% data file
z = z0:dh:(z0+h - dh); % meters - the locations in geometric space that define the bottom of each layer within the cloud'

% define the wavelength vector for integration
wavelength_vec = linspace(wavelength_tau_c, wavelength_tau_c,length(z))';                   % nm - wavelength that defines the cloud optical depth

if strcmp(profile_type,'adiabatic')
    
    % define the droplet profile for adiabatic behavior in geometric
    % coordniates (s.platnick 2000)
    
    b0 = r_bottom^3;
    b1 = r_top^3 - r_bottom^3;
    re_z = (b0 + b1 * (z-z0)./(h-dh)).^(1/3); % droplet profile in geometric coordniate system for an adiabatic cloud
    
    % libRadTran requires liquid water content, but we want to enforce a
    % tau_c_str over the cloud. So we calculate what the total liquid water
    % content using the tau that was input by the user
    
    % ----- I have to fudge the density of water ------
    % This is the only way I can get my estimates of optical depth to match
    % LibRadTrans estimates
    
    rho_liquid_water = 859900;              % grams/m^3 - density of liquid water
    %rho_liquid_water = 1e6;                 % grams/cm^3 - density of liquid water at 0 C
    
    % when we assume an adiabatic cloud, lwc varies linearly with geometric
    % height, z (lwc ~ z). We set a boundary condition at the base of the
    % cloud. Ideally this would be zero, but the boundary of clouds are
    % ambigious (s.platnick 2000)
    
    % ----------------------- ASSUMPTION ----------------------
    % We assume the number concentration is constant with height
    % ---------------------------------------------------------
    
    % First lets compute the extinction efficient
    % Lets assume a monodispersed distribution
    yq = interp_mie_computed_tables([wavelength_vec,re_z'],'mono');
    
    Qext = yq(:,5);                                                                          % Extinction efficiency
    Nc = tau_c/(pi*trapz((z-z0),Qext'.*(re_z*1e-6).^2));                                     % m^(-3) - number concentration
    
    lwc_z = 4/3 * pi * rho_liquid_water * (re_z*1e-6).^3 * Nc;
    
    % create the new file name
    r_bottom_round = round(r_bottom);
    r_top_round = round(r_top);
    tau_c_round = round(tau_c);
    
    r_bottom_str = num2str(r_bottom_round);
    r_top_str = num2str(r_top_round);
    tau_c_str = num2str(tau_c_round);
    
    if r_bottom_round<10 && r_top_round<10 && tau_c_round<10
        newFile = ['WC_profile_',profile_type,'_rbot0',r_bottom_str,'_rtop0',r_top_str,'_T0',tau_c_str,'.DAT'];
    elseif r_bottom_round<10 && r_top_round<10 && tau_c_round>=10
        newFile = ['WC_profile_',profile_type,'_rbot0',r_bottom_str,'_rtop0',r_top_str,'_T',tau_c_str,'.DAT'];
    elseif r_bottom_round<10 && r_top_round>=10 && tau_c_round<10
        newFile = ['WC_profile_',profile_type,'_rbot0',r_bottom_str,'_rtop',r_top_str,'_T0',tau_c_str,'.DAT'];
    elseif r_bottom_round<10 && r_top_round>=10 && tau_c_round>=10
        newFile = ['WC_profile_',profile_type,'_rbot0',r_bottom_str,'_rtop',r_top_str,'_T',tau_c_str,'.DAT'];
    elseif r_bottom_round>=10 && r_top_round<10 && tau_c_round<10
        newFile = ['WC_profile_',profile_type,'_rbot',r_bottom_str,'_rtop0',r_top_str,'_T0',tau_c_str,'.DAT'];
    elseif r_bottom_round>=10 && r_top_round<10 && tau_c_round>=10
        newFile = ['WC_profile_',profile_type,'_rbot',r_bottom_str,'_rtop0',r_top_str,'_T',tau_c_str,'.DAT'];
    elseif r_bottom_round>=10 && r_top_round>=10 && tau_c_round<10
        newFile = ['WC_profile_',profile_type,'_rbot',r_bottom_str,'_rtop',r_top_str,'_T0',tau_c_str,'.DAT'];
    elseif r_bottom_round>=10 && r_top_round>=10 && tau_c_round>=10
        newFile = ['WC_profile_',profile_type,'_rbot',r_bottom_str,'_rtop',r_top_str,'_T',tau_c_str,'.DAT'];
        
    end
    
    
    for xx = 1:length(oldExpr)
        % we need to make sure the new string length is identical to the
        % old one
        
        % lets create the LWC string
        % first lets check to see if there is a decimal point in our string
        num = round(lwc_z(xx),5);
        a = floor(num);
        b = num-a; % this tells us if our number is a fraction, or a whole number
        c = str2double(num2str(num)); % special case where num2str rounds up, but a is a whole integer diferrent
        
        if b==0
            newExpr_lwc_z = [num2str(num),'.0']; % if true, then there is no deicaml point in our string
        elseif (c-a)==1
            % if this is true, then we need to add a decimal point because
            % num2str rounded off and removed it
            newExpr_re = [num2str(num),'.0']; % if true, then there is no deicaml point in our string
        else
            
            newExpr_lwc_z = num2str(num);
        end
        
        % for LWC, we get 5 significant digits, or a string length of 6
        
        % --- check the length of the LWC string ---
        
        
        if length(newExpr_lwc_z)>1
            
            while length(newExpr_lwc_z)<6
                newExpr_lwc_z = [newExpr_lwc_z,'0'];
            end
            
            while length(newExpr_lwc_z)>6
                newExpr_lwc_z = newExpr_lwc_z(1:end-1);
            end
            
        elseif isempty(newExpr_lwc_z) == true
            error('String you are trying to write doesnt exist')
        end
        
        % lets check to make sure the lwc_z is within an acceptable range
        if str2double(newExpr_lwc_z)<0 || str2double(newExpr_lwc_z)>10
            
            disp([newline,'lwc_string value: ',newExpr_lwc_z,'. Layer: ',num2str(xx)])
            error('The value for lwc is either negative or far too large! Something happened when converting to a string')
        end
        
        % lets create the droplet size string
        % first lets check to see if there is a decimal point in our string
        num = round(re_z(xx),5);
        a = floor(num);
        b = num-a; % this tells us if our number is a fraction, or a whole number
        c = str2double(num2str(num)); % special case where num2str rounds up, but a is a whole integer diferrent
        if b==0
            newExpr_re = [num2str(num),'.0']; % if true, then there is no deicaml point in our string
        elseif (c-a)==1
            % if this is true, then we need to add a decimal point because
            % num2str rounded off and removed it
            newExpr_re = [num2str(num),'.0']; % if true, then there is no deicaml point in our string
            
        else
            
            newExpr_re = num2str(num);
        end
        
        
        
        % --- check the length of the droplet size string ---
        
        if length(newExpr_re)>1
            
            while length(newExpr_re)<6
                newExpr_re = [newExpr_re,'0'];
            end
            
            while length(newExpr_re)>6
                newExpr_re = newExpr_re(1:end-1);
            end
            
        elseif isempty(newExpr_re) == true
            error('String you are trying to write doesnt exist')
        end
        
        % quickly check the value, and make sure it is within the bounds of
        % the acceptable droplet sizes for the Hu and Stamnes
        % parameterization
        
        
        if str2double(newExpr_re)<2.5 || str2double(newExpr_re)>60
            
            disp([newline,'re_string value: ',newExpr_re,'. Layer: ',num2str(xx)])
            error('The value for r_e is either too small or too large! Something happened when converting to a string')
        end
        
        
        newExpr{xx} = [oldExpr{xx}(1:14),newExpr_lwc_z,'  ',newExpr_re];
        
    end
    
    % now that we have allocated all of the new values for LWC and re, we
    % can create a new water cloud file
    
    edit_INP_DAT_files(oldFolder,newFolder,oldFile,newFile,oldExpr,newExpr);
    
    
    
    
else
    
    error('I dont understand what kind of droplet profile you want')
    
end










end

