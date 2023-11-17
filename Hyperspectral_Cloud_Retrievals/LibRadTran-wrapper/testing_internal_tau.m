%% Trying to match the internal optical depth calculation with my input

% Andrew John Buggee
%% TESTING TO SEE IF USING THE NETCDF FILE CAN MATCH INTERNAL CALC
% Load files
wc_folder = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/',...
    'Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/data/wc/mie/'];
mie_file = 'wc.gamma_007.0.mie.cdf';

% read the wavelength, reff, and K
K = ncread([wc_folder, mie_file],'ext');        % km^(-1)/(g/m^3) - extinction coefficient per LWC
reff = ncread([wc_folder, mie_file],'reff');        % microns - effective radius
wl = ncread([wc_folder, mie_file],'wavelen');        % microns - wavelength

R = repmat(reff,1,length(wl));                      % meshgrid for reff
WL = repmat(wl',length(reff),1);                    % meshgrid for wwavelength

rho = 1e6;                                          % g/m^3 - density of water

Q_mat = 4*K.*(R*1e-9).*rho/3;

% estimate the extinction efficiency at the wavelength and effective radius
% of my choosing

r = 5;          % microns
wavelength = 0.655;     % microns

Qe = interp2(WL,R,Q_mat, wavelength,r);
K_i = interp2(WL,R,K,wavelength,r);

% estimate the LWC

H = 500;               % meters - depth of cloud layer
Tau = 20;               % optical depth

LWC = 4*(r*1e-6)*rho*Tau/(3*Qe*H)
T = K_i * LWC * H


%% TESTING TO SEE IF SAME CALC FOR TXT AND NETCDF PRODUCE SAME TAU

n_layers = 1;
re = linspace(5,5,n_layers);
tau_c = 23;
dist = 'mono';
H = 0.5;                            % km -  cloud geometric thickness
z0 = 1.5;                           % km - base of cloud
z_topBottom = [z0+H, z0];           % km - altitude above the ground
lambda = 645;                       % nm - wavelength that defines the optical depth

% create a water cloud file

write_wc_file(re,tau_c,z_topBottom, H, lambda, dist);

%% CREATE WATER CLOUD FILE WITH SPECIFIC DROPLET PROFILE

re_TB = [12,5];                         % microns - effective radius at cloud top and bottom
H = 0.5;                                % km - geometric thickness of cloud
n_layers = 5;                          % number of layers to model within cloud

z0 = 1;                                 % km - base height of cloud
z = linspace(z0, z0+H,n_layers);        % km - altitude above ground vector
indVar = 'altitude';                    % string that tells the code which independent variable we used
constraint = 'adiabatic';               % string that tells the code which physical constraint to use

re = create_droplet_profile2(re_TB, z, indVar, constraint);     % microns - effective radius vector

tau_c = 35;                             % optical depth
dist = 'gamma';                         % droplet distribution
lambda = 525;                           % nm - wavelength that defines the optical depth
z_topBottom = [z(end), z(1)];           % km - boundaries of the altitude vector. 

write_wc_file(re, tau_c, z_topBottom, H, lambda, dist);

%% Looking at the Mie calculation data

folder_path = '/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/Mie_Calculations/';

filename = 'Mie_calcs_1nm_sampling_gamma_7.OUT';          % This file computes re=1:1:100 and lambda=100:1:3000

format_spec = '%f %f %f %f %f %f %f %f';        % 8 columns of data

file_id = fopen([folder_path,filename]);

data_raw = textscan(file_id, format_spec, 'CommentStyle','#', 'Delimiter',...
    {'\r', ' '}, 'MultipleDelimsAsOne',1);

% Grab the Qe matrix
wl = unique(data_raw{1});
re = unique(data_raw{2});
Qe = reshape(data_raw{5}, length(unique(data_raw{2})), []);


% Lets make a plot of Qe for 10 different radii 


figure; semilogy(wl, Qe(1:10:end, :))
xlabel('Wavelength (nm)'); ylabel('Q_{ext}')
legend(strcat(string(re(1:10:end)), ' \mu m'), 'Location', 'best');
grid on; grid minor;





