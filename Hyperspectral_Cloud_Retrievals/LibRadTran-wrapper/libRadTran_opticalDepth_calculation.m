%% Determining how LibRadTran estimates optical depth at each layer

% By Andrew John Buggee

%%

% wc.sol.mol.cdf is a matrix of extinction coefficients with the following
% range of independent variables:
%   r_eff: 1 - 25 microns
%      wl: 0.250 - 2.2 microns

% Load the extinction coefficient file wc.sol.mol.cdf
info = ncinfo(['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
    'LibRadTran/libRadtran-2.0.4/data/wc/mie/wc.sol.mie.cdf']);
ext = ncread(['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
    'LibRadTran/libRadtran-2.0.4/data/wc/mie/wc.sol.mie.cdf'], 'ext');

% ext is a matrix that varies by effective radius and wavelength
% load these two independent variables
reff = ncread(['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
    'LibRadTran/libRadtran-2.0.4/data/wc/mie/wc.sol.mie.cdf'], 'reff');
wl = ncread(['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
    'LibRadTran/libRadtran-2.0.4/data/wc/mie/wc.sol.mie.cdf'], 'wavelen');


% plot the extinction coefficient 
figure; imagesc(wl,reff,ext);colorbar; xlabel('wavelength (microns)');ylabel('r_e (microns)')

%% Testing to see if this table matches values in the error file

% create matrix values that define the (wl,reff) grid
[wl_mat,reff_mat] = meshgrid(wl,reff);

% define the quarry points where interpolation occurs
r_q = 12.12;                   % microns
wl_q = 2130;                % nm


% Both variables have to be in microns
% interpolate to get the extinction coefficient at the points of interest
ext_q = interp2(wl_mat,reff_mat,ext,wl_q/1e3,r_q,'linear');                 % km^(-1)/(g/m^3)

% Tau is the product of the ext coefficient, the LWC and the geometric
% depth of the layer 
lwc = 0.3003;               % g/m^3
H = 0.5;                    % km

% compute tau
tau_LRT = ext_q *lwc*H;                 % unitless 

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

% Computing tau with my own precomputed mie table
justQ = false;
mie = interp_mie_computed_tables([wl_q, r_q], 'mono',justQ);

% define the variance of the distribution
dist_var = 7;
% define the index of refraction
index_of_refraction = 'water';
% define the size distribution
size_distribution = 'gamma';

% integrate over a size distribution to get an average
[ssa_avg, Qe_avg, g_avg] = average_mie_over_size_distribution(mie(6), mie(7), mie(5), r_q, dist_var, wl_q,...
    index_of_refraction, size_distribution);

% define the density of liquid water 
rho = 1e6;                  % grams/m^3 - density of liquid water

% Guess at optical depth using constant assumption over layer
tau_myGuess = 3 * lwc * Qe_avg * (H*1e3)/(4 * (r_q*1e-6) * rho);


%% 

% wc.trm.mol.cdf is a matrix of extinction coefficients with the following
% range of independent variables:
%   r_eff: 1 - 25 microns
%      wl: 2.21 - 100 microns

clear variables

% Load the extinction coefficient file wc.trm.mol.cdf
info = ncinfo(['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
    'LibRadTran/libRadtran-2.0.4/data/wc/mie/wc.trm.mie.cdf']);
ext = ncread(['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
    'LibRadTran/libRadtran-2.0.4/data/wc/mie/wc.trm.mie.cdf'], 'ext');

% ext is a matrix that varies by effective radius and wavelength
% load these two independent variables
reff = ncread(['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
    'LibRadTran/libRadtran-2.0.4/data/wc/mie/wc.trm.mie.cdf'], 'reff');
wl = ncread(['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
    'LibRadTran/libRadtran-2.0.4/data/wc/mie/wc.trm.mie.cdf'], 'wavelen');

% Plot the extinction coefficient
figure; imagesc(wl,reff,ext);colorbar; xlabel('wavelength (microns)');ylabel('r_e (microns)')


%%  Plot the results comparing the input tau and libRadTrans computed tau

figure; 
for nn=1:12
    plot(tau_c(nn), tau_lrt(nn),'.','markersize',35); hold on;
    legend_str{nn} = ['$\lambda = $',num2str(wavelength(nn)),' $nm$'];
end

grid on; grid minor;
one2one = linspace(0, 1.1*max([tau_c,tau_lrt],[],'all'),100);
plot(one2one,one2one,'k-','Linewidth',1)
legend(legend_str, 'Interpreter','latex','Location','best')
xlabel('Starting $\tau_c$','Interpreter','latex');
ylabel('LibRadTran $\tau_c$', 'Interpreter', 'latex')
ylim([0, 1.1*max([tau_c,tau_lrt],[],'all')])
xlim([0, 1.1*max([tau_c,tau_lrt],[],'all')])
set(gcf,'position', [0 0 700 700])



