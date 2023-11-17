%% ----- Compute the extinction cross section for a liquid water droplet -----

% wavelength has to be in microns, and is a vector
% r_eff is a vector of the same length as wavelength, units of microns

% By Andrew J. Buggee
%%


function [ext_coeff] = compute_ext_coeff(r_eff,wavelength)

% ensure that both inputs are positive
if sum(r_eff<0)>0
    error('the input, r_eff, is below 0, which is not allowed')
    
elseif sum(wavelength<0)>0
    error('the input, wavelength, is below 0, which is not allowed')
end


% ----- What computer are you using? -----

computer_name = whatComputer;

if strcmp(computer_name,'anbu8374')
    mie_fileName = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/wc/mie/wc.sol.mie.cdf';
elseif strcmp(computer_name,'andrewbuggee')
    error('I dont know where the file is!')
end

% ----- Read in the netcdf precomputed Mie tables from uvspec -----

ext_lwc = ncread(mie_fileName,'ext'); % km^(-1)/(g/m^3) - extinction coefficient
ssa = ncread(mie_fileName,'ssa');
reff = ncread(mie_fileName,'reff');
wavelen = ncread(mie_fileName,'wavelen');

% ------------------------------------------------------------------



rho_l = 1e6; % g/m^3 - liquid water content
lwc_single_drop = 4/3*pi*(r_eff*10^(-6)).^3 * rho_l;   % g/m^3 - lwc for a single droplet

% ----- create mesh grids for wavelenth and effective radius grid for the lookup table ------

[W,R] = meshgrid(wavelen,reff);


% ----- If the user inputs are within the bounds of the pre-computed
% tables, we interpolate -----

index_r = r_eff>=min(reff) & r_eff<=max(reff);
index_w = wavelength>=min(wavelen) & wavelength<=max(wavelen);


% ----- Calculate each extinction coefficient -----

ext_coeff = zeros(1,length(r_eff));  % create a zero vector

for ii = 1:length(lwc_single_drop)
    
    
    
    % change units to m^(-1)
    ext = ext_lwc.*lwc_single_drop(ii) .*(1/1000); % m^(-1) - extinction coefficient

    
    % ----- create mesh grids for wavelenth and effective radius designated by the user ------

    [Wq,Rq] = meshgrid(wavelength(ii),r_eff(ii));


    if index_r(ii)==true && index_w(ii)==true
        % if ture, then we interpolate
        ext_coeff(ii) = interp2(W,R,ext,Wq,Rq);
    
    else
        % if these conditions aren't true, we extrapolate
        ext_coeff(ii) = nan;
    
    end

end





end