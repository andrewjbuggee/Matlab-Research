%% Compute the Voigt spectral line shape using hitran data



% By Andrew John Buggee

%%

function voigt = voigt_lineShape_for_hitran(hitran_lines, wavelength_grid, T, P, P_self)


% Determine which computer you're using
computer_name = whatComputer;

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(computer_name,'andrewbuggee')==true

    hitran_folder = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Radiative_Transfer_Physics/HiTran_data/';

elseif strcmp(computer_name,'anbu8374')==true

    hitran_folder = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Radiative_Transfer_Physics/HiTran_data/';


end


%% We need the Doppler Broadened line shape for each energy transition

f_doppler = gaussian_lineShape_for_hitran(hitran_lines, wavelength_grid, T);

%% We also need the Pressure Broadened line shape for each energy transition

f_pressure = lorentz_lineShape_for_hitran(hitran_lines, wavelength_grid, T, P, P_self);

%% The Voigt Lineshape is the convolution of the doppler broadened and pressure broadened line shapes

% Convert the wavelength grid to a wavenumber grid

% Convert this linearly spaced wavelength grid into a wavenumber grid
% Make sure the wavelength vector is in microns
wavenumber_master_centers = 10^4 ./ (wavelength_grid./1e3);        % cm^(-1)

% Find the line centers closest to each value of the wavenumber grid
[~, w_index] = min(abs(wavenumber_master_centers - hitran_lines.transitionWavenumber'), [], 2);

% Keep only the unique values
w_index = unique(w_index, 'stable');

% grab the the line strength centers
line_center = hitran_lines.transitionWavenumber(w_index);       % cm^(-1)


pressure_induced_line_shift = hitran_lines.pressureShift(w_index);

% compute the Doppler line shape for each line transition

% ------ for the same length convolution -------------
% voigt.shape = zeros(size(f_pressure.shape,1), size(f_pressure.shape,2));
% voigt.wavenum = zeros(size(f_pressure.wavenum,1), size(f_pressure.shape,2));

% ------ for the full length convolution -------------
% voigt.shape = zeros(size(f_pressure.shape,1), size(f_pressure.shape,2)+size(f_doppler.shape,2)-1);
% voigt.wavenum = zeros(size(f_pressure.wavenum,1), size(f_pressure.shape,2)+size(f_doppler.shape,2)-1);

% ------ for calculating the voigt profile directly -------------
% voigt = struct('wavenum', cell(size(f_pressure.wavenum,1),7*size(f_pressure.shape,2)),...
%     'shape', cell(size(f_pressure.wavenum,1),7*size(f_pressure.shape,2)));
voigt = struct('wavenum', cell(size(f_pressure.wavenum,1), 1),...
    'shape', cell(size(f_pressure.wavenum,1),1));

%voigt.shape = zeros(size(voigt.wavenum));

% -------------------------------------------------------------
% --------- Define variables outside the parfor loop ----------

% compute the wavenumber displacement
d_nu = (f_pressure.wavenum(:,2)-f_pressure.wavenum(:,1));     % cm^(-1)

% compute the pressure shifted line centers
pressure_shifted_line_center = line_center + pressure_induced_line_shift*P;   % cm^(-1)

% we need to first and last wavenumber values of pressure line shape wavenumber vectors
wavenumbers_pressureLines_boundaries = [f_pressure.wavenum(:,1), f_pressure.wavenum(:,end)];

% length of wavenumber vector for the Lorentz lineshape
num_wavenumbers_lorentz = size(f_pressure.shape,2);

% length of wavenumber vector for the voigt lineshape
num_wavenumbers_voigt = 7*size(f_pressure.shape,2);

% Store the doppler half-width-half-max
doppler_hwhm = f_doppler.hwhm;

% store the y-constant
y = sqrt(log(2)) * f_pressure.hwhm./doppler_hwhm;

% a finer spacing, dt, produces a smoother voigt profile
% but finer spacing leads to a larger vector which takes a long time to
% integrate.
t = linspace(-100, 100, 10000);

parfor vv = 1:size(voigt,1)


    % ----------------------------------------------------------
    % --- COMPUTE VOIGT FUNCTION USING MATLABS CONV FUNCTION ---
    % ----------------------------------------------------------

% ------ for the same length convolution -------------
%     voigt.wavenum(vv,:) = f_pressure.wavenum(vv,:);     % cm^(-1) - wavenumbers

    % ------ for the full length convolution -------------
%    d_nu = (f_pressure.wavenum(vv,2)-f_pressure.wavenum(vv,1));     % cm^(-1)
%     voigt.wavenum(vv,:) = linspace(f_pressure.wavenum(vv,1) - (d_nu* (size(f_pressure.shape,2)/2)),...
%         f_pressure.wavenum(vv,end) + (d_nu* (size(f_pressure.shape,2)/2 - 1)),...
%         size(f_pressure.shape,2)+size(f_doppler.shape,2)-1);

%    voigt.wavenum(vv,:) = f_pressure.wavenum(vv,1) - (d_nu * 3*size(f_pressure.shape,2)):d_nu:...
%        f_pressure.wavenum(vv,end) + (d_nu * 3*size(f_pressure.shape,2));
    
    
    % Create a lineshape that is the same size as the pressure broadened
    % line shape
    
    % *** The convultion function in matlab is just the convolution sum. To
    % compare with the convolution integral, multiply the result by the
    % discrete step inbetween value of the independent variable ***
%    voigt.shape(vv,:) = d_nu .* conv(f_pressure.shape(vv,:), f_doppler.shape(vv,:));  % (cm^(-1))^(-1) - this is a PDF

    % ----------------------------------------------------------
    % ----------------------------------------------------------



    % ----------------------------------------------
    % ---------- CODE BELOW DOESN'T WORK -----------
    % ----------------------------------------------
    % Testing my own convolution
%     for kk = 1:length(f_pressure.wavenum(vv,:))
%         for jj = 1:length(f_doppler.wavenum(vv,:))
%             while kk-jj+1 > 0
%                 
%                 shape_sum = f_doppler.shape(vv,jj) * f_pressure.shape(vv,jj-kk+1);
%                 
%                 voigt.shape(vv,kk) = (f_pressure.wavenum(1,end)-f_pressure.wavenum(2,1));
%                 
%             end
%         end
%     end
    % ----------------------------------------------
    % ----------------------------------------------


    % ----------------------------------------------
    % ---------- COMPUTE VOIGT DIRECTLY ------------
    % ----------------------------------------------
    % Usuing the formulas outlined in Gharavi and Buckley 2004

    
    
    % Define the wavenumber domain of the voigt line shape.
    % *** The larger the range of this vector, the closer the voigt line
    % shape becomes to a true PDF, i.e. integrates to 1 ***

%     voigt.wavenum(vv,:) = f_pressure.wavenum(vv,1) - (d_nu * 3*size(f_pressure.shape,2)):d_nu:...
%         f_pressure.wavenum(vv,end) + (d_nu * (3*size(f_pressure.shape,2) +1));
    voigt(vv).wavenum = linspace(wavenumbers_pressureLines_boundaries(vv,1) - (d_nu(vv) * 3*num_wavenumbers_lorentz),...
        wavenumbers_pressureLines_boundaries(vv,2) + (d_nu(vv) * 3*num_wavenumbers_lorentz), num_wavenumbers_voigt);

    x = sqrt(log(2))/doppler_hwhm(vv) .*...
        (voigt(vv).wavenum - pressure_shifted_line_center(vv));


    % w = x-t
    w = repmat(x,length(t), 1) - repmat(t', 1, length(x));
    % numerator is exp(-t^2)
    numerator = repmat(exp(-t.^2)',1, length(x));

    K = y(vv)./pi * trapz(t, numerator./(y(vv).^2 + (w).^2), 1);

    voigt(vv).shape = sqrt(log(2))/(sqrt(pi) * doppler_hwhm(vv)) .* K;

    % ----------------------------------------------
    % ----------------------------------------------


end


% ---- Can I do this all outside the for loop? ----
% ---- NOPE! ---


% pressure_shifted_line_center = transition_wavenumber + pressure_induced_line_shift*P;   % cm^(-1)
% 
% d_nu = (f_pressure.wavenum(1,2)-f_pressure.wavenum(1,1));     % cm^(-1)
% 
% for vv = 1:size(f_pressure.wavenum,1)
% 
%     voigt.wavenum(vv,:) = linspace(f_pressure.wavenum(vv,1) - (d_nu * 3*size(f_pressure.shape,2)),...
%         f_pressure.wavenum(vv,end) + (d_nu * 3*size(f_pressure.shape,2)), size(voigt.wavenum,2));
% end
% 
% x = sqrt(log(2))./f_doppler.hwhm .*...
%     (voigt.wavenum - repmat(pressure_shifted_line_center, 1, size(voigt.wavenum,2)));
% 
% y = sqrt(log(2)) .* (f_pressure.hwhm./f_doppler.hwhm);
% 
% % a finer spacing, dt, produces a smoother voigt profile
% t = linspace(-1000, 1000, 40000);
% 
% % w = x-t
% w = repmat(x,1,1, length(t)) - repmat(t', 1, length(x));
% % numerator is exp(-t^2)
% numerator = repmat(exp(-t.^2)',1, length(x));
% 
% K = y./pi * trapz(t, numerator./(y.^2 + (w).^2), 1);
% 
% voigt.shape(vv,:) = sqrt(log(2))/(sqrt(pi) * f_doppler.hwhm(vv)) .* K;


end