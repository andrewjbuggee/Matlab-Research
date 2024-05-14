%% Compute the Voigt spectral line shape using hitran data



% By Andrew John Buggee

%%

function voigt = voigt_lineShape_for_hitran(hitran_lines, T, P, P_self)


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

f_doppler = gaussian_lineShape_for_hitran(hitran_lines, T);

%% We also need the Pressure Broadened line shape for each energy transition

f_pressure = lorentz_lineShape_for_hitran(hitran_lines, T, P, P_self);

%% The Voigt Lineshape is the convolution of the doppler broadened and pressure broadened line shapes

% compute the Doppler line shape for each line transition
voigt.shape = zeros(size(f_pressure.shape,1), size(f_pressure.shape,2));
voigt.wavenum = zeros(size(f_pressure.wavenum,1), size(f_pressure.shape,2));

%idx_half_width = floor(median(1:size(f_doppler.shape,2)));

for vv = 1:size(f_doppler,1)

%     half_width_wavenum = f_pressure.wavenum(vv,idx_median) - f_pressure.wavenum(vv,1);       % cm^(-1)
%     voigt.wavenumbers(vv,:) = [f_pressure.wavenum(vv,1:idx_half_width)-half_width_wavenum,...
%         f_pressure.wavenum(vv,1:idx_half_width+1)+half_width_wavenum];     % cm^(-1) - wavenumbers

    voigt.wavenum(vv,:) = f_pressure.wavenum(vv,:);     % cm^(-1) - wavenumbers
    % Create a lineshape that is the same size as the pressure broadened
    % line shape
    voigt.shape(vv,:) = conv(f_pressure.shape(vv,:), f_doppler.shape(vv,:), "same");  % (cm^(-1))^(-1) - this is a PDF

    % Testing my own convolution
    for kk = 1:length(f_pressure.wavenum(vv,:))
        for jj = 1:length(f_doppler.wavenum(vv,:))
            while kk-jj+1 > 0
                
                shape_sum = f_doppler.shape(vv,jj) * f_pressure.shape(vv,jj-kk+1);

            end
        end
    end
                voigt.shape(vv,kk) = (f_pressure.wavenum(1,end)-f_pressure.wavenum(2,1))

end


end