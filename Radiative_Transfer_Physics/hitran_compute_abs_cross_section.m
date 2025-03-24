%% Compute absorption cross sections from hitran data

% using .par file output from the hitran database, compute the absorption
% cross section [cm^2/molecule]

% ----- INPUTS -----

% (1) T - Temperature of the gas (K)

% (2) P - Pressure of the gas - (atm)

% (3) P_self - partial pressure of the gas in question - (atm)

% (4) wavelength_grid - (nm) - wavelength grid that the absorption cross
% section is computed for.

% (5) solution_type - this is a string with two options (Gharavi and Buckley 2004):
%       (a) 'full_integral' - this approach solves for the Voigt function
%       by integrating over a large and fine grid. This takes a while

%       (b) 'whitting' - this is an approximation that solves for the
%       absorption cross section with a Voigt line shape.

% ----- OUTPUTS -----


% (1) absoprtion - this is a strucutre containing...
%       absorption.cross_section - (cm^2) - absorption cross section
%       absorption.bulk_coefficient - (1/cm) - bulk absorption coefficient
%       absorption.wavenumbers - (cm^(-1)) - wavenumber grid


% By Andrew John Buggee

%%


function [absorption] = hitran_compute_abs_cross_section(hitran_file, T, P, P_self, wavelength_grid_output, solution_type)



% Determine which computer you're using
computer_name = whatComputer;

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(computer_name,'andrewbuggee')==true

    hitran_folder = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Radiative_Transfer_Physics/HiTran_data/';

elseif strcmp(computer_name,'anbu8374')==true

    hitran_folder = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Radiative_Transfer_Physics/HiTran_data/';


end


%% Check to see if the hitran file is a .par

if strcmp(hitran_file(end-2:end), 'par')==true

    % Then we will import the .par file using the function 'importhitran'
    % and save it as a .mat file

    lines = importhitran([hitran_folder, hitran_file]);
    save([hitran_folder, hitran_file(1:end-3), 'mat'], "lines");

else

    % if its a .mat file, simply load the saved data into the workspave
    load([hitran_folder, hitran_file])

end

%% Convert the wavelength grid to a wavenumber grid

% Convert the linearly spaced wavelength grid into a wavenumber grid
% Make sure the wavelength vector is in microns
% This is the wavenumber grid that the cross section and bulk absorption
% will be computed on
wavenumber_grid_output = 10^4 ./ (wavelength_grid_output./1e3);        % cm^(-1)

% Make sure wavenumber_master_centers is a column vector
if size(wavenumber_grid_output,2)>1 && size(wavenumber_grid_output,1)==1

    % take the transpose
    wavenumber_grid_output = wavenumber_grid_output';

end

% ------------------------------------------------------------------------
% ------ Finding the Line Intensity closest to our output master grid ----
% ------------------------------------------------------------------------
% % Find the line centers closest to each value of the wavenumber grid
% [~, w_index] = min(abs(wavenumber_master_grid_output - lines.transitionWavenumber'), [], 2);
% 
% % w_index represents the closest hitran line center to each value on our
% % wavenumber grid
% % Keep only the unique values
% w_index = unique(w_index, 'stable');
% ------------------------------------------------------------------------

% ------------------------------------------------------------------------
% Find all line center transitions within the range of the desired output grid
% ------------------------------------------------------------------------
w_index = lines.transitionWavenumber>=wavenumber_grid_output(end) &...
    lines.transitionWavenumber<=wavenumber_grid_output(1);

% ------------------------------------------------------------------------

%% Define some constants and read the Total Internal Partition Sums

% define the reference temperature
T_ref = 296;            % kelvin

% define the radiometric constants
c2 = 1.4387769;      % cm*K

% Read in the total internal partition Q at the reference temperature for
% the molecule and isotopologue in question
Q = read_total_internal_partition_sum(lines.moleculeNumber(1), lines.isotopologueNumber(1), T);

%% Compute the line intensity

% the line intensity is defined at a reference temperature of 296K
S0 = lines.lineIntensity(w_index);                 % cm^(-1)/(cm^2 molecule)

% grab the wavenumbers of the desired wavelength range
wavenumber_line_center = lines.transitionWavenumber(w_index);      % cm^-1

% do I need to pressure shift the wavenumbers for the line strength
% calculation?
%wavenumber_shifted = wavenumber + lines.pressureShift * P;      % cm^-1

% Grab the lower energy energy state of the transition
E_lower = lines.lowerStateEnergy(w_index);         % cm^-1

% Compute the line intensity at some new temperature
S = S0 .* (Q.ref * exp(-c2 * E_lower./T) .* (1 - exp(-c2*wavenumber_line_center./T))) ./...
    (Q.T .* exp(-c2 * E_lower./T_ref) .* (1 - exp(-c2*wavenumber_line_center./T_ref)));        % cm^(-1)/(cm^2 molecule)


%% Compute the absorption cross section using one of the defined solution types

% If the solution type is 'full_integral', solve the Voigt function

if strcmp(solution_type, 'full_integral')==true

    % Create a Voigt Lineshape
    tic
    voigt = voigt_lineShape_for_hitran(lines, wavelength_grid_output, T, P, P_self);
    toc

    % Compute the Absorption cross-section

    % the aborption cross section is the product of the line strength and the
    % voigt line shape
    abs_cross_sec = zeros(length(S), length(voigt(1).shape));

    % We want to sum over all voigt profiles to have one continuous absorption
    % cross section as a function of wavelength.

    % To do this we will discretize the absorption cross section matrix into
    % wavenumber bin and add the cross sections in the same bin

    % ---------------------- Don't need now -----------------------------------
    %     % Define a master wavenumber vector that can be applied to the entire data
    %     % sat
    %     % Define the steps in terms of nanometers
    %     % Convert 1 nanometer difference to inverse centimeters
    %     d_nu = 10^4 / ((1400.0)/1e3) - 10^4 / ((1400.1)/1e3);        % cm^(-1)
    %
    %     wavenumber_master_edges = (voigt(1).wavenum(1)-(voigt(1).wavenum(1)/1e3)):...
    %         d_nu:(voigt(end).wavenum(end)+(voigt(end).wavenum(end)/1e3));       % cm^(-1)
    %     wavenumber_master_centers = wavenumber_master_edges(1:end-1) + d_nu/2;      % cm^(-1)
    % -------------------------------------------------------------------------

    wavenumber_master_edges = [wavenumber_grid_output(1:end-1) - diff(wavenumber_grid_output)./2,...
        wavenumber_grid_output(end) + (wavenumber_grid_output(end)-wavenumber_grid_output(end-1))];   % cm^-1

    % create a master absorption cross section vector
    abs_cross_sec_master_internal = zeros(size(wavenumber_master_edges));      % cm^2


    for vv = 1:length(S)

        % discretize each voigt profile using the wavenumber_master vector
        D = discretize(voigt(vv).wavenum, wavenumber_master_edges);

        abs_cross_sec(vv,:) = S(vv) .* voigt(vv).shape;         % cm^(2)/molecule

        abs_cross_sec_master_internal(D) = abs_cross_sec_master_internal(D) + abs_cross_sec(vv,:);


    end



elseif strcmp(solution_type, 'whitting')==true

    % Use the Whitting approximation. According to Ghavari and Buckley
    % 2004: The Voigt integral is difficult to evaluate, and simplied
    % expressions are typically used. The expression given by Whitting
    % (E.Whitting 1968) is often used for approximating the Voigt function.

    % ------------------- Need a master grid to sum over -----------------------
    % Create a master wavelength grid to compute the absorption cross
    % section on
%     d_lambda = 1;     % nm
%     wavelength_master_grid_internal = wavelength_grid(1): d_lambda : wavelength_grid(end);    % nm
% 
%     % Convert this linearly spaced wavelength grid into a wavenumber grid
%     % Make sure the wavelength vector is in microns
%     wavenumber_master_grid_internal = 10^4 ./ (wavelength_master_grid_internal./1e3);        % cm^(-1)
    % -------------------------------------------------------------------------

    % ---------------- Create an internal wavenumber grid ----------------
    % The Voigt absorption cross section is computed on this grid
    % The values of the absorption cross section defined by the user input,
    % the desired wavelengths, are determined from this grid
    d_lambda = 0.01;     % nm
    wavelength_master_grid_internal = (wavelength_grid_output(1) - d_lambda): d_lambda :...
        (wavelength_grid_output(end) + d_lambda);    % nm

    % Convert this linearly spaced wavelength grid into a wavenumber grid
    % Make sure the wavelength vector is in microns
    wavenumber_master_grid_internal = 10^4 ./ (wavelength_master_grid_internal./1e3);        % cm^(-1)



    % We need the doppler and pressure line broadening Half-Width-at-Half-Max

    % --------------------------------------------------------------
    % ----------- First compute the doppler shifted HWHM -----------
    % --------------------------------------------------------------
    % We need the molar mass of the isotopologue
    molar_mass = read_isotopologue_molar_mass_hitran(lines);     % g/mol
    % convert molar_mass to kg/mol
    molar_mass = molar_mass/1000;                           % kg/mol
    
    % load physical constants
    con = physical_constants();
    
    % compute the doppler profile half-width at half-max
    % use SI units for everything, except leave the spectral dimension in
    % wavenumbers with units of inverse centimeters
    f_doppler.hwhm = wavenumber_line_center./con.c .* sqrt((2*con.N_A * con.k_B *...
        T * log(2))/molar_mass);  % cm^(-1)
    
    % --------------------------------------------------------------
    % ----------- Next compute the Pressure shifted HWHM -----------
    % --------------------------------------------------------------

    % define the hitran reference temperature
    T_ref = 296;        % K
    
    % compute the pressure broadened half-width at half-max
    f_pressure.hwhm = (T_ref./T).^lines.temperatureDependence(w_index) .*...
        (lines.airBroadenedWidth(w_index) .* (P - P_self) +...
        lines.selfBroadenedWidth(w_index) .* P_self);



    % Now we compute the Voigt profile HWHM
    voigt_hwhm = 0.5346*f_pressure.hwhm + sqrt(0.2166*f_pressure.hwhm.^2 + f_doppler.hwhm.^2);

    % Compute the ratio of the lorentz HWHM and the Voigt HWHM
    x = f_pressure.hwhm./voigt_hwhm;

    % Compute the Voigt absorption cross section at the line center
    abs_cross_sec_center = S./(2*voigt_hwhm .*(1.065 + 0.447*x + 0.058*x.^2));      % cm^2 / molec

    % Now compute the full absorption cross section based on the Voigt profile
    % This is computed on the internal wavenumber grid
    abs_cross_sec_master_internal = zeros(size(wavenumber_master_grid_internal));

    parfor vv = 1:length(wavenumber_line_center)

        % Compute the Voigt HWHM scaled line transition shift
        y = abs(wavenumber_master_grid_internal - wavenumber_line_center(vv))./voigt_hwhm(vv);

        % Compute the absorption cross section on the master grid
        abs_cross_sec_master_internal = abs_cross_sec_master_internal + (abs_cross_sec_center(vv) .* ( (1-x(vv)) .* exp(-0.693*y.^2) +...
            x(vv)./(1+y.^2) + 0.016*(1-x(vv)).*x(vv) .* (exp(-0.0841 * y.^(2.25)) -...
            1./(1+0.0210*y.^(2.25))) ));       % cm^2 / molec

    end




else

    error([newline, 'solution_type intput not recognized', newline])

end



%% Mutiply by the number density and compute the bulk absorption coefficient

% We need to compute the homogenous number density for the given Pressure
% and Tempertaure

% read in physical constants
if exist('con','var')==false
    con = physical_constants();
end

avogadros_number = con.N_A;
pascals_per_atm = con.atm;
gas_constant_universal = con.R_uni;

Nd = (P_self*pascals_per_atm) * avogadros_number / (gas_constant_universal * T);    % molecules/m^3

% ******** IMPORTANT **********
% We don't need to compute the number density in terms of molecule
% /cm^3. That's because, according to the Hitran website, "The
% (line) intensity is defined here for a single molecule, per
% unit volume." Assuming SI, we multiply the absorption cross section by the
% number of molecules in cubic meter.

% % Convert the number density to molecules/cm^3
% Nd_cm = Nd / 1e6;       % molecules/cm^3
% ================================

% ------- Do I need to do this?? ---------
% % Convert the absorption cross section from cm^2/molecule to just cm^2
% abs_cross_sec_master = abs_cross_sec_master .* Nd_cm;       % cm^2
% ----------------------------------------

% Compute the bulk absorption coefficient
k_bulk_internal = abs_cross_sec_master_internal .* Nd;         % cm^(-1)



%% Corral the parameters of the absorption structure

% *** Read only the values desired by the user defined wavelength grid ***

% Perform linear interpolation to find the values desired by the user input

absorption.bulk_coefficient = interp1(wavenumber_master_grid_internal, k_bulk_internal,...
    wavenumber_grid_output, 'linear');
absorption.cross_section = interp1(wavenumber_master_grid_internal, abs_cross_sec_master_internal,...
    wavenumber_grid_output, 'linear');
absorption.wavenumber = wavenumber_grid_output;
absorption.wavelength = wavelength_grid_output;

end