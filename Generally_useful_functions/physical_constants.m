%% Frequently Used Constants!

% Values taken from Astropy: https://docs.astropy.org/en/stable/constants/index.html

% By Andrew John Buggee
%%
function constant = physical_constants()

% universal gravitational constant
constant.G = 6.6743e-11;                % m^3/kg/s^2

% acceleration of gravity at the surface of Earth
constant.g0 = 9.80665;                % m/s^2



% ---- Chemistry stuff ---
% Avogadro's constant
constant.N_A = 6.02214076e+23;                % molecules/mol

% universal gas constant
constant.R_uni = 8.31446262;                % J/K/mol

% Dry air gas constant
constant.R_air_dry = 287.052874;               % J/K/Kg

% water vapor gas constant
constant.R_waterVapor = 461.52;               % J/(kg K) 

% molar mass of dry air
constant.Mol_mass_air_dry = 0.2897;                % kg/mol

% molar mass of water vapor
constant.Mol_mass_h2o_vap = 0.018015;                % kg/mol








% heat capacity of air at constant pressure
constant.cp_air = 1005;                 % J/kg/K

% heat capacity of liquid water at constant pressure - doesn't change much from 0
% to 100 C (engineering toolbox)
constant.cp_h2o = 4220;                 % J/kg/K

% heat of vaporization of water at 100 C (Engineering toolbox)
constant.H_h2o = 2256400;                 % J/kg

% Density of liquid water at 0 C (engineering toolbox)
constant.density_h2o_liquid = 999.87;             % kg/m^3

% adiabatic lapse rate
% constant.L_a = constant.g0/constant.cp_air;                % J/K/mol

% pressure of 1 standard atmosphere
constant.atm = 101325;                % 1 atm = Pa = N/m^2

% Weins Displacement constant
constant.b_wein = 0.00289777196;                % K m

% speed of light in a vacuum
constant.c = 299792458;                % m/s

% Planck's constant
constant.h = 6.62607015e-34;                % J s

% Boltzmann constant
constant.k_B = 1.380649e-23;                % J/K

% proton mass
constant.m_p = 1.67262192e-27;                % kg

% Stefan-Boltzmann constant
constant.sigma_sb = 5.67037442e-08;                % W/K^4/m^2

% 1 Atomic Mass Unit
constant.amu = 1.66053907e-27;                % kg/atom




% ----- Planet Mass ----
% Nominal Earth Mass
constant.M_earth = 5.97216787e+24;                % kg

% Nominal Jupiter Mass
constant.M_jup = 1.8981246e+27;                % kg

% Nominal Sun Mass
constant.M_sun = 1.98840987e+30;                % kg


% --- RADII ---
% Nominal Earth Radius
constant.R_earth = 6378100;                % m

% Nominal Earth Radius
constant.R_moon = 1737500;                % m

% Nominal Jupiter Radius
constant.R_jup = 71492000;                % kg

% Nominal Sun Radius
constant.R_sun = 695700000;                % kg



% 1 Astronomical unit - avg distance between the Earth and the Sun
constant.au = 1.49597871e+11;             % m

% avg distance between the Earth and Moon
constant.d_moon = 384990e3;             % m


% Effective temperature of the Sun
constant.T_sun = 5780;                     % K

% Solar Constant
constant.S0 = 1361;                        % W/m^2




end


