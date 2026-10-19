function c = kohler_constants()
% KOHLER_CONSTANTS  Constants for kappa-Kohler theory, as in Jeramy Dedrick's kappakohler.m.
%
%   BASED ON JERAMY L. DEDRICK'S GITHUB CODE: the constants of his kappakohler.m
%   (https://github.com/jdedrick95/kappakohler, 2024).
%   Please cite his code and Dedrick et al. (2026, J. Atmos. Sci., doi:10.1175/JAS-D-26-0004.1) if you use
%   it. His repository has no licence: ask him before sharing this function outside the group.
%
%   c = kohler_constants()
%
%   Values from Seinfeld & Pandis (2016), as used in kappakohler.m
%   (github.com/jdedrick95/kappakohler). The temperature is the CCN counter's
%   (298.15 K), because these functions describe activation in the instrument
%   and in kappa-Kohler retrievals that use 298 K.
%
%   Sources: the surface tension of water is about 7.5e-2 N m-1 at meteorological
%   temperatures (Rogers & Yau 1989, p. 85). Their linear fit sigma = -1.55e-4 T + 0.118
%   N m-1 (T in K; p. 98, problem 6.9) gives 0.0718 N m-1 at 298 K; Petters &
%   Kreidenweis (2007, p. 1964) used 0.072 J m-2; kappakohler.m uses 0.0728. The
%   differences change A by about 1.4 % and s_c by about 2 % (my calculation). A is
%   the DIAMETER form of the Kelvin coefficient used by Petters & Kreidenweis (2007,
%   p. 1962, eq. 6; A defined on p. 1964 with eq. 10). The RADIUS form,
%   a = 2 sigma / (rho_L R_v T) = A / 2, is the one in Rogers & Yau (1989, p. 88, eq. 6.6).
%
%   Output struct (SI units):
%     c.sigma_w_Nm   surface tension of pure water, N m-1 (kappakohler.m's comment says N/m^2; the value is N m-1)
%     c.M_w_kgmol    molar mass of water, kg mol-1
%     c.R_Jmol_K     universal gas constant, J mol-1 K-1
%     c.rho_w_kgm3   density of liquid water, kg m-3
%     c.T_K          temperature used in the Kelvin term, K
%     c.A_m          Kelvin parameter A = 4 sigma M_w / (R T rho_w), m (~2.1e-9 m = 2.1 nm)
%                    (diameter form: the Kelvin term is exp(A / D) with D the wet DIAMETER)

c.sigma_w_Nm = 0.0728;
c.M_w_kgmol  = 0.018;
c.R_Jmol_K   = 8.314;
c.rho_w_kgm3 = 1000;
c.T_K        = 298.15;
c.A_m        = 4 * c.sigma_w_Nm * c.M_w_kgmol / (c.rho_w_kgm3 * c.R_Jmol_K * c.T_K);

end
