%% Simple plane-parallel model

% Andrew John Buggee
%% Start with the simple two-layer model

clear variables

con = physical_constants;

% define the solar constant
S0 = con.sigma_sb * con.T_sun^4 * (con.R_sun/con.au)^2;   % W/m^2 - higher than the measured value (1368 > 1361)
S0 = con.S0;   % W/m^2 - 1361 measured value

% define the average Earth Albedo
A = 0.3;

% define the average emissivity of the atmosphere in the longwave region
% (primarily water vapor and Oxygen?)
eps = 0.2; 

% compute the surface temperature for a two layer model
% T_0 = @(eps) ( S0 * (1-A) ./ (4 * con.sigma_sb * (3/2 - eps)) ).^(1/4);
T_0 = @(eps) ( S0 * (1-A) ./ (4 * con.sigma_sb * (1 - eps./2)) ).^(1/4);    % Correct equation

eps1 = linspace(0, 1, 100);

figure; plot(eps1, T_0(eps1));
xlabel('$\epsilon_{LW}$','interpreter','latex','FontSize',35);
ylabel('$T_{surface}$ $(K)$','Interpreter','latex')
grid on; grid minor; hold on;
% set figure size
set(gcf,'Position',[0 0 700 810])
title('Simple two-layer model with surface and atm','Interpreter','latex',...
    'FontSize', 30)
