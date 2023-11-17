%% Plot Multiple Planck Emitters in their True RGB Color

% By Andrew John Buggee
%% 
clear variables
close all

scriptPlotting_blk

% Grab the color matching functions
[lambda, x_func, y_func, z_func] = colorMatchFcn('1931_full');


% Lets introduce constants
con = physical_constants();


T_eff = 500:500:12000;                 % K - effective temperature of a star

% f1 = figure(1);

% plot each solar disk on the same plot
f2 = figure(2);
subplot(4,6,1)
pos = [2 4 2 2];

for ii = 1:length(T_eff)
    
    L = plancks_function(lambda, T_eff(ii), 'nanometers');        % W/m^2/nm/sr
    
    % Lets convert radiance to irradiance at the top of the atmosphere for
    % Earth
    F = L* (con.R_sun/con.au)^2 *pi;                              % W/m^2/nm - spectral flux
    
    % Lets look at the spectrum
%     figure(1); hold on;
%     semilogy(lambda, F); grid on; grid minor
    
    % Lets determine the CIE-XYZ color values
    X = trapz(lambda, x_func.*F);
    Y = trapz(lambda, y_func.*F);
    Z = trapz(lambda, z_func.*F);


    % We know normalize these values
    xyz = [X/(X+Y+Z), Y/(X+Y+Z), Z/(X+Y+Z)];


    % now we convert CIE XYZ values to sRBG values
    % Some whitepoint references are:
    %   'a' - CIE standard illuminant A - simulates tungsten filament
    %   'd65' - CIE standard illuminant - simulates noon daylight
    %   'icc' - prgile connection space illumiant
    whitePoint = 'd65';
    %whitePoint = [0.5,0.5,0.5];
    
    % -------------------------------------------------------
    % ---- Calculate RGB Using MATLAB Internal function -----
    %RGB = xyz2rgb(xyz, 'OutputType','uint16', 'ColorSpace', 'adobe-rgb-1998', 'WhitePoint',whitePoint);
    
    % -------------------------------------------------------
    % ----- Calculate RGB Using Inverted Color Matrix -------
    color_mat = [0.41847, -0.15866, -0.082835;...
                 -0.091169, 0.25243, 0.015708;...
                 0.00092090, -0.0025498, 0.17860];
    RGB = color_mat*xyz';
    % Convert to 8-bit image
    RGB = uint8(255*(RGB./max(RGB)));
    % -------------------------------------------------------
    
    % Lets plot the solar disk
    figure(2); subplot(4,6,ii)
    rectangle('Position',pos,'Curvature',[1 1], 'FaceColor',RGB, 'EdgeColor',[0,0,0])
    set(gca,'XColor','k')
    set(gca,'YColor','k')
    axis equal
    title([num2str(T_eff(ii)), ' K'],'color','w');

    
end

% Add labels to figure 1

% figure(1)
% xlabel('Wavelength (nm)'); ylabel('Radiance (W/m^2/nm/sr)')
% title('Planck Spectrum')
% legend(string(T_eff),'Location','best')


% Add labels to figure 2
figure(2)




