%% Plot Solar Disks at Multiple Zenith Angles

% By Andrew John Buggee
%%
scriptPlotting_blk
clear variables
close all

% ----- Introduce the CIE Color Matching Functions -----

% x is the red (long) cell response
% y is the green (medium) cell response
% z is the blue (short) cell response

[lambda, x_func, y_func, z_func] = colorMatchFcn('CIE_1931');



%% Lets read in Earths transmission data
fid = fopen('earth_transmission_std_atm_clear_sky.txt','r');

% define the format spec for the headers
format_spec_headers = '%s %s';

% just read the first row
h = textscan(fid, format_spec_headers, 1, 'CommentStyle', '#');

% define format_spec for two floating point column data
format_spec = '%f %f';

d = textscan(fid, format_spec, 'CommentStyle','#');

% all done! Close the file
fclose(fid);

% lets read in the solar source file

filename = 'kurudz_0.1nm.dat';

% define the wavelength vector
% I only want to look at values within the boundaries of the CIE color
% matching functions

wl = d{1}(d{1}<=max(lambda)&d{1}>=min(lambda))';                      % nm

T = d{2}(d{1}<=max(lambda)&d{1}>=min(lambda))';



% Lets interpolate the XYZ color functions so the have the same resolution
% as our transmission spectrum
x_funcT = interp1(lambda,x_func,wl);
y_funcT = interp1(lambda,y_func,wl);
z_funcT = interp1(lambda,z_func,wl);


%% Define the solar zenith angle

theta_sun_deg = [0:10:80,81:88,88.5,89,89.1:0.2:89.9];       % deg - cannot use 90 degrees. 
theta_sun_rad = theta_sun_deg*(pi/180);                 % rad

%f1 = figure(1);

% plot each solar disk on the same plot
f2 = figure(2);
subplot(4,6,1)
pos = [2 4 2 2];

for ii = 1:length(theta_sun_rad)
    
    % compute the transmission spectrum at several difference solar zenith
    % angles
    T_theta = T.^(1./cos(theta_sun_rad(ii)));
    
    % Lets look at the spectrum
    
    % Lets determine the CIE-XYZ color values
    X = trapz(wl, x_funcT.*T_theta);
    Y = trapz(wl, y_funcT.*T_theta);
    Z = trapz(wl, z_funcT.*T_theta);


    % We know normalize these values
    xyz = [X/(X+Y+Z), Y/(X+Y+Z), Z/(X+Y+Z)];

    % now we convert CIE XYZ values to sRBG values
    %-- Try different color spaces---
    %       'srgb'
    %       'adobe-rgb-1998'

    %-- Try different reference white points ---
    %       'd50' - color temp of 5003 K
    %       'd55' - color temp of 5500 K
    %       'd65' - color temp of 6504 K
    
    
    % -------------------------------------------------------
    % ---- Calculate RGB Using MATLAB Internal function -----
    %RGB = xyz2rgb(xyz, 'OutputType','uint16', 'ColorSpace', 'adobe-rgb-1998', 'WhitePoint',whitePoint);
    
    % -------------------------------------------------------
    % ----- Calculate RGB Using Inverted Color Matrix -------
    % How do I alter the reference white point?
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
    xlabel([num2str(theta_sun_deg(ii)), '$$ ^{\circ} $$'], 'interpreter','latex','color','w');

    
end

% insert Text box

dim = [0.12 0.7 0.3 0.3];
str = {'Color of the Solar Disk at Different Zenith Angles'};
annotation('textbox',dim,'String',str,'FitBoxToText','on', 'EdgeColor','none',...
    'Color',[1,1,1],'interpreter','latex','FontSize',30);


