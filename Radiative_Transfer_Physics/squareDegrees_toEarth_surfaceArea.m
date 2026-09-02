%% Convert square degrees to area on the surface of the Earth



function squareDegrees_toEarth_surfaceArea(lat, dTheta, dPhi)

con = physical_constants;

R_e = con.R_earth;   % meters

theta = linspace(0, dTheta, 100);
phi = linspace(lat, lat+dPhi, 100);

omega = R_e^2 * dTheta * (pi/180) * (cosd(lat) - cosd(lat + dPhi)) * (pi/180);  % steradians



end



dTheta = 0.25; dPhi = 0.25;

lat = linspace(0, 90 - dPhi, 200);
% lat = linspace(66, 90 - dPhi, 100);

A = (con.R_earth^2 / 1e6) .* (dTheta * pi/180) .* (sind(lat + dPhi) - sind(lat));

figure; semilogy(lat, A)
xlabel('Latitude', 'Interpreter','latex'); 
ylabel(['ERA5 Grid Cell Area ($km^2$)'], 'Interpreter','latex')
xline(66.6, 'Color', 'w', 'LineWidth', 2, 'DisplayName', 'Arctic Circle', ...
    'LabelHorizontalAlignment','left', 'LabelVerticalAlignment', 'bottom',...
    'Label', 'Arctic Circle', 'FontSize', 18, 'LineStyle',':')
grid on; grid minor