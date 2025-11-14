% plot the radiosonde water vapor concentration along with 4 different
% atmospheric profiles


% Andrew John Buggee

function [] = plot_radiosonde_wvConcentration_with_US_STD_ATM(radiosonde_wv_concentration, radiosonde_alt)

% first plot the radiosonde data
figure; 
semilogx(radiosonde_wv_concentration, radiosonde_alt./1e3,  'Color', mySavedColors(61, 'fixed'));



% load the standard atm profile
atm_std = read_atmos_prof_data('afglus.dat');
hold on; 
plot(atm_std.H2O, atm_std.altitude, 'Color', mySavedColors(62, 'fixed'));

% load the standard atm profile
atm_trop = read_atmos_prof_data('afglt.dat');
hold on; 
plot(atm_trop.H2O, atm_trop.altitude, 'Color', mySavedColors(63, 'fixed'));

% load the standard atm profile
atm_mlSum = read_atmos_prof_data('afglms.dat');
hold on; 
plot(atm_mlSum.H2O, atm_mlSum.altitude, 'Color', mySavedColors(64, 'fixed'));

% load the standard atm profile
atm_mlWin = read_atmos_prof_data('afglmw.dat');
hold on; 
plot(atm_mlWin.H2O, atm_mlWin.altitude, 'Color', mySavedColors(65, 'fixed'));


legend('radisonde', 'US standard atm', 'tropical atm', 'mid-latitude summer',...
    'mid-latitude winter', 'location',...
    'best','Interpreter','latex', 'Location','best', 'FontSize', 20,...
             'Color', 'white', 'TextColor', 'k')

title('Comparison between radiosonde and model water vapor concentrations', ...
        'FontSize', 20, 'Interpreter', 'latex')

ylabel('Altitude $(km)$', 'FontSize', 20, 'Interpreter', 'latex')
xlabel('Water Vapor Concentration $(cm^{-3})$', 'FontSize', 20, 'Interpreter', 'latex')

grid on; grid minor






end
