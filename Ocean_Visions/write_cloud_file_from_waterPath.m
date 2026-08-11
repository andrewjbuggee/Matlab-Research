%% Write a libRadtran 1D cloud file (wc_file / ic_file) from a column water path
%
% libRadtran cloud files are specified in terms of water content (g/m^3) and
% effective radius (micron) on a level grid -- NOT optical depth. For a single
% homogeneous layer the conversion from column water path is exact arithmetic:
%
%       WC [g/m^3] = WP [g/m^2] / H [m]
%
% so no Mie table, extinction efficiency, or reference wavelength is required.
% This is what makes a liquid-vs-ice comparison at fixed column water mass
% clean: both phases get *identical* condensed mass in an identical geometric
% layer, and the only difference is the single-scattering physics libRadtran
% applies (wc_properties vs ic_properties).
%
% *** File format (libRadtran manual, sec. 3.3.5 / wc_file) ***
% Three columns: level altitude (km), water content (g/m^3), r_eff (micron).
% Since libRadtran 1.4 these are LAYER quantities, and the values on a given
% row apply to the layer between that level and the level ABOVE it. The
% topmost row therefore carries zeros and exists only to close off the cloud
% top; without it the cloud would extend to the top of the atmosphere.
%
%
% INPUTS:
%   (1) water_path_gPerMeterSquared - column liquid or ice water path (g/m^2).
%       Single scalar; one cloud file per call.
%
%   (2) r_e_um - effective radius of the droplets / ice particles (micron).
%       Constant through the layer (vertically homogeneous cloud).
%
%   (3) z_topBottom_km - two-element vector [z_cloudTop_km, z_cloudBottom_km],
%       altitude above sea level. Follows the convention of write_wc_file.m in
%       the LibRadTran-wrapper folder: cloud top FIRST.
%
%   (4) phase_str - 'liquid' or 'ice'. Only affects the header comment and the
%       default file name; the numbers written are identical in form.
%
%   (5) folder_path - folder the .DAT file is written into. Created if absent.
%
%   (6) fileName - name of the .DAT file to write (including the extension).
%
%
% OUTPUTS:
%   (1) fileName - the name of the file written (echoed back for convenience)
%   (2) wc_gPerMeterCubed - the water content actually written (g/m^3)
%
%
% By Andrew John Buggee
%%

function [fileName, wc_gPerMeterCubed] = write_cloud_file_from_waterPath(water_path_gPerMeterSquared,...
    r_e_um, z_topBottom_km, phase_str, folder_path, fileName)

% ------------------------------------------------------------
% ---------------------- CHECK INPUTS ------------------------
% ------------------------------------------------------------

if nargin~=6
    error([newline, 'Need 6 inputs: column water path (g/m^2), effective radius (micron),',...
        ' [z_top, z_bottom] (km), phase string, folder path, and file name.', newline])
end

if numel(water_path_gPerMeterSquared)~=1 || water_path_gPerMeterSquared<=0
    error([newline, 'The column water path must be a single positive value (g/m^2).', newline])
end

if numel(r_e_um)~=1 || r_e_um<=0
    error([newline, 'The effective radius must be a single positive value (micron).', newline])
end

if numel(z_topBottom_km)~=2
    error([newline, 'z_topBottom_km must contain two values: [cloud top, cloud bottom] in km.', newline])
end

% enforce the [top, bottom] ordering used by write_wc_file.m
z_topBottom_km = reshape(z_topBottom_km, 1, []);
if z_topBottom_km(1) <= z_topBottom_km(2)
    error([newline, 'z_topBottom_km must be ordered [cloud top, cloud bottom], with the top higher.', newline])
end

if strcmp(phase_str, 'liquid')==false && strcmp(phase_str, 'ice')==false
    error([newline, 'phase_str must be either ''liquid'' or ''ice''.', newline])
end

if exist(folder_path, 'dir')==0
    mkdir(folder_path)
end

% make sure the folder path ends in a file separator
if folder_path(end)~=filesep
    folder_path = [folder_path, filesep];
end


% ------------------------------------------------------------
% ------------ CONVERT WATER PATH TO WATER CONTENT -----------
% ------------------------------------------------------------

H_m = (z_topBottom_km(1) - z_topBottom_km(2)) * 1e3;                 % m - geometric cloud thickness

wc_gPerMeterCubed = water_path_gPerMeterSquared / H_m;               % g/m^3 - liquid or ice water content


% ------------------------------------------------------------
% ------------------- WRITE THE .DAT FILE --------------------
% ------------------------------------------------------------

if strcmp(phase_str, 'liquid')==true
    content_header = 'LWC';
else
    content_header = 'IWC';
end

fileID = fopen([folder_path, fileName], 'w');
if fileID<0
    error([newline, 'Could not open the cloud file for writing: ', folder_path, fileName, newline])
end

fprintf(fileID, '%s %10s %9s %9s \n', '#', 'z', content_header, 'R_eff');
fprintf(fileID, '%s %10s %9s %9s \n', '#', '(km)', '(g/m^3)', '(micron)');

% Row 1 closes the cloud top: zeros here mean "no condensate above this level"
fprintf(fileID, '%12.4f %9.6f %9.3f \n', z_topBottom_km(1), 0, 0);

% Row 2 carries the properties of the single layer between cloud bottom and
% cloud top (values apply to the layer above the listed level)
fprintf(fileID, '%12.4f %9.6f %9.3f \n', z_topBottom_km(2), wc_gPerMeterCubed, r_e_um);

fclose(fileID);


end
