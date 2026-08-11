%% Read a uvspec output file produced by a thermal twostr run with output_process per_band
%
% Reads the four-column output written by
%
%       output_user lambda edir edn eup
%       output_process per_band
%
% and returns the band-resolved irradiances plus their broadband sums.
%
% NOTE ON THE EXISTING WRAPPER: readUVSPEC_ver2.m in the LibRadTran-wrapper
% folder cannot be used here. Its 'twostr' branch (line ~612) operates on the
% variable 'data', but 'data' is only ever populated inside the 'disort' branch,
% so it is still the empty array [] when the twostr code runs and the indexing
% data(:,2) errors out. This reader is deliberately minimal and specific to the
% column layout requested above.
%
%
% INPUTS:
%   (1) folder_path - folder holding the output file
%   (2) fileName - name of the output file, INCLUDING the .OUT extension
%
%
% OUTPUTS:
%   (1) out - structure with fields:
%           .wavelength_nm      - band centre wavelength (nm), one row per band
%           .edir_WPerMeterSq   - direct irradiance per band (W/m^2). Identically
%                                 zero for a thermal-only run.
%           .edn_WPerMeterSq    - downwelling diffuse irradiance per band (W/m^2)
%           .eup_WPerMeterSq    - upwelling diffuse irradiance per band (W/m^2)
%           .LW_down_WPerMeterSq - broadband downwelling longwave irradiance (W/m^2)
%           .LW_up_WPerMeterSq   - broadband upwelling longwave irradiance (W/m^2)
%
%
% By Andrew John Buggee
%%

function [out] = read_uvspec_thermal_perBand(folder_path, fileName)

if nargin~=2
    error([newline, 'Need 2 inputs: the folder path and the output file name.', newline])
end

if folder_path(end)~=filesep
    folder_path = [folder_path, filesep];
end

full_path = [folder_path, fileName];

if exist(full_path, 'file')==0
    error([newline, 'Could not find the uvspec output file: ', full_path, newline])
end

file_char = fileread(full_path);

if isempty(strtrim(file_char))==true
    error([newline, 'The uvspec output file is empty. Check errMsg.txt in: ', folder_path, newline])
end

data = cell2mat(textscan(file_char, '%f %f %f %f'));

if size(data,2)~=4
    error([newline, 'Expected 4 columns (lambda, edir, edn, eup) but found ',...
        num2str(size(data,2)), '.', newline])
end

out.wavelength_nm       = data(:,1);            % nm - reptran band centre
out.edir_WPerMeterSq    = data(:,2);            % W/m^2 per band (zero for source thermal)
out.edn_WPerMeterSq     = data(:,3);            % W/m^2 per band
out.eup_WPerMeterSq     = data(:,4);            % W/m^2 per band

% band-integrated quantities are already per band, so the broadband flux is a
% plain sum over bands
out.LW_down_WPerMeterSq = sum(data(:,2) + data(:,3));      % W/m^2
out.LW_up_WPerMeterSq   = sum(data(:,4));                  % W/m^2


end
