%% Read and Interpolate Refractive Index of Water

% This function will read and interpolate experimentally determined
% refractive indicies of water.

% INPUTS:

%   (1) wl - wavelength (microns). This is a vector that defines the
%   wavelength query points. The file contains index of refraction values
%   between 0.2 and 200 microns


% OUTPUTS:
%   (1) yq - the values of the refractive index for water, broken up into
%   the real and imaginary parts.
%       (a) refrac_real - real part of the refractive index
%       (b) refrac_imag - imaginary part of the refractive index


% All look up tables were computed using the opensoure radiative transfer
% code, LibRadTran

% By Andrew John Buggee

% ------------------------------------------------------------------------
%%

function [yq] = interp_refractive_index_H2O(wl)

% ----------------------------------------
% ------------- Check Inputs! ------------
% ----------------------------------------

% Check to make sure there are two inputs


if nargin~=1
    error([newline,'Not enough inputs. Need 1: Enter a vector of wavelengths in units',...
        ' of microns.', newline])
end

% Check to make sure the input is within the bounds of the Hale and Query
% (1973) data

if any(wl<0.2) || any(wl>200)
    error([newline, 'Wavelength query points are outside the bounds of the data set.',...
        ' Make sure all entries are greater than 0.2 and less than 200 microns.', newline])
end



% Determine which computer you're using
computer_name = whatComputer;

% use the proper path
if strcmp(computer_name,'anbu8374')==true
    error([newline, 'No Path defined yet!', newline])
    folder_path = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/Mie_Calculations/';
    
elseif strcmp(computer_name,'andrewbuggee')==true
    folder_path = '/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Radiative_Transfer_Physics/';
    
end

%%
% ------------------------------------------
% ------------- Function Stuff! ------------
% ------------------------------------------


% Load the data file
filename = 'refractiveIndex_H2O.txt';
format_spec = '%f %f %f';        % 1 column of data



% ----- READ IN DATA USING TEXTSCAN() ---------------

file_id = fopen([folder_path,filename]);

data = textscan(file_id, format_spec, 'CommentStyle','#', 'Delimiter',...
    {'\r', ' '}, 'MultipleDelimsAsOne',1);

% then we will interpoalte
% Lets grab all of the values we need in the data set

index_real = interp1(data{1}, data{2}, wl);
index_imag = interp1(data{1}, data{3}, wl);

% Lets include the wavelength and effective radius in the
% interpolated data cube

yq = [reshape(index_real,[],1), reshape(index_imag,[],1)];



end
