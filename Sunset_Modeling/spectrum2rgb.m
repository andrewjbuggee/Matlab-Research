%% This function converts a sprectrum to rgb true colors

% This function will read and interpolate precompute Mie calculations for
% water droplets of varrying radii.

% INPUTS:
%   (1) L - radiance spectrum (W/m^2/nm/sr) - spectrum to be converted

%   (2) lambda - wavelength vector (nm) - this is the wavelength vector the
%   defines the spectral location of the radiance input.

%   (3) color_depth - color depth of the output - this can eitehr be 8-bit
%   or 16-bit. The input should be a string
%       (a) '8bit'
%       (b) '16bit'


% OUTPUTS:
%   (1) RGB - true RGB color spanning [0,255] or [0, 65535]


% By Andrew John Buggee

%%

function [RGB] = spectrum2rgb(L, lambda, color_depth)

% ------------------------------------------------------------
% ---------------------- CHECK INPUTS ------------------------
% ------------------------------------------------------------

% Check to make sure there are 3 inputs, droplet radius, cloud optical
% depth, and the altitude vector associated with this cloud


if nargin<3
    error([newline,'Not enough inputs. Need 2: radiance specrum and color depth.', newline])
elseif nargin>3
    error([newline,'Too many inputs. Need 2: radiance specrum and color depth.', newline])

end

% Check to make sure the radiance and wavelength vectors are the same
% length

if length(L)~=length(lambda)
    
    error([newline,'The radiance vector must be the same length as the wavelength vector.', newline])
    
end



% Check to make sure the color depth string is one of two options

if strcmp(color_depth, '8bit')==false && strcmp(color_depth, '16bit')==false
    
    error([newline,'I dont recognize the color_depth string. Must be either "8bit" or "16bit"', newline])
end


% Clip the wavelength vector outside the range of the color matching
% functions and issue a warning
lambda_min = 360;                   % nm
lambda_max = 830;                   % nm


if any(lambda<lambda_min) || any(lambda>lambda_max)
    
    index = lambda<lambda_min | lambda>lambda_max;
    lambda(index) = [];
    
    % We also need to clip the radiance vector
    L(index) = [];
    
    warning([newline, 'Some wavelength values are outside the range of human vision (about 360 to 830 nm).',...
        'Weve computed the RGB color values for the input spectrum within this range.', newline]);
end



%% Read Color matching functions and interpolate

% read in color matching functions
% Grab the color matching functions
[wl, x_func, y_func, z_func] = colorMatchFcn('1931_full');

% Lets interpolate the XYZ color functions so the have the same resolution
% as our transmission spectrum
x_func2 = interp1(wl,x_func,lambda);
y_func2 = interp1(wl,y_func,lambda);
z_func2 = interp1(wl,z_func,lambda);


% Lets determine the CIE-XYZ color values
% we do this by integrating the product of the radiance with the XYZ color
% functions over the wavelength 

X = trapz(lambda, x_func2.*L);
Y = trapz(lambda, y_func2.*L);
Z = trapz(lambda, z_func2.*L);


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

% now we convert CIE XYZ values to sRBG values

if strcmp(color_depth, '8bit')==true
    
    RGB = xyz2rgb(xyz, 'OutputType','uint8', 'ColorSpace', 'srgb', 'WhitePoint','d65');
    
elseif strcmp(color_depth, '16bit')==true
    
    RGB = xyz2rgb(xyz, 'OutputType','uint16', 'ColorSpace', 'srgb', 'WhitePoint','d65');
    
end







end

