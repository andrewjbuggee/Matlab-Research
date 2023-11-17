% Create a lognoraml particle size distribution

% INPUTS:
%   (1) radius_modal - modal radius (microns) - this is a single value
%   that defines the radius that would be most commonly observed if we
%   randomly sampled our gamma particle size distribution a large number of
%   times

%   (2) std_dev - the standard deviation. It defines the width of the 
%   distribution. This has to be greater than 0.

%   (3) N0 - total droplet concentration (cm^(-3)) - this is the total
%   droplet number concentration including all sizes. If n(r) is integrated
%   over r, we get N0.

% OUTPUTS:
%   (1) n(r) - the number concentration of droplets for a given radius r -
%   this output is a vector of the same length as r, which is a hard-coded
%   vector

%   (2) r - the independent variable to defines our distribution n(r). 


% By Andrew John Buggee

%%

function [n_r,r] = lognormal_size_distribution_kokhanovsky(radius_modal, std_dev, N0)

% ------------------------------------------------------------
% ---------------------- CHECK INPUTS ------------------------
% ------------------------------------------------------------

% Check to make sure there are 3 inputs, droplet radius, cloud optical
% depth, and the altitude vector associated with this cloud


if nargin~=3
    error([newline,'Not enough inputs. Need 3: droplet effective radius, total number concentration, and the alpha parameter value.', newline])
end


%Check to make sure the standard deviation is greater than 0

if std_dev<=0

    error([newline,'The standard deviation must be greater than 0', newline])           % says who?
end




if N0<=0
    
    error([newline,'The number concentration must be greater than 0!', newline])
end

if radius_modal<=0
    
    error([newline,'The effective droplet radius must be greater than 0!',newline])
end

%% Define the gamma distribution
% This definition follows the definition of A. Kokhonovsky in 'Cloud
% Optics' on page 5, equation 1.21

% For some modal radius, we have define a droplet size distribution over a
% range of radii values

r = linspace(0.01*radius_modal, 8*radius_modal, 200);                  % microns - vector based on C.Emde (2016)


N = 1/(sqrt(2*pi) * std_dev);  % normalization constant

n_r = N0 * N * 1./r .* exp(-log(r./radius_modal).^2 ./(2*std_dev^2));       % lognormal droplet distribution



end