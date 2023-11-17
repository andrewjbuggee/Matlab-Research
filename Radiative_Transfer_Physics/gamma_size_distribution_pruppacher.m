% Create a gamma particle size distribution

% INPUTS:
%   (1) r_modal - modal radius (microns) - this is a single value
%   that defines the radius that would be most commonly observed if we
%   randomly sampled our gamma particle size distribution a large number of
%   times

%   (2) mu - paramter that characterizes the width of the distribution. 
%   Mu must be a number greater than or equal to 1. The larger mu is,
%   the more narrow the distribution becomes. (Does it neeed to be a
%   integer?)

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

function [n_r,r] = gamma_size_distribution_pruppacher(r_modal, mu, N0)

% ------------------------------------------------------------
% ---------------------- CHECK INPUTS ------------------------
% ------------------------------------------------------------

% Check to make sure there are 3 inputs, droplet radius, cloud optical
% depth, and the altitude vector associated with this cloud


if nargin~=3
    error([newline,'Not enough inputs. Need 3: droplet effective radius, total number concentration, and the alpha parameter value.', newline])
end


%Check to make sure mu is greater than 1

% if mu<1
% 
%     error([newline,'The mu parameter cannot be less than 1', newline])           % says who?
% end
% 
% %Check to make sure that mu is an integer
% if (mu - floor(mu))>0
% 
%     error([newline, 'The mu parameter must be an integer', newline])
% 
% end



if N0<0
    
    error([newline,'The number concentration must be greater than 0!', newline])
end

if r_modal<0
    
    error([newline,'The effective droplet radius must be greater than 0!',newline])
end

%% Define the gamma distribution
% This definition follows the definition of Pruppacher and Klett in
% 'Microphysics of Clouds and Precipitation' on page 483

% For some modal radius, we have define a droplet size distribution over a
% range of radii values

r = linspace(0.01*r_modal, 8*r_modal, 200);                  % microns - vector based on C.Emde (2016)


N = (mu+1)^(mu+1)/(gamma(mu+1));  % normalization constant

n_r = N0 * N * (r.^3./r_modal.^3).^mu .* exp(-(mu+1)*(r.^3./r_modal.^3));                            % gamma droplet distribution



end