% Create a log normal particle size distribution

% A random variable X is said to be log-normally distribution if Y = log(X)
% is a normal distribution. 


% INPUTS:
%   (1) r_modal - modal radius (microns) - this is a single value
%   that defines the radius that would be most commonly observed if we
%   randomly sampled our gamma particle size distribution a large number of
%   times

%   (2) sigma - This parameter defines the spread of our distribution. A
%   value of 0.34 decently replicates the spread of droplets observed in a
%   liquid water cloud

%   (2) N0 - total droplet concentration (cm^(-3)) - this is the total
%   droplet number concentration including all sizes. If n(r) is integrated
%   over r, we get N0.

% OUTPUTS:
%   (1) n(r) - the number concentration of droplets for a given radius r -
%   this output is a vector of the same length as r, which is a hard-coded
%   vector

%   (2) r - the independent variable to defines our distribution n(r). 


% By Andrew John Buggee

%%

function [n_r,r] = lognormal_size_distribution(radius_modal, std_dev, N0)

% ------------------------------------------------------------
% ---------------------- CHECK INPUTS ------------------------
% ------------------------------------------------------------

% Check to make sure there are 3 inputs, droplet radius, cloud optical
% depth, and the altitude vector associated with this cloud


if nargin~=3
    error([newline,'Not enough inputs. Need 3: droplet effective radius, total number concentration, and the alpha parameter value.', newline])
end

% Check to make sure re is the same length as the altitude vector

if std_dev<=0

    error([newline,'The alpha parameter cannot be less than 1', newline])
end

if N0<0
    
    error([newline,'The number concentration must be greater than 0!', newline])
end

if radius_modal<0
    
    error([newline,'The effective droplet radius must be greater than 0!',newline])
end

%%



% Compute the distribution

% -------------------------------------------------------------------------------------
% formula according to https://www.itl.nist.gov/div898/handbook/eda/section3/eda3669.htm
% in this formulation, r_modal is the mean of the log of the
% distribution and must be atleast 0. So if you want the mode to be x, then
% the input to this distribution should be r_modal = log(x)

% For some modal radius, we have define a droplet size distribution over a
% range of radii values

% if radius_modal==0
%     r = linspace(0.001,10,1000);
% elseif radius_modal>0
%     range_factor = 172.7*exp(-4.5*std_dev);         % factor that defines the total range
%     num_points = 23540*exp(-4.28*std_dev);
%     r = linspace(0.001*radius_modal, range_factor*radius_modal, round(num_points));              % microns - vector based on C.Emde (2016)
%     %r = logspace(floor(log10(r_modal/100)), floor(log10(r_modal*100)), 100);              % microns - vector based on C.Emde (2016)
% elseif radius_modal<0
%     error([newline, 'I dont think r_modal can be less than 0...',newline])
% end


% N = 1/(std_dev * sqrt(2*pi));                              % normalization constant
% n_r = N0* N./r .* exp(-(log(r) - radius_modal).^2 ./(2*std_dev.^2));                            % gamma droplet distribution
% -------------------------------------------------------------------------------------


% -------------------------------------------------------------------------------------
%N = 1/(sigma * sqrt(2*pi));                              % normalization constant
% formula according to Cloud Optics by Kokhanovsky

% For some modal radius, we have define a droplet size distribution over a
% range of radii values

% if radius_modal==0
%     r = linspace(0.001,10,1000);
% elseif radius_modal>0
%     range_factor = 172.7*exp(-4.5*std_dev);         % factor that defines the total range
%     num_points = 23540*exp(-4.28*std_dev);
%     r = linspace(0.001*radius_modal, range_factor*radius_modal, round(num_points));              % microns - vector based on C.Emde (2016)
%     %r = logspace(floor(log10(r_modal/100)), floor(log10(r_modal*100)), 100);              % microns - vector based on C.Emde (2016)
% elseif radius_modal<0
%     error([newline, 'I dont think r_modal can be less than 0...',newline])
% end


%n_r = N0* N./r .* exp(-log(r./r_modal).^2 ./(2*sigma.^2));                            % gamma droplet distribution
% -------------------------------------------------------------------------------------



% -------------------------------------------------------------------------------------
% formula according to LibRadTran user's guide page 130
% in this case, r_modal cannot be less than 1
% sigma is>0 and <1 ?
% ***** I still don't understand this version well *****
% The modal radius I select doesn't line up with the peak of the
% distribution. Should it?

% For some modal radius, we have define a droplet size distribution over a
% range of radii values

if radius_modal==0
    r = linspace(0.001,10,1000);
elseif radius_modal>0
    range_factor = 172.7*exp(-4.5*std_dev);         % factor that defines the total range
    num_points = 23540*exp(-4.28*std_dev);
    r = linspace(0.001*radius_modal, range_factor*radius_modal, round(num_points));              % microns - vector based on C.Emde (2016)
    %r = logspace(floor(log10(r_modal/100)), floor(log10(r_modal*100)), 100);              % microns - vector based on C.Emde (2016)
elseif radius_modal<0
    error([newline, 'I dont think r_modal can be less than 0...',newline])
end



N = N0 * (sqrt(1/log(std_dev)^2)/sqrt(2*pi));                              % normalization constant
n_r = N./r .* exp(-(1/2) * ((log(r) - log(radius_modal))./log(std_dev)).^2);                            % gamma droplet distribution
% -------------------------------------------------------------------------------------


end
