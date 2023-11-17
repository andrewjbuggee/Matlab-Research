% Cloud Droplet Physics

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up the droplet distribution
mu = 5e-6; % mean droplet radius
sig = 0.1*mu; % standard deviation of the mean
numDrops = 10;

D = sig*randn(numDrops,1) + mu;

% creat a space for the droplets to exist in
% Let's define a box that is N*mu X M*mu in dimensions
N = 100;
M = 100;
X = N*mu*randn(numDrops,1);
Y = M*mu*randn(numDrops,1);
