%%  Transmission calculations of squares and particles

clear variables
% By Andrew John Buggee
%% --- Loop through random number of trails ---

numTrials = 100;

numParticles = 100;
edgeLength = 2;
plotFlag = false;

particleArea = zeros(1,numParticles);

for ii = 1:numTrials
    
    pgon = squaresOrParticles(edgeLength,numParticles,plotFlag);
    pgonUnion = union(pgon);
    
    particleArea(ii) = area(pgonUnion);
    
end

%% --- Calculations and Plots ---

boxArea = edgeLength^2;


% now calculate the amount of light that makes it through our box
Absorbed = particleArea/boxArea; % amount of light absorbed by our particles
Transmitted = 1 - Absorbed;

