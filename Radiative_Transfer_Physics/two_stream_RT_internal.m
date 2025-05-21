%% Compute Two-Stream Reflectance, Transmittance and Absorption
% This is for a plane-parallel homogenous medium composed of a single
% substance


% INPUTS

%   (1): tau_vector - a vector spanning from 0 to some tau max. Used to
%   compute the upward and downward fluxes within the cloud

%   (2): ssa - single scattering albdeo of the particles that make up the
%   homogenous medium

%   (3): g - asymmetry parameter of the particles that make up the
%   homogenous medium

%   (4): albedo_maxTau - surface albedo below the layer. If there is
%   absorption and the medium isn't infintely thick, we must account for
%   the albedo below the layer. Are photons bouncing back up through the
%   medium, or are they absorbed? This is a number between 0 and 1


% OUTPUTS:



% By Andrew John Buggee

%%

function [photon_fraction_up, photon_fraction_down] = two_stream_RT_internal(tau_vector, ssa, g, albedo_maxTau)





%%


photon_fraction_up = zeros(length(tau_vector), 1);
photon_fraction_down = zeros(length(tau_vector), 1);

% grab the maximum tau value, the boundary value
tau_max = max(tau_vector);



% Check to see if there is absorption

if ssa<1


    % Next, check to see if our layer is infinitely thick, or has a finite
    % thickness

    if tau_max == inf

        % K is defined in Bohren and Clothiaux (eq. 5.70)
        K = sqrt((1 - ssa) .* (1 - g.*ssa));
        R_inf = (sqrt(1-ssa*g) - sqrt(1 - ssa))/(sqrt(1-ssa*g) + sqrt(1 - ssa));


        photon_fraction_up = R_inf * exp(-K*tau_vector);
        photon_fraction_down = exp(-K*tau_vector);






    elseif tau_max>0 && tau_max~=inf




        % K is defined in Bohren and Clothiaux (eq. 5.70)
        K = sqrt((1 - ssa)*(1 - g*ssa));
        % Define the reflectivity at the top of a semi-infinite cloud
        % R_inf = F_up(tau=0)/F_0
        % This is the fraction of light that scatters out the cloud top
        R_inf = (sqrt(1-ssa*g) - sqrt(1 - ssa))/(sqrt(1-ssa*g) + sqrt(1 - ssa));

        % Define the constants
        A = (R_inf - albedo_maxTau)*exp(-K*tau_max)/...
            (R_inf*(R_inf - albedo_maxTau)*exp(-K*tau_max) - (1 - albedo_maxTau*R_inf)*exp(K*tau_max));

        B = -(1 - R_inf*albedo_maxTau)*exp(K*tau_max)/...
            (R_inf*(R_inf - albedo_maxTau)*exp(-K*tau_max) - (1 - albedo_maxTau*R_inf)*exp(K*tau_max));


        photon_fraction_up = A*exp(K*tau_vector) + B*R_inf*exp(-K*tau_vector);


        photon_fraction_down = A*R_inf*exp(K*tau_vector) + B*exp(-K*tau_vector);






    end



    % Is our medium purely scattering incident light?

elseif ssa==1

    % Next, check to see if our layer is infinitely thick, or has a finite
    % thickness

    error([newline,'I dont know what to do for conservative scattering.',newline])


    if tau_max == inf




    elseif tau_max>0 && tau_max<inf

        % For a non-absorbing layer of finite thickness, the analytical
        % solutions to the two stream radiative transfer equations are...









    end




end







end
