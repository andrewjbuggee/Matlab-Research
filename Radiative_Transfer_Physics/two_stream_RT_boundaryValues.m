%% Compute Two-Stream Reflectance, Transmittance and Absorption
% This is for a plane-parallel homogenous medium composed of a single
% substance


% INPUTS

%   (1): tau_layer - optical depth of the layer being modeled

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

function [reflectance, transmittance, absorptance] = two_stream_RT_boundaryValues(tau_layer, ssa, g, albedo_maxTau)


% check to see if the optical depth, single scattering albedo and asymmetry
% parameter are the same length

if length(tau_layer)~=length(ssa) || length(tau_layer)~=length(g) || length(ssa)~=length(g)

    error([newline, 'The first three inputs must have the same number of values', newline])

end


%%


reflectance = zeros(1, length(tau_layer));
transmittance = zeros(1, length(tau_layer));
absorptance = zeros(1, length(tau_layer));


for nn = 1:length(tau_layer)


    % Check to see if there is absorption

    if ssa(nn)<1


        % Next, check to see if our layer is infinitely thick, or has a finite
        % thickness

        if tau_layer(nn) == inf

            % K is defined in Bohren and Clothiaux (eq. 5.70)
            K = sqrt((1 - ssa(nn)) .* (1 - g(nn).*ssa(nn)));
            reflectance(nn) = (sqrt(1-ssa(nn)*g(nn)) - sqrt(1 - ssa(nn)))/(sqrt(1-ssa(nn)*g(nn)) + sqrt(1 - ssa(nn)));

            transmittance(nn) = 0;      % nothing transmits in a medium that is infinitely thick

            absorptance(nn) = 1 - reflectance(nn);


            photon_fraction_up = @(tau) reflectance(nn) * exp(-K*tau);
            photon_fraction_down = @(tau) exp(-K*tau);






        elseif tau_layer(nn)>0 && tau_layer(nn)~=inf

            % Tau_layer is the maximum tau.


            % K is defined in Bohren and Clothiaux (eq. 5.70)
            K = sqrt((1 - ssa(nn))*(1 - g(nn)*ssa(nn)));
            % Define the reflectivity at the top of a semi-infinite cloud
            % R_inf = F_up(tau=0)/F_0
            % This is the fraction of light that scatters out the cloud top
            R_inf = (sqrt(1-ssa(nn)*g(nn)) - sqrt(1 - ssa(nn)))/(sqrt(1-ssa(nn)*g(nn)) + sqrt(1 - ssa(nn)));

            % Define the constants
            A = (R_inf - albedo_maxTau)*exp(-K*tau_layer(nn))/...
                (R_inf*(R_inf - albedo_maxTau)*exp(-K*tau_layer(nn)) - (1 - albedo_maxTau*R_inf)*exp(K*tau_layer(nn)));

            B = -(1 - R_inf*albedo_maxTau)*exp(K*tau_layer(nn))/...
                (R_inf*(R_inf - albedo_maxTau)*exp(-K*tau_layer(nn)) - (1 - albedo_maxTau*R_inf)*exp(K*tau_layer(nn)));


            photon_fraction_up = @(tau) A*exp(K*tau) + B*R_inf*exp(-K*tau);


            photon_fraction_down = @(tau) A*R_inf*exp(K*tau) + B*exp(-K*tau);


            reflectance(nn) = photon_fraction_up(0);      % fraction of incident energy that is scattered out the cloud top

            transmittance(nn) = photon_fraction_down(tau_layer(nn));            % fraction of incident energy that is transmitted through the later

            absorptance(nn) = 1 - reflectance(nn) - transmittance(nn);






        end



        % Is our medium purely scattering incident light?

    elseif ssa(nn)==1

        % Next, check to see if our layer is infinitely thick, or has a finite
        % thickness

        if tau_layer(nn) == inf


            error([newline,'I dont know what to do with a layer of infinte thickness and conservative scattering.',newline])


        elseif tau_layer(nn)>0 && tau_layer(nn)<inf

            % For a non-absorbing layer of finite thickness, the analytical
            % solutions to the two stream radiative transfer equations are...

            reflectance(nn) = (tau_layer(nn) * (1 - g(nn))/2)/(1 + tau_layer(nn)*(1 - g(nn))/2);

            transmittance(nn) = 1/(1 + tau_layer(nn)*(1 - g(nn))/2);

            absorptance(nn) = 0;        % no absorption if ssa=1




        end




    end



end



end
