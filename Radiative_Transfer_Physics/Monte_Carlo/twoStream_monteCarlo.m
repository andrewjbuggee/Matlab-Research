% 2-stream radiative transfer using Monte Carlo methods
% by Andrew John Buggee

function [F_norm, final_state, photon_tracking, inputs] = twoStream_monteCarlo(inputs)

% ---------------------------------------
% ------- Unpack input structure --------
% ---------------------------------------

% Define the boundaries of the medium
tau_upper_limit = inputs.tau_upper_limit;
tau_lower_limit = inputs.tau_lower_limit;

% Define the number of layers within the medium that differ
N_layers = inputs.N_layers;

% Layers can differ by radius, by material, or both. If the layers differ
% by radius, create a vector describing each layer radius from top to
% bottom (tau = 0 to tau = tau0). If the layers differ by index of
% refraction, create a vector describing each layer indec of refraction
% from top to bottom. If both are true, create a cell array for both the
% changing radii and the changing index of refraction

layerRadii = inputs.layerRadii;      % radius change between two layers

% Define the layer boundaries given the number of layers and the boundaries
% of the entire medium
layerBoundaries = inputs.layerBoundaries;

% Define the albedo of the lower tau boundary
albedo_maxTau = inputs.albedo_maxTau;

% Define the number of photons to inject into the medium
N_photons = inputs.N_photons;

% define the wavelength
if inputs.mie.wavelength(3)>0
    wavelength = inputs.mie.wavelength(1):inputs.mie.wavelength(3):inputs.mie.wavelength(2);           % nanometers
else
    wavelength = inputs.mie.wavelength(1);           % nanometers
end


% Define the size of the scatterer and its scattering properties
% Assuming a pure homogenous medium composed of a single substance
% Define the single scattering albedo of the substance in question
if inputs.mie.radius(3)>0
    r = inputs.mie.radius(1):inputs.mie.radius(3):inputs.mie.radius(2);                     % microns
else
    r = inputs.mie.radius(1);
end

% ------------------------------------------------------------------------
% Organize the scattering properties in the layer order from top to bottom
% ------------------------------------------------------------------------
g = zeros(length(wavelength), length(r));
ssa = zeros(length(wavelength), length(r));
for LL = 1:N_layers
    % define the asymmetry parameter
    [~,index] = min(abs(r-layerRadii(LL)));
    g(:,LL) = inputs.g(:,index);

    % define the single scattering albedo
    ssa(:,LL) = inputs.ssa(:,index);

end

% Redefine the new order of g and ssa
inputs.g = g;
inputs.ssa = ssa;

% -----------------------------------------------------------
% *** Define the probability of travelling some depth tau ***
% -----------------------------------------------------------

% Let's define T as the random variable that describes the probability
% of a photon travelling an optical depth T. P(T) is that probability

% Let's define x to be the uniform PDF where x is between [0,1].

% define the relationship between x and tau

tau = @(x) -log(1 - x);


% reset the random number generator
rng('default');

% -------------------------------------------------------------
% *** Create a bunch of empty arrays to keep track of stuff ***
% -------------------------------------------------------------

% we have to keep track of how our photon travels within the medium
depth_travelled = cell(1,N_photons);

% Save all of the direction changes so we can easily compute upwards and
% downwards irradiance
direction = cell(1, N_photons);

% Store the values of the maximum penetration depth for each photon
max_depth = zeros(1, N_photons);

% Store the value of the depth at which a photon is absorbed
depth_absorbed = [];

% Store the values of the number of scattering events for each photon
number_of_scattering_events = zeros(1, N_photons);

% Let's track the final state of each photon
% do they end up leaving the medium? Do they scatter out the top of the
% bottom? Or are they absorbed along their journey?
scatter_out_top = 0;
scatter_out_bottom = 0;
absorbed = 0;


% I also want to know which photons end up where. So we will keep track of
% the index of each final state. We can then trace back the path of each
% photon
scatter_out_top_index = [];
scatter_out_bottom_index = [];
absorbed_index = [];


% ----------------------------------------------------------------------
% *** Inject photons into the medium and keep track of what happens! ***
% ----------------------------------------------------------------------




for nn = 1:N_photons



    % Each photon inject has a probability of travelling a certain distance
    % before some event occurs. This event can be either scattering or
    % absorption.

    % when photons are injected into our medium, they start with a downward
    % direction

    % travelling downward is defined as positive increase in tau
    % travelling upwards is defined as a negative increase in tau

    % Start with two values, the first 1 tells us the photon was moving
    % down at t = 0, the second says it was movng down at until the tau
    % that is sampled
    direction{nn} = [1,1];


    % --------------------------------------------------------
    % *****----- Roll dice and move photon forward -----******
    % --------------------------------------------------------

    % We start by drawing a value from our uniform distribution and
    % inserting it into the distribution that describes the probability of
    % travelling a distance tau (dT/dx relationship)
    tau_sample = tau(rand(1,1));

    % This tells us how far our photon travels before something happens.

    % we need to keep track of where the photon is within our medium
    % Start by giving each vector a 0 to show the starting point at time
    % t=0
    depth_travelled{nn} = [0, tau_sample];


    % -----------------------------------------
    % **** Has the photon left out medium? ****
    % -----------------------------------------

    % There is a chance we draw a number that causes the photon to transmit
    % through without any scattering events!

    if depth_travelled{nn}(end)>tau_upper_limit

        % The photon has transmitted right through our medium with any
        % scattering or absorption events.

        % ------------------------------------------------------
        % **** Did our photon reflect back into our medium? ****
        % ------------------------------------------------------

        % When our photon reaches the lower boudnary, we will check to see
        % if the photon is absorbed or scattered back into the medium by
        % comparing a random uniform number drawing with the albedo

        scat_or_abs = rand(1,1);

        if scat_or_abs<=albedo_maxTau
            % The photon is scattered back into our medium!!

            % First, set the distance travelled to be the boundary
            % condition value, since we are treating this as a hard
            % boundary
            depth_travelled{nn}(end) = tau_upper_limit;

            % Photon switches direction and is back scattered
            direction{nn}(end+1) = -1*direction{nn}(end);       % switches direction (back scattered)

            % --------------------------------------------------------
            % *****----- Roll dice and move photon forward -----******
            % --------------------------------------------------------

            % Determine how far the photon travels after being scattered
            % back into our medium
            tau_sample = tau(rand(1,1));

            % we need to keep track of where the photon is within our medium
            % Start by giving each vector a 0 to show the starting point at time
            % t=0
            depth_travelled{nn}(end+1) = direction{nn}(end)*tau_sample + depth_travelled{nn}(end);

            % ------------------------------------------------------
            % ----- Check to see what layer is our photon in!! -----
            % ------------------------------------------------------
            if N_layers>1
                I_lessThan = find(layerBoundaries<depth_travelled{nn}(end));
                I_greaterThan = find(layerBoundaries>depth_travelled{nn}(end));

                % set the index describing which layer the photon is in
                index_layer = I_lessThan(end);

            else

                % the index layer is 1, since there is only 1
                index_layer = 1;

            end
            % ------------------------------------------------------


        else

            % Our photon was absorbed


            % let's set the final depth_travelled value to be the
            % tau limit
            depth_travelled{nn}(end) = tau_upper_limit;
            % Let's record which boundary the photon surpassed
            scatter_out_bottom = scatter_out_bottom + 1;

            % record the index
            scatter_out_bottom_index = [scatter_out_bottom_index, nn];



        end



    else

        % If the photon does not travel through the medium on the first tau
        % draw, or if it is scattered back into our medium due to the bottom
        % boundary, then we redraw a tau and continue its life.

        % ------------------------------------------------------
        % ----- Check to see what layer is our photon in!! -----
        % ------------------------------------------------------
        if N_layers>1
            I_lessThan = find(layerBoundaries<depth_travelled{nn}(end));

            % set the index describing which layer the photon is in
            index_layer = I_lessThan(end);

        else

            % the index layer is 1, since there is only 1
            index_layer = 1;

        end
        % ------------------------------------------------------


        % --------------------------------------------------------
        % ***---- Roll dice. Did photon scatter or absorb? ---****
        % --------------------------------------------------------

        % To determine if the photon scattered or absorbed, draw a random
        % number from a uniform distribution. If this number is less than or
        % equal to the single scattering albedo, the photon scattered.
        % Otherwise we assume it absorbed

        scat_or_abs = rand(1,1);

        while scat_or_abs<=ssa(index_layer)

            % If this is true, our photon scattered! We need a way of keeping
            % track of the direction our photon is moving because we have
            % boundaries.

            % ---------------------------------------------------------
            % *** Determine which direction the photon scattered in ***
            % ---------------------------------------------------------

            if direction{nn}(end)==1
                % photon is moving down
                % draw a random number
                g_sample = rand(1,1);
                % is this less than or equal to the probability of forward
                % scattering?
                if g_sample <= ((1+g(index_layer))/2)
                    % the the photon scattered in the forward direction
                    direction{nn}(end+1) = direction{nn}(end);       % Continues along the same direction

                else
                    % Photon switches direction and is back scattered
                    direction{nn}(end+1) = -1*direction{nn}(end);       % switches direction (back scattered)

                end


            else

                % photon is moving up
                % draw a random number
                g_sample = rand(1,1);
                % is this less than or equal to the probability of forward
                % scattering?
                if g_sample <= ((1+g(index_layer))/2)
                    % the the photon scattered in the forward direction
                    direction{nn}(end+1) = direction{nn}(end);       % Continues along the same direction

                else
                    % Photon switches direction and is back scattered
                    direction{nn}(end+1) = -1*direction{nn}(end);       % switches direction (back scattered)

                end


            end



            % -----------------------------------------------------------
            % ** Determine how far the photon travels after scattering **
            % -----------------------------------------------------------
            % draw another random variable and determine how far the photon
            % travels
            tau_sample = tau(rand(1,1));

            % record where the photon is
            depth_travelled{nn}(end+1) = direction{nn}(end)*tau_sample + depth_travelled{nn}(end);




            % -------------------------------------------------
            % **** Has the photon breached any boundaries? ****
            % -------------------------------------------------

            if depth_travelled{nn}(end)<tau_lower_limit

                % The photon has scattered out the top of our medium
                % if so, let's set the final depth_travelled value to be 0
                depth_travelled{nn}(end) = 0;

                % Let's record which boundary the photon surpassed
                scatter_out_top = scatter_out_top + 1;

                % record the index
                scatter_out_top_index = [scatter_out_top_index, nn];


                break

            elseif depth_travelled{nn}(end)>tau_upper_limit

                % The photon has scattered out the bottom of our medium


                % ------------------------------------------------------
                % **** Did our photon reflect back into our medium? ****
                % ------------------------------------------------------

                % When our photon reaches the lower boudnary, we will check to see
                % if the photon is absorbed or scattered back into the medium by
                % comparing a random uniform number drawing with the albedo

                scat_or_abs = rand(1,1);

                if scat_or_abs<=albedo_maxTau
                    % The photon is scattered back into our medium!!

                    % First, set the distance travelled to be the boundary
                    % condition value, since we are treating this as a hard
                    % boundary
                    depth_travelled{nn}(end) = tau_upper_limit;

                    % Photon switches direction and is back scattered
                    direction{nn}(end+1) = -1*direction{nn}(end);       % switches direction (back scattered)

                    % --------------------------------------------------------
                    % *****----- Roll dice and move photon forward -----******
                    % --------------------------------------------------------

                    % Determine how far the photon travels after being scattered
                    % back into our medium
                    tau_sample = tau(rand(1,1));

                    % we need to keep track of where the photon is within our medium
                    % Start by giving each vector a 0 to show the starting point at time
                    % t=0
                    depth_travelled{nn}(end+1) = direction{nn}(end)*tau_sample + depth_travelled{nn}(end);


                    % ------------------------------------------------------
                    % ----- Check to see what layer is our photon in!! -----
                    % ------------------------------------------------------
                    if N_layers>1
                        I_lessThan = find(layerBoundaries<depth_travelled{nn}(end));

                        % set the index describing which layer the photon is in
                        index_layer = I_lessThan(end);

                    else

                        % the index layer is 1, since there is only 1
                        index_layer = 1;

                    end
                    % ------------------------------------------------------

                else

                    % Our photon was absorbed by the lower boundary


                    % let's set the final depth_travelled value to be the
                    % tau limit
                    depth_travelled{nn}(end) = tau_upper_limit;
                    % Let's record which boundary the photon surpassed
                    scatter_out_bottom = scatter_out_bottom + 1;

                    % record the index
                    scatter_out_bottom_index = [scatter_out_bottom_index, nn];


                    break



                end


            end


            % ------------------------------------------------------
            % ----- Check to see what layer is our photon in!! -----
            % ------------------------------------------------------
            if N_layers>1
                I_lessThan = find(layerBoundaries<depth_travelled{nn}(end));

                % set the index describing which layer the photon is in
                index_layer = I_lessThan(end);

            else

                % the index layer is 1, since there is only 1
                index_layer = 1;

            end
            % ------------------------------------------------------


            % draw a new random number! Did our photon scatter or absorb?
            scat_or_abs = rand(1,1);


        end


    end



    % --------------------------------------------------
    % **** Did the photon perish due to absorption? ****
    % --------------------------------------------------
    % If the photon ended in absorption, record it!

    if scat_or_abs>ssa(index_layer) && depth_travelled{nn}(end)>tau_lower_limit && depth_travelled{nn}(end)<tau_upper_limit

        % Make sure all three conditions are met!
        % The random number, scat_or_abs is greater than the ssa
        % the photon is still within the boundaries of our defined medium

        absorbed = absorbed + 1;

        % record the index
        absorbed_index = [absorbed_index, nn];

        % record the depth of absorption
        depth_absorbed = [depth_absorbed, depth_travelled{nn}(end)];

    end



    % Save the maximum depth reached.
    max_depth(nn) = max(depth_travelled{nn});

    % Save the number of scattering events
    % (1) If the final event was absorption, we substract two from the
    % depth_travelled vector
    % (2) If the photon was not absorbed, it left our medium. We ignore the
    % first entry in depth_travlled, since it is a boundary condition.
    % We also must ignore the last entry, since we are no longer
    % tracking scattering events outside the medium
    number_of_scattering_events(nn) = length(depth_travelled{nn})-2;



end


%%

% ----------------------------------------------------------------------------------
% *** Divide the tau space into a grid and compute reflectance and transmittance ***
% ----------------------------------------------------------------------------------

% Set up the tau grid first

N_bins = 100;

if tau_upper_limit==inf
    binEdges = logspace(-3, ceil(log10(max(max_depth))),N_bins+1);

else
    binEdges = linspace(tau_lower_limit, tau_upper_limit,N_bins+1);
end


% Create the bins that will keep tally
N_counts_moving_up = zeros(1, N_bins);
N_counts_moving_down = zeros(1, N_bins);

% We want to tally the direction a photon moves through each bin defined
% above

% Check to see if we switched directions. If we did, we want to
% make sure we tally that a photon was moving in both
% directions in the bin where it turned around


% assign each x to one of our bins
for nn=1:N_photons



    for tt = 1:length(depth_travelled{nn})-1

        % Check to see if the photon is moving down or up and check to see
        % if its continuing along the same direction, or it it's switched
        % directions

        % ----------------------------------------
        % **** Photon continuing to move down ****
        % ----------------------------------------
        if direction{nn}(tt+1)==1 && direction{nn}(tt)==1
            % *** the photon is moving down ***
            % Check the see which tau bins the photon moves through. Tally each
            % bin that is found


            % If the photon is moving down and was already heading
            % down, we don't need to account for an additional bin

            % If the photon is moving down, then when it cross bin-edge 3,
            % its in bucket 3.
            bins_photon_moves_through = binEdges>=depth_travelled{nn}(tt) & ...
                binEdges<depth_travelled{nn}(tt+1);

            N_counts_moving_down(bins_photon_moves_through) = N_counts_moving_down(bins_photon_moves_through) +1;


            % --------------------------------------------
            % **** Photon moving down after moving up ****
            % --------------------------------------------
        elseif direction{nn}(tt+1)==1 && direction{nn}(tt)==-1
            % If the photon is moving down, then when it cross bin-edge 3,
            % its in bucket 3.
            bins_photon_moves_through = binEdges>=depth_travelled{nn}(tt) & ...
                binEdges<depth_travelled{nn}(tt+1);

            % Just make sure the first bin is not 1! We can't have a 0
            % index in matlab
            if find(bins_photon_moves_through,1)~=1
                bins_photon_moves_through(find(bins_photon_moves_through,1)-1) = 1;
            end

            N_counts_moving_down(bins_photon_moves_through) = N_counts_moving_down(bins_photon_moves_through) +1;



            % --------------------------------------
            % **** Photon continuing to move up ****
            % --------------------------------------
        elseif direction{nn}(tt+1)==-1 && direction{nn}(tt)==-1
            % *** the photon is moving up ***
            % When this is true, depth_travlled{nn}(tt+1) is always less
            % than depth_travelled{nn}(tt) so we need to determine the tau
            % bins that the photon passes through between
            % depth_travlled{nn}(tt) and depth_travlled{nn}(tt+1)


            % If the photon is continuing up then we only have to count
            % the bins it passes through
            bins_photon_moves_through = binEdges<=depth_travelled{nn}(tt) & ...
                binEdges>=depth_travelled{nn}(tt+1);

            bins_photon_moves_through = find(bins_photon_moves_through)-1;


            % Just make sure the first bin is not 1! We can't have a 0
            % index in matlab
            bins_photon_moves_through(bins_photon_moves_through<1)=[];


            N_counts_moving_up(bins_photon_moves_through) = N_counts_moving_up(bins_photon_moves_through) +1;




            % --------------------------------------------
            % **** Photon moving up after moving down ****
            % --------------------------------------------
        elseif direction{nn}(tt+1)==-1 && direction{nn}(tt)==1


            % If the photon switched direction, we have to account for
            % the final bin it ends up in

            bins_photon_moves_through = binEdges<=depth_travelled{nn}(tt) & ...
                binEdges>=depth_travelled{nn}(tt+1);

            % Just make sure the first bin is not 1! We can't have a 0
            % index in matlab
            if find(bins_photon_moves_through,1)~=1
                bins_photon_moves_through(find(bins_photon_moves_through,1)-1) = 1;
            end

            % Some photons in the category will reflect off the bottom
            % boudnary. If this happened, we simply need to remove the
            % logical true value for the binEdge equal to tau_upper_limit.
            % If we don't we get an error, and all we need to keep track of
            % is whether or not the photon passed through this bin
            if bins_photon_moves_through(end)==true && (depth_travelled{nn}(tt)==tau_upper_limit || depth_travelled{nn}(tt+1)==tau_upper_limit)
                bins_photon_moves_through(end) = 0;
            end


            N_counts_moving_up(bins_photon_moves_through) = N_counts_moving_up(bins_photon_moves_through) +1;



        else


            error([newline,'Im not sure what the photon is doing',newline])

        end


    end

end





% ----------------------------------------------
% *** Let's collect all the output variables ***
% ----------------------------------------------

% compute the normalized upward moving irradiance
F_norm.up = N_counts_moving_up./N_photons;


% compute the normalized downward moving irradiance
F_norm.down = N_counts_moving_down./N_photons;

% save the bin edges
F_norm.binEdges = binEdges;


% output all final states of each photon
final_state.scatter_out_top = scatter_out_top;
final_state.scatter_out_bottom = scatter_out_bottom;
final_state.absorbed = absorbed;

% output the indices for each final state
final_state.scatter_out_top_INDEX = scatter_out_top_index;
final_state.scatter_out_bottom_INDEX = scatter_out_bottom_index;
final_state.absorbed_INDEX = absorbed_index;


% output other photon tracking variables that I've computed

% keep track of the maximum penetration depth of each photon
photon_tracking.maxDepth = max_depth;

% Count how many times a scatter/absorption event occurs during a photons
% lifetime. A value of 1 means it was absorbed right away or transmitted
% through
photon_tracking.number_of_scattering_events = number_of_scattering_events;

% Record the depths at which absorption occured
photon_tracking.absorbed_depth = depth_absorbed;







end