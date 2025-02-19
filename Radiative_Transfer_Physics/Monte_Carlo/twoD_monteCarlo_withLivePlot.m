% 2-stream radiative transfer using Monte Carlo methods
% by Andrew John Buggee

function [F_norm, final_state, photon_tracking, inputs] = twoD_monteCarlo_withLivePlot(inputs)

% ---------------------------------------
% ------- Unpack input structure --------
% ---------------------------------------

% Define the boundaries of the medium
tau_y_upper_limit = inputs.tau_y_upper_limit;
tau_y_lower_limit = inputs.tau_y_lower_limit;

% Define the solar zenith angle
solar_zenith_angle = inputs.solar_zenith_angle;

% Define the number of layers within the medium that differ
N_layers = inputs.N_layers;

% Layers can differ by radius, by material, or both. If the layers differ
% by radius, create a vector describing each layer radius from top to
% bottom (tau = 0 to tau = tau0). If the layers differ by index of
% refraction, create a vector describing each layer indec of refraction
% from top to bottom. If both are true, create a cell array for both the
% changing radii and the changing index of refraction

% Define the layer boundaries given the number of layers and the boundaries
% of the entire medium
layerBoundaries = inputs.layerBoundaries;

% Define the albedo of the lower tau boundary
albedo_maxTau = inputs.albedo_maxTau;

% Define the number of photons to inject into the medium
N_photons = inputs.N_photons;

% unpack the mie scattering properties!

% asymmetry parameter for each layer
if inputs.mie.integrate_over_size_distribution==false
    
    % then we use the values directly computed by the mie program
    
    % define the asymmetry parameter
    g = inputs.g;

    % define the single scattering albedo
    ssa = inputs.ssa;

else

    % then we use the values that have been integrated over some size
    % distribution

    % define the average asymmetry parameter
    g = inputs.g_avg;

    % define the average single scattering albedo
    ssa = inputs.ssa_avg;

end




% -----------------------------------------------------------
% *** Define the probability of travelling some depth tau ***
% -----------------------------------------------------------

% Let's define T as the random variable that describes the probability
% of a photon travelling an optical depth T. P(T) is that probability

% Let's define x to be the uniform PDF where x is between [0,1].

% define the relationship between x and tau

tau = @(x) -log(1 - x);



% ------------------------------------------------------------------------
% *** Define the probability of scattering into a particular direction ***
% ------------------------------------------------------------------------

% The Henyey-Greenstein function describes the probability of a photon
% scattering into a direction mu for some asymmetry parameter g (mu =
% cos(theta). Therefore, to use this in a monte carlo model, we have to
% convert P(mu) into P(x(mu)) using the uniformly distributed random
% variable x. P(mu)*dmu/dx * dmu = P(x) dx. Solving this, we get a function
% mu(x):

mu = @(g, x) 1/(2*g) * (1 + g^2 - ((1 - g^2)./(1 - g + 2*g*x)).^2);

% Values of mu between (0,1] are scattered in the forward direction.
% Values of mu between [-1,0) are scattered in the backwards direction


% ---------------------------------------------------------------------------------------
% *** Define the clockwise and counter-clockwise rotation transformation matrix in 2D ***
% ---------------------------------------------------------------------------------------

% There are from the perspective of the most recent tau vector

% clockwise rotation
M_clockwise = @(Mu) [Mu, sqrt(1 - Mu^2);...
                    -sqrt(1 - Mu^2), Mu];

% counter-clockwise rotation
M_counterClockwise = @(Mu) [Mu, -sqrt(1 - Mu^2);...
                            sqrt(1 - Mu^2), Mu];




% ---------------------------------------------------------------------------------------------
% *** Define the function that computes the x and y components in the new coordinate system ***
% ---------------------------------------------------------------------------------------------
XandY = @(tau,Mu) tau.*[rand_plusORminus_one(1,1)*sqrt(1 - Mu^2), Mu];



% reset the random number generator
rng('default');


% -------------------------------------------------------------
% *** Create a bunch of empty arrays to keep track of stuff ***
% -------------------------------------------------------------

% we have to keep track of how our photon travels within the medium
photon_tau_position = cell(1,N_photons);


% Store the values of the maximum penetration depth for each photon
max_depth = zeros(1, N_photons);
max_horizontal_position = zeros(1, N_photons);

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
% photon. We set these to be 0 to start in case the vector is empty when it
% is probed in a logical question. We then delete the zero value in the
% end, because matlab indices start with 1.
scatter_out_top_index = 0;
scatter_out_bottom_index = 0;
absorbed_index = 0;


% ----------------------------------------------------------------------
% ************************* START LIVE PLOT ****************************
% ----------------------------------------------------------------------


figure;
xline(0,'LineWidth',0.5,'Color','black')
hold on;
yline(0,'LineWidth',1.5,'Color','black')
hold on
yline(tau_y_upper_limit,'LineWidth',1.5,'Color','black')
hold on;
grid on; grid minor

% Set figure position and size
set(gcf, 'Position',[0 0 1200 700])

ylim([-1,tau_y_upper_limit+1])
xlim([-round(3*tau_y_upper_limit), round(3*tau_y_upper_limit)])
title('Photon position','Interpreter','latex')
xlabel('$\tau_x$','Interpreter','latex', 'FontSize',32);
ylabel('$\tau_y$','Interpreter','latex', 'FontSize',32);

set(gca,'YDir','reverse')
set(gca,'FontSize',20)


% Set the text box location
dim = [0.7 0.875 0.236083824157715 0.0985912584612054];

drawnow


% ----------------------------------------------------------------------
% *** Inject photons into the medium and keep track of what happens! ***
% ----------------------------------------------------------------------




for nn = 1:N_photons



    % Each photon injected has a probability of travelling a certain distance
    % before some event occurs. This event can be either scattering or
    % absorption.

    % when photons are injected into our medium, they start with a downward
    % direction

    % travelling downward is defined as positive increase in tau
    % travelling upwards is defined as a negative increase in tau


    % Set the transformation matrix to be the identity matrix for the first
    % photon. After that, for each scattering event we will multiply the
    % transformation matrix used in the previous transformation with the 
    % transformation of the lateset scattering event.
    M_transformation = eye(2);
    
    % We need to keep track of the current mu sample and the previous one. The
    % new coordinate system depends on the current mu and the previous mu. 
    mu_vector = cosd(solar_zenith_angle);
    %mu_vector = 0;   



    % --------------------------------------------------------
    % *****----- Roll dice and move photon forward -----******
    % --------------------------------------------------------

    % We start by drawing a value from our uniform distribution and
    % inserting it into the distribution that describes the probability of
    % travelling a distance tau (dT/dx relationship)
    tau_sample = tau(rand(1,1));

    % This tells us how far our photon travels before something happens.

    % Each new tau and mu define a distance travelled and a direction. And the
    % coordinate system for each new photon is defined by the previous tau and
    % mu. We need to keep track of the photon position with respect to the
    % previous vector so we know what kind of transformation to make. We
    % only need to keep track of the current and previous position. The
    % first column is the x position in tau space, the second column is the
    % y position in tau space.
    photon_position_in_new_coordniates = XandY(tau_sample, mu_vector(end));
    %photon_position_in_new_coordniates = [0, tau_sample];
    
    % we need to keep track of where the photon is within our medium
    % Start by giving each vector a 0 to show the starting point at time
    % t=0. The first column is the x position in tau space. The second
    % column is the y position in tau space.
    photon_tau_position{nn} = [0, 0;...
                               photon_position_in_new_coordniates];


    % -----------------------------------------
    % **** Has the photon left our medium? ****
    % -----------------------------------------

    % There is a chance we draw a number that causes the photon to transmit
    % through without any scattering events!

    if photon_tau_position{nn}(end,2)>tau_y_upper_limit

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
            % ----------------------------------------------
            % The photon is scattered back into our medium!!
            % ----------------------------------------------

            % First, set the distance travelled to be the boundary
            % condition value, since we are treating this as a hard
            % boundary
            photon_tau_position{nn}(end,2) = tau_y_upper_limit;

            % ------------------------------------------------------
            % ----- Check to see what layer is our photon in!! -----
            % ------------------------------------------------------
            if N_layers>1
                I_lessThan = find(layerBoundaries<photon_tau_position{nn}(end,2));

                % set the index describing which layer the photon is in
                index_layer = I_lessThan(end);

            else

                % the index layer is 1, since there is only 1
                index_layer = 1;

            end
            % ------------------------------------------------------



            % ---------------------------------------------------------
            % *** Determine which direction the photon scattered in ***
            % ---------------------------------------------------------

            % The photon has scattered back into the medium, because the
            % albedo is defined as the fraction of irradiance (light from
            % all directions) that is reflected. So we draw a mu value
            % until we have one that is between [-1,0] since we know our
            % photon back scatters into the medium

            % draw a random number and compute a random mu value but make
            % sure it is between -1 and 0
            mu_vector(end+1) = mu(g(index_layer), rand(1,1));

            while mu_vector(end)>=0
                mu_vector(end) = mu(g(index_layer), rand(1,1));
            end

            % Now that we have a value that represents a backscattering
            % event, record the new Y position

            % Is our photon continuing to move down, or is it moving
            % up, despite the angle between the direction of motion and
            % the plane parallel layers?


            % --------------------------------------------------------
            % *****----- Roll dice and move photon forward -----******
            % --------------------------------------------------------

            % Determine how far the photon travels after being scattered
            % back into our medium
            tau_sample = tau(rand(1,1));



            % ---------------------------------------------------------------------------------
            % ** Determine if the previous transformation was clockwise or coutner-clockwise **
            % ---------------------------------------------------------------------------------
            if photon_position_in_new_coordniates(end-1,1)<0
                % If the x position is less than 0, we perform a counter-clockwise transformation
                M_transformation = M_transformation * M_counterClockwise(mu_vector(end-1));
            else
                % If the x position is greater than 0, we perform a clockwise transformation
                M_transformation = M_transformation * M_clockwise(mu_vector(end-1));
            end


            % record where the photon along the y axis and along the x axis
            photon_tau_position{nn}(end+1,:) = photon_tau_position{nn}(end,:) + (M_transformation * photon_position_in_new_coordniates(end,:)')';

            % ------------------------------------------------------
            % ----- Check to see what layer is our photon in!! -----
            % ------------------------------------------------------
            if N_layers>1
                I_lessThan = find(layerBoundaries<photon_tau_position{nn}(end,2));

                % set the index describing which layer the photon is in
                index_layer = I_lessThan(end);

            else

                % the index layer is 1, since there is only 1
                index_layer = 1;

            end
            % ------------------------------------------------------


        else

            % Our photon was absorbed - but what we really mean is that our
            % photon has scattered out of the bottom of our medium and will
            % not return. We consider this event 'scattered out the bottom'


            % let's set the final depth_travelled value to be the
            % tau limit
            photon_tau_position{nn}(end,2) = tau_y_upper_limit;
            % Let's record which boundary the photon surpassed
            scatter_out_bottom = scatter_out_bottom + 1;

            % record the index
            scatter_out_bottom_index = [scatter_out_bottom_index, nn];



        end



    else

        % If the photon does not travel through the medium on the first tau
        % draw, or if it was scattered back into our medium due to the bottom
        % boundary, then we redraw a tau and continue its life.

        % ------------------------------------------------------
        % ----- Check to see what layer is our photon in!! -----
        % ------------------------------------------------------
        if N_layers>1
            I_lessThan = find(layerBoundaries<photon_tau_position{nn}(end,2));

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

            % draw a random mu value
            mu_vector(end+1) = mu(g(index_layer), rand(1,1)); 

            % Is our photon continuing to move down, or is it moving
            % up, despite the angle between the direction of motion and
            % the plane parallel layers?


            % -----------------------------------------------------------
            % ** Determine how far the photon travels after scattering **
            % -----------------------------------------------------------
            % draw another random variable and determine how far the photon
            % travels
            tau_sample = tau(rand(1,1));

            photon_position_in_new_coordniates(end+1,:) = XandY(tau_sample, mu_vector(end)); 


            % ---------------------------------------------------------------------------------
            % ** Determine if the previous transformation was clockwise or coutner-clockwise **
            % ---------------------------------------------------------------------------------
            if photon_position_in_new_coordniates(end-1,1)<0
                % If the x position is less than 0, we perform a counter-clockwise transformation
                M_transformation = M_transformation * M_counterClockwise(mu_vector(end-1));
            else
                % If the x position is greater than 0, we perform a clockwise transformation
                M_transformation = M_transformation * M_clockwise(mu_vector(end-1));
            end

            
            % record where the photon along the y axis and along the x axis
            photon_tau_position{nn}(end+1,:) = photon_tau_position{nn}(end,:) + (M_transformation * photon_position_in_new_coordniates(end,:)')';




            % -------------------------------------------------
            % **** Has the photon breached any boundaries? ****
            % -------------------------------------------------

            if photon_tau_position{nn}(end,2)<tau_y_lower_limit

                % The photon has scattered out the top of our medium
                % if so, let's set the final depth_travelled value to be 0
                photon_tau_position{nn}(end,2) = 0;

                % Let's record which boundary the photon surpassed
                scatter_out_top = scatter_out_top + 1;

                % record the index
                scatter_out_top_index = [scatter_out_top_index, nn];


                break

            elseif photon_tau_position{nn}(end,2)>tau_y_upper_limit

                % The photon has scattered out the bottom of our medium


                % ------------------------------------------------------
                % **** Did our photon reflect back into our medium? ****
                % ------------------------------------------------------

                % When our photon reaches the lower boudnary, we will check to see
                % if the photon is absorbed or scattered back into the medium by
                % comparing a random uniform number drawing with the albedo

                scat_or_abs = rand(1,1);

                if scat_or_abs<=albedo_maxTau
                    % The photon is scattered back into the medium!!

                    % First, set the distance travelled to be the boundary
                    % condition value, since we are treating this as a hard
                    % boundary
                    photon_tau_position{nn}(end,2) = tau_y_upper_limit;

                    % ---------------------------------------------------------
                    % *** Determine which direction the photon scattered in ***
                    % ---------------------------------------------------------

                    % The photon has scattered back into the medium, because the
                    % albedo is defined as the fraction of irradiance (light from
                    % all directions) that is reflected. So we draw a mu value
                    % until we have one that is between [-1,0] since we know our
                    % photon back scatters into the medium

                    % draw a random number and compute a random mu value but make
                    % sure it is between -1 and 0
                    mu_vector(end+1) = mu(g(index_layer), rand(1,1));

                    while mu_vector(end)>0
                        mu_vector(end) = mu(g(index_layer), rand(1,1));
                    end

                    % Now that we have a value that represents a backscattering
                    % event, record the new Y position

                    % Is our photon continuing to move down, or is it moving
                    % up, despite the angle between the direction of motion and
                    % the plane parallel layers?


                    % --------------------------------------------------------
                    % *****----- Roll dice and move photon forward -----******
                    % --------------------------------------------------------

                    % Determine how far the photon travels after being scattered
                    % back into our medium
                    tau_sample = tau(rand(1,1));

                    % Compute the new position
                    photon_position_in_new_coordniates(end+1,:) = XandY(tau_sample, mu_vector(end));

                    % We have to make sure this next photon path length
                    % turns in the positive x direction in the new
                    % coordniate space in order to keep in in the
                    % medium
                    photon_position_in_new_coordniates(end,1) = abs(photon_position_in_new_coordniates(end,1));

                    % ---------------------------------------------------------------------------------
                    % ** Determine if the previous transformation was clockwise or coutner-clockwise **
                    % ---------------------------------------------------------------------------------
                    if photon_position_in_new_coordniates(end-1,1)<0
                        % If the x position is less than 0, we perform a counter-clockwise transformation
                        M_transformation = M_transformation * M_counterClockwise(mu_vector(end-1));

                        
                    else
                        % If the x position is greater than 0, we perform a clockwise transformation
                        M_transformation = M_transformation * M_clockwise(mu_vector(end-1));
                    end


                    % record where the photon along the y axis and along the x axis
                    photon_tau_position{nn}(end+1,:) = photon_tau_position{nn}(end,:) + (M_transformation * photon_position_in_new_coordniates(end,:)')';


                   

                    % ------------------------------------------------------
                    % ----- Check to see what layer is our photon in!! -----
                    % ------------------------------------------------------
                    if N_layers>1
                        I_lessThan = find(layerBoundaries<photon_tau_position{nn}(end,2));

                        % set the index describing which layer the photon is in
                        index_layer = I_lessThan(end);

                    else

                        % the index layer is 1, since there is only 1
                        index_layer = 1;

                    end
                    % ------------------------------------------------------



                else
                    % ---------------------------------------------
                    % Our photon was absorbed by the lower boundary
                    % ---------------------------------------------

                    % let's set the final depth_travelled value to be the
                    % tau limit
                    photon_tau_position{nn}(end,2) = tau_y_upper_limit;
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
                I_lessThan = find(layerBoundaries<photon_tau_position{nn}(end,2));

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

    if photon_tau_position{nn}(end,2)>tau_y_lower_limit && photon_tau_position{nn}(end,2)<tau_y_upper_limit

        % Make sure both conditions are met!
        % The random number, scat_or_abs is greater than the ssa
        % the photon is still within the boundaries of our defined medium

        absorbed = absorbed + 1;

        % record the index
        absorbed_index = [absorbed_index, nn];

        % record the depth of absorption
        depth_absorbed = [depth_absorbed, photon_tau_position{nn}(end,2)];

    end



    % Save the maximum y position reached.
    max_depth(nn) = max(photon_tau_position{nn}(:,2));
    max_horizontal_position(nn) = max(photon_tau_position{nn}(:,1));

    % Save the number of scattering events
    % (1) If the final event was absorption, we substract two from the
    % depth_travelled vector
    % (2) If the photon was not absorbed, it left our medium. We ignore the
    % first entry in depth_travlled, since it is a boundary condition.
    % We also must ignore the last entry, since we are no longer
    % tracking scattering events outside the medium
    number_of_scattering_events(nn) = size(photon_tau_position{nn},1)-2;



    % ---------------------------
    % ***** Make live plot! *****
    % ---------------------------
    % Create plot object
    p = zeros(1, size(photon_tau_position{nn},1));

    for ff = 1:size(photon_tau_position{nn},1)-1

        % Plot photon path
        p(ff) = plot([photon_tau_position{nn}(ff,1), photon_tau_position{nn}(ff+1,1)],...
            [photon_tau_position{nn}(ff,2), photon_tau_position{nn}(ff+1,2)], ...
            'LineWidth',6);
        hold on


        drawnow
        pause(0.2)



    end

    % plot latest tally of absorbed and scattered photons
    if exist('t', 'var')
        delete(t)
    end

    texBox_str = {['Num Absorbed = ',num2str(absorbed)],...
        ['Num Scattered Out Top = ',num2str(scatter_out_top )],...
        ['Num Scattered Out Bottom = ',num2str(scatter_out_bottom)]};

    t = annotation('textbox',dim,'string',texBox_str,'Interpreter','latex');
    t.Color = 'black';
    t.FontSize = 15;
    t.FontWeight = 'bold';
    t.EdgeColor = 'black';
    t.FitBoxToText = 'on';

    drawnow

    % Clear figure and draw a new photon path
    delete(p)





end


%%

% ----------------------------------------------------------------------------------
% *** Divide the tau space into a grid and compute reflectance and transmittance ***
% ----------------------------------------------------------------------------------

% Set up the tau grid first

N_bins = 500;

if tau_y_upper_limit==inf
    binEdges = logspace(-3, ceil(log10(max(max_depth))),N_bins+1);

else
    binEdges = linspace(tau_y_lower_limit, tau_y_upper_limit,N_bins+1);
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

    

    % Compute whether or not the photon was moving in the positive y
    % direction or the negative y direction. A value of 1 tells the code
    % that the photon was increasing in tau along the y axis, and a value
    % of negative 1 tells the code the photon was decreasing in tau along
    % the y axis. The first value is always 1 since the photons starts its
    % journey by moving down, increasing its value of tau along the y axis
    y_direction = [1; sign(diff(photon_tau_position{nn}(:,2)))];


    for tt = 1:size(photon_tau_position{nn},1)-1

        % Check to see if the photon is moving down or up and check to see
        % if its continuing along the same direction, or if it's switched
        % directions
        
        %[nn,tt]

        % ----------------------------------------
        % **** Photon continuing to move down ****
        % ----------------------------------------
        if y_direction(tt+1)==1 && y_direction(tt)==1
            % *** the photon is moving down ***
            % Check the see which tau bins the photon moves through. Tally each
            % bin that is found


            % If the photon is moving down and was already heading
            % down, we don't need to account for an additional bin

            % If the photon is moving down, then when it cross bin-edge 3,
            % its in bucket 3.
            bins_photon_moves_through = binEdges>=photon_tau_position{nn}(tt,2) & ...
                binEdges<photon_tau_position{nn}(tt+1,2);

            N_counts_moving_down(bins_photon_moves_through) = N_counts_moving_down(bins_photon_moves_through) +1;


            % --------------------------------------------
            % **** Photon moving down after moving up ****
            % --------------------------------------------
        elseif y_direction(tt+1)==1 && y_direction(tt)==-1
            % In this case the photon has turned around in the bin it was
            % left off in. Therefore we have to count this as a photon
            % moving down in this bin.
            
            % Photon is moving down, so the starting bin will be smaller
            % than, or equal to, the final bin
            bin_edges_less_than_start_position = find(binEdges<=photon_tau_position{nn}(tt,2));
            
            % This bin should be greater than, or equal to, the starting
            % bin
            bin_edges_less_than_end_position = find(binEdges<photon_tau_position{nn}(tt+1,2));

            bins_photon_moves_through = zeros(1,N_bins+1);
            bins_photon_moves_through(bin_edges_less_than_start_position(end):bin_edges_less_than_end_position(end)) = 1;
            bins_photon_moves_through = logical(bins_photon_moves_through);
            

            N_counts_moving_down(bins_photon_moves_through) = N_counts_moving_down(bins_photon_moves_through) +1;



            % --------------------------------------
            % **** Photon continuing to move up ****
            % --------------------------------------
        elseif y_direction(tt+1)==-1 && y_direction(tt)==-1
            % *** the photon is moving up ***
            % When this is true, depth_travlled{nn}(tt+1) is always less
            % than depth_travelled{nn}(tt) so we need to determine the tau
            % bins that the photon passes through between
            % depth_travlled{nn}(tt) and depth_travlled{nn}(tt+1)


            % If the photon is continuing up then we only have to count
            % the bins it passes through
            
            % Photon is moving up, so the starting bin will be larger
            % than, or equal to, the final bin
            bin_edges_less_than_start_position = find(binEdges<=photon_tau_position{nn}(tt,2));
            
            % This bin should be less than, or equal to, the starting
            % bin
            bin_edges_less_than_end_position = find(binEdges<=photon_tau_position{nn}(tt+1,2));

            bins_photon_moves_through = zeros(1,N_bins+1);

            % let's check to see if the photon is still in the same bin. If
            % it is, we've already accounted for it's upward motion, so we
            % don't want to double count this. 
            if bin_edges_less_than_end_position(end)==bin_edges_less_than_start_position(end)
                % If this is true, the photon is continuing to move upward
                % in the same bin, and we don't count it
            else

                % If the photon is continuing to move upwards, that means
                % the previous loop already accounted for the bin it's
                % currently in, either because the photon turned around and
                % ended up in the current bin, or because the photon
                % continued to move upwards into the current bin. The logic
                % above will count the bin its currently in, so we need to
                % remove this value. We do this by subtracting 1 from the
                % star position.
                bins_photon_moves_through(bin_edges_less_than_end_position(end):bin_edges_less_than_start_position(end-1)) = 1;           

            end

            bins_photon_moves_through = logical(bins_photon_moves_through);

            N_counts_moving_up(bins_photon_moves_through) = N_counts_moving_up(bins_photon_moves_through) +1;




            % --------------------------------------------
            % **** Photon moving up after moving down ****
            % --------------------------------------------
        elseif y_direction(tt+1)==-1 && y_direction(tt)==1

            % If the photon switched direction, we have to account for
            % the final bin it ends up in

            % Photon is moving up, so the starting bin will be larger
            % than, or equal to, the final bin
            bin_edges_less_than_start_position = find(binEdges<=photon_tau_position{nn}(tt,2));

            % The final bin will be smaller than, or equal to, to starting
            % bin
            bin_edges_less_than_end_position = find(binEdges<=photon_tau_position{nn}(tt+1,2));

            bins_photon_moves_through = zeros(1,N_bins+1);
            bins_photon_moves_through(bin_edges_less_than_end_position(end):bin_edges_less_than_start_position(end)) = 1;
            bins_photon_moves_through = logical(bins_photon_moves_through);
            

            % Some photons in this category will reflect off the bottom
            % boudnary. If this happened, we simply need to remove the
            % logical true value for the binEdge equal to tau_upper_limit.
            % If we don't we get an error, and all we need to keep track of
            % is whether or not the photon passed through this bin
            if bins_photon_moves_through(end)==true && (photon_tau_position{nn}(tt,2)==tau_y_upper_limit || photon_tau_position{nn}(tt+1,2)==tau_y_upper_limit)
                bins_photon_moves_through(end) = 0;
            end




            N_counts_moving_up(bins_photon_moves_through) = N_counts_moving_up(bins_photon_moves_through) +1;



        else


            error([newline,'Im not sure what the photon is doing',newline])

        end


    end

end


%%



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
% *** delete the zero indices, they were simply place holders ***
scatter_out_top_index(1) = [];
scatter_out_bottom_index(1) = [];
absorbed_index(1) = [];
final_state.scatter_out_top_INDEX = scatter_out_top_index;
final_state.scatter_out_bottom_INDEX = scatter_out_bottom_index;
final_state.absorbed_INDEX = absorbed_index;


% output other photon tracking variables that I've computed

% keep track of the maximum penetration depth of each photon
photon_tracking.maxDepth = max_depth;
photon_tracking.max_horizontal_position = max_horizontal_position;

% Count how many times a scatter/absorption event occurs during a photons
% lifetime. A value of 1 means it was absorbed right away or transmitted
% through
photon_tracking.number_of_scattering_events = number_of_scattering_events;

% Record the depths at which absorption occured
photon_tracking.absorbed_depth = depth_absorbed;







end