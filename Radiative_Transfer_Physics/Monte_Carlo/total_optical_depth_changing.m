% For loop for different total optical depths


%%
% We start by injecting a single photon into our medium.
clear variables
% Let's define T as the random variable that describes the probability
% of a photon travelling an optical depth T. P(T) is that probability

% Let's define x to be the uniform PDF where x is between [0,1].

% define the relationship between x and tau

tau = @(x) -log(1 - x);

tau_lower_limit = 0;
tau_upper_limit = 1:100;


% We start by injecting a single photon into our medium.

% define the number of photons that will run through our simulation
N_photons = 5000;


% ***** Custom ssa and g *****
ssa = 0.9;                % between 0 and 1
g = 0.85;               % between -1 and 1
% ----------------------------

% reset the random number generator
rng('default');

% we have to keep track of how our photon travels within the medium
depth_travelled = cell(1,N_photons);

% Save all of the direction changes so we can easily compute upwards and
% downwards irradiance
direction = cell(1, N_photons);

% Store the values of the maximum penetration depth for each photon
maxDepth = zeros(1, N_photons);

% Store the values of the number of scattering events for each photon
number_of_scattering_events = zeros(1, N_photons);

% Let's track the final state of each photon
% do they end up leaving the medium? Do they scatter out the top of the
% bottom? Or are they absorbed along their journey?

final_state.scatter_out_top = zeros(1, length(tau_upper_limit));
final_state.scatter_out_bottom = zeros(1, length(tau_upper_limit));
final_state.absorbed = zeros(1, length(tau_upper_limit));

%% Run the Monte Carlo Calculations in Parallel

for tt = 1:length(tau_upper_limit)

    for nn = 1:N_photons



        % Each photon injected has a probability of travelling a certain distance
        % before some event occurs. This event can be either scattering or
        % absorption.

        % when photons are injected into our medium, they start with a downward
        % direction

        % travelling downward is defined as positive increase in tau
        % travelling upwards is defined as a negative increase in tau

        % Start with two values, the first 1 tells us the photon was moving
        % down at t = 0, the second says it was moving down until the tau
        % that is sampled
        direction{nn} = [1,1];


        % --------------------------------------------------------
        % *****----- Roll dice and move photon forward -----******
        % --------------------------------------------------------

        % We start by drawing a value from our uniform distribution
        x_draw = rand(1,1);
        % Next we plug this in to our dT/dx relationship
        tau_sample = tau(x_draw);

        % This tells us how far our photon travels before something happens.

        % we need to keep track of where the photon is within our medium
        % Start by giving each vector a 0 to show the starting point at time
        % t=0
        depth_travelled{nn} = [0, tau_sample];

        % -------------------------------------------------
        % **** Has the photon breached any boundaries? ****
        % -------------------------------------------------

        % There is a chance we draw a number that causes the photon to transmit
        % through without any scattering events!

        if depth_travelled{nn}(end)>tau_upper_limit(tt)

            % The photon has scattered out the bottom of our medium with a
            % single scattering or absorption event!

            % let's set the final depth_travelled value to be the
            % tau limit
            depth_travelled{nn}(end) = tau_upper_limit(tt);
            % Let's record which boundary the photon surpassed
            final_state.scatter_out_bottom(tt) = final_state.scatter_out_bottom(tt) + 1;



        else

            % If the photon does not travel through the medium on the first tau
            % draw, then it is either scattered or absorbed


            % --------------------------------------------------------
            % ***---- Roll dice. Did photon scatter or absorb? ---****
            % --------------------------------------------------------

            % To determine if the photon scattered or absorbed, draw a random
            % number from a uniform distribution. If this number is less than or
            % equal to the single scattering albedo, the photon scattered.
            % Otherwise we assume it absorbed

            scat_or_abs = rand(1,1);

            while scat_or_abs<=ssa

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
                    if g_sample <= ((1+g)/2)
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
                    if g_sample <= ((1+g)/2)
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
                    final_state.scatter_out_top(tt) = final_state.scatter_out_top(tt) + 1;


                    break

                elseif depth_travelled{nn}(end)>tau_upper_limit(tt)

                    % The photon has scattered out the bottom of our medium
                    % If so, let's set the final depth_travelled value to be the
                    % tau limit
                    depth_travelled{nn}(end) = tau_upper_limit(tt);
                    % Let's record which boundary the photon surpassed
                    final_state.scatter_out_bottom(tt) = final_state.scatter_out_bottom(tt) + 1;


                    break

                end


                % draw a new random number! Did our photon scatter or absorb?
                scat_or_abs = rand(1,1);


            end


        end



        % --------------------------------------------------
        % **** Did the photon perish due to absorption? ****
        % --------------------------------------------------
        % If the photon ended in absorption, record it!

%         if scat_or_abs>ssa && depth_travelled{nn}(end)>tau_lower_limit && depth_travelled{nn}(end)<tau_upper_limit(tt)
% 
%             % Make sure all three conditions are met!
%             % The random number, scat_or_abs is greater than the ssa
%             % the photon is still within the boundaries of our defined medium
% 
%             final_state.absorbed(tt) = final_state.absorbed(tt) + 1;
% 
%         end



        % Save the maximum depth reached.
        maxDepth(nn) = max(depth_travelled{nn});

        % Save the number of scattering events
        number_of_scattering_events(nn) = length(depth_travelled{nn})-1;


    end

end

%%

figure; 
plot(tau_upper_limit, final_state.scatter_out_bottom./N_photons)
hold on
ylabel('Transmission')
yyaxis right
hold on
plot(tau_upper_limit, final_state.scatter_out_top./N_photons)
xlabel('Total Optical Depth')
ylabel('Reflectivity')
grid on; grid minor
title(['$\tilde{\omega} = $', num2str(ssa), '$\;\;\;\; g = $',...
    num2str(g), '$\;\;\;\; N_{photons} = $', num2str(N_photons)], ...
    'Interpreter','latex', 'FontSize', 25)

