%% This Monte Carlo code compute 2-stream radiative transfer

% Both scattering and Absorption are included


% By Andrew John Buggee

%% For an Infintely Thick Medium...
% We start by injecting a single photon into our medium.
clear variables
% Let's define T as the random variable that describes the probability
% of a photon travelling an optical depth T. P(T) is that probability

% Let's define x to be the uniform PDF where x is between [0,1].

% define the relationship between x and tau

tau = @(x) -log(1 - x);

tau_lower_limit = 0;
tau_upper_limit = 7;


% We start by injecting a single photon into our medium.

% define the number of photons that will run through our simulation
N_photons = 10000;

% define the wavelength
%wavelength = 1950;            % nanometers - ssa ~ 0.9
wavelength = 500;            % nanometers - conservative scattering

% -----------------------------------------------------------
% ************* Define Scattering Properties ****************
% -----------------------------------------------------------

% Assuming a pure homogenous medium composed of a single substance
% Define the single scattering albedo of the substance in question
% assuming water droplets, what is the radius?
r = 10;                     % microns
mie_properties_water = interp_mie_computed_tables([wavelength, r],'mono',false);
ssa = mie_properties_water(6);
g = mie_properties_water(7);

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
final_state.scatter_out_top = 0;
final_state.scatter_out_bottom = 0;
final_state.absorbed = 0;

%% Run the Monte Carlo Calculations in Parallel


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

    if depth_travelled{nn}(end)>tau_upper_limit

        % The photon has scattered out the bottom of our medium with a
        % single scattering or absorption event!

        % let's set the final depth_travelled value to be the
        % tau limit
        depth_travelled{nn}(end) = tau_upper_limit;
        % Let's record which boundary the photon surpassed
        final_state.scatter_out_bottom = final_state.scatter_out_bottom + 1;



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
                final_state.scatter_out_top = final_state.scatter_out_top + 1;


                break

            elseif depth_travelled{nn}(end)>tau_upper_limit

                % The photon has scattered out the bottom of our medium
                % If so, let's set the final depth_travelled value to be the
                % tau limit
                depth_travelled{nn}(end) = tau_upper_limit;
                % Let's record which boundary the photon surpassed
                final_state.scatter_out_bottom = final_state.scatter_out_bottom + 1;


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

    if scat_or_abs>ssa && depth_travelled{nn}(end)>tau_lower_limit && depth_travelled{nn}(end)<tau_upper_limit

        % Make sure all three conditions are met!
        % The random number, scat_or_abs is greater than the ssa
        % the photon is still within the boundaries of our defined medium

        final_state.absorbed = final_state.absorbed + 1;

    end



    % Save the maximum depth reached.
    maxDepth(nn) = max(depth_travelled{nn});

    % Save the number of scattering events
    number_of_scattering_events(nn) = length(depth_travelled{nn})-1;


end


%% Let's Plot the number of scattering events for each photon

% The number of scattering events is equal to the length of the
% depth_travelled vector - 1

nBins = 75;
% scattering events can only take on integer values
binEdges = unique(floor(logspace(0, ceil(log10(max(number_of_scattering_events))),nBins+1)));

figure; subplot(2,1,1)
h = histogram(number_of_scattering_events,nBins,"BinEdges",binEdges);
set(gca,'YScale','log')
set(gca, 'XScale', 'log')
%set(gca, 'YDir','reverse')
grid on; grid minor
xlabel('Number of Scattering Events','Interpreter','latex');
ylabel('Counts','Interpreter','latex')
title([num2str(N_photons),' photons,   ',...
    '$\lambda$ = ',num2str(wavelength), ' $nm$',...
    ', $\tilde{\omega}$ = ', num2str(ssa), ...
    ', $g$ = ', num2str(g)], 'Interpreter','latex')
set(gcf, 'Position',[0 0 1000 630])

% Let's plot the maximum depth of each photon

nBins = 50;
binEdges = logspace(-3, ceil(log10(max(maxDepth))),nBins+1);

subplot(2,1,2)
h = histogram(maxDepth,nBins,"BinEdges",binEdges,'Orientation','horizontal');
set(gca,'YScale','log')
set(gca, 'XScale', 'log')
set(gca, 'YDir','reverse')
grid on; grid minor
xlabel('Counts','Interpreter','latex');
ylabel('$\tau$','Interpreter','latex')
title('Max depth reached by each photon', 'Interpreter','latex')
ylim([10^(-3), 10^(5)])

% Make a text box with the parameters used
% dim = [0 0.5 0 0];
%
% texBox_str = ['$\lambda$ = ',num2str(wavelength), ' $nm$',...
%               '$, \tilde{\omega}$ = ', num2str(ssa), ...
%               ', $g$ = ', num2str(g)];
% t = annotation('textbox',dim,'string',texBox_str,'Interpreter','latex');
% t.Color = 'white';
% t.FontSize = 12;
% t.EdgeColor = 'white';
% t.FitBoxToText = 'on';

%% What were the end states of each photon?


figure;
X = categorical({'Scatter out top','Scatter out Bottom','Absorbed'});
b = bar(X,[final_state.scatter_out_top, final_state.scatter_out_bottom, final_state.absorbed]);
set(gca,'TickLabelInterpreter', 'latex','FontSize',25);
set(gca, 'TitleFontSizeMultiplier',1.3)
bar_label_xLocation = b.XEndPoints;
bar_label_yLocation = b.YEndPoints;
bar_labels = string(b.YData./N_photons);
text(bar_label_xLocation,bar_label_yLocation,bar_labels,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',20,'FontWeight','bold','Interpreter','latex')
grid on; grid minor
ylabel('Counts','Interpreter','latex')
title([num2str(N_photons),' photons,   ',...
    '$\lambda$ = ',num2str(wavelength), ' $nm$',...
    ', $\tilde{\omega}$ = ', num2str(ssa), ...
    ', $g$ = ', num2str(g)], 'Interpreter','latex')
set(gcf, 'Position',[0 0 1000 630])



%% What if I look at ALL logged movements, both up and down, for each photon?

% Set up the tau grid first

N_bins = 100;
edge_numbers = 1:N_bins+1;

if tau_upper_limit==inf
    binEdges = logspace(-3, ceil(log10(max(maxDepth))),N_bins+1);

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
            % down, we don't want to additionally count that bin

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

            N_counts_moving_up(bins_photon_moves_through) = N_counts_moving_up(bins_photon_moves_through) +1;



        else


            error([newline,'Im not sure what the photon is doing',newline])

        end


    end

end




% Let's plot the histogram of photons moving up and down

figure;
subplot(2,1,1)
h = histogram("BinEdges",binEdges,"BinCounts",N_counts_moving_up,...
    'Orientation','horizontal');
set(gca,'YScale','log')
set(gca, 'XScale', 'log')
set(gca, 'YDir','reverse')
grid on; grid minor
xlabel('Counts','Interpreter','latex');
ylabel('$\tau$','Interpreter','latex')
title('Photons moving up', 'Interpreter','latex')
ylim([10^(-3), 10^(5)])
set(gcf, 'Position',[0 0 1000 630])


subplot(2,1,2)
h = histogram("BinEdges",binEdges,"BinCounts",N_counts_moving_down,...
    'Orientation','horizontal');
set(gca,'YScale','log')
set(gca, 'XScale', 'log')
set(gca, 'YDir','reverse')
grid on; grid minor
xlabel('Counts','Interpreter','latex');
ylabel('$\tau$','Interpreter','latex')
title('Photons moving down', 'Interpreter','latex')
ylim([10^(-3), 10^(5)])


%% How does this compare with the analytical solutions?


% Check to see if there is absorption

if ssa<1

    % Next, check to see if our layer is infinitely thick, or has a finite
    % thickness

    if tau_upper_limit == inf

        K = sqrt((1 - ssa)*(1 - g*ssa));
        R_inf = (sqrt(1-ssa*g) - sqrt(1 - ssa))/(sqrt(1-ssa*g) + sqrt(1 - ssa));

        photon_fraction_up = @(tau) R_inf * exp(-K*tau);
        photon_fraction_down = @(tau) exp(-K*tau);

        % lets plot some range of tau
        tau = logspace(-2,2,100);

        C1 = [8.492724501845559e-01     5.437108503990062e-02     9.681090252965144e-01];
        C2 = [4.843239544631898e-01     7.049687413641152e-01     6.568033076805693e-01];

        figure; plot(tau, photon_fraction_up(tau),'Color',C1)
        hold on;
        plot(tau, photon_fraction_down(tau), 'Color',C2)
        plot(binEdges(1:end-1)+diff(binEdges)/2, N_counts_moving_up./N_photons,'Color',C1,'LineStyle','--')
        plot(binEdges(1:end-1)+diff(binEdges)/2, N_counts_moving_down./N_photons,'Color',C2,'LineStyle','--')
        grid on; grid minor
        ylabel('$F/F_0$','Interpreter','latex');
        xlabel('$\tau$','Interpreter','latex');
        title({'Comparing analytical 2 stream with Monte Carlo',...
            'for an infinitely thick absorbing medium'}, 'Interpreter','latex')
        legend('$F_{\uparrow}/F_0$ analytical','$F_{\downarrow}/F_0$ analytical',...
            '$F_{\uparrow}/F_0$ Monte Carlo','$F_{\downarrow}/F_0$ Monte Carlo','Interpreter','latex',...
            'Location','best')
        set(gcf, 'Position',[0 0 1000 630])

        % Make a text box with the parameters used
        dim = [0.65 0.65 0 0];

        texBox_str = ['$R_{\infty} = \frac{\sqrt{1 - g \tilde{\omega}}\, - \,\sqrt{1 - \tilde{\omega}}}',...
            '{\sqrt{1 - g \tilde{\omega}}\, + \,\sqrt{1 - \tilde{\omega}}}$'];
        t = annotation('textbox',dim,'string',texBox_str,'Interpreter','latex');
        t.Color = 'white';
        t.FontSize = 25;
        t.FontWeight = 'bold';
        t.EdgeColor = 'white';
        t.FitBoxToText = 'on';


    elseif tau_upper_limit>0 && tau_upper_limit<inf
        
        % Define the albedo of the lower boundary in our model
        albedo_lowerBoundary = 0;
        
        % K is defined in Bohren and Clothiaux (eq. 5.70)
        K = sqrt((1 - ssa)*(1 - g*ssa));
        % Define the reflectivity at the top of our layer, the photons that
        % scatter out the cloud top
        R_inf = (sqrt(1-ssa*g) - sqrt(1 - ssa))/(sqrt(1-ssa*g) + sqrt(1 - ssa));

        % Define the constants
        A = (R_inf - albedo_lowerBoundary)*exp(-K*tau_upper_limit)/...
            (R_inf*(R_inf - albedo_lowerBoundary)*exp(-K*tau_upper_limit) - (1 - albedo_lowerBoundary*R_inf)*exp(K*tau_upper_limit));

        B = -(1 - R_inf*albedo_lowerBoundary)*exp(K*tau_upper_limit)/...
            (R_inf*(R_inf - albedo_lowerBoundary)*exp(-K*tau_upper_limit) - (1 - albedo_lowerBoundary*R_inf)*exp(K*tau_upper_limit));

        photon_fraction_up = @(tau) A*exp(K*tau) + B*R_inf*exp(-K*tau);

                
        photon_fraction_down = @(tau) A*R_inf*exp(K*tau) + B*exp(-K*tau);


        % lets plot some range of tau
        tau = linspace(tau_lower_limit,tau_upper_limit,100);

        C1 = [8.492724501845559e-01     5.437108503990062e-02     9.681090252965144e-01];
        C2 = [4.843239544631898e-01     7.049687413641152e-01     6.568033076805693e-01];

        figure; plot(tau, photon_fraction_up(tau),'Color',C1)
        hold on;
        plot(tau, photon_fraction_down(tau), 'Color',C2)
        plot(binEdges(1:end-1)+diff(binEdges)/2, N_counts_moving_up./N_photons,'Color',C1,'LineStyle',':')
        plot(binEdges(1:end-1)+diff(binEdges)/2, N_counts_moving_down./N_photons,'Color',C2,'LineStyle',':')
        grid on; grid minor
        ylabel('$F/F_0$','Interpreter','latex');
        xlabel('$\tau$','Interpreter','latex');
        title({'Comparing analytical 2 stream with Monte Carlo',...
            'for an absorbing medium of finite thickness'}, 'Interpreter','latex')
        legend('$F_{\uparrow}/F_0$ analytical','$F_{\downarrow}/F_0$ analytical',...
            '$F_{\uparrow}/F_0$ Monte Carlo','$F_{\downarrow}/F_0$ Monte Carlo','Interpreter','latex',...
            'Location','best')
        set(gcf, 'Position',[0 0 1000 630])


    end


elseif ssa==1

    % Next, check to see if our layer is infinitely thick, or has a finite
    % thickness

    if tau_upper_limit == inf


        error([newline,'I dont know what to do with a layer of finte thickness and conservative scattering.',newline])

    elseif tau_upper_limit>0 && tau_upper_limit<inf

        % For a non-absorbing layer of finite thickness, the analytical
        % solutions to the two stream radiative transfer equations are...

        R_atTop = (tau_upper_limit * (1 - g)/2)/(1 + tau_upper_limit*(1 - g)/2);

        T_atBottom = 1/(1 + tau_upper_limit*(1 - g)/2);

        % lets plot some range of tau
        tau = linspace(tau_lower_limit,tau_upper_limit,100);

        C1 = [8.492724501845559e-01     5.437108503990062e-02     9.681090252965144e-01];
        C2 = [4.843239544631898e-01     7.049687413641152e-01     6.568033076805693e-01];

        figure; xline(R_atTop,'Color',C1,'LineWidth',4)
        hold on;
        xline(T_atBottom,'Color',C2,'LineWidth',4)
        plot(N_counts_moving_up./N_photons,binEdges(1:end-1)+diff(binEdges)/2,'Color',C1,'LineStyle','--')
        plot(N_counts_moving_down./N_photons,binEdges(1:end-1)+diff(binEdges)/2, 'Color',C2,'LineStyle','--')
        grid on; grid minor
        xlabel('$F/F_0$','Interpreter','latex');
        ylabel('$\tau$','Interpreter','latex');
        set(gca,'YDir','reverse')
        title({'Comparing analytical 2 stream with Monte Carlo',...
            'Conservative scattering for a finite layer'},...
            'Interpreter','latex')
        legend('$R_{\infty}$ analytical','$T(\bar{\tau})$ analytical',...
            '$F_{\uparrow}/F_0$ Monte Carlo','$F_{\downarrow}/F_0$ Monte Carlo','Interpreter','latex',...
            'Location','best')
        set(gcf, 'Position',[0 0 1000 630])

        % Make a text box with the parameters used
        dim = [0.41 0.85 0 0];

        texBox_str = {'$R(\tau = 0) = \frac{\bar{\tau}(1 - g)/2}{1 + \bar{\tau}(1 - g)/2}$',...
            '$T(\tau = \bar{\tau}) = \frac{1}{1 + \bar{\tau}(1 - g)/2}$'};
        t = annotation('textbox',dim,'string',texBox_str,'Interpreter','latex');
        t.Color = 'white';
        t.FontSize = 25;
        t.FontWeight = 'bold';
        t.EdgeColor = 'white';
        t.FitBoxToText = 'on';


    end




end



%% How does the analytical reflectivity change with wavelength?

wavelength = 100:2300;           % nanometers -

% Assuming a pure homogenous medium composed of a single substance
% Define the single scattering albedo of the substance in question
% assuming water droplets, what is the radius?
r = 10;                     % microns
mie_properties_water = interp_mie_computed_tables([wavelength',...
    linspace(r,r,length(wavelength))'],'mono',false);
ssa = mie_properties_water(:,6);
g = mie_properties_water(:,7);

R_inf = (sqrt(1-ssa.*g) - sqrt(1 - ssa))./(sqrt(1-ssa.*g) + sqrt(1 - ssa));

figure; plot(wavelength, R_inf)
grid on; grid minor
ylabel('$\frac{F_{\uparrow}(\tau=0)}{F_0}$','Interpreter','latex');
xlabel('Wavelength $(nm)$','Interpreter','latex');
title('Reflectivity - Theoretical 2 stream - infinitely absorbing medium', 'Interpreter','latex')
set(gcf, 'Position',[0 0 1000 630])


