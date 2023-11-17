%% Plot the Liquid Water Content and the Effective Radius versus Altitude 


% By Andrew John Buggee
%%

function [] = plot_LWC_and_re_vs_opticalDepth(vert_profiles, indices, normalize_opticalDepth)

% define the denisty of liquid water
density = 1e6;          % g/m^3

% Plot the number of curves within the vert_profiles structure

N_cuvres = length(indices);

figure;
for nn = 1:N_cuvres

    tau = zeros(1, length(vert_profiles.altitude{indices(nn)}));

    % Profiles measured while the plane was descending will start with values
    % at the cloud top
    
    if (vert_profiles.altitude{indices(nn)}(end)-vert_profiles.altitude{indices(nn)}(1))>0
        % This profile is ascending, meaning the first data points are at
        % the cloud bottom, when tau is largest, since tau is defined from
        % top to bottom
        starting_idx = length(tau)+1;

        % compute optical depth
        % Let's assume the extinction coefficient is 2, meaning the radius is
        % much larger than the incident light
        for ii = length(vert_profiles.altitude{indices(nn)}):-1:2
            tau(starting_idx-ii) = 1/2 * 1/density * trapz(vert_profiles.altitude{indices(nn)}(1:ii),...
                            vert_profiles.lwc{indices(nn)}(1:ii)./(vert_profiles.re{indices(nn)}(1:ii)*1e-6));
        end

    elseif (vert_profiles.altitude{indices(nn)}(end)-vert_profiles.altitude{indices(nn)}(1))<0
        % This profile is descending, meaning the first data points are
        % measured at cloud top
        
        % compute optical depth
        % Let's assume the extinction coefficient is 2, meaning the radius is
        % much larger than the incident light
        for ii = 2:length(vert_profiles.altitude{indices(nn)})
            tau(ii) = -1/2 * 1/density * trapz(vert_profiles.altitude{indices(nn)}(1:ii),...
                            vert_profiles.lwc{indices(nn)}(1:ii)./(vert_profiles.re{indices(nn)}(1:ii)*1e-6));
        end
    end

    
    % Because the data is oriented by values at the cloud bottom first
    % the tau vector goes from the largest values to 0
    
    if normalize_opticalDepth==false

        subplot(1,2,1); plot(vert_profiles.lwc{indices(nn)}, tau)
        hold on

        % next plot the effective radius
        subplot(1,2,2); plot(vert_profiles.re{indices(nn)}, tau);
        hold on

    else

        subplot(1,2,1); plot(vert_profiles.lwc{indices(nn)}, tau./max(tau))
        hold on

        % next plot the effective radius
        subplot(1,2,2); plot(vert_profiles.re{indices(nn)}, tau./max(tau));
        hold on

    end



end

% Make each subplot pretty
subplot(1,2,1)
set(gca, 'ydir', 'reverse')
grid on; grid minor; 
xlabel('LWC ($g/m^3$)', 'Interpreter','latex');
ylabel('$\tau$', 'Interpreter','latex');


subplot(1,2,2)
set(gca, 'ydir', 'reverse')
grid on; grid minor; 
xlabel('$r_{e}$ ($\mu m$)', 'Interpreter','latex')

% set plot size
set(gcf, 'Position', [0 0 1000 500])





end