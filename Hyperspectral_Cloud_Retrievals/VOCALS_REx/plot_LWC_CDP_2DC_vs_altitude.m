%% Plot the Liquid Water Content measured by the CDP Probe and the 2DC probe

% only plot the vertical profiles corresponding to the input: indices

% noramlize altitude is a true/false input. If true, the y-axis will be
% normalized so that all profiles have a vertical depth between 0 and 1


% By Andrew John Buggee
%%

function [] = plot_LWC_CDP_2DC_vs_altitude(vert_profiles, indices, normalize_altitude)

% Check to make sure there are 3 inputs


if nargin~=3
    error([newline,'Wrong number of inputs. Need 3: vertical profiles, indices to plot, normalization flag', newline])
end



% Plot the number of curves within the vert_profiles structure

N_plots = length(indices);


figure;
for nn = 1:N_plots

    % if normalize altitude is true, all altitude vectors will be
    % normalized between 0 and 1

    % if there is more than 1 index, we create a subplot
    if N_plots>1
        ax{nn} = subplot(1,N_plots,nn);

    end

    if normalize_altitude==true

        norm_alt = (vert_profiles.altitude{indices(nn)} - min(vert_profiles.altitude{indices(nn)}))./...
            (max(vert_profiles.altitude{indices(nn)}) - min(vert_profiles.altitude{indices(nn)}));

        
        % Plot the CDP LWC first
        plot(vert_profiles.lwc_CDP{indices(nn)}, norm_alt, 'Color', mySavedColors(1, 'fixed'));
        hold on
        % Plot the 2DC LWC next
        plot(vert_profiles.lwc_2DC{indices(nn)}, norm_alt, 'Color', mySavedColors(2, 'fixed'));



    else

        % Plot the CDP LWC first
        plot(vert_profiles.lwc_CDP{indices(nn)}, vert_profiles.altitude{indices(nn)},...
            'Color', mySavedColors(1, 'fixed'));
        hold on
        % Plot the 2DC LWC next
        plot(vert_profiles.lwc_2DC{indices(nn)}, vert_profiles.altitude{indices(nn)},...
            'Color', mySavedColors(2, 'fixed'));



    end

    % Provide a title that is the index
    title(['idx = ', num2str(indices(nn))], 'Interpreter','latex')

    % define plot labels
    if nn ==1

        % Define the x and y labels
        grid on; grid minor;
        xlabel('LWC ($g/m^3$)', 'Interpreter','latex');
        ylabel('Altitude ($m$)', 'Interpreter','latex');
        
        % Define the legend for the two different instruments
        legend({'CDP', '2DC'}, 'Interpreter','latex', 'Location','best')


    else

        grid on; grid minor

    end


end


% set plot size
set(gcf, 'Position', [0 0 1200 625])

if N_plots>1
    % link the yaxes so that they all have the same bounds
    linkaxes([ax{:}],'y')
end




end






