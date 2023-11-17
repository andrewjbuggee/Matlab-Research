%% Plot Droplet Effective Radius measured by the CDP Probe and the 2DC probe
% FOR HORIZONTAL PROFILES ONLY

% only plot the vertical profiles corresponding to the input: indices

% noramlize altitude is a true/false input. If true, the x-axis will be
% normalized so that all profiles have a horizontal distance between 0 and 1


% By Andrew John Buggee
%%

function [] = plot_re_CDP_2DC_vs_horizontal_distance(horz_profiles, indices, normalize_dist)

% Check to make sure there are 3 inputs


if nargin~=3
    error([newline,'Wrong number of inputs. Need 3: horizontal profiles, indices to plot, normalization flag', newline])
end



% Plot the number of curves within the horz_profiles structure

N_plots = length(indices);


figure;
for nn = 1:N_plots

    % if normalize altitude is true, all altitude vectors will be
    % normalized between 0 and 1

    % if there is more than 1 index, we create a subplot
    if N_plots>1
        ax{nn} = subplot(1,N_plots,nn);

    end

    if normalize_dist==true

        norm_alt = (horz_profiles.horz_dist{indices(nn)} - min(horz_profiles.horz_dist{indices(nn)}))./...
            (max(horz_profiles.horz_dist{indices(nn)}) - min(horz_profiles.horz_dist{indices(nn)}));


        % Plot the CDP re first
        plot(norm_alt./1e3, horz_profiles.re_CDP{indices(nn)}, 'Color', mySavedColors(1, 'fixed'));
        hold on
        % Plot the 2DC re next
        % Sometimes VOCALS-REx doesn't publish the number concentration
        % data or an effective radius estimate for the 2DC data. In those
        % cases they usualy report a mean value, which appears to be the
        % first moment. Plot this if the 2DC effective radius isn't
        % available
        if isfield(horz_profiles, 're_2DC')==true
            plot(norm_alt./1e3, horz_profiles.re_2DC{indices(nn)}, 'Color', mySavedColors(2, 'fixed'));

        elseif isfield(horz_profiles, 'mean_r_2DC')==true
            plot(norm_alt./1e3, horz_profiles.mean_r_2DC{indices(nn)}, 'Color', mySavedColors(2, 'fixed'));

        end


    else

        % Plot the CDP re first
        plot(horz_profiles.horz_dist{indices(nn)}./1e3, horz_profiles.re_CDP{indices(nn)},...
            'Color', mySavedColors(1, 'fixed'));
        hold on
        % Plot the 2DC re next
        % Sometimes VOCALS-REx doesn't publish the number concentration
        % data or an effective radius estimate for the 2DC data. In those
        % cases they usualy report a mean value, which appears to be the
        % first moment. Plot this if the 2DC effective radius isn't
        % available
        if isfield(horz_profiles, 're_2DC')==true
            plot(horz_profiles.horz_dist{indices(nn)}./1e3, horz_profiles.re_2DC{indices(nn)}, 'Color', mySavedColors(2, 'fixed'));

        elseif isfield(horz_profiles, 'mean_r_2DC')==true
            plot(horz_profiles.horz_dist{indices(nn)}./1e3, horz_profiles.mean_r_2DC{indices(nn)}, 'Color', mySavedColors(2, 'fixed'));

        end



    end

    % Provide a title that is the index
    title(['idx = ', num2str(indices(nn))], 'Interpreter','latex')

    % define plot labels
    if nn ==1

        % Define the x and y labels
        grid on; grid minor;
        ylabel('$r_{e}$ ($\mu m$)', 'Interpreter','latex');
        xlabel('Horizontal Distance Travelled ($km$)', 'Interpreter','latex');

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
    linkaxes([ax{:}],'x')
end




end






