%% Plot the droplet distribution of some vertical profile as a histogram at a specific altitude

% INPUTS:
% -------
%       (1) vert_profs - vertical profiles cell array found from the
%       find_verticalProfiles_VOCALS_REx function

%       (2) indexes2plot - the indexes associated with the vertical
%       profiles you wish to plot

%       (3) index_altitude - this is the index within the vertical profile
%       defining the altitude to plot. ** This must be the same length as
%       indexes2plot. There should be a unique index_altitude for each
%       index2plot

%       (3) radius_limits - the range of radii you wish to display on this
%       histrogram. This will trim the x-axis. Input as [r_min, r_max] in
%       units of microns


% By Andrew John Buggee
%%

function [] = plot_dropDistribution_hist(vert_profs, indexes2plot, index_altitude, radius_limits)



% Define the min and max radius values to plot
% r_min = 0;      % microns
% r_max = 30;     % microns

r_min = radius_limits(1);
r_max = radius_limits(2);


for nn = 1:length(indexes2plot)


    f1 = figure;

    lgnd_str = cell(1, 2*length(index_altitude));

    for ll = 1:length(index_altitude)

        hold on

        % define the index for cloud top and cloud bottom
        time2plot = index_altitude(ll);

        % Compute the effective radius for the two distributions and plot it as a solid vertical line
        re = vert_profs(indexes2plot(nn)).re(time2plot);

        % Plot the distribution at cloud bottom first
        h = histogram('BinEdges',vert_profs(indexes2plot(nn)).drop_radius_bin_edges ,'BinCounts',...
            vert_profs(indexes2plot(nn)).Nc(:,(time2plot))');
        h.FaceColor = mySavedColors(61+(ll-1), 'fixed');
        h.FaceAlpha = 0.7;
        h.EdgeAlpha = 1;
        hold on
        xline(re,':', 'LineWidth',4, 'Color',mySavedColors(61+(ll-1), 'fixed'))

        lgnd_str{2*ll -1} = ['$r_e(z(', num2str(round(vert_profs(indexes2plot(nn)).altitude(time2plot))/1e3),...
            ')=$ ',num2str(round(re,1)), ' $\mu m$'];
        lgnd_str{2*ll} = '';


    end


    % what profile are we plotting?
    title(['Vertical Profile ', num2str(indexes2plot(nn))], 'Interpreter','latex')

    % set axes limits and labels
    xlabel('Droplet Radius ($\mu m$)', 'Interpreter','latex', 'FontSize',32);
    ylabel('$n(r)$ ($cm^{-3} \, \mu m^{-1}$)', 'Interpreter','latex', 'FontSize',32);
    grid on; grid minor; hold on;
    xlim([r_min, r_max])
    %ylim([10^(-2) 10^2])
    set(gca, 'YScale', 'log')
    set(gcf, 'Position',[0 0 1000, 600])


    legend(lgnd_str,'Location','best', 'Interpreter','latex')





end











end