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

    % define the index for cloud top and cloud bottom
    time2plot = index_altitude(nn);

    % Compute the effective radius for the two distributions and plot it as a solid vertical line
    re = vert_profs.re{indexes2plot(nn)}(time2plot);

    % Plot the distribution at cloud bottom first
    h1 = histogram('BinEdges',vert_profs.drop_radius_bin_edges ,'BinCounts',...
        vert_profs.nr{indexes2plot(nn)}(:,(time2plot)));
    h1.FaceColor = mySavedColors(1, 'fixed');
    h1.FaceAlpha = 0.7;
    h1.EdgeAlpha = 1;
    hold on
    xline(re,':', 'LineWidth',4, 'Color',mySavedColors(1, 'fixed'))




    % what profile are we plotting?
    title(['Vertical Profile ', num2str(indexes2plot(nn)), ...
        '  Altitude ', num2str(round(vert_profs.altitude{indexes2plot(nn)}(index_altitude(nn)))),...
        ' m'], 'Interpreter','latex')

    % set axes limits and labels
    xlabel('Droplet Radius ($\mu m$)', 'Interpreter','latex', 'FontSize',32);
    ylabel('$n(r)$ ($cm^{-3} \, \mu m^{-1}$)', 'Interpreter','latex', 'FontSize',32);
    grid on; grid minor; hold on;
    xlim([r_min, r_max])
    %ylim([10^(-2) 10^2])
    set(gca, 'YScale', 'log')
    set(gcf, 'Position',[0 0 1000, 600])


    legend('', ['$r_e$ = ',num2str(round(re(1),2)), ' $\mu m$'],'Location','best', 'Interpreter','latex')



end











end