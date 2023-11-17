%% Plot the droplet distribution as a histogram for multiple any number of heights within a cloud

% INPUTS:
% -------
%       (1) vert_profs - vertical profiles cell array found from the
%       find_verticalProfiles_VOCALS_REx function

%       (2) indexes2plot - the indexes associated with the vertical
%       profiles you wish to plot

%       (3) radius_limits - the range of radii you wish to display on this
%       histrogram. This will trim the x-axis. Input as [r_min, r_max]


% By Andrew John Buggee
%%

function [] = plot_dropDistribution_topAndBottom_hist(vert_profs, indexes2plot, radius_limits)



% Define the min and max radius values to plot
% r_min = 0;      % microns
% r_max = 30;     % microns

r_min = radius_limits(1);
r_max = radius_limits(2);


for nn = 1:length(indexes2plot)


    f1 = figure;

    % define the index for cloud top and cloud bottom
    times2plot = [1, numel(vert_profs.Nc{indexes2plot(nn)})];
    
    % set the effective radius vector
    re = zeros(1, 2);

    for tt = 1:length(times2plot)


        % Compute the effective radius for the two distributions and plot it as a solid vertical line
        re(tt) = vert_profs.re{indexes2plot(nn)}(times2plot(tt));

        % Plot the distribution at cloud bottom first
        h1 = histogram('BinEdges',vert_profs.drop_radius_bin_edges ,'BinCounts',...
            vert_profs.nr{indexes2plot(nn)}(:,(times2plot(tt))));
        h1.FaceColor = mySavedColors(tt, 'fixed');
        h1.FaceAlpha = 0.7;
        h1.EdgeAlpha = 1;
        hold on
        xline(re(tt),':', 'LineWidth',4, 'Color',mySavedColors(tt, 'fixed'))

    end
    

    % what profile are we plotting?
    title(['Vertical Profile ', num2str(indexes2plot(nn))], 'Interpreter','latex')

    % set axes limits and labels
    xlabel('Droplet Radius ($\mu m$)', 'Interpreter','latex', 'FontSize',32);
    ylabel('$n(r)$ ($cm^{-3} \, \mu m^{-1}$)', 'Interpreter','latex', 'FontSize',32);
    grid on; grid minor; hold on;
    xlim([r_min, r_max])
    ylim([10^(-2) 10^2])
    set(gca, 'YScale', 'log')
    set(gcf, 'Position',[0 0 1000, 600])
    
    
    legend('Cloud Bottom', ['$r_e$ = ',num2str(round(re(1),2)), ' $\mu m$'], 'Cloud Top',...
            ['$r_e$ = ',num2str(round(re(2),2)), ' $\mu m$'], 'Location','best', 'Interpreter','latex')



end











end