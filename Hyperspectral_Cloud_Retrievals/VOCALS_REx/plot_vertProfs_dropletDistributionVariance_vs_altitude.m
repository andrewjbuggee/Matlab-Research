%% Plot the droplet number distribution as an image showing dependence on altitude and droplet size

% by Andrew John Buggee
%%

function [] = plot_vertProfs_dropletDistributionVariance_vs_altitude(vocalsRex, vert_profs, index2plot, dist_name)

% how many vertical profiles do you wish to plot?

n_profiles = length(index2plot);


% What is the minimum value you wish to keep in the number concentration
% data
min_Nc = 10^(-2);




% let's make a subplot of all the profiles.
% if the total number of profiles is less that 11, we will use two rows.
% if the total number is greater than or equal to 11, we will use 3 rows

if n_profiles<=3

    % only create 1 row of subplots

    % imagesc assumes linearly spaced vectors along the x and y axes. Use p
    % color for non-linear sampling


    figure;

    for nn = 1:n_profiles
        

      
        % define the subplot
        subplot(1,n_profiles, nn)

        for zz = 1:length(vert_profs.time{index2plot(nn)})

            % what if we truncate data below 10^{-2}
            data = vocalsRex.Nc(:, vert_profs.time{index2plot(nn)}(zz));
            % truncate the data
            data = data(data>=min_Nc);



            % fit distribution type to the data provided
            pd = fitdist(data, dist_name);

        end

        % plot the altitude on the y-axis
        Y_data = vocalsRex.altitude(vert_profs.time{index2plot(nn)});


        [X,Y] = meshgrid(X_data, Y_data);

        s = pcolor(X,Y,vocalsRex.Nc(:, vert_profs.time{index2plot(nn)})');

        % reduce the transparancy of the bin edges
        s.EdgeAlpha = 0;

        xlabel('$r$ $(\mu m)$', 'Interpreter','latex')
        % only define the ylabel for the first (leftmost) subplot
        if nn==1
            ylabel('Altitude ($m$)', 'Interpreter','latex')
        end

        c = colorbar;

        xlim([r_min r_max])
        set(gca, 'YDir', 'normal')
        

        % only define the colorbar label for the last (right-most) subplot
        if nn==n_profiles
            ylabel(c,'$n(r)$ $(m^{-3}$)','FontSize',25, 'interpreter', 'latex');
        end

        
        % set the colormap to vary logarithmically
        if log_color==true
            set(gca,'ColorScale','log')
        end


        set(gca,'CLim', [10^(-1) , 10^2])

    end



elseif n_profiles>3 && n_profiles<=10

    % create a subplot with two rows

    % determine if the number of profiles is even or not
    iseven = rem(n_profiles,2)==0;

    if iseven==true
        n_columns = n_profiles/2;

    else
        n_columns = (n_profiles+1)/2;
    end



    % imagesc assumes linearly spaced vectors along the x and y axes. Use p
    % color for non-linear sampling

    % plot the droplet sizes on the x-axis
    X_data = vocalsRex.drop_radius_bin_center;

    figure;

    for nn = 1:n_profiles

        % define the subplot
        subplot(2,n_columns, nn)

        % plot the altitude on the y-axis
        Y_data = vocalsRex.altitude(vert_profs.time{index2plot(nn)});


        [X,Y] = meshgrid(X_data, Y_data);

        s = pcolor(X,Y,vocalsRex.Nc(:, vert_profs.time{index2plot(nn)})');

        % reduce the transparancy of the bin edges
        s.EdgeAlpha = 0;

        xlabel('$r$ $(\mu m)$', 'Interpreter','latex')
        % only define the ylabel for the first (leftmost) subplot
        if nn==1
            ylabel('Altitude ($m$)', 'Interpreter','latex')
        end

        c = colorbar;

        xlim([r_min r_max])
        set(gca, 'YDir', 'normal')
        

        % only define the colorbar label for the last (right-most) subplot
        if nn==n_profiles
            ylabel(c,'$n(r)$ $(m^{-3}$)','FontSize',25, 'interpreter', 'latex');
        end

        
        % set the colormap to vary logarithmically
        if log_color==true
            set(gca,'ColorScale','log')
        end


        set(gca,'CLim', [10^(-1) , 10^2])

    end




elseif n_profiles>10

    % create a subplot with 3 rows


end



% Define teh desired figure properties

% set the figure size
set(gcf, 'Position', [0 0 1200 600])









end