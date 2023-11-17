%% Plot the droplet number distribution as an image showing dependence on altitude and droplet size

% by Andrew John Buggee
%%

function [] = plot_vocalsRex_vertical_profiles_pcolor_vs_optical_depth(vocalsRex, vert_profs, index2plot, log_color)

% how many vertical profiles do you wish to plot?

n_profiles = length(index2plot);

% define the denisty of liquid water
density = 1e6;          % g/m^3


% define radius limits to plot
r_min = 0;          % microns
r_max = 30;         % microns


% let's make a subplot of all the profiles.
% if the total number of profiles is less that 11, we will use two rows.
% if the total number is greater than or equal to 11, we will use 3 rows

if n_profiles<=3

    % only create 1 row of subplots

    % imagesc assumes linearly spaced vectors along the x and y axes. Use p
    % color for non-linear sampling

    % plot the droplet sizes on the x-axis
    X_data = vocalsRex.drop_radius_bin_center;

    figure;

    for nn = 1:n_profiles

        tau = zeros(1, length(vert_profs.altitude{index2plot(nn)}));

        % Profiles measured while the plane was descending will start with values
        % at the cloud top

        if (vert_profs.altitude{index2plot(nn)}(end)-vert_profs.altitude{index2plot(nn)}(1))>0
            % This profile is ascending, meaning the first data points are at
            % the cloud bottom, when tau is largest, since tau is defined from
            % top to bottom
            starting_idx = length(tau)+1;

            % compute optical depth
            % Let's assume the extinction coefficient is 2, meaning the radius is
            % much larger than the incident light
            for ii = length(vert_profs.altitude{index2plot(nn)}):-1:2
                tau(starting_idx-ii) = 1/2 * 1/density * trapz(vert_profs.altitude{index2plot(nn)}(1:ii),...
                    vert_profs.lwc{index2plot(nn)}(1:ii)./(vert_profs.re{index2plot(nn)}(1:ii)*1e-6));
            end

        elseif (vert_profs.altitude{index2plot(nn)}(end)-vert_profs.altitude{index2plot(nn)}(1))<0
            % This profile is descending, meaning the first data points are
            % measured at cloud top

            % compute optical depth
            % Let's assume the extinction coefficient is 2, meaning the radius is
            % much larger than the incident light
            for ii = 2:length(vert_profs.altitude{index2plot(nn)})
                tau(ii) = -1/2 * 1/density * trapz(vert_profs.altitude{index2plot(nn)}(1:ii),...
                    vert_profs.lwc{index2plot(nn)}(1:ii)./(vert_profs.re{index2plot(nn)}(1:ii)*1e-6));
            end
        end


        % plot the optical depth on the y-axis
        Y_data = tau;


        % define the subplot
        subplot(1,n_profiles, nn)


        [X,Y] = meshgrid(X_data, Y_data);

        s = pcolor(X,Y,vocalsRex.Nc(:, vert_profs.time{index2plot(nn)})');

        % reduce the transparancy of the bin edges
        s.EdgeAlpha = 0;

        xlabel('$r$ $(\mu m)$', 'Interpreter','latex')
        % only define the ylabel for the first (leftmost) subplot
        if nn==1
            ylabel('$\tau$', 'Interpreter','latex')
        end

        % set the title as the UTC time at the middle of the vertical
        % profile
        title(['UTC Time: ',num2str(median(vert_profs.time_utc{index2plot(nn)}))], 'Interpreter','latex', ...
            'FontSize', 20)

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


        % Set y-direction to be in reverse
        set(gca, 'ydir', 'reverse')

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

        tau = zeros(1, length(vert_profs.altitude{index2plot(nn)}));

        % Profiles measured while the plane was descending will start with values
        % at the cloud top

        if (vert_profs.altitude{index2plot(nn)}(end)-vert_profs.altitude{index2plot(nn)}(1))>0
            % This profile is ascending, meaning the first data points are at
            % the cloud bottom, when tau is largest, since tau is defined from
            % top to bottom
            starting_idx = length(tau)+1;

            % compute optical depth
            % Let's assume the extinction coefficient is 2, meaning the radius is
            % much larger than the incident light
            for ii = length(vert_profs.altitude{index2plot(nn)}):-1:2
                tau(starting_idx-ii) = 1/2 * 1/density * trapz(vert_profs.altitude{index2plot(nn)}(1:ii),...
                    vert_profs.lwc{index2plot(nn)}(1:ii)./(vert_profs.re{index2plot(nn)}(1:ii)*1e-6));
            end

        elseif (vert_profs.altitude{index2plot(nn)}(end)-vert_profs.altitude{index2plot(nn)}(1))<0
            % This profile is descending, meaning the first data points are
            % measured at cloud top

            % compute optical depth
            % Let's assume the extinction coefficient is 2, meaning the radius is
            % much larger than the incident light
            for ii = 2:length(vert_profs.altitude{index2plot(nn)})
                tau(ii) = -1/2 * 1/density * trapz(vert_profs.altitude{index2plot(nn)}(1:ii),...
                    vert_profs.lwc{index2plot(nn)}(1:ii)./(vert_profs.re{index2plot(nn)}(1:ii)*1e-6));
            end
        end



        % plot the optical depth on the y-axis
        Y_data = tau;


        % define the subplot
        subplot(2,n_columns, nn)


        [X,Y] = meshgrid(X_data, Y_data);

        s = pcolor(X,Y,vocalsRex.Nc(:, vert_profs.time{index2plot(nn)})');

        % reduce the transparancy of the bin edges
        s.EdgeAlpha = 0;
       
        
        % only define the xlabel ylabel for the first (leftmost) subplot
        if nn==1
            ylabel('Optical Depth', 'Interpreter','latex')
            xlabel('$r$ $(\mu m)$', 'Interpreter','latex')
        end

        % set the title as the UTC time at the middle of the vertical
        % profile
        title(['UTC Time: ',num2str(median(vert_profs.time_utc{index2plot(nn)}))], 'Interpreter','latex', ...
            'FontSize', 20)

        c = colorbar;

        xlim([r_min r_max])

        % flip the yaxis
        set(gca, 'YDir', 'reverse')


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


    % determine if the number of profiles is divisible by 3
    is_divisible = rem(n_profiles,3)==0;

    if is_divisible==true
        n_columns = n_profiles/3;

    else

        x = n_profiles+1;

        while rem(x,3)~=0
            x = x+1;

        end

        n_columns = x/3;

    end



    % imagesc assumes linearly spaced vectors along the x and y axes. Use p
    % color for non-linear sampling

    % plot the droplet sizes on the x-axis
    X_data = vocalsRex.drop_radius_bin_center;

    figure;

    for nn = 1:n_profiles


        tau = zeros(1, length(vert_profs.altitude{index2plot(nn)}));

        % Profiles measured while the plane was descending will start with values
        % at the cloud top

        if (vert_profs.altitude{index2plot(nn)}(end)-vert_profs.altitude{index2plot(nn)}(1))>0
            % This profile is ascending, meaning the first data points are at
            % the cloud bottom, when tau is largest, since tau is defined from
            % top to bottom
            starting_idx = length(tau)+1;

            % compute optical depth
            % Let's assume the extinction coefficient is 2, meaning the radius is
            % much larger than the incident light
            for ii = length(vert_profs.altitude{index2plot(nn)}):-1:2
                tau(starting_idx-ii) = 1/2 * 1/density * trapz(vert_profs.altitude{index2plot(nn)}(1:ii),...
                    vert_profs.lwc{index2plot(nn)}(1:ii)./(vert_profs.re{index2plot(nn)}(1:ii)*1e-6));
            end

        elseif (vert_profs.altitude{index2plot(nn)}(end)-vert_profs.altitude{index2plot(nn)}(1))<0
            % This profile is descending, meaning the first data points are
            % measured at cloud top

            % compute optical depth
            % Let's assume the extinction coefficient is 2, meaning the radius is
            % much larger than the incident light
            for ii = 2:length(vert_profs.altitude{index2plot(nn)})
                tau(ii) = -1/2 * 1/density * trapz(vert_profs.altitude{index2plot(nn)}(1:ii),...
                    vert_profs.lwc{index2plot(nn)}(1:ii)./(vert_profs.re{index2plot(nn)}(1:ii)*1e-6));
            end
        end



        % plot the optical depth on the y-axis
        Y_data = tau;


        % define the subplot
        subplot(3,n_columns, nn)


        [X,Y] = meshgrid(X_data, Y_data);

        s = pcolor(X,Y,vocalsRex.Nc(:, vert_profs.time{index2plot(nn)})');

        % reduce the transparancy of the bin edges
        s.EdgeAlpha = 0;

        
        % only define the xlabel and ylabel for the first (leftmost) subplot
        if nn==1
            ylabel('Optical Depth', 'Interpreter','latex')
            xlabel('$r$ $(\mu m)$', 'Interpreter','latex')
        end

        c = colorbar;

        xlim([r_min r_max])
        set(gca, 'YDir', 'normal')


        % only define the colorbar label for the last (right-most) subplot
        if nn==n_profiles
            ylabel(c,'$n(r)$ $(m^{-3}$)','FontSize',25, 'interpreter', 'latex');
        end

        % set the title as the UTC time at the middle of the vertical
        % profile
        title(['UTC Time: ',num2str(median(vert_profs.time_utc{index2plot(nn)}))], 'Interpreter','latex',...
            'FontSize', 20)


        % set the colormap to vary logarithmically
        if log_color==true
            set(gca,'ColorScale','log')
        end


        set(gca,'CLim', [10^(-1) , 10^2])

    end



end



% Define the desired figure properties

% set the figure size
set(gcf, 'Position', [0 0 1200 600])









end