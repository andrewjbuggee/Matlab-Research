%% Plot the medain profile for LWC, re, and Nc

% only use the profiles defined by the logical index provided

% the changing variable is a string that defines the different indexed
% ensembles

% legend_str defines each indexed ensemble and is plotted in a legend

function plot_multiple_median_ensemble_LWC_re_NC_vs_norm_optical_depth(ensemble_profiles, index, changing_variable, legend_str)


% Plot multiple indexed ensembles on top of one another

figure;

for ii = 1:size(index, 1)

    % Grab the porfiles you wish to plot



    if sum(index(ii,:))<length(ensemble_profiles.re)

        profiles2plot.altitude = ensemble_profiles.altitude(index(ii,:));
        profiles2plot.tau = ensemble_profiles.tau(index(ii, :));
        profiles2plot.re = ensemble_profiles.re(index(ii, :));
        profiles2plot.lwc = ensemble_profiles.lwc(index(ii, :));
        profiles2plot.Nc = ensemble_profiles.Nc(index(ii, :));
        profiles2plot.time = ensemble_profiles.time(index(ii, :));

    end





    %%  Segment re, LWC, Nc into N bins along optical depth

    % In order to compute a median vertical profile, we have to first normalize
    % the vertical extent so that all profiles lie between values [0,1]. Then
    % we break up the vertical component in n discrete bins. Within each bin we
    % can compute the mean, median and standard deviation

    n_bins = 30; % number of segments the noramlized vertical component is broken up into

    bin_edges = 0:1/n_bins:1;

    % set up an empty cell array for all the values of each variable of interest
    % within each segment boundaries. Let's do this for droplet size, total
    % number concentration and liquid water content
    vertically_segmented_attributes = cell(n_bins, 3);


    normalized_tau = cell(1, length(profiles2plot.altitude));


    for nn = 1:length(profiles2plot.altitude)

        % first we need to normalize the vertical component of all profiles
        normalized_tau{nn} = profiles2plot.tau{nn}./profiles2plot.tau{nn}(end);

        % the data is stored in altitude space. If we wish to have the data
        % oriented in optical depth space, we need to reverse the order,
        % because optical depth is meausred from cloud top to cloud bottom
        % so we need to check if the plane is ascending or descending so we
        % know whether the data starts at cloud top or cloud bottom
        % We want all data oriented in tau space. So look for ascending
        % profiles, and flip the vector. If the plane is descending, we don't
        % need to do anything
        if mean(diff(profiles2plot.altitude{nn})./diff(profiles2plot.time{nn}))>0
            % if true, then the plane is ascending, grab and flip the variables
            % of interest
            re = flipud(profiles2plot.re{nn});
            lwc = flipud(profiles2plot.lwc{nn});
            Nc = flipud(profiles2plot.Nc{nn});

        else
            % if false, then the plane is descending, just grab the variables
            % of interest
            re = profiles2plot.re{nn};
            lwc = profiles2plot.lwc{nn};
            Nc = profiles2plot.Nc{nn};

        end





        % for each profile, we need to segment the variables of interest into n
        % bins.

        for bb = 1:length(bin_edges)-1

            % grab all re, LWC, and Nc values within each bin. Segment them
            % accordingly
            if bb==1
                index_segment = normalized_tau{nn}>=bin_edges(bb) & normalized_tau{nn}<=bin_edges(bb+1);

            else
                index_segment = normalized_tau{nn}>bin_edges(bb) & normalized_tau{nn}<=bin_edges(bb+1);
            end

            % store the effective radius values
            vertically_segmented_attributes{bb, 1} = [vertically_segmented_attributes{bb, 1}; re(index_segment)];

            % store the liquid water content values
            vertically_segmented_attributes{bb, 2} = [vertically_segmented_attributes{bb, 2}; lwc(index_segment)];

            % store the total droplet number concentration values
            vertically_segmented_attributes{bb, 3} = [vertically_segmented_attributes{bb, 3}; Nc(index_segment)];



        end



    end


    %%  Compute the Median LWC, re, and Nc of each layer

    re_median = zeros(n_bins, 1);
    re_std = zeros(n_bins, 1);
    re_avg_deviation_from_median = zeros(n_bins, 1);
    re_min = zeros(n_bins, 1);

    lwc_median = zeros(n_bins, 1);
    lwc_std = zeros(n_bins, 1);

    Nc_median = zeros(n_bins, 1);
    Nc_std = zeros(n_bins, 1);

    bin_center = zeros(n_bins, 1);



    for bb = 1:n_bins

        % compute the median value for the current bin
        re_median(bb) = median(vertically_segmented_attributes{bb,1});       % microns - mean droplet radius value

        % compute the standard deviation of the current bin
        re_std(bb) = std(vertically_segmented_attributes{bb,1});         % microns - standard deviation

        % The std is the average deviation from the mean value. Compute the
        % average deviation from the median value
        re_avg_deviation_from_median(bb) = sqrt(mean((vertically_segmented_attributes{bb,1} - re_median(bb)).^2));

        % find the smallest droplet size value at each bin
        re_min(bb) = min(vertically_segmented_attributes{bb,1});

        % compute the median value for the current bin
        lwc_median(bb) = median(vertically_segmented_attributes{bb,2});       % microns - mean droplet radius value

        % compute the standard deviation of the current bin
        lwc_std(bb) = std(vertically_segmented_attributes{bb,2});         % microns - standard deviation

        % compute the median value for the current bin
        Nc_median(bb) = median(vertically_segmented_attributes{bb,3});       % microns - mean droplet radius value

        % compute the standard deviation of the current bin
        Nc_std(bb) = std(vertically_segmented_attributes{bb,3});         % microns - standard deviation


        % compute the bin center which is the tau location of the mean data
        bin_center(bb) = (bin_edges(bb+1) - bin_edges(bb))/2 + bin_edges(bb);

    end


    %% Make a subplot of all 3 median profiles



    % plot the median effective radius
    subplot(1,3,1)

    hold on

    % plot the median droplet profile
    plot(re_median, bin_center, 'Color', mySavedColors(ii, 'fixed'))

    hold on

    if ii == size(index,1)
        
        set(gca, 'YDir', 'reverse')
    
        grid on; grid minor
        xlabel('$<r_e(z)>$ $(\mu m)$', 'Interpreter','latex')
        ylabel('Normalized Optical Depth', 'Interpreter', 'latex')

    end



    % plot the median liquid water content profile
    subplot(1,3,2)

    hold on

    % plot the median droplet profile
    plot(lwc_median, bin_center, 'Color', mySavedColors(ii, 'fixed'))
    
    if ii == size(index,1)
        set(gca, 'YDir', 'reverse')
    
        grid on; grid minor
        xlabel('$<LWC(z)>$ $(g/m^{3})$', 'Interpreter','latex')
        ylabel('Normalized Optical Depth', 'Interpreter','latex')

            % set the figure title
        title(['Median Profiles of different ', changing_variable, ' regimes'],...
            'Interpreter','latex')

    end





    % plot the median droplet number concentration
    subplot(1,3,3)

    hold on

    % plot the median droplet profile
    plot(Nc_median, bin_center, 'Color', mySavedColors(ii, 'fixed'))

    if ii == size(index, 1)
       
        set(gca, 'YDir', 'reverse')
    
        grid on; grid minor
        xlabel('$<N_c(z)>$ $(cm^{-3})$', 'Interpreter','latex')
        ylabel('Normalized Optical Depth', 'Interpreter', 'latex')

        legend(legend_str, 'Interpreter','latex', 'Location', 'best')

    end



    % Clear profiles to plot and prepare for the next index
    clear profiles2plot


end

% set the size of the figure
set(gcf, 'Position', [0 0 1200 625])

% Create textbox
annotation('textbox',[0.87 0.4352 0.07567 0.08],'String', ...
    {['$LWC \geq$', num2str(ensemble_profiles.inputs.LWC_threshold), ' $g/m^{2}$', ...
    newline,'   $N_c \geq$', num2str(ensemble_profiles.inputs.Nc_threshold), ' $cm^{-3}$']},...
    'LineWidth',2,...
    'FitBoxToText','on',...
    'FontSize',17,...
    'Interpreter','latex');








end

