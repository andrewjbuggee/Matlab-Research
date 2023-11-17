%% Plot distribution of effective radii for the cloud top and bottom

% -------------- INPUTS -----------------

% (1) ensemble_profiles - this is an output from the script
% ensemble_statistics_cloud_top_and_bottom. It is a collection of vertical
% profiles found over all VOCALS-REx data based on a LWC and Nc threshold

% (2) averaging_length - a number that tells the code to average the cloud
% bottom radius over the bottom averaging_length measurements. The same
% goes for cloud top. If this value is greater than 1, it will be used to
% take an average.

% (3) normalization - the normalization inputs accepted by MATLAB's
% histogram function
%   (a) 'count' - (default) Count or frequency of observations.
%   (b) 'probability' - Relative probability.
%   (c) 'percentage' - Relative percentage.
%   (d) 'countdensity' - Count or frequency scaled by width of bin.
%   (e) 'cumcount' - Cumulative count, or the number of observations in each bin and all previous bins.
%   (f) 'pdf' - Probability density function estimate.
%   (g) 'cdf' - Cumulative distribution function estimate.



% By Andrew John Buggee

%%


function plot_ensemble_re_distribution_for_cloud_top_and_bottom(ensemble_profiles, averaging_length, normalization)




% Create a histogram of the ensemble of effective radii at some cloud level
r_top = zeros(1, length(ensemble_profiles.re));
r_bot = zeros(1, length(ensemble_profiles.re));


% grab all effective radii at cloud bottom and cloud top
for nn = 1:length(ensemble_profiles.re)

    % first we need to check if the plane is ascending or descending. That
    % tells us whether the measurement was made at cloud top or cloud
    % bottom

    dz_dt = diff(reshape(ensemble_profiles.altitude{nn}, [],1))./diff(reshape(ensemble_profiles.time{nn}, [], 1));

    if mean(dz_dt)>0
        % then the plane is ascending!

        if averaging_length==1

            % The first value is taken at cloud bottom
            r_bot(nn) = ensemble_profiles.re{nn}(1);            % microns

            % The last value is take at cloud top
            r_top(nn) = ensemble_profiles.re{nn}(end);          % microns

        elseif averaging_length>1
            % Then we average over several values
            % The first value is taken at cloud bottom
            r_bot(nn) = mean(ensemble_profiles.re{nn}(1:averaging_length));            % microns

            % The last value is take at cloud top
            r_top(nn) = mean(ensemble_profiles.re{nn}(end-averaging_length+1:end));          % microns

        else

            error([newline, 'The averaging length must be an integer greater than or equal to 1', newline])



        end



    elseif mean(dz_dt<0)
        % Then the plane is descending!

        if averaging_length==1

            % The first value is taken at cloud top
            r_top(nn) = ensemble_profiles.re{nn}(1);            % microns

            % The last value is take at cloud bottom
            r_bot(nn) = ensemble_profiles.re{nn}(end);          % microns

        elseif averaging_length>1
            % Then we average over several values
            % The first value is taken at cloud top
            r_top(nn) = mean(ensemble_profiles.re{nn}(1:averaging_length));            % microns

            % The last value is take at cloud bottom
            r_bot(nn) = mean(ensemble_profiles.re{nn}(end-averaging_length+1:end));          % microns

        end

    else

        error([newline, 'The averaging length must be an integer greater than or equal to 1', newline])


    end


end


% ----- Create Histogram -----
n_bins = 20;

figure;

% ----- Plot r_top first -----
subplot(1,2,1)

if strcmp(normalization, 'count')==true

    histogram(r_top, n_bins, 'Normalization',normalization)
    ylabel('Counts', 'Interpreter','latex')

elseif strcmp(normalization, 'pdf')==true

    histogram(r_top, 'Normalization',normalization)
    ylabel('PDF', 'Interpreter','latex')

end

xlabel('$r_{e}(\tau = 0)$', 'Interpreter','latex')


grid on; grid minor

% Plot the median value as a vertical line
hold on
xline(median(r_top), 'LineWidth',2, 'Color','k', 'Label','median',...
    'LineStyle','--', 'FontSize',15)



% ----- Plot r_bot next -----
subplot(1,2,2)

if strcmp(normalization, 'count')==true

    histogram(r_bot, n_bins, 'Normalization',normalization)
    ylabel('Counts', 'Interpreter','latex')

elseif strcmp(normalization, 'pdf')==true

    histogram(r_bot, 'Normalization',normalization)
    ylabel('PDF', 'Interpreter','latex')

end

xlabel('$r_{e}(\tau = \tau_c)$', 'Interpreter','latex')


grid on; grid minor

% Plot the median value as a vertical line
hold on
xline(median(r_bot), 'LineWidth',2, 'Color','k', 'Label','median',...
    'LineStyle','--', 'FontSize',15)


% set figure size
set(gcf, 'Position', [0 0 1200 650])







end