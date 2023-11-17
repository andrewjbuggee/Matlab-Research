%% Plot the Liquid Water Content and the Effective Radius versus Altitude


%   (3) noramlize - a logic entry (true of false) that tells the code, if
%   true, to normalize all optical depths to be between 0 and 1

% By Andrew John Buggee
%%

function [] = plot_LWC_and_re_and_Nc_vs_opticalDepth(vert_profiles, indices, normalize_opticalDepth)

% Check to make sure there are 3 inputs

if nargin~=3
    error([newline,'Wrong number of inputs. Need 3: vertical profiles, indices to plot, normalization flag', newline])
end

% define the denisty of liquid water
density = 1e6;          % g/m^3

% Plot the number of curves within the vert_profiles structure

N_curves = length(indices);
legend_str = cell(1, N_curves);

figure;

if N_curves<4

    for nn = 1:N_curves

        legend_str{nn} = ['idx = ', num2str(indices(nn))];
        % Profiles measured while the plane was descending will start with values
        % at the cloud top

        if (vert_profiles.altitude{indices(nn)}(end)-vert_profiles.altitude{indices(nn)}(1))>0
            % This profile is ascending, meaning the first data points are at
            % the cloud bottom, when tau is largest, since tau is defined from
            % top to bottom

            % all stored tau vectors start from 0. So if the plane is
            % ascending, we need to flip the tau vector
            tau = flipud(reshape(vert_profiles.tau{indices(nn)}, [], 1));



        elseif (vert_profiles.altitude{indices(nn)}(end)-vert_profiles.altitude{indices(nn)}(1))<0
            % This profile is descending, meaning the first data points are
            % measured at cloud top

            % in this case we don't flip the tau vector
            tau = reshape(vert_profiles.tau{indices(nn)}, [], 1);

        end


        % Because the data is oriented by values at the cloud bottom first
        % the tau vector goes from the largest values to 0



        % if normalize optical depth is true, all altitude vectors will be
        % normalized between 0 and 1

        if normalize_opticalDepth==true

            norm_tau = tau./max(tau);

            ax1 = subplot(1,3,1); plot(vert_profiles.lwc{indices(nn)}, norm_tau)
            hold on

            % next plot the effective radius
            % if the 2DC data is compliant, plot the effective radius computed
            % using both instruments
            if vert_profiles.flag_2DC_data_is_conforming==true
                ax2 = subplot(1,3,2); plot(vert_profiles.re{indices(nn)}, norm_tau);
            else
                % if the 2DC data is non-conforming, use only the CDP data and
                % make a note of it
                ax2 = subplot(1,3,2); plot(vert_profiles.re_CDP{indices(nn)}, norm_tau);
            end
            hold on

            % Finally, plot the total droplet number concentration
            ax3 = subplot(1,3,3); plot(vert_profiles.Nc{indices(nn)}, norm_tau);
            hold on


        else



            ax1 = subplot(1,3,1); plot(vert_profiles.lwc{indices(nn)}, tau)
            hold on

            % next plot the effective radius
            % if the 2DC data is compliant, plot the effective radius computed
            % using both instruments
            if vert_profiles.flag_2DC_data_is_conforming==true
                ax2 = subplot(1,3,2); plot(vert_profiles.re{indices(nn)}, tau);
            else
                % if the 2DC data is non-conforming, use only the CDP data and
                % make a note of it
                ax2 = subplot(1,3,2); plot(vert_profiles.re_CDP{indices(nn)}, tau);
            end
            hold on

            % Finally, plot the total droplet number concentration
            ax3 = subplot(1,3,3); plot(vert_profiles.Nc{indices(nn)}, tau);
            hold on

        end



    end


else
    % if there is more than 3 curves, let's make them transparent
    % Make all of the curves semi transparent
    for nn = 1:N_curves

        tau = zeros(1, length(vert_profiles.altitude{indices(nn)}));

        legend_str{nn} = ['idx = ', num2str(indices(nn))];

        % Profiles measured while the plane was descending will start with values
        % at the cloud top

        if (vert_profiles.altitude{indices(nn)}(end)-vert_profiles.altitude{indices(nn)}(1))>0
            % This profile is ascending, meaning the first data points are at
            % the cloud bottom, when tau is largest, since tau is defined from
            % top to bottom
            % all stored tau vectors start from 0. So if the plane is
            % ascending, we need to flip the tau vector
            tau = flipud(reshape(vert_profiles.tau{nn}, [], 1));

        elseif (vert_profiles.altitude{indices(nn)}(end)-vert_profiles.altitude{indices(nn)}(1))<0
            % This profile is descending, meaning the first data points are
            % measured at cloud top

            % in this case we don't flip the tau vector
            tau = reshape(vert_profiles.tau{nn}, [], 1);
        end


        % Because the data is oriented by values at the cloud bottom first
        % the tau vector goes from the largest values to 0

        % if normalize optical depth is true, all altitude vectors will be
        % normalized between 0 and 1

        if normalize_opticalDepth==true

            norm_tau = tau./tau(1);

            ax1 = subplot(1,3,1); l = plot(vert_profiles.lwc{indices(nn)}, norm_tau);
            % Set the transparency to 50%
            %l.Color(4) = 0.5;
            hold on

            % next plot the effective radius
            % if the 2DC data is compliant, plot the effective radius computed
            % using both instruments
            if vert_profiles.flag_2DC_data_is_conforming==true
                ax2 = subplot(1,3,2); plot(vert_profiles.re{indices(nn)}, norm_tau);
            else
                % if the 2DC data is non-conforming, use only the CDP data and
                % make a note of it
                ax2 = subplot(1,3,2); plot(vert_profiles.re_CDP{indices(nn)}, norm_tau);
            end
            % Set the transparency to 50%
            %l.Color(4) = 0.5;
            hold on

            % Finally, plot the total droplet number concentration
            ax3 = subplot(1,3,3); l = plot(vert_profiles.Nc{indices(nn)}, norm_tau);
            % Set the transparency to 50%
            %l.Color(4) = 0.5;
            hold on



        else


            ax1 = subplot(1,3,1); l = plot(vert_profiles.lwc{indices(nn)}, tau);
            % Set the transparency to 50%
            %l.Color(4) = 0.5;
            hold on

            % next plot the effective radius
            % if the 2DC data is compliant, plot the effective radius computed
            % using both instruments
            if vert_profiles.flag_2DC_data_is_conforming==true
                ax2 = subplot(1,3,2); plot(vert_profiles.re{indices(nn)}, tau);
            else
                % if the 2DC data is non-conforming, use only the CDP data and
                % make a note of it
                ax2 = subplot(1,3,2); plot(vert_profiles.re_CDP{indices(nn)}, tau);
            end
            % Set the transparency to 50%
            %l.Color(4) = 0.5;
            hold on

            % Finally, plot the total droplet number concentration
            ax3 = subplot(1,3,3); l = plot(vert_profiles.Nc{indices(nn)}, tau);
            % Set the transparency to 50%
            %l.Color(4) = 0.5;
            hold on

        end


    end

end


% Make each subplot pretty
subplot(1,3,1)
set(gca, 'ydir', 'reverse')
grid on; grid minor;
xlabel('LWC ($g/m^3$)', 'Interpreter','latex');
ylabel('$\tau$', 'Interpreter','latex');


subplot(1,3,2)
set(gca, 'ydir', 'reverse')
grid on; grid minor;
xlabel('$r_e$ ($\mu m$)', 'Interpreter','latex')
% include a title in the middle plot
if isfield(vert_profiles, 'LWC_threshold')==true
    title(['LWC $\geq$ ', num2str(vert_profiles.LWC_threshold),' $g/m^{3}$'], 'interpreter', 'latex')

elseif isfield(vert_profiles.inputs, 'LWC_threshold')==true
    title(['LWC $\geq$ ', num2str(vert_profiles.inputs.LWC_threshold),' $g/m^{3}$'], 'interpreter', 'latex')

end

% define a legend
legend(legend_str, 'Interpreter','latex', 'FontSize', 18, 'Location','southeast')


subplot(1,3,3)
set(gca, 'ydir', 'reverse')
grid on; grid minor;
xlabel('$N_c$ ($cm^{-3}$)', 'Interpreter','latex')

% set plot size
set(gcf, 'Position', [0 0 1200 650])

% Link the vertical axes together so that when you zoom, all subplots have
% the same limits
linkaxes([ax1, ax2, ax3], 'y')





end