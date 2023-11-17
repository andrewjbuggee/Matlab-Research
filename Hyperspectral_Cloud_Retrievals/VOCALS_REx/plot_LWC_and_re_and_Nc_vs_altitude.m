%% Plot the Liquid Water Content and the Effective Radius versus Altitude

% only plot the vertical profiles corresponding to the input: indices

% noramlize altitude is a true/false input. If true, the y-axis will be
% normalized so that all profiles have a vertical depth between 0 and 1


% By Andrew John Buggee
%%

function [] = plot_LWC_and_re_and_Nc_vs_altitude(vert_profiles, indices, normalize_altitude)

% Check to make sure there are 3 inputs


if nargin~=3
    error([newline,'Wrong number of inputs. Need 3: vertical profiles, indices to plot, normalization flag', newline])
end



% Plot the number of curves within the vert_profiles structure

N_cuvres = length(indices);

legend_str = cell(1, N_cuvres);

figure;
for nn = 1:N_cuvres

    % if normalize altitude is true, all altitude vectors will be
    % normalized between 0 and 1

    if normalize_altitude==true

        norm_alt = (vert_profiles.altitude{indices(nn)} - min(vert_profiles.altitude{indices(nn)}))./...
            (max(vert_profiles.altitude{indices(nn)}) - min(vert_profiles.altitude{indices(nn)}));

        % First plot the LWC
        ax1 = subplot(1,3,1); plot(vert_profiles.lwc{indices(nn)}, norm_alt);
        hold on

        % next plot the effective radius
        % if the 2DC data is compliant, plot the effective radius computed
        % using both instruments
        if vert_profiles.flag_2DC_data_is_conforming==true
            ax2 = subplot(1,3,2); plot(vert_profiles.re{indices(nn)}, norm_alt);
        else
            % if the 2DC data is non-conforming, use only the CDP data and
            % make a note of it
            ax2 = subplot(1,3,2); plot(vert_profiles.re_CDP{indices(nn)}, norm_alt);
        end
        hold on

        % Lastly, plot the total droplet number concentration
        ax3 = subplot(1,3,3); plot(vert_profiles.Nc{indices(nn)}, norm_alt);
        hold on

    else

        % First plot the LWC
        ax1 = subplot(1,3,1); plot(vert_profiles.lwc{indices(nn)}, vert_profiles.altitude{indices(nn)});
        hold on

        % next plot the effective radius
        % if the 2DC data is compliant, plot the effective radius computed
        % using both instruments
        if vert_profiles.flag_2DC_data_is_conforming==true
            ax2 = subplot(1,3,2); plot(vert_profiles.re{indices(nn)}, vert_profiles.altitude{indices(nn)});
        else
            % if the 2DC data is non-conforming, use only the CDP data and
            % make a note of it
            ax2 = subplot(1,3,2); plot(vert_profiles.re_CDP{indices(nn)}, vert_profiles.altitude{indices(nn)});
        end
        hold on

        % Lastly, plot the total droplet number concentration
        ax3 = subplot(1,3,3); plot(vert_profiles.Nc{indices(nn)}, vert_profiles.altitude{indices(nn)});
        hold on

    end

    legend_str{nn} = ['idx = ', num2str(indices(nn))];


end

% Make each subplot pretty
subplot(1,3,1)
grid on; grid minor;
xlabel('LWC ($g/m^3$)', 'Interpreter','latex');
ylabel('Altitude ($m$)', 'Interpreter','latex');



subplot(1,3,2)
grid on; grid minor;
% if the 2DC data is compliant, plot the effective radius computed
        % using both instruments
if vert_profiles.flag_2DC_data_is_conforming==true
    xlabel('$r_e$ ($\mu m$)', 'Interpreter','latex')
else
    % if the 2DC data is non-conforming, use only the CDP data and
    % make a note of it
    xlabel('$r_e$ ($\mu m$) - (CDP only)', 'Interpreter','latex')
end

% include a title in the middle plot
if isfield(vert_profiles, 'LWC_threshold')==true
    title(['$LWC \geq$ ', num2str(vert_profiles.LWC_threshold),' $g/m^{3}$',...
        '   $N_c \geq$ ', num2str(vert_profiles.Nc_threshold), ' $cm^{-3}$'], 'interpreter', 'latex')

elseif isfield(vert_profiles.inputs, 'LWC_threshold')==true
    title(['$LWC \geq$ ', num2str(vert_profiles.inputs.LWC_threshold),' $g/m^{3}$',...
        '   $N_c \geq$ ', num2str(vert_profiles.inputs.Nc_threshold), ' $cm^{-3}$'], 'interpreter', 'latex')

end




subplot(1,3,3)
grid on; grid minor;
xlabel('$N_c$ ($cm^{-3}$)', 'Interpreter','latex')

% in the third subplot, define the indices being plotted
legend(legend_str, 'Interpreter','latex', 'Location','best')

% set plot size
set(gcf, 'Position', [0 0 1200 625])

% link the yaxes so that they all have the same bounds
linkaxes([ax1 ax2 ax3],'y')




end