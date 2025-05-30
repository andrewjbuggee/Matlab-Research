%% Plot horizontal profiles of Liquid Water Content and the Effective Radius versus ground distance travelled

% only plot the horizontal profiles corresponding to the input: indices

% noramlize altitude is a true/false input. If true, the x-axis will be
% normalized so that all profiles have a horizontal distance travelled of between 0 and 1


% By Andrew John Buggee
%%

function [] = plot_horiztonal_profiles_LWC_and_re_and_Nc(horz_profiles, indices, normalize_distance, line_colors)

% Check to make sure there are 3 inputs


if nargin~=4
    error([newline,'Wrong number of inputs. Need 4: horizontal profiles, indices to plot, normalization flag and colors for each line', newline])
end



% Plot the number of curves within the horz_profiles structure

N_cuvres = length(indices);

legend_str = cell(1, N_cuvres);

figure;
for nn = 1:N_cuvres

    % if normalize distance is true, all distance vectors will be
    % normalized between 0 and 1

    if normalize_distance==true

        norm_dist = (horz_profiles.horz_dist{indices(nn)} - min(horz_profiles.horz_dist{indices(nn)}))./...
            (max(horz_profiles.horz_dist{indices(nn)}) - min(horz_profiles.horz_dist{indices(nn)}));

        % First plot the LWC
        if isempty(line_colors)==true
            ax1 = subplot(3,1,1); plot(norm_dist, horz_profiles.lwc{indices(nn)});
        else
            ax1 = subplot(3,1,1); plot(norm_dist, horz_profiles.lwc{indices(nn)}, 'Color', line_colors(nn,:));
        end

        hold on

        % next plot the effective radius
        % if the 2DC data is compliant, plot the effective radius computed
        % using both instruments
        if horz_profiles.flag_2DC_data_is_conforming==true
            if isempty(line_colors)==true
                ax2 = subplot(3,1,2); plot(norm_dist, horz_profiles.re{indices(nn)});
            else
                 ax2 = subplot(3,1,2); plot(norm_dist, horz_profiles.re{indices(nn)}, 'Color', line_colors(nn,:));
            end
        else
            % if the 2DC data is non-conforming, use only the CDP data and
            % make a note of it
            if isempty(line_colors)==true
                ax2 = subplot(3,1,2); plot(norm_dist, horz_profiles.re_CDP{indices(nn)});
            else
                ax2 = subplot(3,1,2); plot(norm_dist, horz_profiles.re_CDP{indices(nn)}, 'Color', line_colors(nn,:));
            end
        end
        hold on

        % Lastly, plot the total droplet number concentration
        if isempty(line_colors)==true
            ax3 = subplot(3,1,3); plot(norm_dist, horz_profiles.Nc{indices(nn)});
        else
            ax3 = subplot(3,1,3); plot(norm_dist, horz_profiles.Nc{indices(nn)}, 'Color', line_colors(nn,:));
        end


        hold on

    else

        % --- DATA IS IN METERS - PLOT IN KILOMETERS ---
        % First plot the LWC
        if isempty(line_colors)==true
            ax1 = subplot(3,1,1); plot(horz_profiles.horz_dist{indices(nn)}./1e3, horz_profiles.lwc{indices(nn)});
        else
             ax1 = subplot(3,1,1); plot(horz_profiles.horz_dist{indices(nn)}./1e3, horz_profiles.lwc{indices(nn)}, 'Color', line_colors(nn,:));
        end
        hold on

        % next plot the effective radius
        % if the 2DC data is compliant, plot the effective radius computed
        % using both instruments

        if horz_profiles.flag_2DC_data_is_conforming==true
            if isempty(line_colors)==true
                ax2 = subplot(3,1,2); plot(horz_profiles.horz_dist{indices(nn)}./1e3, horz_profiles.re{indices(nn)});
            else
                ax2 = subplot(3,1,2); plot(horz_profiles.horz_dist{indices(nn)}./1e3, horz_profiles.re{indices(nn)}, 'Color', line_colors(nn,:));
            end

        else
            % if the 2DC data is non-conforming, use only the CDP data and
            % make a note of it
            if isempty(line_colors)==true
                ax2 = subplot(3,1,2); plot(horz_profiles.horz_dist{indices(nn)}./1e3, horz_profiles.re_CDP{indices(nn)});
            else
                ax2 = subplot(3,1,2); plot(horz_profiles.horz_dist{indices(nn)}./1e3, horz_profiles.re_CDP{indices(nn)}, 'Color', line_colors(nn,:));
            end

        end
        hold on

        % Lastly, plot the total droplet number concentration
        if isempty(line_colors)==true
            ax3 = subplot(3,1,3); plot(horz_profiles.horz_dist{indices(nn)}./1e3, horz_profiles.Nc{indices(nn)});
        else
            ax3 = subplot(3,1,3); plot(horz_profiles.horz_dist{indices(nn)}./1e3, horz_profiles.Nc{indices(nn)}, 'Color', line_colors(nn,:));
        end

        hold on

    end

    legend_str{nn} = ['idx = ', num2str(indices(nn))];


end

% Make each subplot pretty
subplot(3,1,1)
grid on; grid minor;
ylabel('LWC ($g/m^3$)', 'Interpreter','latex');




subplot(3,1,2)
grid on; grid minor;
% if the 2DC data is compliant, plot the effective radius computed
        % using both instruments
if horz_profiles.flag_2DC_data_is_conforming==true
    ylabel('$r_e$ ($\mu m$)', 'Interpreter','latex')
else
    % if the 2DC data is non-conforming, use only the CDP data and
    % make a note of it
    ylabel('$r_e$ ($\mu m$) - (CDP only)', 'Interpreter','latex')
end

% include a title in the middle plot
if isfield(horz_profiles, 'LWC_threshold')==true
    title(['$LWC \geq$ ', num2str(horz_profiles.LWC_threshold),' $g/m^{3}$',...
        '     $N_c \geq$ ', num2str(horz_profiles.Nc_threshold), ' $cm^{-3}$',...
        '     Max vert displacement: ', num2str(horz_profiles.max_vert_displacement), ' $m$'], 'interpreter', 'latex')

elseif isfield(horz_profiles.inputs, 'LWC_threshold')==true
    title(['$LWC \geq$ ', num2str(horz_profiles.inputs.LWC_threshold),' $g/m^{3}$',...
        '     $N_c \geq$ ', num2str(horz_profiles.inputs.Nc_threshold), ' $cm^{-3}$',...
        '     Max vert displacement: ', num2str(horz_profiles.inputs.max_vert_displacement), ' $m$'], 'interpreter', 'latex')

end






subplot(3,1,3)
grid on; grid minor;
ylabel('$N_c$ ($cm^{-3}$)', 'Interpreter','latex')

% in the third subplot, define the indices being plotted
legend(legend_str, 'Interpreter','latex', 'Location','best')

% set plot size
set(gcf, 'Position', [0 0 1200 625])

% Include an x axis label on the middle plot
if normalize_distance==true

    xlabel('Normalized Horizontal Distance Travelled', 'Interpreter','latex');
else

    xlabel('Horizontal Distance Travelled ($km$)', 'Interpreter','latex');
end

% link the yaxes so that they all have the same bounds
linkaxes([ax1 ax2 ax3],'x')




end