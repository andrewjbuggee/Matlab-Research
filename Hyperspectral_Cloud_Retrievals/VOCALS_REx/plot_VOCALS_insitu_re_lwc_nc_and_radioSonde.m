%% Plot a subplot with VOCALS-REx CDP measured data and radiosonde measured data


function [] = plot_VOCALS_insitu_re_lwc_nc_and_radioSonde(cloudDropletprofiles, radiosondeProfiles, indices, normalize_altitude)


% Check to make sure there are 4 inputs


if nargin~=4
    error([newline,'Wrong number of inputs. Need 4: vertical profiles, radioSonde_profiles, indices to plot', newline])
end



% Plot the number of curves within the cloudDropletprofiles structure

N_cuvres = length(indices);

legend_str = cell(1, N_cuvres);

figure;

clr_start = 61;

ttl_fnt = 20;
ax_fnt = 30;
lgnd_fnt = 18;






if iscell(cloudDropletprofiles)

    for nn = 1:N_cuvres

        % if normalize altitude is true, all altitude vectors will be
        % normalized between 0 and 1

        if normalize_altitude==true

            norm_alt = (cloudDropletprofiles{indices(nn)}.altitude - min(cloudDropletprofiles{indices(nn)}.altitude))./...
                (max(cloudDropletprofiles{indices(nn)}.altitude) - min(cloudDropletprofiles{indices(nn)}.altitude));

            % First plot the LWC
            ax1 = subplot(1,3,1); plot(cloudDropletprofiles{indices(nn)}.lwc, norm_alt, ...
                'Color', mySavedColors(clr_start + (nn-1), 'fixed'));
            hold on

            % next plot the effective radius
            % if the 2DC data is compliant, plot the effective radius computed
            % using both instruments
            if cloudDropletprofiles{indices(nn)}.flag_2DC_data_is_conforming==true
                ax2 = subplot(1,3,2); plot(cloudDropletprofiles{indices(nn)}.re, norm_alt, ...
                    'Color', mySavedColors(clr_start + (nn-1), 'fixed'));
            else
                % if the 2DC data is non-conforming, use only the CDP data and
                % make a note of it
                ax2 = subplot(1,3,2); plot(cloudDropletprofiles{indices(nn)}.re_CDP, norm_alt, ...
                    'Color', mySavedColors(clr_start + (nn-1), 'fixed'));
            end
            hold on

            % Lastly, plot the total droplet number concentration
            ax3 = subplot(1,3,3); plot(cloudDropletprofiles{indices(nn)}.total_Nc, norm_alt, ...
                'Color', mySavedColors(clr_start + (nn-1), 'fixed'));
            hold on


        else

            % First plot the LWC
            ax1 = subplot(1,4,1); plot(cloudDropletprofiles{indices(nn)}.lwc, cloudDropletprofiles{indices(nn)}.altitude, ...
                'Color', mySavedColors(clr_start + (nn-1), 'fixed'));
            hold on

            % next plot the effective radius
            % if the 2DC data is compliant, plot the effective radius computed
            % using both instruments
            if cloudDropletprofiles{indices(nn)}.flag_2DC_data_is_conforming==true
                ax2 = subplot(1,4,2); plot(cloudDropletprofiles{indices(nn)}.re, cloudDropletprofiles{indices(nn)}.altitude, ...
                    'Color', mySavedColors(clr_start + (nn-1), 'fixed'));
            else
                % if the 2DC data is non-conforming, use only the CDP data and
                % make a note of it
                ax2 = subplot(1,4,2); plot(cloudDropletprofiles{indices(nn)}.re_CDP, cloudDropletprofiles{indices(nn)}.altitude, ...
                    'Color', mySavedColors(clr_start + (nn-1), 'fixed'));
            end
            hold on

            % plot the total droplet number concentration
            ax3 = subplot(1,4,3); plot(cloudDropletprofiles{indices(nn)}.total_Nc, cloudDropletprofiles{indices(nn)}.altitude, ...
                'Color', mySavedColors(clr_start + (nn-1), 'fixed'));
            hold on

            % finally, plot the radisonde RH profiles, along with the
            % radiosonde determined cloud top height and base
            ax4 = subplot(1,4,4); plot(radiosondeProfiles.watVap_prof{indices(nn)}, radiosondeProfiles.altitude{indices(nn)}./1e3, ...
                'Color', mySavedColors(clr_start + (nn-1), 'fixed'));
            hold on
            
            % add horizontal a horizontal line depicting cloud top height
            cloudTopHeight = radiosondeProfiles.cloud_topHeight(indices(nn)) / 1e3;
            % yline(cloudTopHeight, '--', 'Color', 'k', "LineWidth", 2,...
            %     'LabelHorizontalAlignment','left', 'LabelVerticalAlignment','top',...
            %     'DisplayName', ['CTH = ', num2str(cloudTopHeight*1e3), ' m'],...
            %     'FontSize', 16);
            yline(cloudTopHeight, ':', 'Color', 'k', "LineWidth", 2)
            


        end

        legend_str{nn} = ['idx = ', num2str(indices(nn))];
        legend_str{nn+1} = ['CTH = ', num2str(cloudTopHeight*1e3), ' m'];



    end


    % Make each subplot pretty
    subplot(1,4,1)
    grid on; grid minor;
    xlabel('LWC ($g/m^3$)', 'Interpreter','latex', 'FontSize', ax_fnt)
    ylabel('Altitude ($m$)', 'Interpreter','latex', 'FontSize', ax_fnt)



    subplot(1,4,2)
    grid on; grid minor;
    % if the 2DC data is compliant, plot the effective radius computed
    % using both instruments
    if cloudDropletprofiles{indices(nn)}.flag_2DC_data_is_conforming==true
        xlabel('$r_e$ ($\mu m$)', 'Interpreter','latex', 'FontSize', ax_fnt)
    else
        % if the 2DC data is non-conforming, use only the CDP data and
        % make a note of it
        xlabel('$r_e$ ($\mu m$) - (CDP only)', 'Interpreter','latex')
    end

    % % include a title in the middle plot
    % if isfield(cloudDropletprofiles{indices(nn)}, 'LWC_threshold')==true
    %     title(['$LWC \geq$ ', num2str(cloudDropletprofiles{indices(nn)}.LWC_threshold),' $g/m^{3}$',...
    %         '   $N_c \geq$ ', num2str(cloudDropletprofiles{indices(nn)}.Nc_threshold), ' $cm^{-3}$'],...
    %         'interpreter', 'latex', 'FontSize', ttl_fnt)
    % 
    % elseif isfield(cloudDropletprofiles{indices(nn)}.inputs, 'LWC_threshold')==true
    %     title(['$LWC \geq$ ', num2str(cloudDropletprofiles{indices(nn)}.inputs.LWC_threshold),' $g/m^{3}$',...
    %         '   $N_c \geq$ ', num2str(cloudDropletprofiles{indices(nn)}.inputs.Nc_threshold), ' $cm^{-3}$'],...
    %         'interpreter', 'latex', 'FontSize', ttl_fnt)
    % 
    % end




    subplot(1,4,3)
    grid on; grid minor;
    xlabel('$N_c$ ($cm^{-3}$)', 'Interpreter','latex', 'FontSize', ax_fnt)




    subplot(1,4,4)
    grid on; grid minor;
    xlabel('$RH$ ($\%$)', 'Interpreter','latex', 'FontSize', ax_fnt)
    ylabel('Altitude ($km$)', 'Interpreter','latex', 'FontSize', ax_fnt)



    % in the third subplot, define the indices being plotted
    legend(legend_str, 'Interpreter','latex', 'Location','best',...
        'Color', 'white', 'TextColor', 'k', 'FontSize', lgnd_fnt)

    % set plot size
    set(gcf, 'Position', [0 0 1300 850])

    % link the yaxes so that they all have the same bounds
    linkaxes([ax1 ax2 ax3],'y')











elseif isstruct(cloudDropletprofiles)


    for nn = 1:N_cuvres

        % if normalize altitude is true, all altitude vectors will be
        % normalized between 0 and 1

        if normalize_altitude==true

            norm_alt = (cloudDropletprofiles(indices(nn)).altitude - min(cloudDropletprofiles(indices(nn)).altitude))./...
                (max(cloudDropletprofiles(indices(nn)).altitude) - min(cloudDropletprofiles(indices(nn)).altitude));

            % First plot the LWC
            ax1 = subplot(1,3,1); plot(cloudDropletprofiles(indices(nn)).lwc, norm_alt, ...
                'Color', mySavedColors(clr_start + (nn-1), 'fixed'));
            hold on

            % next plot the effective radius
            % if the 2DC data is compliant, plot the effective radius computed
            % using both instruments
            if cloudDropletprofiles(indices(nn)).flag_2DC_data_is_conforming==true
                ax2 = subplot(1,3,2); plot(cloudDropletprofiles(indices(nn)).re, norm_alt, ...
                    'Color', mySavedColors(clr_start + (nn-1), 'fixed'));
            else
                % if the 2DC data is non-conforming, use only the CDP data and
                % make a note of it
                ax2 = subplot(1,3,2); plot(cloudDropletprofiles(indices(nn)).re_CDP, norm_alt, ...
                    'Color', mySavedColors(clr_start + (nn-1), 'fixed'));
            end
            hold on

            % Lastly, plot the total droplet number concentration
            ax3 = subplot(1,3,3); plot(cloudDropletprofiles(indices(nn)).total_Nc, norm_alt, ...
                'Color', mySavedColors(clr_start + (nn-1), 'fixed'));
            hold on

        else

            % First plot the LWC
            ax1 = subplot(1,3,1); plot(cloudDropletprofiles(indices(nn)).lwc, cloudDropletprofiles(indices(nn)).altitude, ...
                'Color', mySavedColors(clr_start + (nn-1), 'fixed'));
            hold on

            % next plot the effective radius
            % if the 2DC data is compliant, plot the effective radius computed
            % using both instruments
            if cloudDropletprofiles(indices(nn)).flag_2DC_data_is_conforming==true
                ax2 = subplot(1,3,2); plot(cloudDropletprofiles(indices(nn)).re, cloudDropletprofiles(indices(nn)).altitude, ...
                    'Color', mySavedColors(clr_start + (nn-1), 'fixed'));
            else
                % if the 2DC data is non-conforming, use only the CDP data and
                % make a note of it
                ax2 = subplot(1,3,2); plot(cloudDropletprofiles(indices(nn)).re_CDP, cloudDropletprofiles(indices(nn)).altitude, ...
                    'Color', mySavedColors(clr_start + (nn-1), 'fixed'));
            end
            hold on

            % Lastly, plot the total droplet number concentration
            ax3 = subplot(1,3,3); plot(cloudDropletprofiles(indices(nn)).total_Nc, cloudDropletprofiles(indices(nn)).altitude, ...
                'Color', mySavedColors(clr_start + (nn-1), 'fixed'));
            hold on

        end

        legend_str{nn} = ['idx = ', num2str(indices(nn))];


    end



    % Make each subplot pretty
    subplot(1,3,1)
    grid on; grid minor;
    xlabel('LWC ($g/m^3$)', 'Interpreter','latex', 'FontSize', ax_fnt)
    ylabel('Altitude ($m$)', 'Interpreter','latex', 'FontSize', ax_fnt)



    subplot(1,3,2)
    grid on; grid minor;
    % if the 2DC data is compliant, plot the effective radius computed
    % using both instruments
    if cloudDropletprofiles(indices(nn)).flag_2DC_data_is_conforming==true
        xlabel('$r_e$ ($\mu m$)', 'Interpreter','latex', 'FontSize', ax_fnt)
    else
        % if the 2DC data is non-conforming, use only the CDP data and
        % make a note of it
        xlabel('$r_e$ ($\mu m$) - (CDP only)', 'Interpreter','latex', 'FontSize', ax_fnt)
    end

    % include a title in the middle plot
    if isfield(cloudDropletprofiles, 'LWC_threshold')==true
        title(['$LWC \geq$ ', num2str(cloudDropletprofiles(indices(nn)).LWC_threshold),' $g/m^{3}$',...
            '   $N_c \geq$ ', num2str(cloudDropletprofiles(indices(nn)).Nc_threshold), ' $cm^{-3}$'],...
            'interpreter', 'latex', 'FontSize', ttl_fnt)

    elseif isfield(cloudDropletprofiles(indices(nn)).inputs, 'LWC_threshold')==true
        title(['$LWC \geq$ ', num2str(cloudDropletprofiles(indices(nn)).inputs.LWC_threshold),' $g/m^{3}$',...
            '   $N_c \geq$ ', num2str(cloudDropletprofiles(indices(nn)).inputs.Nc_threshold), ' $cm^{-3}$'],...
            'interpreter', 'latex', 'FontSize', ttl_fnt)

    end




    subplot(1,3,3)
    grid on; grid minor;
    xlabel('$N_c$ ($cm^{-3}$)', 'Interpreter','latex', 'FontSize', ax_fnt)

    % in the third subplot, define the indices being plotted
    legend(legend_str, 'Interpreter','latex', 'Location','best',...
        'Color', 'white', 'TextColor', 'k', 'FontSize', lgnd_fnt)

    % set plot size
    set(gcf, 'Position', [0 0 1300 850])

    % link the yaxes so that they all have the same bounds
    linkaxes([ax1 ax2 ax3],'y')





end









end