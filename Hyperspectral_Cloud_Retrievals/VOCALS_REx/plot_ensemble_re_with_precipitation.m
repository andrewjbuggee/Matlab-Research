%% Plot ensemble vertical profiles with precipitation

%% Simple plot, check to see what direction the plane was moving when the data was collected

srt_idx = 1;
end_idx = 15;
legend_str = [];

figure;

for nn = srt_idx:end_idx

    % Determine whether or not the plane ascended or descended while
    % samplign the cloud
    dz = ensemble_profiles.altitude{nn}(end) - ensemble_profiles.altitude{nn}(1);

    
    if dz>0
        % if true, then the plane was ascending.

        plot(ensemble_profiles.re{nn}, ...
            (ensemble_profiles.altitude{nn} - ensemble_profiles.altitude{nn}(1))./...
            (ensemble_profiles.altitude{nn}(end) - ensemble_profiles.altitude{nn}(1)))
        hold on


    elseif dz<0
        % if true then the plane was descending.

        plot(ensemble_profiles.re{nn}, ...
            (ensemble_profiles.altitude{nn} - ensemble_profiles.altitude{nn}(end))./...
            (ensemble_profiles.altitude{nn}(1) - ensemble_profiles.altitude{nn}(end)))
        hold on

    end

    legend_str{nn-srt_idx+1} = ['$LWP_{2DC} = $', num2str(abs(ensemble_profiles.lwp_2dc(nn))), ' ($g/m^2$)'];


end

grid on; grid minor
xlabel('$r_e$ ($\mu m$)', 'Interpreter','latex');
ylabel('Altitude ($m$)', 'Interpreter','latex');

legend(legend_str, 'Location','best', 'Interpreter','latex')


% set plot size
set(gcf, 'Position', [0 0 1000 550])


