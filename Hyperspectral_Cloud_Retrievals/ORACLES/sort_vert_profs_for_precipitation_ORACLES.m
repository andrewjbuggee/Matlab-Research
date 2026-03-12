%% Sort ORACLES vertical profiles by presence of drizzle/precipitation

% Analogous to sort_vert_profs_for_precipitation.m (VOCALS-REx), adapted
% for ORACLES profiles where precipitation is estimated from the 2DS+HVPS-3
% liquid water path (lwp_2DS_HVPS), which covers droplets with D > 50 µm.

% INPUTS:
% -------
%   vert_profs                   - array of vertical profile structures from
%                                  find_verticalProfiles_ORACLES.m
%   precipitation_drizzle_threshold - (g/m^2) rain water path threshold above
%                                     which a profile is considered precipitating

% OUTPUTS:
% --------
%   index_precip_drizzle - indices into vert_profs of profiles that exceed
%                          the precipitation threshold

% By Andrew John Buggee

%%

function [index_precip_drizzle] = sort_vert_profs_for_precipitation_ORACLES(vert_profs, precipitation_drizzle_threshold)


index_precip_drizzle = [];

for nn = 1:length(vert_profs)

    % Check that the 2DS+HVPS liquid water path is non-negative
    if vert_profs(nn).lwp_2DS_HVPS < 0
        error([newline, 'lwp_2DS_HVPS is negative for profile ', num2str(nn), ...
            '! Check the liquid water path calculation.', newline])
    end

    % If the rain water path exceeds the threshold, flag as precipitating
    if vert_profs(nn).lwp_2DS_HVPS > precipitation_drizzle_threshold
        index_precip_drizzle = [index_precip_drizzle, nn];
    end

end


end
