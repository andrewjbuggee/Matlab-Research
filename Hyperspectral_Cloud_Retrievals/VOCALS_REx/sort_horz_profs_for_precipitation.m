%% Search through the vertical profiles and split them into precipitating
% and non-precipitating cases

% --------------------- INPUTS ------------------------



% By Andrew John Buggee

%%

function [index_precip_drizzle] = sort_horz_profs_for_precipitation(horz_profs, precipitation_drizzle_threshold)


% Step through each horizontal profile found
% If the total Liquid water Path calcualted by the 2DC instruemnt, which
% measures droplet sizes from 50 to 1500 microns in diameter, we can safely
% assume some drizzle and/or precipitation

% define a 2DC LWP threshold that will define the presence of precipitation
%precipitation_drizzle_threshold = 10;       % g/m^2

index_precip_drizzle = [];

for nn = 1:length(horz_profs.lwc_2DC)

    % integrate the lwc content over the distance travelled
    LWP_2DC = trapz(horz_profs.horz_dist{nn}, horz_profs.lwc_2DC{nn});      % g/m^2

    if LWP_2DC>precipitation_drizzle_threshold
        index_precip_drizzle = [index_precip_drizzle, nn];
    end

end



% precipitation is defined if the 2DC LWP is larger than the CDP LWP

% index_precipitation = [];
% 
% for nn = 1:length(vert_profs.lwp_2DC)
% 
%     if vert_profs.lwp_2DC{nn}>vert_profs.lwp_CDP{ii}
%         index_precipitation = [index_precipitation, nn];
%     end
% 
% end







end