%% Computing LWC assuming total number concentration is constant with height

% pick an index
idx = 8;

% compute the mean difference between each vertical sample point
dH = mean(diff(vert_profs.altitude{idx} - min(vert_profs.altitude{idx})));

% define the density of water in grams per cubic meter
density  = 10^6;        % g/m^3

% grab the optical depth of the cloud
tau_c = vert_profs.tau{idx}(end);

% convert the effective radius to meters
re_meters = (vert_profs.re{idx}*1e-6);            % meters

test_lwc = (2*tau_c * density * re_meters.^3)./(3 * dH * sum(re_meters.^2));

