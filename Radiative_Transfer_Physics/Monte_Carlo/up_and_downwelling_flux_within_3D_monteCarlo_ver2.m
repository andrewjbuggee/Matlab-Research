% ----------------------------------------------------------------------------------
% *** Divide the tau space into a grid and compute reflectance and transmittance ***
% ----------------------------------------------------------------------------------

% Set up the tau grid first

N_bins = 500;

if tau_z_upper_limit==inf
    binEdges = logspace(-3, ceil(log10(max(max_position_reached))),N_bins+1);

else
    binEdges = linspace(tau_z_lower_limit, tau_z_upper_limit, N_bins+1);
end


% Initialize final count arrays before parfor
N_counts_moving_down = zeros(1, N_bins);
N_counts_moving_up = zeros(1, N_bins);

parfor nn = 1:N_photons

    % Extract photon data safely for each iteration
    photon_data = photon_tau_absolute_position{nn};

    % Local tallies for this specific photon (memory-efficient)
    local_counts_down = zeros(1, N_bins);
    local_counts_up = zeros(1, N_bins);

    % Compute z-direction movement
    z_direction = [1; sign(diff(photon_data(:,3)))];

    for tt = 1:size(photon_data,1)-1

        % ----------------------------------------
        % **** Photon continuing to move down ****
        % ----------------------------------------
        if z_direction(tt+1)==1 && z_direction(tt)==1
            % *** the photon is moving down ***
            % Check the see which tau bins the photon moves through. Tally each
            % bin that is found

            bins_photon_moves_through = find(binEdges>=photon_data(tt,3) & ...
                                             binEdges<photon_data(tt+1,3));
            local_counts_down(bins_photon_moves_through) = local_counts_down(bins_photon_moves_through) + 1;

            % --------------------------------------------
            % **** Photon moving down after moving up ****
            % --------------------------------------------
        elseif z_direction(tt+1)==1 && z_direction(tt)==-1
            % In this case the photon has turned around in the bin it was
            % left off in. Therefore we have to count this as a photon
            % moving down in this bin.

            bins_photon_moves_through = find(binEdges<=photon_data(tt,3) & binEdges<photon_data(tt+1,3));
            local_counts_down(bins_photon_moves_through) = local_counts_down(bins_photon_moves_through) + 1;

            % --------------------------------------
            % **** Photon continuing to move up ****
            % --------------------------------------
        elseif z_direction(tt+1)==-1 && z_direction(tt)==-1
            % *** the photon is moving up ***
            % When this is true, depth_travlled{nn}(tt+1) is always less
            % than depth_travelled{nn}(tt) so we need to determine the tau
            % bins that the photon passes through between
            % depth_travlled{nn}(tt) and depth_travlled{nn}(tt+1)

            bins_photon_moves_through = find(binEdges<=photon_data(tt,3) & binEdges<=photon_data(tt+1,3));
            local_counts_up(bins_photon_moves_through) = local_counts_up(bins_photon_moves_through) + 1;

            % --------------------------------------------
            % **** Photon moving up after moving down ****
            % --------------------------------------------
        elseif z_direction(tt+1)==-1 && z_direction(tt)==1
            % If the photon switched direction, we have to account for
            % the final bin it ends up in
            
            bins_photon_moves_through = find(binEdges<=photon_data(tt,3) & binEdges<=photon_data(tt+1,3));
            local_counts_up(bins_photon_moves_through) = local_counts_up(bins_photon_moves_through) + 1;
        else
            error("Unexpected photon behavior detected.");
        end
    end

    % Accumulate results using MATLAB's Reduction feature
    % This ensures no large temporary arrays are created
    N_counts_moving_down = N_counts_moving_down + local_counts_down;
    N_counts_moving_up = N_counts_moving_up + local_counts_up;
end
