% function that tallys photons moving along the positive or negative z
% axis. Computes the downwelling and upwelling flux for 3D Monte Carlo
% Model

% By Andrew John Buggee

%%

function [] = up_and_downwelling_flux_within_3D_monteCarlo()

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


% Create the bins that will keep tally
N_counts_moving_up = zeros(1, N_bins);
N_counts_moving_down = zeros(1, N_bins);

% We want to tally the direction a photon moves through each bin defined
% above

% Check to see if we switched directions. If we did, we want to
% make sure we tally that a photon was moving in both
% directions in the bin where it turned around


% assign each x to one of our bins
for nn=1:N_photons



    % Compute whether or not the photon was moving in the positive z
    % direction or the negative z direction. A value of 1 tells the code
    % that the photon was increasing in tau along the z axis, and a value
    % of negative 1 tells the code the photon was decreasing in tau along
    % the y axis. The first value is always 1 since the photons starts its
    % journey by moving down, increasing its value of tau along the z axis
    z_direction = [1; sign(diff(photon_tau_absolute_position{nn}(:,3)))];


    for tt = 1:size(photon_tau_absolute_position{nn},1)-1

        % Check to see if the photon is moving down or up and check to see
        % if its continuing along the same direction, or if it's switched
        % directions

        %[nn,tt]

        % ----------------------------------------
        % **** Photon continuing to move down ****
        % ----------------------------------------
        if z_direction(tt+1)==1 && z_direction(tt)==1
            % *** the photon is moving down ***
            % Check the see which tau bins the photon moves through. Tally each
            % bin that is found


            % If the photon is moving down and was already heading
            % down, we don't need to account for an additional bin

            % If the photon is moving down, then when it cross bin-edge 3,
            % its in bucket 3.
            bins_photon_moves_through = binEdges>=photon_tau_absolute_position{nn}(tt,3) & ...
                binEdges<photon_tau_absolute_position{nn}(tt+1,3);

            N_counts_moving_down(bins_photon_moves_through) = N_counts_moving_down(bins_photon_moves_through) +1;


            % --------------------------------------------
            % **** Photon moving down after moving up ****
            % --------------------------------------------
        elseif z_direction(tt+1)==1 && z_direction(tt)==-1
            % In this case the photon has turned around in the bin it was
            % left off in. Therefore we have to count this as a photon
            % moving down in this bin.

            % Photon is moving down, so the starting bin will be smaller
            % than, or equal to, the final bin
            bin_edges_less_than_start_position = find(binEdges<=photon_tau_absolute_position{nn}(tt,3));

            % This bin should be greater than, or equal to, the starting
            % bin
            bin_edges_less_than_end_position = find(binEdges<photon_tau_absolute_position{nn}(tt+1,3));

            bins_photon_moves_through = zeros(1,N_bins+1);
            bins_photon_moves_through(bin_edges_less_than_start_position(end):bin_edges_less_than_end_position(end)) = 1;
            bins_photon_moves_through = logical(bins_photon_moves_through);


            N_counts_moving_down(bins_photon_moves_through) = N_counts_moving_down(bins_photon_moves_through) +1;



            % --------------------------------------
            % **** Photon continuing to move up ****
            % --------------------------------------
        elseif z_direction(tt+1)==-1 && z_direction(tt)==-1
            % *** the photon is moving up ***
            % When this is true, depth_travlled{nn}(tt+1) is always less
            % than depth_travelled{nn}(tt) so we need to determine the tau
            % bins that the photon passes through between
            % depth_travlled{nn}(tt) and depth_travlled{nn}(tt+1)


            % If the photon is continuing up then we only have to count
            % the bins it passes through

            % Photon is moving up, so the starting bin will be larger
            % than, or equal to, the final bin
            bin_edges_less_than_start_position = find(binEdges<=photon_tau_absolute_position{nn}(tt,3));

            % This bin should be less than, or equal to, the starting
            % bin
            bin_edges_less_than_end_position = find(binEdges<=photon_tau_absolute_position{nn}(tt+1,3));

            bins_photon_moves_through = zeros(1,N_bins+1);

            % let's check to see if the photon is still in the same bin. If
            % it is, we've already accounted for it's upward motion, so we
            % don't want to double count this.
            if bin_edges_less_than_end_position(end)==bin_edges_less_than_start_position(end)
                % If this is true, the photon is continuing to move upward
                % in the same bin, and we don't count it
            else

                % If the photon is continuing to move upwards, that means
                % the previous loop already accounted for the bin it's
                % currently in, either because the photon turned around and
                % ended up in the current bin, or because the photon
                % continued to move upwards into the current bin. The logic
                % above will count the bin its currently in, so we need to
                % remove this value. We do this by subtracting 1 from the
                % star position.
                bins_photon_moves_through(bin_edges_less_than_end_position(end):bin_edges_less_than_start_position(end-1)) = 1;

            end

            bins_photon_moves_through = logical(bins_photon_moves_through);

            N_counts_moving_up(bins_photon_moves_through) = N_counts_moving_up(bins_photon_moves_through) +1;




            % --------------------------------------------
            % **** Photon moving up after moving down ****
            % --------------------------------------------
        elseif z_direction(tt+1)==-1 && z_direction(tt)==1

            % If the photon switched direction, we have to account for
            % the final bin it ends up in

            % Photon is moving up, so the starting bin will be larger
            % than, or equal to, the final bin
            bin_edges_less_than_start_position = find(binEdges<=photon_tau_absolute_position{nn}(tt,3));

            % The final bin will be smaller than, or equal to, to starting
            % bin
            bin_edges_less_than_end_position = find(binEdges<=photon_tau_absolute_position{nn}(tt+1,3));

            bins_photon_moves_through = zeros(1,N_bins+1);
            bins_photon_moves_through(bin_edges_less_than_end_position(end):bin_edges_less_than_start_position(end)) = 1;
            bins_photon_moves_through = logical(bins_photon_moves_through);


            % Some photons in this category will reflect off the bottom
            % boudnary. If this happened, we simply need to remove the
            % logical true value for the binEdge equal to tau_upper_limit.
            % If we don't we get an error, and all we need to keep track of
            % is whether or not the photon passed through this bin
            if bins_photon_moves_through(end)==true && (photon_tau_absolute_position{nn}(tt,3)==tau_z_upper_limit || photon_tau_absolute_position{nn}(tt+1,3)==tau_z_upper_limit)
                bins_photon_moves_through(end) = 0;
            end




            N_counts_moving_up(bins_photon_moves_through) = N_counts_moving_up(bins_photon_moves_through) +1;



        else


            error([newline,'Im not sure what the photon is doing',newline])

        end


    end

end


end
