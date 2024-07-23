%% Script for profiling the rms residual calcualtions. Figure out how to run faster!



% By Andrew John Buggee

%% 3D Interpolate the radiance calculations on a finer grid and compare with EMIT measurement


% Meshgrid is defined on x,y,z space, not row, column, depth space
% In 3D space, z = row, x = column, y = depth
[R_bot, R_top, Tau_c] = meshgrid(r_bot, r_top, tau_c);

% Create the new fine grid to interpolate on
% define the discrete step length of each variable
d_r_top = 0.1;      % microns
d_r_bot = 0.1;      % microns
d_tau_c = 0.05;

r_top_fine = r_top(1):d_r_top:r_top(end);
r_bot_fine = r_bot(1):d_r_bot:r_bot(end);
tau_c_fine = tau_c(1):d_tau_c:tau_c(end);

[R_bot_fine, R_top_fine, Tau_c_fine] = meshgrid(r_bot_fine, r_top_fine, tau_c_fine);

Refl_model_fine = zeros(length(r_top_fine), length(r_bot_fine), length(tau_c_fine), size(Refl_model,4));


tic
for wl = 1:size(Refl_model,4)

    Refl_model_fine(:,:,:,wl) = interp3(R_bot, R_top, Tau_c, Refl_model(:, :, :, wl),...
        R_bot_fine, R_top_fine, Tau_c_fine);


end
toc


% Using the new fine grid, calculate how many sets of measurements are
% within the EMIT measurement and it's uncertainty


redundant_states = [];
rms_residual = zeros(length(r_top_fine), length(r_bot_fine), length(tau_c_fine));


tic
for rt = 1:size(Refl_model_fine,1)


    for rb = 1:size(Refl_model_fine,2)


        parfor tc = 1:size(Refl_model_fine,3)

            % Check to see if the radiance computed by the model is
            % within the listed uncertainty for EMIT
            %redundant_states(rt,rb,tc) = all(abs(R_emit - reshape(R_model_fine(rt,rb,tc,:), 1, [])) <= R_emit_uncert);
            redundant_states = [redundant_states, abs(Refl_emit - reshape(Refl_model_fine(rt,rb,tc,:), [], 1)) <= Refl_emit_uncertainty];
            rms_residual(rt, rb, tc) = sqrt(mean( (Refl_emit - reshape(Refl_model_fine(rt,rb,tc,:), [], 1)).^2) );


        end
    end
end
toc

% Save Refl_model_file and the rms_residual, because these calculations
% take a while!
save(filename,"Refl_model_fine", "rms_residual", "-append");

% Find the number states that lead to modeled measurements within the EMIT
% measurement uncertainty
num_states = sum(all(redundant_states, 1));

% print the percentage of redundant states
disp([newline, num2str(100* (num_states/(numel(r_top_fine)*numel(r_bot_fine)*numel(tau_c_fine)))),...
    '% of modeled states that produce measurements within the EMIT uncertainty', newline])