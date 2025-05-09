%% 3D interpolate HySICS simulated reflectance for each wavelength to create a lookup table

% By Andrew John Buggee

%%

function [Refl_model_fine, r_top_fine, r_bot_fine, tau_c_fine] = lookUp_table_3D_interpolate_HySICS(forward_model_calcs)


% I could start with a coarse grid and cascade down to finer scales. This way I
% don't have to run the 3D interpolation on the entire look-up table at a
% fine resolution

% Meshgrid is defined on x,y,z space, not row, column, depth space
% In 3D space, z = row, x = column, y = depth
[R_bot, R_top, Tau_c] = meshgrid(forward_model_calcs.inputs.RT.r_bot, forward_model_calcs.inputs.RT.r_top,...
    forward_model_calcs.inputs.RT.tau_c);

% Create the new fine grid to interpolate on
% define the discrete step length of each variable
d_r_top = 0.05;      % microns
d_r_bot = 0.05;      % microns
d_tau_c = 0.01;

r_top_fine = forward_model_calcs.inputs.RT.r_top(1):d_r_top:forward_model_calcs.inputs.RT.r_top(end);
r_bot_fine = forward_model_calcs.inputs.RT.r_bot(1):d_r_bot:forward_model_calcs.inputs.RT.r_bot(end);
tau_c_fine = forward_model_calcs.inputs.RT.tau_c(1):d_tau_c:forward_model_calcs.inputs.RT.tau_c(end);

[R_bot_fine, R_top_fine, Tau_c_fine] = meshgrid(r_bot_fine, r_top_fine, tau_c_fine);

Refl_model_fine = zeros(size(forward_model_calcs.Refl_model,1), length(r_top_fine),...
    length(r_bot_fine), length(tau_c_fine));


tic
for wl = 1:size(forward_model_calcs.Refl_model,1)

    Refl_model_fine(wl,:,:,:) = interp3(R_bot, R_top, Tau_c, reshape(forward_model_calcs.Refl_model(wl, :, :, :),...
        length(forward_model_calcs.inputs.RT.r_top), length(forward_model_calcs.inputs.RT.r_bot), []),...
        R_bot_fine, R_top_fine, Tau_c_fine);


end
toc



end
