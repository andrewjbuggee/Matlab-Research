% -----------------------------------------------------------------
% ------ Compute INDIVIDUAL PARAMETER information content --------
% -----------------------------------------------------------------

% Assume you have defined parameter groups (modify indices as needed)
% Example parameter groupings - adjust these to match your state vector
temp_indices = 1:66;        % Temperature profile indices
water_indices = 67:132;     % Water vapor profile indices  
surface_indices = 133;      % Surface temperature index
% Add more parameter groups as needed

% Create parameter group definitions
param_groups = struct();
param_groups.temperature = temp_indices;
param_groups.water_vapor = water_indices;
param_groups.surface = surface_indices;
param_names = fieldnames(param_groups);

% Initialize storage for individual parameter information content
info_content_by_param = zeros(num_channels, length(param_names));
info_content_by_param_sequential = zeros(num_channels, length(param_names));

S_prior = GN_inputs.model.covariance;

%%
% ----------------------------------------------------------------------
% ------ Individual parameter information content (relative to prior) --
% ----------------------------------------------------------------------

for nn = 1:num_channels
    
    % Get the weighting function for channel nn
    k_i = GN_outputs.Jacobian_final(nn,:)';
    
    % Measurement error variance
    sigma_i_squared = GN_inputs.measurement.covariance(nn,nn);
    
    % Compute total posterior covariance for this channel
    numerator = (S_prior * k_i) * (S_prior * k_i)';
    denominator = sigma_i_squared + k_i' * S_prior * k_i;
    S_post_total = S_prior - numerator / denominator;
    
    % Compute information content for each parameter group
    for pp = 1:length(param_names)
        param_indices = param_groups.(param_names{pp});
        
        % Extract submatrices for this parameter group
        S_prior_param = S_prior(param_indices, param_indices);
        S_post_param = S_post_total(param_indices, param_indices);
        
        % Information content for this parameter group
        % H = -log2(det(S_post)/det(S_prior))
        info_content_by_param(nn, pp) = log2(det(S_prior_param) / det(S_post_param));
    end
end

%%
% ----------------------------------------------------------------------
% ------ Sequential parameter information content ---------------------
% ----------------------------------------------------------------------

% Initialize sequential posterior covariance
posterior_cov_param_sequential = zeros(num_parameters, num_parameters, num_channels+1);
posterior_cov_param_sequential(:,:,1) = S_prior;

for nn = 1:num_channels
    
    % Get the weighting function for channel nn
    k_i = GN_outputs.Jacobian_final(nn,:)';
    
    % Measurement error variance
    sigma_i_squared = GN_inputs.measurement.covariance(nn,nn);
    
    % Update posterior covariance sequentially
    S_prev = posterior_cov_param_sequential(:,:,nn);
    numerator = (S_prev * k_i) * (S_prev * k_i)';
    denominator = sigma_i_squared + k_i' * S_prev * k_i;
    
    S_current = S_prev - numerator / denominator;
    posterior_cov_param_sequential(:,:,nn+1) = S_current;
    
    % Compute sequential information content for each parameter group
    for pp = 1:length(param_names)
        param_indices = param_groups.(param_names{pp});
        
        % Extract submatrices for this parameter group
        S_prev_param = S_prev(param_indices, param_indices);
        S_current_param = S_current(param_indices, param_indices);
        
        % Sequential information content change for this parameter group
        info_content_by_param_sequential(nn, pp) = log2(det(S_prev_param) / det(S_current_param));
    end
end

%%
% -------------------------------------------------------------
% ----------------------- Enhanced Plotting ------------------
% -------------------------------------------------------------

fntSze1 = 20;
fntSze2 = 35;

% Colors for different parameters
colors = lines(length(param_names));

% Plot 1: Total and individual parameter information content
figure;
subplot(2,2,1)
plot(1:num_channels, info_content_individual, 'k.-', 'MarkerSize', 8, 'LineWidth', 2)
hold on
for pp = 1:length(param_names)
    plot(1:num_channels, info_content_by_param(:,pp), '.-', 'Color', colors(pp,:), ...
         'MarkerSize', 6, 'LineWidth', 1.5, 'DisplayName', param_names{pp})
end
grid on; grid minor
xlabel('Channel Number', 'FontSize', fntSze1, 'Interpreter','latex')
ylabel('Information Content (bits)', 'FontSize', fntSze1, 'Interpreter','latex')
title('Individual Parameter Information Content', 'FontSize', fntSze1, 'Interpreter','latex')
legend('Total', param_names{:}, 'Location', 'best')

% Plot 2: Sequential parameter information content
subplot(2,2,2)
plot(1:num_channels, info_content_sequential, 'k.-', 'MarkerSize', 8, 'LineWidth', 2)
hold on
for pp = 1:length(param_names)
    plot(1:num_channels, info_content_by_param_sequential(:,pp), '.-', 'Color', colors(pp,:), ...
         'MarkerSize', 6, 'LineWidth', 1.5, 'DisplayName', param_names{pp})
end
grid on; grid minor
xlabel('Channel Number', 'FontSize', fntSze1, 'Interpreter','latex')
ylabel('Information Content (bits)', 'FontSize', fntSze1, 'Interpreter','latex')
title('Sequential Parameter Information Content', 'FontSize', fntSze1, 'Interpreter','latex')
legend('Total', param_names{:}, 'Location', 'best')

% Plot 3: Cumulative parameter information content
subplot(2,2,3)
plot(1:num_channels, cumulative_info_content, 'k.-', 'MarkerSize', 8, 'LineWidth', 2)
hold on
for pp = 1:length(param_names)
    cumulative_param = cumsum(info_content_by_param_sequential(:,pp));
    plot(1:num_channels, cumulative_param, '.-', 'Color', colors(pp,:), ...
         'MarkerSize', 6, 'LineWidth', 1.5, 'DisplayName', param_names{pp})
end
grid on; grid minor
xlabel('Channel Number', 'FontSize', fntSze1, 'Interpreter','latex')
ylabel('Cumulative Information Content (bits)', 'FontSize', fntSze1, 'Interpreter','latex')
title('Cumulative Parameter Information Content', 'FontSize', fntSze1, 'Interpreter','latex')
legend('Total', param_names{:}, 'Location', 'best')

% Plot 4: Parameter information spectrum vs wavelength
subplot(2,2,4)
for pp = 1:length(param_names)
    plot(mean(GN_inputs.RT.wavelengths2run, 2), info_content_by_param(:,pp), '.-', ...
         'Color', colors(pp,:), 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', param_names{pp})
    hold on
end
grid on; grid minor
xlabel('Wavelength', 'FontSize', fntSze1, 'Interpreter','latex')
ylabel('Information Content (bits)', 'FontSize', fntSze1, 'Interpreter','latex')
title('Parameter Information Spectrum', 'FontSize', fntSze1, 'Interpreter','latex')
legend(param_names{:}, 'Location', 'best')

% Enhanced information spectrum plot (like Figure 1 in the paper)
figure;
for pp = 1:length(param_names)
    plot(mean(GN_inputs.RT.wavelengths2run, 2), info_content_by_param(:,pp), '.-', ...
         'Color', colors(pp,:), 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', param_names{pp})
    hold on
end
% Add total information content
plot(mean(GN_inputs.RT.wavelengths2run, 2), info_content_individual, 'k.-', ...
     'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'Total')
grid on; grid minor
xlabel('Wavelength', 'FontSize', fntSze2, 'Interpreter','latex')
ylabel('Information Content (bits)', 'FontSize', fntSze2, 'Interpreter','latex')
title('Information Spectrum by Parameter', 'FontSize', fntSze2, 'Interpreter','latex')
legend('Location', 'best')

% Optional: Create a stacked area plot to show parameter contributions
figure;
area_data = info_content_by_param;
area(mean(GN_inputs.RT.wavelengths2run, 2), area_data);
grid on; grid minor
xlabel('Wavelength', 'FontSize', fntSze2, 'Interpreter','latex')
ylabel('Information Content (bits)', 'FontSize', fntSze2, 'Interpreter','latex')
title('Stacked Parameter Information Spectrum', 'FontSize', fntSze2, 'Interpreter','latex')
legend(param_names{:}, 'Location', 'best')
colormap(colors)

% Print summary statistics
fprintf('\n=== PARAMETER INFORMATION CONTENT SUMMARY ===\n');
fprintf('Total information content: %.2f bits\n', sum(info_content_individual));
for pp = 1:length(param_names)
    param_total = sum(info_content_by_param(:,pp));
    param_percent = 100 * param_total / sum(info_content_individual);
    fprintf('%s: %.2f bits (%.1f%%)\n', param_names{pp}, param_total, param_percent);
end