% -----------------------------------------------------------------
% ------ Compute SEQUENTIAL information content (Rodgers 1996) ----
% ------------------------------------------------------------------

num_channels = numel(GN_inputs.bands2run);
num_parameters = GN_inputs.num_model_parameters;

% Initialize arrays
posterior_cov_sequential = zeros(num_parameters, num_parameters, num_channels+1);
info_content_sequential = zeros(num_channels, 1);
dfs_sequential = zeros(num_channels, 1);

% Start with prior covariance
posterior_cov_sequential(:,:,1) = GN_inputs.model.covariance;

% Sequential channel addition
for nn = 1:num_channels
    
    % Get the weighting function (Jacobian) for channel nn
    k_i = GN_outputs.Jacobian_final(nn,:)';  % Make it a column vector
    
    % Assumed measurement error variance
    sigma_i_squared = GN_inputs.measurement.covariance(nn,nn);  % Assuming diagonal measurement covariance
    
    % Update posterior covariance using equation (11) from Rodgers
    S_prev = posterior_cov_sequential(:,:,nn);
    numerator = (S_prev * k_i) * (S_prev * k_i)';
    denominator = sigma_i_squared + k_i' * S_prev * k_i;
    
    posterior_cov_sequential(:,:,nn+1) = S_prev - numerator / denominator;
    
    % Calculate information content change using equation (13)
    % δS_i = ln(1 + k_i^T * S_(i-1) * k_i / σ_i^2)
    info_content_sequential(nn) = log2(1 + k_i' * S_prev * k_i / sigma_i_squared);
    
    % Calculate degrees of freedom change using equation (15)
    % δd_s = (S_(i-1) * k_i)^T * (S_(i-1) * k_i) / (σ_i^2 + k_i^T * S_(i-1) * k_i)
    S_k = S_prev * k_i;
    dfs_sequential(nn) = (S_k' * S_k) / (sigma_i_squared + k_i' * S_prev * k_i);
    
end

% Cumulative information content
cumulative_info_content = cumsum(info_content_sequential);

%%
% ----------------------------------------------------------------------
% ------ Compute INDIVIDUAL information content (relative to prior) ----
% ----------------------------------------------------------------------

S_prior = GN_inputs.model.covariance;
info_content_individual = zeros(num_channels, 1);
dfs_individual = zeros(num_channels, 1);

% Initialize storage for individual parameter information content
info_content_by_param = zeros(num_channels, num_parameters);
info_content_by_param_sequential = zeros(num_channels, num_parameters);

for nn = 1:num_channels
    
    % Get the weighting function for channel nn
    k_i = GN_outputs.Jacobian_final(nn,:)';
    
    % Measurement error variance
    sigma_i_squared = GN_inputs.measurement.covariance(nn,nn);
    
    % Calculate posterior covariance using only this channel
    numerator = (S_prior * k_i) * (S_prior * k_i)';
    denominator = sigma_i_squared + k_i' * S_prior * k_i;
    
    S_post = S_prior - numerator / denominator;
    
    % Information content relative to prior
    info_content_individual(nn) = log2(1 + k_i' * S_prior * k_i / sigma_i_squared);
    
    % Degrees of freedom relative to prior
    S_k = S_prior * k_i;
    dfs_individual(nn) = (S_k' * S_k) / (sigma_i_squared + k_i' * S_prior * k_i);


    % Compute information content for each parameter group
    for pp = 1:num_parameters

        % Extract submatrices for this parameter group
        S_prior_param = S_prior(pp, pp);
        S_post_param = S_post_total(pp, pp);

        % Information content for this parameter group
        % H = -log2(det(S_post)/det(S_prior))
        info_content_by_param(nn, pp) = log2(det(S_prior_param) / det(S_post_param));

    end

end

%%
% -------------------------------------------------------------
% ----------------------- Plotting ----------------------------
% -------------------------------------------------------------

fntSze1 = 20;
fntSze2 = 35;

figure;
subplot(2,2,1)
plot(1:num_channels, info_content_sequential, 'b.-', 'MarkerSize', 10, 'LineWidth', 1.5)
hold on
plot(1:num_channels, info_content_individual, 'r.-', 'MarkerSize', 10, 'LineWidth', 1.5)
grid on; grid minor
xlabel('Channel Number', 'FontSize', fntSze1, 'Interpreter','latex')
ylabel('Information Content (bits)', 'FontSize', fntSze1, 'Interpreter','latex')
title('Sequential vs Individual Information Content', 'FontSize', fntSze1, 'Interpreter','latex')
legend('Sequential', 'Individual', 'Location', 'best')

subplot(2,2,2)
plot(1:num_channels, cumulative_info_content, 'g.-', 'MarkerSize', 10, 'LineWidth', 1.5)
grid on; grid minor
xlabel('Channel Number', 'FontSize', fntSze1, 'Interpreter','latex')
ylabel('Cumulative Information Content (bits)', 'FontSize', fntSze1, 'Interpreter','latex')
title('Cumulative Information Content', 'FontSize', fntSze1, 'Interpreter','latex')

subplot(2,2,3)
plot(1:num_channels, dfs_sequential, 'b.-', 'MarkerSize', 10, 'LineWidth', 1.5)
hold on
plot(1:num_channels, dfs_individual, 'r.-', 'MarkerSize', 10, 'LineWidth', 1.5)
grid on; grid minor
xlabel('Channel Number', 'FontSize', fntSze1, 'Interpreter','latex')
ylabel('Degrees of Freedom for Signal', 'FontSize', fntSze1, 'Interpreter','latex')
title('Sequential vs Individual DFS', 'FontSize', fntSze1, 'Interpreter','latex')
legend('Sequential', 'Individual', 'Location', 'best')

subplot(2,2,4)
plot(1:num_channels, cumsum(dfs_sequential), 'g.-', 'MarkerSize', 10, 'LineWidth', 1.5)
grid on; grid minor
xlabel('Channel Number', 'FontSize', fntSze1, 'Interpreter','latex')
ylabel('Cumulative DFS', 'FontSize', fntSze1, 'Interpreter','latex')
title('Cumulative Degrees of Freedom', 'FontSize', fntSze1, 'Interpreter','latex')

% -------------------------------------------------------------
% ------ Information spectrum (like Figure 1 in the paper) ----
% -------------------------------------------------------------

figure;
plot(mean(GN_inputs.RT.wavelengths2run, 2), info_content_individual, '.-', 'MarkerSize', 15, 'LineWidth', 1.5)
grid on; grid minor
xlabel('Wavelength', 'FontSize', fntSze2, 'Interpreter','latex')
ylabel('Information Content (bits)', 'FontSize', fntSze2, 'Interpreter','latex')
title('Information Spectrum - Individual Channel Contributions', 'FontSize', fntSze2, 'Interpreter','latex')


% If you want to see the effect of selecting channels sequentially
% (like the lower curve in Figure 1), you would need to implement
% the channel selection algorithm from Section 3