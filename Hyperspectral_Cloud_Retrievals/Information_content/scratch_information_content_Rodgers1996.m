% -------------------------------------------------------------
% ------ Compute the retrieval covariance for each channel ----
posterior_cov_perChannel = zeros(num_parameters, num_parameters, num_bands+1);
posterior_cov_perChannel(:,:,1) = model_cov;

dH = zeros(num_parameters, num_parameters, num_bands);

H = zeros(num_parameters, num_bands);

for nn = 2:(num_bands+1)

    
    % The retrieval covariance for channels 1:nn
    posterior_cov_perChannel(:,:, nn) = posterior_cov_perChannel(:,:, nn-1) -...
        (posterior_cov_perChannel(:,:, nn-1) * Jacobian(nn-1,:)')*(posterior_cov_perChannel(:,:, nn-1) * Jacobian(nn-1,:)')'/...
        (1 + (posterior_cov_perChannel(:,:, nn-1) * Jacobian(nn-1,:)')' * Jacobian(nn-1,:)');


    % The change in information content 
    dH(:,:, nn-1) = 1/2 * (log2(posterior_cov_perChannel(:,:,nn-1)) - log2(posterior_cov_perChannel(:,:,nn)));

    % Let's grab just the main diagonal components and take the square root
    H(:, nn-1) = sqrt(diag(dH(:,:,nn-1)));

end


% plot the change in bits for each channel and each parameter
figure; plot(mean(GN_inputs.RT.wavelengths2run, 2), H, '.-', 'MarkerSize', 20, 'LineWidth', 1.5)
grid on; grid minor
xlabel('Spectral Channel')
ylabel('Change in bits')
title('Sequential information change in bits w/ increasing channel', 'FontSize',20)

legend('r-top', 'r-bot', 'tau_c', 'cwv', 'Location', 'best')


%% Which channels have the highest information content above that of the a priori?



H_above_aPriori = zeros(num_parameters, num_bands);

for nn = 1:num_bands

    posterior_cov_perChannel_above_apriori = [];

    % The retrieval covariance using only the model covaraince
    posterior_cov_perChannel_above_apriori = model_cov - (model_cov * Jacobian(nn,:)')*(model_cov * Jacobian(nn,:)')' /...
        (1 + (model_cov * Jacobian(nn,:)')' * Jacobian(nn,:)');

    dH = [];

    % The change in information content from the model a priori for each channel 
    dH = 1/2 * (log2(model_cov) - log2(posterior_cov_perChannel_above_apriori));

    % Let's grab just the main diagonal components and take the square root
    H_above_aPriori(:, nn) = sqrt(diag(dH));

end


% Plot the gain in information from each channel for each parameter
figure; plot(mean(GN_inputs.RT.wavelengths2run, 2), H_above_aPriori, '.-', 'MarkerSize', 20, 'LineWidth', 1.5)
grid on; grid minor
xlabel('Spectral Channel')
ylabel('Change in Bits')
title('Information change compared to the a priori')
legend('r-top', 'r-bot', 'tau_c', 'cwv', 'Location', 'best')



%% What if I transform the a priori covariance matrix so that it is a unit matrix?

%% Transform the state and measurement space!

x_new = model_cov^(-0.5) * current_guess;

y_new = measurement_cov^(-0.5) * measurements;

% now let's calculate the jacobian? 