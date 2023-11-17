
function [GN_output, GN_inputs] = calc_retrieval_gauss_newton_4modis(GN_inputs,modis,modisInputs, pixels2use)

% ----- unpack inputs -----

model_apriori = GN_inputs.model.apriori'; % a priori expected values for the model parameters
model_cov = GN_inputs.model.covariance; % model parameter covariance matrix
measurement_cov = GN_inputs.measurement.covariance; % measurement covaraince matrix
initialGuess = GN_inputs.model.initialGuess';      % Initial guess to start the Gauss-Newton iteration

% Retrieve the convergence limit
convergence_limit = GN_inputs.convergence_limit;

% retrieve the percent limit change between successive iterations
percent_change_limit = GN_inputs.percent_change_limit;

% Create the measurement vectors for each pixel!
% Each column is associated with a specific pixel, and each row represents
% the reflectance measurement at a specific modis band
measurements = create_measurement_vector(modis,GN_inputs, pixels2use); % each column represents one pixel

% % The data from 11-11-2008 at 18:50 measured erroneous values in the 1.6
% % micron channel. If using this data, lets ignore this measurement
% if strcmp(modisInputs.modisDataFolder(96:end), '/2008_11_11_1850/')==true
%
%     % get rid of the 6th MODIS channel
%     measurements(6,:) = [];
%
% else
%
% end

% -----------------------------------------------------------------------
% --------------------- PLOT JACOBIAN BAR PLOT?!?! ----------------------
% -----------------------------------------------------------------------

jacobian_barPlot_flag = false;



% -----------------------------------------------------------------------
% --------------- Define the spectral response function -----------------
% -----------------------------------------------------------------------

% using the input 'bands2run', this defines the spectral bands that will be
% written into INP files.

% check to see if the MODIS instrument is aboard Terra or Aqua
if strcmp(modisInputs.L1B_filename(1:3), 'MOD')==true
    % Then read in the spectral response functions for the terra instrument
    GN_inputs.spec_response = modis_terra_specResponse_func(GN_inputs.bands2use, GN_inputs.RT.sourceFile_resolution);

elseif strcmp(modisInputs.L1B_filename(1:3), 'MYD')==true
    % Then read in the spectral response functions for the Aqua instrument
    GN_inputs.spec_response = modis_aqua_specResponse_func(GN_inputs.bands2use, GN_inputs.RT.sourceFile_resolution);
end


% ------------------------------------------------------------------------





num_iterations = GN_inputs.GN_iterations; % number of iterations to preform
num_parameters = GN_inputs.num_model_parameters; % number of parameters to solve for

% The number of pixels to solve for using the Gauss-Newton method comes
% from the Gauss-Newton input structure
num_pixels = GN_inputs.numPixels2Calculate;

% ----- define number of spectral bands to use -----

num_bands = length(GN_inputs.bands2use);



% --- Create iterative Gauss-Newton Solver ----
% the retrieval matrix should be composed from a cell arary
% this is necessary because different pixels converge with a differing
% number of iterations
retrieval = cell(1, num_pixels);
residual = cell(1, num_pixels);
rms_residual = cell(1, num_pixels);
diff_guess_prior = cell(1, num_pixels);
jacobian_diff_guess_prior = cell(1, num_pixels);
posterior_cov = zeros(num_parameters,num_parameters,num_pixels); % my posterior covariance matrix

for pp = 1:num_pixels
    

    % define a matrix of zeros for the retrival matrix
    retrieval{pp} = zeros(num_parameters,num_iterations+1); % we include the starting point, which is outside the number of iterations
    
    % define the residual matrix
    residual{pp} = zeros(num_bands,num_iterations); % we include the starting point, which is outside the number of iterations
    
    % define the rms residual vector between the true measurements and the
    % measurement estimate
    rms_residual{pp} = zeros(1, num_iterations);      % RMS of the residual across all bands

    % create a matrix for the difference between the guess and the prior
    diff_guess_prior{pp} = zeros(num_parameters,num_iterations); % we include the starting point, which is outside the number of iterations

    % create a matrix for the product of the Jacobian and the difference
    % between the state vector guess and the state vector prior
    jacobian_diff_guess_prior{pp} = zeros(num_bands,num_iterations);      % this is an expression that multiplies the Jacobian with the difference between the current iteration and the a priori



    % --- define an initial guess ----

    % we've set the effective radius at the top and bottom to be the same
    % value. This is our initial guess. By setting this to be the a priori
    % guess we can compute the inforamtion gain between the bi-spectral
    % approach and the hyperspectal approach

    % define the initial guess
    % Here we define it to be the same re retireved for r_top and r_bottom
    retrieval{pp}(:,1) = initialGuess(:,pp);

    % -----------------------------------------------
    % ----- USING King and Vaughn (2012) Method -----
    % -----------------------------------------------
    % we set the a priori guess to be our initial guess. That way, any
    % change in our posterior tells us the information gained over the
    % bi-spectral method


    % ---- define which pixel to use -----
    pixel_row = pixels2use.res1km.row(pp);
    pixel_col = pixels2use.res1km.col(pp);



    for ii = 1:num_iterations

        disp(['(iteration,pixel): (',num2str(ii),',',num2str(pp),')']) % this tells the user where we are in the process

        % at each iteration I need to compute the forward model at my current
        % state vector estimate


        current_guess = retrieval{pp}(:,ii);

        % we compute the forward model at our previous estimate of the state vector
        % Therefore, we ask, 'what is the reflectance of a cloud with our
        % current state vector guess?'
        measurement_estimate = compute_forward_model_4modis(modis,current_guess,GN_inputs,pixel_row,pixel_col,modisInputs, pp)';

        Jacobian = compute_jacobian_4modis(modis,current_guess,measurement_estimate,GN_inputs,modisInputs, pixel_row,pixel_col, pp, jacobian_barPlot_flag);




        residual{pp}(:,ii) = measurements(:,pp) - measurement_estimate;
        diff_guess_prior{pp}(:,ii) = current_guess - model_apriori(:,pp);
        jacobian_diff_guess_prior{pp}(:,ii) = Jacobian*diff_guess_prior{pp}(:,ii);

        % -----------------------------------------------------------------
        % -----------------------------------------------------------------
        % new_guess using the previous iteration
        %new_guess = current_guess + (model_cov(:,:,pp)^(-1) + jacobian' * measurement_cov^(-1) *jacobian)^(-1) * (jacobian' *  measurement_cov(:,:,pp)^(-1) * residual(:,ii,pp) - model_cov(:,:,pp)^(-1) *diff_guess_prior(:,ii,pp));

        % new_guess using the model prior mean value
        new_guess = model_apriori(:,pp) + model_cov(:,:,pp) * Jacobian' * (Jacobian * model_cov(:,:,pp) * Jacobian' + measurement_cov(:,:,pp))^(-1) * (residual{pp}(:,ii) + jacobian_diff_guess_prior{pp}(:,ii));
        % -----------------------------------------------------------------
        % -----------------------------------------------------------------


        % If the new guess is outside the bounds of the pre-computed mie
        % table, then we must reset the value.
        if new_guess(1)>25
            disp([newline,'r_top = ',num2str(new_guess(1)),'. Set to 20 \mum'])
            new_guess(1) = 20; % microns - this may just bump back up to 60, but maybe not. The model prior should help with that
        elseif new_guess(1)<3.5
            disp([newline,'r_top = ',num2str(new_guess(1)),'. Set to 3.5 \mum'])
            new_guess(1) = 3.5; % microns
        end

        if new_guess(2)>25
            disp([newline,'r_bottom = ',num2str(new_guess(2)),'. Set to 20 \mum'])
            new_guess(2) = 20; % microns - this may just bump back up to 60, but maybe not. The model prior should help with that
        elseif new_guess(2)<3.5
            disp([newline,'r_bottom = ',num2str(new_guess(2)),'. Set to 3.5 \mum'])
            new_guess(2) = 3.5; % microns
        end


        % store the latest guess
            retrieval{pp}(:,ii+1) = new_guess;

        % If the residual is below a certain threshold as defined in the
        % GN_inputs strucute, break the for loop. We've converged
        rms_residual{pp}(ii) = sqrt(sum(residual{pp}(:,ii).^2)/size(residual{pp},1));

        if rms_residual{pp}(ii)<convergence_limit

            disp([newline, 'Convergence reached in ', num2str(ii),' iterations.', newline,...
                'RMS = ', num2str(rms_residual{pp}(ii))])

            % Clear the rest of the zeros that are place holders for later
            % iterations
            retrieval{pp}(:,ii+2:end) = [];
            rms_residual{pp}(ii+1:end) = [];
            residual{pp}(:,ii+1:end) = [];
            diff_guess_prior{pp}(:,ii+1:end) = [];
            jacobian_diff_guess_prior{pp}(:,ii+1:end) = [];

            break

        end

        % if the rms residual starts to increase, break the loop
        if ii>1 && ii<5

            if rms_residual{pp}(ii)>rms_residual{pp}(ii-1)

                disp([newline, 'RMS residual started to increase. Lowest value was: ',...
                    'RMS = ', num2str(rms_residual{pp}(ii-1)), newline])

                % Clear the rest of the zeros that are place holders for later
                % iterations
                retrieval{pp}(:,ii+2:end) = [];
                rms_residual{pp}(ii+1:end) = [];
                residual{pp}(:,ii+1:end) = [];
                diff_guess_prior{pp}(:,ii+1:end) = [];
                jacobian_diff_guess_prior{pp}(:,ii+1:end) = [];

                break

            end

        end


        % if the rms residual changes by less than 3% of the previous iteration, break the loop.
        % You're not going to do any better
        if ii>1 && ii<5

            if abs(rms_residual{pp}(ii) - rms_residual{pp}(ii-1))/rms_residual{pp}(ii-1)<percent_change_limit

                disp([newline, 'RMS residual has plataued. The current value differs from the previous value by less than 1 percent', newline,...
                    'Lowest value was: ','RMS = ', num2str(rms_residual{pp}(ii))])

                % Clear the rest of the zeros that are place holders for later
                % iterations
                retrieval{pp}(:,ii+2:end) = [];
                rms_residual{pp}(ii+1:end) = [];
                residual{pp}(:,ii+1:end) = [];
                diff_guess_prior{pp}(:,ii+1:end) = [];
                jacobian_diff_guess_prior{pp}(:,ii+1:end) = [];

                break

            end

        end



    end

        
        
    % ----------------- COMPUTE THE POSTERIOR COVARIANCE ------------------
    % once convergence has occured, we can compute the posterior covariance
    % matrix
    % First compute the latest measurement estimate
    measurement_estimate = compute_forward_model_4modis(modis, retrieval{pp}(:,end), GN_inputs, pixel_row, pixel_col, modisInputs, pp)';
    % we need to compute the jacobian using the solution state
    Jacobian = compute_jacobian_4modis(modis, retrieval{pp}(:,end), measurement_estimate, GN_inputs, modisInputs, pixel_row, pixel_col, pp, jacobian_barPlot_flag);

    posterior_cov(:,:,pp) = (Jacobian' * measurement_cov(:,:,pp)^(-1) * Jacobian + model_cov(:,:,pp)^(-1))^(-1);


    
    % ---------------- COMPUTE LIQUID WATER PATH ------------------
    % Compute the retireved Liquid water path with the final profile

    % define the altitude vector
    z = linspace(GN_inputs.RT.cloudTop_height(pp) - GN_inputs.RT.cloudDepth,...
        GN_inputs.RT.cloudTop_height(pp), GN_inputs.RT.cloud_layers)*1e3;        % m - altitude above ground vector
    
    
    
    % --------------- assuming geometric optics limit ------------------
    % For now I'm assuming that the extinction efficiency is equal to 2,
    % which is a good approximation the scattering particle is larger than
    % the incident wavelength (r > wl). I can make this more exact by 
    % computing Qe for each r value in my profile later on

    density_liquid_water = 10^6;                % g/m^3
    re_profile = create_droplet_profile2([retrieval{pp}(1,end), retrieval{pp}(2,end)],...
                    z, 'altitude', GN_inputs.model.profile.type);                               % microns          

    GN_output.LWP(pp) = 2/3 * density_liquid_water * retrieval{pp}(3,end) * trapz(z, (re_profile*1e-6).^3)/...
                                                                          trapz(z, (re_profile*1e-6).^2);           %g/m^2
    % -------------------------------------------------------------------


    
    % -------------------------------------------------------------------
    % Get rid of the Nan values and create droplet profiles with the
    % retrieval
    idx_nans = find(isnan(retrieval{pp}(3,:)));
    
    if isempty(idx_nans)~=true
       
        GN_output.tau_vector(:, pp) = linspace(0, retrieval{pp}(3,idx_nans(1)-1), 100);
        
        GN_output.re_profile(:, pp) = create_droplet_profile2([retrieval{pp}(1,idx_nans(1)-1), retrieval{pp}(2,idx_nans(1)-1)],...
                                    GN_output.tau_vector(:, pp), 'optical_depth', GN_inputs.model.profile.type);

    else
        
        GN_output.tau_vector(:, pp) = linspace(0, retrieval{pp}(3, end), 100);
        
        GN_output.re_profile(:, pp) = create_droplet_profile2([retrieval{pp}(1, end), retrieval{pp}(2, end)],...
                                    GN_output.tau_vector(:, pp), 'optical_depth', GN_inputs.model.profile.type);


    end
    % -------------------------------------------------------------------



end


GN_output.retrieval = retrieval;
GN_output.residual = residual;
GN_output.rms_residual = rms_residual;
GN_output.diff_guess_prior = diff_guess_prior;
GN_output.jacobian_diff_guess_prior = jacobian_diff_guess_prior;
GN_output.posterior_cov = posterior_cov;






end















