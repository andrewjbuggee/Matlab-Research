


function [retrieval_output, GN_inputs] = calc_retrieval_gauss_newton_4EMIT_top_bottom(GN_inputs, emit, spec_response, folder_paths)


% ----- unpack inputs -----

model_apriori = GN_inputs.model.apriori'; % a priori expected values for the model parameters
model_cov = GN_inputs.model.covariance; % model parameter covariance matrix
measurement_cov = GN_inputs.measurement.covariance; % measurement covaraince matrix
initialGuess = GN_inputs.model.initialGuess';      % Initial guess to start the Gauss-Newton iteration

% Retrieve the convergence limit
absolute_convergence_limit = GN_inputs.convergence_limit;

% retrieve the percent limit change between successive iterations
percent_change_limit = GN_inputs.percent_change_limit;


% Create the measurement vectors for each pixel!
% Each column is associated with a specific pixel, and each row represents
% the reflectance measurement at a specific modis band
measurements = emit.reflectance.value(GN_inputs.bands2run,:); % each column represents one pixel



% -----------------------------------------------------------------------
% --------------------- PLOT JACOBIAN BAR PLOT?!?! ----------------------
jacobian_barPlot_flag = false;
% -----------------------------------------------------------------------


num_iterations = GN_inputs.GN_iterations; % number of iterations to preform
num_parameters = GN_inputs.num_model_parameters; % number of parameters to solve for


% ----- define number of spectral bands to use -----
% If a measurement vector has a nan value, ignore this spectral channel
num_bands = length(measurements);



% --- Create iterative Gauss-Newton Solver ----
% the retrieval matrix should be composed from a cell arary
% this is necessary because different pixels converge with a differing
% number of iterations
% define a matrix of zeros for the retrival matrix
retrieval = zeros(num_parameters,num_iterations+1); % we include the starting point, which is outside the number of iterations

% define the residual matrix
residual = zeros(num_bands, num_iterations); % we include the starting point, which is outside the number of iterationsrss_residual = cell(1, 1);

% define the rss residual vector between the true measurements and the
% measurement estimate
rss_residual = zeros(1, num_iterations);      % Root-Sum-Square of the residual across all bands

% create a matrix for the difference between the guess and the prior
diff_guess_prior = zeros(num_parameters,num_iterations); % we include the starting point, which is outside the number of iterations

% create a matrix for the product of the Jacobian and the difference
% between the state vector guess and the state vector prior
jacobian_diff_guess_prior = zeros(num_bands,num_iterations);      % this is an expression that multiplies the Jacobian with the difference between the current iteration and the a priori







% --- define an initial guess ----

% we've set the effective radius at the top and bottom to be the same
% value. This is our initial guess. By setting this to be the a priori
% guess we can compute the inforamtion gain between the bi-spectral
% approach and the hyperspectal approach

% define the initial guess
retrieval(:,1) = initialGuess;

% -----------------------------------------------
% ----- USING King and Vaughn (2012) Method -----
% -----------------------------------------------
% we set the a priori guess to be our initial guess. That way, any
% change in our posterior tells us the information gained over the
% bi-spectral method





for ii = 1:num_iterations

    disp(['iteration: ',num2str(ii)]) % this tells the user where we are in the process

    % at each iteration I need to compute the forward model at my current
    % state vector estimate


    current_guess = retrieval(:,ii);


    if ii==1

        % Compute the measurement estimate for the first time.
        % we compute the forward model at our previous estimate of the state vector
        % Therefore, we ask, 'what is the reflectance of a cloud with our
        % current state vector guess?'
        measurement_estimate = compute_forward_model_4EMIT_top_bottom(current_guess, GN_inputs, spec_response.value, folder_paths);

        % compute residual, rms residual, the difference between the
        % iterate and the prior, and the product of the jacobian with
        % the difference between the current guess and the prior
        residual(:,ii) = measurements - measurement_estimate;
        rss_residual(ii) = sqrt(sum(residual(:,ii).^2));

    else

        % We've already calculated the measurement estimate!
        measurement_estimate = new_measurement_estimate;

    end


    % compute the jacobian
    Jacobian = compute_jacobian_4EMIT_top_bottom(emit,current_guess,measurement_estimate,GN_inputs,...
        pixels2use, pp, jacobian_barPlot_flag);


    diff_guess_prior(:,ii) = current_guess - model_apriori;
    jacobian_diff_guess_prior(:,ii) = Jacobian*diff_guess_prior(:,ii);


    % -----------------------------------------------------------------
    % -------------- Compute the new state vector ---------------------
    % -----------------------------------------------------------------

    % ----- using constraint -----
    % new guess using the modified bound-constraint algorithm (Docicu
    % et al 2003)
    % compute the Gauss-Newton direction for each retrevial variable
    new_direction = (model_cov(:,:,pp)^(-1) + Jacobian' * measurement_cov(:,:,pp)^(-1) *Jacobian)^(-1) *...
        (Jacobian' *  measurement_cov(:,:,pp)^(-1) * residual(:,ii) - model_cov(:,:,pp)^(-1) * diff_guess_prior(:,ii));

    % find the maximum non-negative value, a, that satisfies the
    % following: l< current_guess + new_direction <u
    % where each variable is bounded: l<x1<u
    % we want to compute the maximum non-negative feasible step within
    % our bounds
    a = linspace(0, 20, 2000);
    constrained_guesses = current_guess + new_direction*a;

    % let's find the new guesses that satisfy the following
    % constraints: r_bot< r_top + new_direction <inf  and
    % 0< r_bot + new_direction <r_top
    % the first row is r_top. This has to be greater than r_bot which
    % is the value of the second row.
    % find the maximum a where this is satisfied

    % --- the max re value depends on the parameterization used ---
    if strcmp(GN_inputs.RT.parameterization_str, 'mie')==true
        % mie parameterization is valid for effective radii between
        % [1,25] microns
        [max_a, ~] = max(a(constrained_guesses(1,:)>=constrained_guesses(2,:) & ...
            constrained_guesses(1,:) <= 25 & ...
            constrained_guesses(1,:) >= 1 & ...
            constrained_guesses(2,:) >= 1   & ...
            constrained_guesses(3,:) > 0));

    elseif strcmp(GN_inputs.RT.parameterization_str, 'hu')==true
        % Hu and Stamnes  parameterization is valid for effective
        % radii between [3.5, 60] microns
        [max_a, ~] = max(a(constrained_guesses(1,:)>=constrained_guesses(2,:) & ...
            constrained_guesses(1,:) <= 60 & ...
            constrained_guesses(1,:) >= 3.5 & ...
            constrained_guesses(2,:) >= 3.5   & ...
            constrained_guesses(3,:) > 0));
    end


    % if the maximum value of a is 0, then there is no solution space
    % with the current Gauss-Newton direction that will result in r_top
    % being larger than r_bot. If this occurs on the initial iteration
    % when the a priori value for the radius at cloud top is equal to
    % the radius at cloud bottom, we will reset the guess.
    if max_a==0 && current_guess(1)==current_guess(2)

        % In this case, the new guass-newton direction is causing the
        % radius at cloud top to be larger than that at cloud bottom.
        % Reset the initial guess to be the TBLUT value at cloud top,
        % and 72% that value at cloud bottom. This is based on the
        % median results from analyzing non-precipitating clouds from
        % the VOCALS-REx flight campaing. The median vertical profile
        % shows the radius at cloud bottom to be 72% the value at cloud
        % top.

        disp([newline, 'The Gauss-Newton direction causes r_top<r_bot.',newline, ...
            'Trying a different inital guess...', newline])

        new_guess = [current_guess(1), 0.7273*current_guess(2), current_guess(3)];

        disp([newline, 'Initial guess was: r_t = ', num2str(current_guess(1)),...
            '  r_b = ', num2str(current_guess(2)), '  Tau_c = ', num2str(current_guess(3)), newline])
        disp(['New guess is : r_t = ', num2str(new_guess(1)),...
            '  r_b = ', num2str(new_guess(2)), '  Tau_c = ', num2str(new_guess(3)), newline])

        % Use the new guess to compute the rms residual, which is used
        % to detmerine convergence
        new_measurement_estimate = compute_forward_model_4EMIT_top_bottom(emit, new_guess, GN_inputs, pixels2use, pp)';
        residual(:,ii+1) = measurements - new_measurement_estimate;
        rms_residual(ii+1) = sqrt(mean(residual(:,ii+1).^2));


    elseif max_a==0 && ii>1

        % in this case, the only direction the algorithm finds is one
        % that forces the profile to be non-adiabatic. Let's break the
        % loop at this point

        disp([newline, 'The only solution forces r_top<r_bot.',newline, ...
            'Breaking the loop...', newline])


        % Clear the rest of the zeros that are place holders for later
        % iterations
        retrieval(:,ii+1:end) = [];
        rms_residual(ii+1:end) = [];
        residual(:,ii+1:end) = [];
        diff_guess_prior(:,ii+1:end) = [];
        jacobian_diff_guess_prior(:,ii+1:end) = [];

        break


    else

        % We want to make sure the new step is within the feasible
        % range, not at the boundaries. So we only accept a values that
        % are less than the max a value.
        percent_of_maxA = 0.95;
        % Set the a vector to values between 0 and some fraction of the max a
        a = linspace(0, percent_of_maxA * max_a, 10);
        % recompute the constrained guesses
        constrained_guesses = current_guess + new_direction*a;



        % We need to compute the measurement estimate of the constrained
        % solution. We want to determine a value for which the rms
        % difference between the new estimated measurement vector and the true
        % measurement vector is less than the previous RMS difference
        constrained_measurement_estimate = zeros(num_bands, length(a));
        for mm = 1:length(a)
            constrained_measurement_estimate(:,mm) = compute_forward_model_4EMIT_top_bottom(emit, constrained_guesses(:,mm),...
                GN_inputs, pixels2use, pp)';
        end

        % compute the rms_residual for the constrained state vector
        rms_residual_constrained = sqrt(mean((constrained_measurement_estimate - repmat(measurements, 1, length(a))).^2, 1));
        % find the smallest rms residual that is less than the previus
        % itereates rms residual
        [min_val_lessThanPrevious, ~] = min(rms_residual_constrained(rms_residual_constrained < rms_residual(ii)));

        % Check to see if all rms_residuals are greater than the
        % previous iterate
        if isempty(min_val_lessThanPrevious)

            % If no rms_residual is less than the previous iterate,
            % find the minimum and move foward. The algorithm will flag
            % this as find the minimum rms_residual
            [~, min_residual_idx] = min(rms_residual_constrained);


        else

            % find the index for the smallest rms that is less than the
            % previous iterate rms residual
            min_residual_idx = find(min_val_lessThanPrevious == rms_residual_constrained);

        end


        % Select the step length by choosing the a value with the minimumum
        % residual
        new_measurement_estimate = constrained_measurement_estimate(:, min_residual_idx);
        residual(:,ii+1) = measurements - new_measurement_estimate;
        rms_residual(ii+1) = sqrt(mean(residual(:,ii+1).^2));
        new_guess = constrained_guesses(:, min_residual_idx);



    end




    % ----- new_guess using the previous iteration -----
    %new_guess = current_guess + (model_cov(:,:,pp)^(-1) + jacobian' * measurement_cov^(-1) *jacobian)^(-1) * (jacobian' *  measurement_cov(:,:,pp)^(-1) * residual(:,ii,pp) - model_cov(:,:,pp)^(-1) *diff_guess_prior(:,ii,pp));


    % ----- new_guess using the model prior mean value -----
    %new_guess = model_apriori + model_cov(:,:,pp) * Jacobian' * (Jacobian * model_cov(:,:,pp) * Jacobian' + measurement_cov(:,:,pp))^(-1) * (residual(:,ii) + jacobian_diff_guess_prior(:,ii));

    % -----------------------------------------------------------------
    % -----------------------------------------------------------------

    if strcmp(GN_inputs.RT.parameterization_str, 'mie')==true
        % If the new guess is outside the bounds of the pre-computed mie
        % table, then we must reset the value.

        if new_guess(1)>25
            disp([newline,'r_top = ',num2str(new_guess(1)),'. Set to 20 \mum'])
            new_guess(1) = 20; % microns - this may just bump back up to 60, but maybe not. The model prior should help with that
        elseif new_guess(1)<1
            disp([newline,'r_top = ',num2str(new_guess(1)),'. Set to 3.5 \mum'])
            new_guess(1) = 3.5; % microns
        end

        if new_guess(2)>25
            disp([newline,'r_bottom = ',num2str(new_guess(2)),'. Set to 20 \mum'])
            new_guess(2) = 20; % microns - this may just bump back up to 60, but maybe not. The model prior should help with that
        elseif new_guess(2)<1
            disp([newline,'r_bottom = ',num2str(new_guess(2)),'. Set to 3.5 \mum'])
            new_guess(2) = 3.5; % microns
        end



    elseif strcmp(GN_inputs.RT.parameterization_str, 'hu')==true

        % If the new guess is outside the bounds of the Hu and
        % Stamnes parameterization, we must reset the values

        if new_guess(1)>60
            disp([newline,'r_top = ',num2str(new_guess(1)),'. Set to 20 \mum'])
            new_guess(1) = 20; % microns - this may just bump back up to 60, but maybe not. The model prior should help with that
        elseif new_guess(1)<2.5
            disp([newline,'r_top = ',num2str(new_guess(1)),'. Set to 3.5 \mum'])
            new_guess(1) = 3.5; % microns
        end

        if new_guess(2)>60
            disp([newline,'r_bottom = ',num2str(new_guess(2)),'. Set to 20 \mum'])
            new_guess(2) = 20; % microns - this may just bump back up to 60, but maybe not. The model prior should help with that
        elseif new_guess(2)<2.5
            disp([newline,'r_bottom = ',num2str(new_guess(2)),'. Set to 3.5 \mum'])
            new_guess(2) = 3.5; % microns
        end


    end



    % store the latest guess
    retrieval(:,ii+1) = new_guess;

    % -------------------------------------------------------------
    % -------------------- Convergence checks ---------------------
    % -------------------------------------------------------------
    % If the residual is below a certain threshold as defined in the
    % GN_inputs strucute, break the for loop. We've converged


    if rms_residual(ii+1)<absolute_convergence_limit

        disp([newline, 'Convergence reached in ', num2str(ii),' iterations.', newline,...
            'RMS = ', num2str(rms_residual(ii+1))])

        % Clear the rest of the zeros that are place holders for later
        % iterations
        retrieval(:,ii+2:end) = [];
        rms_residual(ii+2:end) = [];
        residual(:,ii+2:end) = [];
        diff_guess_prior(:,ii+1:end) = [];
        jacobian_diff_guess_prior(:,ii+1:end) = [];

        break

    end

    % if the rms residual starts to increase, break the loop
    if ii>1 && ii<5

        if rms_residual(ii+1)>rms_residual(ii)

            disp([newline, 'RMS residual started to increase. Lowest value was: ',...
                'RMS = ', num2str(rms_residual(ii)), newline])

            % Clear the rest of the zeros that are place holders for later
            % iterations
            retrieval(:,ii+2:end) = [];
            rms_residual(ii+2:end) = [];
            residual(:,ii+2:end) = [];
            diff_guess_prior(:,ii+1:end) = [];
            jacobian_diff_guess_prior(:,ii+1:end) = [];

            break

        end

    end


    % if the rms residual changes by less than x% of the previous iteration, break the loop.
    % You're not going to do any better
    if ii>1 && ii<5

        if abs(rms_residual(ii+1) - rms_residual(ii))/rms_residual(ii)==0

            if isempty(min_val_lessThanPrevious)==true
                disp([newline, 'There is no new direction that reduces the RMS uncertainty.', newline,...
                    'Lowest value was: ','RMS = ', num2str(rms_residual(ii+1))])

            else

                disp([newline, 'RMS residual did not change at all from the previous iteration.', newline,...
                    'Lowest value was: ','RMS = ', num2str(rms_residual(ii+1))])

            end

            % Clear the rest of the zeros that are place holders for later
            % iterations
            retrieval(:,ii+2:end) = [];
            rms_residual(ii+2:end) = [];
            residual(:,ii+2:end) = [];
            diff_guess_prior(:,ii+1:end) = [];
            jacobian_diff_guess_prior(:,ii+1:end) = [];

            break

        elseif abs(rms_residual(ii+1) - rms_residual(ii))/rms_residual(ii)<percent_change_limit

            disp([newline, 'RMS residual has plataued. The current value differs from the previous value by ',...
                num2str(abs(rms_residual(ii+1) - rms_residual(ii))/rms_residual(ii) * 100), '%', newline,...
                'Lowest value was: ','RMS = ', num2str(rms_residual(ii+1))])

            % Clear the rest of the zeros that are place holders for later
            % iterations
            retrieval(:,ii+2:end) = [];
            rms_residual(ii+2:end) = [];
            residual(:,ii+2:end) = [];
            diff_guess_prior(:,ii+1:end) = [];
            jacobian_diff_guess_prior(:,ii+1:end) = [];

            break

        end

    end



end



% ----------------- COMPUTE THE POSTERIOR COVARIANCE ------------------
% once convergence has occured, we can compute the posterior covariance
% matrix
% First compute the latest measurement estimate

% we need to compute the jacobian using the solution state
Jacobian = compute_jacobian_4EMIT_top_bottom(emit, retrieval(:,end), new_measurement_estimate, GN_inputs,...
    pixels2use, pp, jacobian_barPlot_flag);

posterior_cov(:,:,pp) = (Jacobian' * measurement_cov(:,:,pp)^(-1) *...
    Jacobian + model_cov(:,:,pp)^(-1))^(-1);



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
re_profile = create_droplet_profile2([retrieval(1,end), retrieval(2,end)],...
    z, 'altitude', GN_inputs.model.profile.type);                               % microns

retrieval_output.LWP(pp) = 2/3 * density_liquid_water * retrieval(3,end) * trapz(z, (re_profile*1e-6).^3)/...
    trapz(z, (re_profile*1e-6).^2);           %g/m^2
% -------------------------------------------------------------------



% -------------------------------------------------------------------
% Get rid of the Nan values and create droplet profiles with the
% retrieval
idx_nans = find(isnan(retrieval(3,:)));

if isempty(idx_nans)~=true

    retrieval_output.tau_vector(:, pp) = linspace(0, retrieval(3,idx_nans(1)-1), 100);

    retrieval_output.re_profile(:, pp) = create_droplet_profile2([retrieval(1,idx_nans(1)-1), retrieval(2,idx_nans(1)-1)],...
        retrieval_output.tau_vector(:, pp), 'optical_depth', GN_inputs.model.profile.type);

else

    retrieval_output.tau_vector(:, pp) = linspace(0, retrieval(3, end), 100);

    retrieval_output.re_profile(:, pp) = create_droplet_profile2([retrieval(1, end), retrieval(2, end)],...
        retrieval_output.tau_vector(:, pp), 'optical_depth', GN_inputs.model.profile.type);


end
% -------------------------------------------------------------------





retrieval_output.variables = retrieval;
retrieval_output.computed_reflectance = new_measurement_estimate;
retrieval_output.residual = residual;
retrieval_output.rms_residual = rms_residual;
retrieval_output.diff_guess_prior = diff_guess_prior;
retrieval_output.jacobian_diff_guess_prior = jacobian_diff_guess_prior;
retrieval_output.posterior_cov = posterior_cov;






end