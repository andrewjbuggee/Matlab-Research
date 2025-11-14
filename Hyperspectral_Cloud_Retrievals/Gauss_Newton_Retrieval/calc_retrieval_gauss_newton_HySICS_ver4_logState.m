% This script retrieves 4 variables: ln(r_top), ln(r_bot), ln(tau_c), and ln(acpw)


function [GN_output, GN_inputs] = calc_retrieval_gauss_newton_HySICS_ver4_logState(GN_inputs, hysics, folder_paths, print_status_updates)


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
measurements = hysics.Refl_model; % column vector of the reflectance measurements



% --------------------------------------------
% ---- set direction of change parameters ----
% --------------------------------------------

% set the maximum scalar value used to create a set of scalar values that
% will be multiplied along the direction of greatest change
a_largestVal = 3;

% if trying to have specific values, like a=1, use a step size
a_stepSize = 0.333;

% define length of initial array used to check which state vectors meet the
% defined constraints
array_length_initialConstraints = 2000;

% define the array of values between 0 and the maximum scalar value
array_length_newMax = 10;

% We want to make sure the new step is within the feasible
% range, not at the boundaries. So we only accept a values that
% are less than the max a value.
percent_of_maxA = 0.95;
% --------------------------------------------



% -----------------------------------------------------------------------
% --------------------- PLOT JACOBIAN BAR PLOT?!?! ----------------------
jacobian_barPlot_flag = false;
% -----------------------------------------------------------------------


num_iterations = GN_inputs.GN_iterations; % number of iterations to preform
num_parameters = GN_inputs.num_model_parameters; % number of parameters to solve for

% ----- define number of spectral bands to use -----
% If a measurement vector has a nan value, ignore this spectral channel
num_bands = length(measurements);


% define the spectral response function
spec_response = hysics.spec_response.value;


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
% Here we define it to be the same re retireved for r_top and r_bottom
retrieval(:,1) = initialGuess;

% -----------------------------------------------
% ----- USING King and Vaughn (2012) Method -----
% -----------------------------------------------
% we set the a priori guess to be our initial guess. That way, any
% change in our posterior tells us the information gained over the
% bi-spectral method



if print_status_updates==true

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

            % For the retrieval of r_top, r_bot, tau_c, cwv
            disp([newline, 'Estimating spectral measurements...', newline])
            measurement_estimate = compute_forward_model_HySICS_ver2(exp(current_guess), GN_inputs, spec_response, folder_paths);


            % compute residual, rss residual, the difference between the
            % iterate and the prior, and the product of the jacobian with
            % the difference between the current guess and the prior
            residual(:,ii) = measurements - measurement_estimate;
            rss_residual(ii) = sqrt(sum(residual(:,ii).^2));

        else

            % We've already calculated the measurement estimate!
            measurement_estimate = new_measurement_estimate;

        end


        % **** compute the jacobian ****
        % For the retrieval of r_top, r_bot, tau_c, cwv
        disp([newline, 'Computing the Jacobian...', newline])
        Jacobian = compute_jacobian_HySICS_ver2(exp(current_guess), measurement_estimate, GN_inputs,...
            hysics.spec_response.value, jacobian_barPlot_flag, folder_paths);



        diff_guess_prior(:,ii) = current_guess - model_apriori;
        jacobian_diff_guess_prior(:,ii) = Jacobian*diff_guess_prior(:,ii);



        % -------------- Compute the new state vector ---------------------
        % -----------------------------------------------------------------

        % ----- using constraint -----
        % new guess using the modified bound-constraint algorithm (Docicu
        % et al 2003)
        % compute the Gauss-Newton direction for each retrevial variable
        new_direction = (model_cov^(-1) + Jacobian' * measurement_cov^(-1) *Jacobian)^(-1) *...
            (Jacobian' *  measurement_cov^(-1) * residual(:,ii) - model_cov^(-1) * diff_guess_prior(:,ii));

        % fine the maximum non-negative value, a, that satisfies the
        % following: l< current_guess + a*new_direction <u
        % where the variable is bounded: l<x1<u
        % we want to compute the maximum non-negative feasible step within
        % our bounds
        % a = linspace(0, 20, 2000);
        a = linspace(0, a_largestVal, array_length_initialConstraints);
        constrained_guesses = current_guess + new_direction*a;

        % let's find the new guesses that satisfy the following
        % constraints: r_bot< r_top + new_direction <inf  and
        % 0< r_bot + new_direction <r_top
        % the first row is r_top. This has to be greater than r_bot which
        % is the value of the second row.
        % find the maximum a where this is satisfied
        [max_a, ~] = max(a(constrained_guesses(1,:)>=constrained_guesses(2,:) & ...
            constrained_guesses(1,:)<=log(25) & ...
            constrained_guesses(2,:)>0   & ...
            constrained_guesses(3,:)>0   & ...
            constrained_guesses(4,:)>0));

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
            
            % Check to see if both guesses are 3.5, the lower limit of our
            % lookup table
            if current_guess(1)==3.5 && current_guess(2)==3.5

                new_guess = [9, 5, current_guess(3), current_guess(4)];

            else

                new_guess = [current_guess(1), 0.7273*current_guess(2), current_guess(3), current_guess(4)];

            end
    
            disp([newline, 'Initial guess was: r_t = ', num2str(current_guess(1)),...
                '  r_b = ', num2str(current_guess(2)), '  Tau_c = ', num2str(current_guess(3)),...
                '  acpw = ', num2str(current_guess(4)), newline])
            disp(['New guess is : r_t = ', num2str(new_guess(1)),...
                '  r_b = ', num2str(new_guess(2)), '  Tau_c = ', num2str(new_guess(3)),...
                '  acpw = ', num2str(new_guess(4)), newline])

            % Use the new guess to compute the rss residual, which is used
            % to detmerine convergence
            disp([newline, 'Estimating spectral measurements...', newline])
            new_measurement_estimate = compute_forward_model_HySICS_ver2(new_guess, GN_inputs, spec_response, folder_paths);
            residual(:,ii+1) = measurements - new_measurement_estimate;
            rss_residual(ii+1) = sqrt(sum(residual(:,ii+1).^2));


        elseif max_a==0 && ii>1

            % in this case, the only direction the algorithm finds is one
            % that forces the profile to be non-adiabatic. Let's break the
            % loop at this point

            disp([newline, 'The only solution forces r_top<r_bot.',newline, ...
                'Breaking the loop...', newline])


            % Clear the rest of the zeros that are place holders for later
            % iterations
            retrieval(:,ii+1:end) = [];
            rss_residual(ii+1:end) = [];
            residual(:,ii+1:end) = [];
            diff_guess_prior(:,ii+1:end) = [];
            jacobian_diff_guess_prior(:,ii+1:end) = [];

            break


        else

            disp([newline, 'Computing new direction using predefined constraints...', newline])

            
            % Set the a vector to values between 0 and some fraction of the max a
            a = linspace(0, percent_of_maxA * max_a, array_length_newMax);
            if max_a>1
                % include a=1
                a = sort([a, 1]);
            end
            % a = 0:a_stepSize:max_a;
            % recompute the constrained guesses
            constrained_guesses = current_guess + new_direction*a;



            % We need to compute the measurement estimate of the constrained
            % solution. We want to determine a value for which the L2 norm of
            % the difference between the new constrained guess and the true
            % measurements is less than the previous guess and the measurements
            constrained_measurement_estimate = zeros(num_bands, length(a));

            % This loop cannot be a parfor loop because the function within the
            % loop also uses parfor! compute_forward_model_HySICS_ver2 uses a
            % for loop. You could rewrite this whole loop so ALL constrained
            % guess can run in a single parfor loop...
            for mm = 1:length(a)

                % some guesses might be out of the appropriate range for
                % the Mie Interpolation function. If so, set the
                % constrained measurement estimates to 0
                if constrained_guesses(1,mm)>0 && constrained_guesses(1,mm)<log(25) && ...
                        constrained_guesses(2,mm)>0 && constrained_guesses(2,mm)<log(25)

                    disp([newline, 'Estimating spectral measurements...', newline])
                    constrained_measurement_estimate(:,mm)= compute_forward_model_HySICS_ver2(exp(constrained_guesses(:,mm)),...
                        GN_inputs, spec_response, folder_paths);

                else

                    constrained_measurement_estimate(:,mm) = 0;
                end

            end

            % compute the rss_residual for the constrained state vector
            rss_residual_constrained = sqrt(sum((constrained_measurement_estimate - repmat(measurements, 1, length(a))).^2, 1));
            % find the smallest rss residual that is less than the previus
            % itereates rss residual
            [min_val_lessThanPrevious, ~] = min(rss_residual_constrained(rss_residual_constrained < rss_residual(ii)));

            % Check to see if all rss_residuals are greater than the
            % previous iteration
            if isempty(min_val_lessThanPrevious)

                % If no rss_residual is less than the previous iterate,
                % find the minimum and move foward. Tha algorithm will flag
                % this as find the minimum rss_residual
                [~, min_residual_idx] = min(rss_residual_constrained);


            else

                % find the index for the smallest rss that is less than the
                % previous iterate rss residual
                min_residual_idx = find(min_val_lessThanPrevious == rss_residual_constrained);

            end


            % Select the step length by choosing the a value with the minimumum
            % residual
            new_measurement_estimate = constrained_measurement_estimate(:, min_residual_idx);
            residual(:,ii+1) = measurements - new_measurement_estimate;
            rss_residual(ii+1) = sqrt(sum(residual(:,ii+1).^2));
            new_guess = constrained_guesses(:, min_residual_idx);



        end




        % ----- new_guess using the previous iteration -----
        %new_guess = current_guess + (model_cov(:,:,pp)^(-1) + jacobian' * measurement_cov^(-1) *jacobian)^(-1) * (jacobian' *  measurement_cov(:,:,pp)^(-1) * residual(:,ii,pp) - model_cov(:,:,pp)^(-1) *diff_guess_prior(:,ii,pp));


        % ----- new_guess using the model prior mean value -----
        %new_guess = model_apriori(:,pp) + model_cov(:,:,pp) * Jacobian' * (Jacobian * model_cov(:,:,pp) * Jacobian' + measurement_cov(:,:,pp))^(-1) * (residual(:,ii) + jacobian_diff_guess_prior(:,ii));

        % -----------------------------------------------------------------
        % -----------------------------------------------------------------


        % If the new guess is outside the bounds of the pre-computed mie
        % table, then we must reset the value.
        if new_guess(1)>log(25)
            disp([newline,'r_top = ',num2str(exp(new_guess(1))),'. Set to 15 \mum'])
            new_guess(1) = log(15); % microns - this may just bump back up to 60, but maybe not. The model prior should help with that
        elseif new_guess(1)<log(3.5)
            disp([newline,'r_top = ',num2str(exp(new_guess(1))),'. Set to 3.5 \mum'])
            new_guess(1) = log(3.5); % microns
        end

        if new_guess(2)>log(25)
            disp([newline,'r_bottom = ',num2str(exp(new_guess(2))),'. Set to 20 \mum'])
            new_guess(2) = log(15); % microns - this may just bump back up to 60, but maybe not. The model prior should help with that
        elseif new_guess(2)<log(3.5)
            disp([newline,'r_bottom = ',num2str(exp(new_guess(2))),'. Set to 3.5 \mum'])
            new_guess(2) = log(3.5); % microns
        end


        % store the latest guess
        retrieval(:,ii+1) = new_guess;

        % If the residual is below a certain threshold as defined in the
        % GN_inputs strucute, break the for loop. We've converged


        if rss_residual(ii+1)<convergence_limit

            disp([newline, 'Convergence reached in ', num2str(ii),' iterations.', newline,...
                'RSS Limit = ', num2str(convergence_limit),newline,...
                'RSS = ', num2str(rss_residual(ii+1))])

            % Clear the rest of the zeros that are place holders for later
            % iterations
            retrieval(:,ii+2:end) = [];
            rss_residual(ii+2:end) = [];
            residual(:,ii+2:end) = [];
            diff_guess_prior(:,ii+1:end) = [];
            jacobian_diff_guess_prior(:,ii+1:end) = [];

            break

        end

        % if the rss residual starts to increase, break the loop
        if ii>1 && ii<5

            if rss_residual(ii+1)>rss_residual(ii)

                disp([newline, 'RSS residual started to increase. Lowest value was: ',...
                    'RSS Limit = ', num2str(convergence_limit),newline,...
                    'RSS = ', num2str(rss_residual(ii)), newline])

                % Clear the rest of the zeros that are place holders for later
                % iterations
                retrieval(:,ii+2:end) = [];
                rss_residual(ii+2:end) = [];
                residual(:,ii+2:end) = [];
                diff_guess_prior(:,ii+1:end) = [];
                jacobian_diff_guess_prior(:,ii+1:end) = [];

                break

            end

        end


        % if the rss residual changes by less than 3% of the previous iteration, break the loop.
        % You're not going to do any better
        if ii>1 && ii<5

            if abs(rss_residual(ii+1) - rss_residual(ii))/rss_residual(ii)<percent_change_limit

                disp([newline, 'RSS residual has plataued. The current value differs from the previous value by less than ',...
                    num2str(100*percent_change_limit), '%', newline,...
                    'RSS Limit = ', num2str(convergence_limit),newline,...
                    'Lowest value was: ','RSS = ', num2str(rss_residual(ii+1))])

                % Clear the rest of the zeros that are place holders for later
                % iterations
                retrieval(:,ii+2:end) = [];
                rss_residual(ii+2:end) = [];
                residual(:,ii+2:end) = [];
                diff_guess_prior(:,ii+1:end) = [];
                jacobian_diff_guess_prior(:,ii+1:end) = [];

                break

            end

        end



    end








else
    
    % ---------------------------------------------
    % Dont print any messages to the command window
    % ---------------------------------------------


    for ii = 1:num_iterations


        % at each iteration I need to compute the forward model at my current
        % state vector estimate


        current_guess = retrieval(:,ii);


        if ii==1

            % Compute the measurement estimate for the first time.
            % we compute the forward model at our previous estimate of the state vector
            % Therefore, we ask, 'what is the reflectance of a cloud with our
            % current state vector guess?'

            % For the retrieval of r_top, r_bot, tau_c, cwv
            measurement_estimate = compute_forward_model_HySICS_ver2(current_guess, GN_inputs, spec_response, folder_paths);


            % compute residual, rss residual, the difference between the
            % iterate and the prior, and the product of the jacobian with
            % the difference between the current guess and the prior
            residual(:,ii) = measurements - measurement_estimate;
            rss_residual(ii) = sqrt(sum(residual(:,ii).^2));

        else

            % We've already calculated the measurement estimate!
            measurement_estimate = new_measurement_estimate;

        end


        % **** compute the jacobian ****
        % For the retrieval of r_top, r_bot, tau_c, cwv
        Jacobian = compute_jacobian_HySICS_ver2(current_guess, measurement_estimate, GN_inputs,...
            hysics.spec_response.value, jacobian_barPlot_flag, folder_paths);



        diff_guess_prior(:,ii) = current_guess - model_apriori;
        jacobian_diff_guess_prior(:,ii) = Jacobian*diff_guess_prior(:,ii);



        % -------------- Compute the new state vector ---------------------
        % -----------------------------------------------------------------

        % ----- using constraint -----
        % new guess using the modified bound-constraint algorithm (Docicu
        % et al 2003)
        % compute the Gauss-Newton direction for each retrevial variable
        new_direction = (model_cov^(-1) + Jacobian' * measurement_cov^(-1) *Jacobian)^(-1) *...
            (Jacobian' *  measurement_cov^(-1) * residual(:,ii) - model_cov^(-1) * diff_guess_prior(:,ii));

        % fine the maximum non-negative value, a, that satisfies the
        % following: l< current_guess + a*new_direction <u
        % where the variable is bounded: l<x1<u
        % we want to compute the maximum non-negative feasible step within
        % our bounds
        a = linspace(0, a_largestVal, array_length_initialConstraints);
        constrained_guesses = current_guess + new_direction*a;

        % let's find the new guesses that satisfy the following
        % constraints: r_bot< r_top + new_direction <inf  and
        % 0< r_bot + new_direction <r_top
        % the first row is r_top. This has to be greater than r_bot which
        % is the value of the second row.
        % find the maximum a where this is satisfied
        [max_a, ~] = max(a(constrained_guesses(1,:)>=constrained_guesses(2,:) & ...
            constrained_guesses(1,:)<=30 & ...
            constrained_guesses(2,:)>0   & ...
            constrained_guesses(3,:)>0   & ...
            constrained_guesses(4,:)>0));

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



            
            % Check to see if both guesses are 3.5, the lower limit of our
            % lookup table
            if current_guess(1)==3.5 && current_guess(2)==3.5

                new_guess = [9, 5, current_guess(3), current_guess(4)];

            else

                new_guess = [current_guess(1), 0.7273*current_guess(2), current_guess(3), current_guess(4)];

            end



            % Use the new guess to compute the rss residual, which is used
            % to detmerine convergence
            new_measurement_estimate = compute_forward_model_HySICS_ver2(new_guess, GN_inputs, spec_response, folder_paths);
            residual(:,ii+1) = measurements - new_measurement_estimate;
            rss_residual(ii+1) = sqrt(sum(residual(:,ii+1).^2));


        elseif max_a==0 && ii>1

            % in this case, the only direction the algorithm finds is one
            % that forces the profile to be non-adiabatic. Let's break the
            % loop at this point




            % Clear the rest of the zeros that are place holders for later
            % iterations
            retrieval(:,ii+1:end) = [];
            rss_residual(ii+1:end) = [];
            residual(:,ii+1:end) = [];
            diff_guess_prior(:,ii+1:end) = [];
            jacobian_diff_guess_prior(:,ii+1:end) = [];

            break


        else



            % Set the a vector to values between 0 and some fraction of the max a
            a = linspace(0, percent_of_maxA * max_a, array_length_newMax);
            if max_a>1
                % include a=1
                a = sort([a, 1]);
            end
            % a = 0:a_stepSize:max_a;
            % recompute the constrained guesses
            constrained_guesses = current_guess + new_direction*a;



            % We need to compute the measurement estimate of the constrained
            % solution. We want to determine a value for which the L2 norm of
            % the difference between the new constrained guess and the true
            % measurements is less than the previous guess and the measurements
            constrained_measurement_estimate = zeros(num_bands, length(a));

            % This loop cannot be a parfor loop because the function within the
            % loop also uses parfor! compute_forward_model_HySICS_ver2 uses a
            % for loop. You could rewrite this whole loop so ALL constrained
            % guess can run in a single parfor loop...
            for mm = 1:length(a)

                % some guesses might be out of the appropriate range for
                % the Mie Interpolation function. If so, set the
                % constrained measurement estimates to 0
                if constrained_guesses(1,mm)>1 && constrained_guesses(1,mm)<25 && ...
                        constrained_guesses(2,mm)>1 && constrained_guesses(2,mm)<25

                    constrained_measurement_estimate(:,mm)= compute_forward_model_HySICS_ver2(constrained_guesses(:,mm),...
                        GN_inputs, spec_response, folder_paths);

                else

                    constrained_measurement_estimate(:,mm) = 0;
                end

            end

            % compute the rss_residual for the constrained state vector
            rss_residual_constrained = sqrt(sum((constrained_measurement_estimate - repmat(measurements, 1, length(a))).^2, 1));
            % find the smallest rss residual that is less than the previus
            % itereates rss residual
            [min_val_lessThanPrevious, ~] = min(rss_residual_constrained(rss_residual_constrained < rss_residual(ii)));

            % Check to see if all rss_residuals are greater than the
            % previous iterate
            if isempty(min_val_lessThanPrevious)

                % If no rss_residual is less than the previous iterate,
                % find the minimum and move foward. Tha algorithm will flag
                % this as find the minimum rss_residual
                [~, min_residual_idx] = min(rss_residual_constrained);


            else

                % find the index for the smallest rss that is less than the
                % previous iterate rss residual
                min_residual_idx = find(min_val_lessThanPrevious == rss_residual_constrained);

            end


            % Select the step length by choosing the a value with the minimumum
            % residual
            new_measurement_estimate = constrained_measurement_estimate(:, min_residual_idx);
            residual(:,ii+1) = measurements - new_measurement_estimate;
            rss_residual(ii+1) = sqrt(sum(residual(:,ii+1).^2));
            new_guess = constrained_guesses(:, min_residual_idx);



        end




        % ----- new_guess using the previous iteration -----
        %new_guess = current_guess + (model_cov(:,:,pp)^(-1) + jacobian' * measurement_cov^(-1) *jacobian)^(-1) * (jacobian' *  measurement_cov(:,:,pp)^(-1) * residual(:,ii,pp) - model_cov(:,:,pp)^(-1) *diff_guess_prior(:,ii,pp));


        % ----- new_guess using the model prior mean value -----
        %new_guess = model_apriori(:,pp) + model_cov(:,:,pp) * Jacobian' * (Jacobian * model_cov(:,:,pp) * Jacobian' + measurement_cov(:,:,pp))^(-1) * (residual(:,ii) + jacobian_diff_guess_prior(:,ii));

        % -----------------------------------------------------------------
        % -----------------------------------------------------------------


        % If the new guess is outside the bounds of the pre-computed mie
        % table, then we must reset the value.
        if new_guess(1)>25
            new_guess(1) = 20; % microns - this may just bump back up to 60, but maybe not. The model prior should help with that
        elseif new_guess(1)<3.5
            new_guess(1) = 3.5; % microns
        end

        if new_guess(2)>25
            new_guess(2) = 20; % microns - this may just bump back up to 60, but maybe not. The model prior should help with that
        elseif new_guess(2)<3.5
            new_guess(2) = 3.5; % microns
        end


        % store the latest guess
        retrieval(:,ii+1) = new_guess;

        % If the residual is below a certain threshold as defined in the
        % GN_inputs strucute, break the for loop. We've converged


        if rss_residual(ii+1)<convergence_limit



            % Clear the rest of the zeros that are place holders for later
            % iterations
            retrieval(:,ii+2:end) = [];
            rss_residual(ii+2:end) = [];
            residual(:,ii+2:end) = [];
            diff_guess_prior(:,ii+1:end) = [];
            jacobian_diff_guess_prior(:,ii+1:end) = [];

            break

        end

        % if the rss residual starts to increase, break the loop
        if ii>1 && ii<5

            if rss_residual(ii+1)>rss_residual(ii)


                % Clear the rest of the zeros that are place holders for later
                % iterations
                retrieval(:,ii+2:end) = [];
                rss_residual(ii+2:end) = [];
                residual(:,ii+2:end) = [];
                diff_guess_prior(:,ii+1:end) = [];
                jacobian_diff_guess_prior(:,ii+1:end) = [];

                break

            end

        end


        % if the rss residual changes by less than 3% of the previous iteration, break the loop.
        % You're not going to do any better
        if ii>1 && ii<5

            if abs(rss_residual(ii+1) - rss_residual(ii))/rss_residual(ii)<percent_change_limit


                % Clear the rest of the zeros that are place holders for later
                % iterations
                retrieval(:,ii+2:end) = [];
                rss_residual(ii+2:end) = [];
                residual(:,ii+2:end) = [];
                diff_guess_prior(:,ii+1:end) = [];
                jacobian_diff_guess_prior(:,ii+1:end) = [];

                break

            end

        end



    end


end







% transform the retrieval back to linear space
retrieval = exp(retrieval);








% ----------------- COMPUTE THE POSTERIOR COVARIANCE ------------------
% once convergence has occured, we can compute the posterior covariance
% matrix

% we need to compute the jacobian using the solution state
Jacobian = compute_jacobian_HySICS_ver2(exp(retrieval(:,end)), new_measurement_estimate, GN_inputs,...
    hysics.spec_response.value, jacobian_barPlot_flag, folder_paths);

posterior_cov = ((Jacobian' * measurement_cov^(-1) * Jacobian) + model_cov^(-1))^(-1);






% ---------------- COMPUTE LIQUID WATER PATH ------------------
% Compute the retireved Liquid water path with the final profile


% --------------- assuming geometric optics limit ------------------
% For now I'm assuming that the extinction efficiency is equal to 2,
% which is a good approximation the scattering particle is larger than
% the incident wavelength (r > wl). I can make this more exact by
% computing Qe for each r value in my profile later on

density_liquid_water = 10^6;                % g/m^3

re_profile = create_droplet_profile2([retrieval(1,end), retrieval(2,end)],...
    GN_inputs.RT.z, 'altitude', GN_inputs.model.profile.type);                               % microns

% compute LWP
GN_output.LWP = 2/3 * density_liquid_water * retrieval(3,end) * trapz(GN_inputs.RT.z, (re_profile*1e-6).^3)/...
    trapz(GN_inputs.RT.z, (re_profile*1e-6).^2);           %g/m^2
% -------------------------------------------------------------------



% -------------------------------------------------------------------
% Get rid of the Nan values and create droplet profiles with the
% retrieval
idx_nans = find(isnan(retrieval(3,:)));

if isempty(idx_nans)~=true

    GN_output.tau_vector = linspace(0, retrieval(3,idx_nans(1)-1), 100);

    GN_output.re_profile = create_droplet_profile2([retrieval(1,idx_nans(1)-1), retrieval(2,idx_nans(1)-1)],...
        GN_output.tau_vector, 'optical_depth', GN_inputs.model.profile.type);

else

    GN_output.tau_vector = linspace(0, retrieval(3, end), 100);

    GN_output.re_profile = create_droplet_profile2([retrieval(1, end), retrieval(2, end)],...
        GN_output.tau_vector, 'optical_depth', GN_inputs.model.profile.type);


end
% -------------------------------------------------------------------







% -------------------------------------------------------------
% ----------- Determine the degree of non-linearity -----------
% -------------------------------------------------------------
% Rodgers (2000) pg. 83
% Define the state vector used to evaluate the non linearity of our problem
% by defining a state vector with values 1 standard deviation away from the
% retrieved state
state_vec_plus_1sigma = retrieval(:,end) + sqrt(diag(posterior_cov));
state_vec_minus_1sigma = retrieval(:,end) - sqrt(diag(posterior_cov));

% compute the measurement estimate for this state vector
meas_est_plus_1sigma = compute_forward_model_HySICS_ver2(state_vec_plus_1sigma,...
                        GN_inputs, spec_response, folder_paths);

% lastly, we need the Jacobian of the difference between the retrieved
% state vector and the state vector 1 standard deviation away
Jacobian_plus_1sigma = compute_jacobian_HySICS_ver2(state_vec_plus_1sigma,...
    meas_est_plus_1sigma, GN_inputs,...
    hysics.spec_response.value, jacobian_barPlot_flag, folder_paths);


Jacobian_plus_1sigma = compute_jacobian_HySICS_ver2(sqrt(diag(posterior_cov)),...
    meas_est_plus_1sigma, GN_inputs,...
    hysics.spec_response.value, jacobian_barPlot_flag, folder_paths);

% nonLin_degree = (new_measurement_estimate - meas_est_plus_1sigma - Jacobian_plus_1sigma);

nonLin_degree = (new_measurement_estimate - meas_est_plus_1sigma - (Jacobian_plus_1sigma - Jacobian_plus_1sigma));


% -------------------------------------------------------------


% -------------------------------------------------------------
% ------ Compute the retrieval covariance for each channel ----
% Which channels have the highest information content above that of the a priori?



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


% -----------------------------------------------------------------------------
% ----- compute percent difference between true and retrieved state vector ----
% -----------------------------------------------------------------------------

percentDiff_abs = 100.* [abs((GN_inputs.measurement.r_top - retrieval(1,end))/GN_inputs.measurement.r_top),...
    abs((GN_inputs.measurement.r_bot - retrieval(2,end))/GN_inputs.measurement.r_bot),...
    abs((GN_inputs.measurement.tau_c - retrieval(3,end))/GN_inputs.measurement.tau_c),...
    abs((GN_inputs.measurement.actpw - retrieval(4,end))/GN_inputs.measurement.actpw)]';



% -----------------------------------------------------------------------------
% ----- compute absolute difference between true and retrieved state vector ----
% -----------------------------------------------------------------------------

absDiff_stateVec = [abs((GN_inputs.measurement.r_top - retrieval(1,end))),...
    abs((GN_inputs.measurement.r_bot - retrieval(2,end))),...
    abs((GN_inputs.measurement.tau_c - retrieval(3,end))),...
    abs((GN_inputs.measurement.actpw - retrieval(4,end)))]';



% ---- Collect all outputs ----

GN_output.retrieval = retrieval;
GN_output.residual = residual;
GN_output.rss_residual = rss_residual;
GN_output.percentDiff_abs = percentDiff_abs;
GN_output.absDiff = absDiff_stateVec;
GN_output.diff_guess_prior = diff_guess_prior;
GN_output.jacobian_diff_guess_prior = jacobian_diff_guess_prior;
GN_output.posterior_cov = posterior_cov;
GN_output.Jacobian_final = Jacobian;
GN_output.H_above_aPriori = H_above_aPriori;






end