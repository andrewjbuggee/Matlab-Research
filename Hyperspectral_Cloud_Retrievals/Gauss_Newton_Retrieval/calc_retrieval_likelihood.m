% -Compute retrieval likelihood for gaussian model and measurement pdf -


% By Andrew J. Buggee
%%

function retrieval = calc_retrieval_likelihood(inputs,data_inputs,modis,K,offset)

% --- parse through inputs

num_parameters = inputs.num_model_parameters;
measurement_covariance_inv = (inputs.measurement.covaraince)^(-1);
model_covariance_inv = (inputs.model.covaraince)^(-1);
model_mean = inputs.model.mean';
model_var = inputs.model.variance;

pixels2use = data_inputs.pixels2use;
% first, collect the measurements into vector format
% reshape measurement data into Rodgers format for gaussian model and
% measurement pdfs
y = create_measurement_vector(modis,data_inputs);

num_pixels = size(y,2);
retrieval = zeros(num_parameters,num_pixels);

% -------------------------------------------------------------
% ---***--- NEED TO ACCOUNT FOR OFFSET IN NEW K MODEL ---***---


retrieval = zeros(num_parameters+1,num_pixels);
model_mean = [1;inputs.model.mean'];
model_var = inputs.model.variance;
model_covariance_inv = (diag([0.001,model_var]))^(-1);


% ----------------------------------------------------------------
% ---------------------------------------------------------------



for pp = 1:num_pixels
    
    % --- **** we need to redefine a new model guess each pixel **** ---
    row = pixels2use.res1km.row(pp);
    col = pixels2use.res1km.col(pp);
    
    % --- Want to perturb your initial guesses with gaussian noise? ---
    mean_r = modis.cloud.effRadius17(row,col);
    var_r = 0.1*mean_r;
    mean_T = modis.cloud.optThickness17(row,col);
    var_T = 0.1*mean_T;
    model_mean(2) = normrnd(mean_r,var_r); % effetive radius expected value
    model_mean(3) = normrnd(mean_T,var_T); % optical depth expected value
    
    % lets find the minimum K by taking the sum of the squares for the
    % derivative of each wavelength, and each model parameter

%     sum_squares = sqrt(sum(sum(K{pp},1),2));
%     [min_val,index_min] = min(sum_squares);
    
%     index_min = randi(size(K{pp},3),1);
%     some_K = K{pp}(:,:,index_min) * diag([normrnd(model_mean(1),model_var(1)),normrnd(model_mean(2),model_var(2))]); % picking a random locaiton in our grid, but this doesn't seem right
     
    retrieval(:,pp) = (K{pp}' * measurement_covariance_inv * K{pp} + model_covariance_inv)^(-1) * ...
                    (K{pp}' * measurement_covariance_inv * y(:,pp) + model_covariance_inv * model_mean);
                
end

    





end