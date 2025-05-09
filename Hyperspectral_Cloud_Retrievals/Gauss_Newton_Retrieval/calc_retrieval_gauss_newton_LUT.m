%% --- Calculate Optical Retrieval Using Gauss Newton Method



% By Andrew J. Buggee
%%

function retrieval = calc_retrieval_gauss_newton_LUT(bayes_inputs,data,data_inputs)

% ----- unpack inputs -----

model_cov = bayes_inputs.model.covaraince; % model parameter covariance matrix
measurement_cov = bayes_inputs.measurement.covariance; % measurement covaraince matrix

%model_mean = bayes_inputs.model.mean; % a priori expected value for the model parameters

measurements = create_measurement_vector(data,data_inputs); % each column represents one pixel
biSpectral_estimate = [data_inputs.truthTable.estR17,data_inputs.truthTable.estT17];

num_iterations = bayes_inputs.GN_iterations; % number of iterations to preform
num_parameters = bayes_inputs.num_model_parameters; % number of parameters to solve for
num_pixels = data_inputs.inputs.pixels.num_2calculate;

% ----- unpack look-up table information -----
num_r = size(data_inputs.R,2); % number of r values ran in the RT model
num_T = size(data_inputs.R,3); % number of tau values ran in the RT model
num_bands = size(data_inputs.R,4); % number of spectral bands ran in the RT model

re = data_inputs.inputs.re; % values of re that were ran through the RT model to create the look-up table
tau_c = data_inputs.inputs.tau_c; % values of tau_c that were ran through the RT model to create the look-up table
[Tau,Re] = meshgrid(tau_c,re);

% --- Create iterative Gauss-Newton Solver ----


retrieval = zeros(num_parameters,num_iterations+1,num_pixels); % we include the starting point, which is outside the number of iterations

for pp = 1:num_pixels
    
    % we've set the effective radius at the top and bottom to be the same
    % value. This is our initial guess. By setting this to be the a priori
    % guess we can compute the inforamtion gain between the bi-spectral
    % approach and the hyperspectal approach
    
    retrieval(:,1,pp) = [biSpectral_estimate(1),biSpectral_estimate(1),biSpectral_estimate(2)];
    model_mean = retrieval(:,1,pp);
    
    look_up_tables = reshape(data_inputs.R(pp,:,:,:),num_r,num_T,num_bands);
    
    for ii = 1:num_iterations
        
        % at each iteration I need to compute the forward model at my current
        % state vector estimate
        for bb = 1:num_bands
            
            forward_model(bb) = interp2(Tau,Re,look_up_tables(:,:,bb),
        
        
        
        
        
        
        
        
        
        
        
        
    end
    
end






end