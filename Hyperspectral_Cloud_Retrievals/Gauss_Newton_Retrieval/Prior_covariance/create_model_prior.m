%% --- Create Model Prior ---

% takes in the input structure created for bayesian inverse problems and
% the bi-spectral calculations


% By Andrew J. Buggee
%%

function inputs = create_model_prior(inputs,data_inputs)




if strcmp(model_prior,'custom')

    % lets create the variance and mean for each model parameter
    bayes_inputs.model.variance = [4,10]; % variance for the effective radius (microns squared) and optical thickness respectively
    bayes_inputs.model.mean = [15,15]; % expected values for the effective radius (microns) and the optical depth
    bayes_inputs.model.covaraince = diag(bayes_inputs.prior.variance);

    error('mvnpdf is not properly set up yet!')

    bayes_inputs.model.prior = mvnpdf(bayes_inputs.model.mean,bayes_inputs.model.covaraince);



elseif strcmp(model_prior,'Standard_Normal')
    % if the standard normal option is chosen, the code assumes the model
    % parameters take on the standard normal gaussian, where the the
    % expected value for each parameter is 0, and the variance is 1

    bayes_inputs.model.varaince = ones(1,bayes_inputs.num_model_parameters);
    bayes_inputs.model.mean = zeros(1,bayes_inputs.num_model_parameters);

    % create the covaraince matrix of the model parameters
    bayes_inputs.model.covaraince = diag(bayes_inputs.prior.variance);

    error('standard normal option is not properly set up yet!')


elseif strcmp(model_prior,'gaussian')
    % if the 'gaussian' option is chosen, a multivariate gaussian pdf with
    % custom variance and mean will be implemented

    % lets create the variance and mean for each model parameter
    bayes_inputs.model.variance = [4,10]; % variance for the effective radius (microns squared) and optical thickness respectively
    bayes_inputs.model.mean = [2,15]; % expected values for the effective radius (microns) and the optical depth
    bayes_inputs.model.covaraince = diag(bayes_inputs.model.variance);


    % since we are assuming the model prior is gaussian, we don't need to
    % keep track of the PDF.

end





end
