%% Create the model covariance matrix and define the a priori guesses

% Use for retrieving droplet size at cloud top and bottom

% By Andrew John Buggee

%%

function [inputs] = create_model_prior_covariance_EMIT_top_bottom(inputs, pixels2use, tblut_retrieval,use_TBLUT_estimates)

% -------------------------------------------------------------
% -------------------------------------------------------------
% THIS CODE ASSUMES THE RETIREVAL VARIABLES ARE ALL INDEPENDENT
% -------------------------------------------------------------
% -------------------------------------------------------------


% Define the number of pixels to run
num_pixels_2run = length(pixels2use.idx);



% define the model variance and mean using the Truth Table found by my
% TBLUT algorithm




if use_TBLUT_estimates==false

    %-----------------------------------------------------
    % ---- use my own TBLUT estimates to define priors ---
    %-----------------------------------------------------


    % the order of the values below: (r_top, r_bottom, tau_c)
    % The model mean is the a priori, which is separate from our initial
    % guess.


    %----------------------------------------------------
    % ----------- Set the a priori value ----------------
    %----------------------------------------------------

    inputs.model.apriori = [15, 5, 10]; % expected values for the effective radius (microns) and the optical depth

    % lets create the variance and mean for each model parameter
    % Using the same values defined by King and Vaughn (2012)
    % King and Vaughn define the standard deviation of each variable
    % The first two values are the standard deviation of the effective
    % radius at the top of the cloud and the bottom of the cloud, measured
    % in microns. The third value is the percentage of the optical depth
    % that defines the standard deviation.
    stdev_variables = [3, 5, 0.5];




    inputs.model.variance = [linspace(stdev_variables(1)^2,stdev_variables(1)^2, num_pixels_2run)',...
        linspace(stdev_variables(2)^2,stdev_variables(2)^2, num_pixels_2run)',...
        stdev_variables(3)^2 ]; % variance for the effective radius (microns squared) and optical thickness respectively



    %----------------------------------------------------
    % ----------- Set the Initial guess  ----------------
    %----------------------------------------------------

    % Define the initial guess as being similar to the a priori except that
    % we define the initial guess as having the same value for the
    % effective radius at cloud top and bottom
    inputs.model.initialGuess = inputs.model.apriori;



    %----------------------------------------------------------
    % ----------- Define the Covariance Matrix ----------------
    %----------------------------------------------------------

    % For now lets claim the desired variables are independent
    for nn = 1:num_pixels_2run
        inputs.model.covariance(:,:,nn) = diag(inputs.model.variance(nn,:));
    end




else

    %--------------------------------------------------------------
    % ----- use TBLUT retrievals for initial guess and priori -----
    %--------------------------------------------------------------


    %----------------------------------------------------
    % ----------- Set the a priori value ----------------
    %----------------------------------------------------

    for nn = 1:num_pixels_2run


        %inputs.model.apriori = [1.5*truthTable.modisR17(1:n), 0.5*truthTable.modisR17(1:n), truthTable.modisT17(1:n)]; % expected values for the effective radius (microns) and the optical depth

        % set the apriori value of cloud bottom radius as some
        % percentage of the value at the top of the cloud.
        % Using in-situ measurements of non-precipitating clouds,
        % the median value of droplet size at cloud bottom was
        %  70% the value at cloud top. This is r_bot = 0.7058*r_top
        inputs.model.apriori(nn,:) = [tblut_retrieval.minRe(nn), 0.7058*tblut_retrieval.minRe(nn),...
            tblut_retrieval.minTau(nn)];


        % lets create the variance and mean for each model parameter
        % Using the same values defined by King and Vaughn (2012)
        % King and Vaughn define the standard deviation of each variable...
        % The first two values are the standard deviation of the effective
        % radius at the top of the cloud and the bottom of the cloud, measured
        % in microns. The third value is the percentage of the optical depth
        % that defines the standard deviation.

        % Using the ensemble results from in-situ measurements of
        % non-precipitating cloud from the VOCALS-REx campaign,
        % the average deviation above the median value at cloud
        % bottom is 58.3% larger. the average deviation below the
        % median value at cloud bottom is 25% smaller.

        % Set the uncertainty of the radius at cloud top to be the
        % retireval uncertainty

        %stdev_variables = [sqrt(3), sqrt(10), sqrt(0.1)];
%         stdev_variables = [inputs.model.apriori(1)*0.1, inputs.model.apriori(2)*0.3,...
%             inputs.model.apriori(3)*0.05];
        stdev_variables = [sqrt(0.1), sqrt(1), sqrt(0.1)];


        % variance for the effective radius (microns squared) and optical thickness respectively
        inputs.model.variance(nn, :) = [stdev_variables(1)^2, stdev_variables(2)^2, stdev_variables(3)^2];




        %----------------------------------------------------
        % ----------- Set the Initial guess  ----------------
        %----------------------------------------------------

        % Define teh initial guess as the a priori value
        inputs.model.initialGuess(nn,:) = inputs.model.apriori(nn,:);



        %----------------------------------------------------------
        % ----------- Define the Covariance Matrix ----------------
        %----------------------------------------------------------

        % For now lets claim the desired variables are independent
        inputs.model.covariance(:,:,nn) = diag(inputs.model.variance(nn,:));



    end






end





end

