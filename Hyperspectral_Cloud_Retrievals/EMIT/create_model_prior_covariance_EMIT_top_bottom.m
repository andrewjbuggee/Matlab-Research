%% Create the model covariance matrix and define the a priori guesses

% Use for retrieving droplet size at cloud top and bottom

% By Andrew John Buggee

%%

function [GN_inputs] = create_model_prior_covariance_EMIT_top_bottom(GN_inputs, tblut_retrieval,use_TBLUT_estimates)

% -------------------------------------------------------------
% -------------------------------------------------------------
% THIS CODE ASSUMES THE RETIREVAL VARIABLES ARE ALL INDEPENDENT
% -------------------------------------------------------------
% -------------------------------------------------------------




% define the model variance and mean using the Truth Table found by my
% TBLUT algorithm




if use_TBLUT_estimates==false

    %-----------------------------------------------------
    % ------ use my own estimates to define priors ------
    %-----------------------------------------------------


    % the order of the values below: (r_top, r_bottom, tau_c)
    % The model mean is the a priori, which is separate from our initial
    % guess.


    %----------------------------------------------------
    % ----------- Set the a priori value ----------------
    %----------------------------------------------------

    GN_inputs.model.apriori = [15, 5, 10]; % expected values for the effective radius (microns) and the optical depth

    % lets create the variance and mean for each model parameter
    % Using the same values defined by King and Vaughn (2012)
    % King and Vaughn define the standard deviation of each variable
    % The first two values are the standard deviation of the effective
    % radius at the top of the cloud and the bottom of the cloud, measured
    % in microns. The third value is the percentage of the optical depth
    % that defines the standard deviation.
    stdev_variables = [3, 5, 0.5];




    GN_inputs.model.variance = [linspace(stdev_variables(1)^2,stdev_variables(1)^2, num_pixels_2run)',...
        linspace(stdev_variables(2)^2,stdev_variables(2)^2, num_pixels_2run)',...
        stdev_variables(3)^2 ]; % variance for the effective radius (microns squared) and optical thickness respectively



    %----------------------------------------------------
    % ----------- Set the Initial guess  ----------------
    %----------------------------------------------------

    % Define the initial guess as being similar to the a priori except that
    % we define the initial guess as having the same value for the
    % effective radius at cloud top and bottom
    GN_inputs.model.initialGuess = GN_inputs.model.apriori;



    %----------------------------------------------------------
    % ----------- Define the Covariance Matrix ----------------
    %----------------------------------------------------------

    % For now lets claim the desired variables are independent
    for nn = 1:num_pixels_2run
        GN_inputs.model.covariance(:,:,nn) = diag(GN_inputs.model.variance(nn,:));
    end




else

    %--------------------------------------------------------------
    % ----- use TBLUT retrievals for initial guess and priori -----
    %--------------------------------------------------------------






    %----------------------------------------------------
    % ----------- Set the a priori value ----------------
    %----------------------------------------------------

    %inputs.model.apriori = [1.5*truthTable.modisR17(1:n), 0.5*truthTable.modisR17(1:n), truthTable.modisT17(1:n)]; % expected values for the effective radius (microns) and the optical depth

    % set the apriori value of cloud bottom radius as some
    % percentage of the value at the top of the cloud.
    % Using in-situ measurements of non-precipitating clouds,
    % the median value of droplet size at cloud bottom was
    %  70% the value at cloud top. This is r_bot = 0.7058*r_top
    %         inputs.model.apriori = [tblut_retrieval.minRe, 0.7058*tblut_retrieval.minRe,...
    %             tblut_retrieval.minTau];

    % *** Test Values for new EMIT retrieval ***
    GN_inputs.model.apriori = [tblut_retrieval.minRe, 0.5*tblut_retrieval.minRe,...
        tblut_retrieval.minTau];



    %---------------------------------------------------------------
    % ------- Set the variable covariance matrix values ------------
    %---------------------------------------------------------------

    % The main diagonal of the covariance matrix for the retrieved
    % variables is the variance. That is, if each variable is gaussian
    % distributed, the squareroot of the main diagonal is the standard
    % deviation of each variables distribution. So the units for r_top
    % and r_bot are microns squared.

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

    % The standard deviation should be in units of microns for the
    % first two variables, and unitless for the last variable

    % [sigma(r_top), sigma(r_bot), sigma(tau_c)]

    %         stdev_variables = [3, 10, 0.1];       % Values used by King and Vaughan (2012)


    %         stdev_variables = [3, 5, 0.25];       % standard deviation of the three retrieved variables

    %         stdev_variables = [sqrt(3), sqrt(10), sqrt(0.3)];       % Values used by King and Vaughan (2012)

    %         stdev_variables = [inputs.model.apriori(1)*0.1, inputs.model.apriori(2)*0.3,...
    %             inputs.model.apriori(3)*0.05];
    %stdev_variables = [sqrt(0.1), sqrt(1), sqrt(0.1)];



    % let's define the uncertainty of the effective radius retrieval as
    % 20%. This is a factor of two larger than the modis retireval uncertanties for
    % liquid water clouds over ocean with an optical thickness of
    % atleast 3
    effRad_uncert = 0.1;

    % let's define the uncertainty of the optical depth retrieval as
    % 5%. This is simular to the modis retireval uncertanties for
    % liquid water clouds over ocean with an optical thickness of
    % atleast 3
    optThick_uncert = 0.1;

    % stdev_variables = [GN_inputs.model.apriori(1) * effRad_uncert ...
    %     GN_inputs.model.apriori(2) * 6*effRad_uncert,...
    %     GN_inputs.model.apriori(3) * optThick_uncert];

    % -------------------------------
    % *** TESTING NEW UNCERTAINTY ***
    % -------------------------------
    % The uncertainty for the a priori at cloud bottom is 50%
    stdev_variables = [GN_inputs.model.apriori(1) * effRad_uncert ...
        GN_inputs.model.apriori(2) * 3*effRad_uncert,...
        GN_inputs.model.apriori(3) * optThick_uncert];




    % variance for the effective radius (microns squared) and optical thickness respectively
    GN_inputs.model.variance = [stdev_variables(1)^2, stdev_variables(2)^2, stdev_variables(3)^2];




    %----------------------------------------------------
    % ----------- Set the Initial guess  ----------------
    %----------------------------------------------------

    % Define the initial guess as the a priori value
    GN_inputs.model.initialGuess = GN_inputs.model.apriori;



    %----------------------------------------------------------
    % ----------- Define the Covariance Matrix ----------------
    %----------------------------------------------------------

    % For now lets claim the desired variables are independent
    GN_inputs.model.covariance = diag(GN_inputs.model.variance);










end





end

