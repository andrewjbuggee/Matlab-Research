%% Create model priori covariance inputs


% Create inputs to retrieve r_top, r_bot, tau_c, cwv



% (3) use_TBLUT_estimate - this is a true or false flag that tells the code
% whether or not to use the TBLUT estimate as the initial guess


% By Andrew John Buggee

%%

function [GN_inputs] = create_model_prior_covariance_HySICS_ver2(GN_inputs, tblut, use_TBLUT_estimates, acpw)

% -------------------------------------------------------------
% -------------------------------------------------------------
% THIS CODE ASSUMES THE RETIREVAL VARIABLES ARE ALL INDEPENDENT
% -------------------------------------------------------------
% -------------------------------------------------------------


% define the model variance and mean using the Truth Table found by my
% TBLUT algorithm. Either use my own retireval estiamtes or the values
% derived by the MODIS Level 6 cloud products




if use_TBLUT_estimates==true

    %-----------------------------------------------------
    % ---- use my own TBLUT estimates to define priors ---
    %-----------------------------------------------------



    % the order of the values below: (r_top, r_bottom, tau_c)
    % The model mean is the a priori, which is separate from our initial
    % guess.





    %----------------------------------------------------
    % ----------- Set the a priori value ----------------
    %----------------------------------------------------

    % set the apriori value of cloud bottom radius as some
    % percentage of the value at the top of the cloud.
    % Using in-situ measurements of non-precipitating clouds,
    % the median value of droplet size at cloud bottom was
    %  70% the value at cloud top. This is r_bot = 0.7058*r_top
    % GN_inputs.model.apriori = [1.5*tblut.minRe, 0.5*tblut.minRe, tblut.minTau]; % expected values for the effective radius (microns) and the optical depth
    % GN_inputs.model.apriori = [tblut.minRe, 0.7058*tblut.minRe, tblut.minTau];

    % *** Test Values for new HySICS retrieval ***
    % Use for now the input used to create the measurements, but scaled to
    % represent a prior
    % for column water vapor, the units are kg/m^2 (where 1 kg/m^2 is about
    % equivelant to 1 mm)
    GN_inputs.model.apriori = [tblut.minRe, 0.5*tblut.minRe, tblut.minTau, acpw.min_interpolated];


    % The first two values are the standard deviation of the effective
    % radius at the top of the cloud and the bottom of the cloud, measured
    % in microns. The third value is the percentage of the optical depth
    % that defines the standard deviation.
    %stdev_variables = [sqrt(3), sqrt(10), sqrt(0.1 *truthTable.modisT17(1:n))];

    % Using the ensemble results from in-situ measurements of
    % non-precipitating cloud from the VOCALS-REx campaign,
    % the average deviation above the median value at cloud
    % bottom is 58.3% larger. the average deviation below the
    % median value at cloud bottom is 25% smaller.

    % Set the uncertainty of the radius at cloud top to be the
    % retireval uncertainty

    % For the a priori uncertainty of the radius at cloud bottom,
    % we scaled the bi-spectral retrieval uncertainty of effective radius
    % using the weighting function for 2.13 𝜇𝑚. over 50% of the measured
    % signal comes from the upper quartile of the cloud. Only 8% of the
    % total signal comes from the lowest quartile. Thus, we adopted a cloud bottom
    % uncertainty of a factor 6 larger than retrieved effective radius uncertainty.

    % let's define the uncertainty of the effective radius retrieval as
    % 10%. This is simular to the modis retireval uncertanties for
    % liquid water clouds over ocean with an optical thickness of
    % atleast 3
    effRad_uncert = 0.1;

    % let's define the uncertainty of the optical depth retrieval as
    % 5%. This is simular to the modis retireval uncertanties for
    % liquid water clouds over ocean with an optical thickness of
    % atleast 3
    optThick_uncert = 0.05;

    % let's define the uncertainty of the above cloud column water vapor
    colWaterVapor_uncert = 0.05;

    % stdev_variables = [GN_inputs.model.apriori(1) * effRad_uncert ...
    %     GN_inputs.model.apriori(2) * 6*effRad_uncert,...
    %     GN_inputs.model.apriori(3) * optThick_uncert];


    % *** TESTING NEW UNCERTAINTY ***
    stdev_variables = [GN_inputs.model.apriori(1) * effRad_uncert ...
        GN_inputs.model.apriori(2) * 3*effRad_uncert,...
        GN_inputs.model.apriori(3) * optThick_uncert,...
        GN_inputs.model.apriori(4) * colWaterVapor_uncert];

    % variance for the effective radius (microns squared) and optical thickness respectively
    GN_inputs.model.variance = [stdev_variables(1)^2, stdev_variables(2)^2, stdev_variables(3)^2, stdev_variables(4)^2];



    %----------------------------------------------------
    % ----------- Set the Initial guess  ----------------
    %----------------------------------------------------

    % Define teh initial guess as the a priori value
    GN_inputs.model.initialGuess = GN_inputs.model.apriori;



    %----------------------------------------------------------
    % ----------- Define the Covariance Matrix ----------------
    %----------------------------------------------------------

    % For now lets claim the desired variables are independent
    GN_inputs.model.covariance = diag(GN_inputs.model.variance);






else

    %--------------------------------------------------------------
    % ----- Use some other guess for the a priori values -----
    %--------------------------------------------------------------

    error([newline, 'No inputs defined for this option.', newline])






end









end










