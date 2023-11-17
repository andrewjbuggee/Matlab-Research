function [GN_inputs] = create_model_prior_covariance_andCloudHeight(GN_inputs, pixels2use, truthTable, use_MODIS_estimates, modis, vocalsRex)

% -------------------------------------------------------------
% -------------------------------------------------------------
% THIS CODE ASSUMES THE RETIREVAL VARIABLES ARE ALL INDEPENDENT
% -------------------------------------------------------------
% -------------------------------------------------------------


% define the model variance and mean using the Truth Table found by my
% TBLUT algorithm


% grab the pixel indices that will be used in the following analysis
indexes2run = pixels2use.res1km.linearIndex;

if use_MODIS_estimates==false

    %-----------------------------------------------------
    % ---- use my own TBLUT estimates to define priors ---
    %-----------------------------------------------------

    if GN_inputs.numPixels2Calculate<=size(truthTable,1)

        % Pick the first n pixels in the table
        n = GN_inputs.numPixels2Calculate;

        % the order of the values below: (r_top, r_bottom, tau_c)
        % The model mean is the a priori, which is separate from our initial
        % guess.





        %----------------------------------------------------
        % ----------- Set the a priori value ----------------
        %----------------------------------------------------

        GN_inputs.model.apriori = [1.5*truthTable.estR17(1:n), 0.5*truthTable.estR17(1:n), truthTable.estT17(1:n)]; % expected values for the effective radius (microns) and the optical depth
        %bayes_inputs.model.apriori = [truthTable.estR17(1:n), truthTable.estR17(1:n), truthTable.estT17(1:n)];

        % lets create the variance and mean for each model parameter
        % Using the same values defined by King and Vaughn (2012)
        % King and Vaughn define the standard deviation of each variable
        % The first two values are the standard deviation of the effective
        % radius at the top of the cloud and the bottom of the cloud, measured
        % in microns. The third value is the percentage of the optical depth
        % that defines the standard deviation.
        stdev_variables = [3, 5, 0.1 * truthTable.estT17(1:n)];




        GN_inputs.model.variance = [linspace(stdev_variables(1)^2,stdev_variables(1)^2,n)',...
            linspace(stdev_variables(2)^2,stdev_variables(2)^2,n)',...
            stdev_variables(3)^2 ]; % variance for the effective radius (microns squared) and optical thickness respectively



        %----------------------------------------------------
        % ----------- Set the Initial guess  ----------------
        %----------------------------------------------------

        % Define the initial guess as being similar to the a priori except that
        % we define the initial guess as having the same value for the
        % effective radius at cloud top and bottom
        %bayes_inputs.model.initialGuess = [1.5*truthTable.estR17(1:n), 0.5*truthTable.estR17(1:n), truthTable.estT17(1:n)];
        GN_inputs.model.initialGuess = [truthTable.estR17(1:n), truthTable.estR17(1:n), truthTable.estT17(1:n)];



        %----------------------------------------------------------
        % ----------- Define the Covariance Matrix ----------------
        %----------------------------------------------------------

        % For now lets claim the desired variables are independent
        for ii = 1:n
            GN_inputs.model.covariance(:,:,ii) = diag(GN_inputs.model.variance(ii,:));
        end







    else

        % We will only use what is in the truth table!
        GN_inputs.numPixels2Calculate = size(truthTable,1);
        n = size(truthTable,1);
        % the order of the values below: (r_top, r_bottom, tau_c)

        %----------------------------------------------------
        % ----------- Set the a priori value ----------------
        %----------------------------------------------------

        % The model mean is the a priori; our first guess
        GN_inputs.model.apriori = [truthTable.estR17, truthTable.estR17, truthTable.estT17]; % expected values for the effective radius (microns) and the optical depth


        % lets create the variance and mean for each model parameter
        % Using the same values defined by King and Vaughn (2012)
        r_top_var = 3^2;
        r_bot_var = 10^2;
        tau_var = (0.1 * truthTable.estT17)^2;

        GN_inputs.model.variance = [linspace(r_top_var,r_top_var,n)',linspace(r_bot_var,r_bot_var,n)', tau_var]; % variance for the effective radius (microns squared) and optical thickness respectively


        % For now lets claim the desired variables are independent
        for ii = 1:n

            GN_inputs.model.covariance(:,:,ii) = diag(GN_inputs.model.variance(ii,:));
        end


    end




else

    %--------------------------------------------------------------
    % ----- use MODIS retrievals for initial guess and priori -----
    %--------------------------------------------------------------


    %----------------------------------------------------
    % ----------- Set the a priori value ----------------
    %----------------------------------------------------

    for nn = 1:numel(indexes2run)


        %bayes_inputs.model.apriori = [1.5*truthTable.modisR17(1:n), 0.5*truthTable.modisR17(1:n), truthTable.modisT17(1:n)]; % expected values for the effective radius (microns) and the optical depth
        GN_inputs.model.apriori(nn,:) = [modis.cloud.effRadius17(indexes2run(nn)), modis.cloud.effRadius17(indexes2run(nn)), modis.cloud.optThickness17(indexes2run(nn))];


        % lets create the variance and mean for each model parameter
        % Using the same values defined by King and Vaughn (2012)
        % King and Vaughn define the standard deviation of each variable...
        % The first two values are the standard deviation of the effective
        % radius at the top of the cloud and the bottom of the cloud, measured
        % in microns. The third value is the percentage of the optical depth
        % that defines the standard deviation.

        %stdev_variables = [sqrt(3), sqrt(10), sqrt(0.1 *truthTable.modisT17(1:n))];
        stdev_variables = [1.5, 7, (0.2 * modis.cloud.optThickness17(indexes2run(nn)))];

        % variance for the effective radius (microns squared) and optical thickness respectively
        GN_inputs.model.variance(nn, :) = [stdev_variables(1)^2, stdev_variables(2)^2, stdev_variables(3)^2];




        %----------------------------------------------------
        % ----------- Set the Initial guess  ----------------
        %----------------------------------------------------

        % Define the initial guess as being similar to the a priori except that
        % we define the initial guess as having the same value for the
        % effective radius at cloud top and bottom

        %bayes_inputs.model.initialGuess = [1.5*truthTable.modisR17(1:n), 0.5*truthTable.modisR17(1:n), truthTable.modisT17(1:n)];
        GN_inputs.model.initialGuess(nn,:) = [modis.cloud.effRadius17(indexes2run(nn)), modis.cloud.effRadius17(indexes2run(nn)), modis.cloud.optThickness17(indexes2run(nn))];



        %----------------------------------------------------------
        % ----------- Define the Covariance Matrix ----------------
        %----------------------------------------------------------

        % For now lets claim the desired variables are independent
        GN_inputs.model.covariance(:,:,nn) = diag(GN_inputs.model.variance(nn,:));




        %----------------------------------------------------------
        % ----------- Define the Cloud Top Height ----------------
        %----------------------------------------------------------

        if GN_inputs.RT.use_MODIS_cloudTopHeight==true
            % Define cloud top height using MODIS data
            GN_inputs.RT.cloudTop_height(nn) = modis.cloud.topHeight(indexes2run(nn))/1e3;        % km

        elseif GN_inputs.RT.use_VOCALS_cloudTopHeight==true

            % Define cloud top height using Vocals Rex
            GN_inputs.RT.cloudTop_height(nn) = vocalsRex.altitude(end)/1e3;        % km

        else

            % Define a fixed cloud top height
            GN_inputs.RT.cloudTop_height(nn) = 6;           % km
            
        end










    end








    %----------------------------------------------------------
    % ----------- Define the Cloud Top Depth ------------------
    %----------------------------------------------------------

    if GN_inputs.RT.use_VOCALS_cloudDepth==true

        % Define cloud depth using Vocals Rex
        GN_inputs.RT.cloudDepth = (vocalsRex.altitude(end) - vocalsRex.altitude(1))/1e3;      % km

    else

        % Define a custom cloud depth
        GN_inputs.RT.cloudDepth = 1;            % km

    end






    % Define number of layers to use in libRadTran when defining
    % vertically inhomogenous clouds
    GN_inputs.RT.cloud_layers = 10;



end












end



















% ------------------ OLD CODE --------------------------


% if bayes_inputs.numPixels2Calculate<=size(truthTable,1)
%
%     % Pick the first n pixels in the table
%     n = bayes_inputs.numPixels2Calculate;
%
%     % the order of the values below: (r_top, r_bottom, tau_c)
%     % The model mean is the a priori, which is separate from our initial
%     % guess.
%
%
%     if use_MODIS_estimates==false
%
%         %--------------------------------------------
%         % use my own TBLUT estimates to define priors
%         %--------------------------------------------
%
%
%
%
%         %----------------------------------------------------
%         % ----------- Set the a priori value ----------------
%         %----------------------------------------------------
%
%         bayes_inputs.model.apriori = [1.5*truthTable.estR17(1:n), 0.5*truthTable.estR17(1:n), truthTable.estT17(1:n)]; % expected values for the effective radius (microns) and the optical depth
%         %bayes_inputs.model.apriori = [truthTable.estR17(1:n), truthTable.estR17(1:n), truthTable.estT17(1:n)];
%
%         % lets create the variance and mean for each model parameter
%         % Using the same values defined by King and Vaughn (2012)
%         % King and Vaughn define the standard deviation of each variable
%         % The first two values are the standard deviation of the effective
%         % radius at the top of the cloud and the bottom of the cloud, measured
%         % in microns. The third value is the percentage of the optical depth
%         % that defines the standard deviation.
%         stdev_variables = [sqrt(3), sqrt(10), sqrt(0.3)];
%
%
%
%
%         bayes_inputs.model.variance = [linspace(stdev_variables(1)^2,stdev_variables(1)^2,n)',...
%             linspace(stdev_variables(2)^2,stdev_variables(2)^2,n)',...
%             stdev_variables(3)^2 *truthTable.estT17(1:n)]; % variance for the effective radius (microns squared) and optical thickness respectively
%
%
%
%         %----------------------------------------------------
%         % ----------- Set the Initial guess  ----------------
%         %----------------------------------------------------
%
%         % Define the initial guess as being similar to the a priori except that
%         % we define the initial guess as having the same value for the
%         % effective radius at cloud top and bottom
%         %bayes_inputs.model.initialGuess = [1.5*truthTable.estR17(1:n), 0.5*truthTable.estR17(1:n), truthTable.estT17(1:n)];
%         bayes_inputs.model.initialGuess = [truthTable.estR17(1:n), truthTable.estR17(1:n), truthTable.estT17(1:n)];
%
%
%
%         %----------------------------------------------------------
%         % ----------- Define the Covariance Matrix ----------------
%         %----------------------------------------------------------
%
%         % For now lets claim the desired variables are independent
%         for ii = 1:n
%             bayes_inputs.model.covariance(:,:,ii) = diag(bayes_inputs.model.variance(ii,:));
%         end
%
%
%
%
%
%
%
%     else
%
%         %---------------------------------------------------
%         % use MODIS retrievals for initial guess and priori
%         %---------------------------------------------------
%
%
%         %----------------------------------------------------
%         % ----------- Set the a priori value ----------------
%         %----------------------------------------------------
%
%         bayes_inputs.model.apriori = [1.5*truthTable.modisR17(1:n), 0.5*truthTable.modisR17(1:n), truthTable.modisT17(1:n)]; % expected values for the effective radius (microns) and the optical depth
%         %bayes_inputs.model.apriori = [truthTable.modisR17(1:n), truthTable.modisR17(1:n), truthTable.modisT17(1:n)];
%
%
%         % lets create the variance and mean for each model parameter
%         % Using the same values defined by King and Vaughn (2012)
%         % King and Vaughn define the standard deviation of each variable...
%         % The first two values are the standard deviation of the effective
%         % radius at the top of the cloud and the bottom of the cloud, measured
%         % in microns. The third value is the percentage of the optical depth
%         % that defines the standard deviation.
%         %stdev_variables = [sqrt(3), sqrt(10), sqrt(0.1 *truthTable.modisT17(1:n))];
%         stdev_variables = [sqrt(0.75), sqrt(2.25), sqrt(0.05 *truthTable.modisT17(1:n))];
%
%         bayes_inputs.model.variance = [linspace(stdev_variables(1)^2,stdev_variables(1)^2,n)',...
%             linspace(stdev_variables(2)^2,stdev_variables(2)^2,n)',...
%             stdev_variables(3).^2 ]; % variance for the effective radius (microns squared) and optical thickness respectively
%
%
%
%
%
%         %----------------------------------------------------
%         % ----------- Set the Initial guess  ----------------
%         %----------------------------------------------------
%
%         % Define the initial guess as being similar to the a priori except that
%         % we define the initial guess as having the same value for the
%         % effective radius at cloud top and bottom
%
%         %bayes_inputs.model.initialGuess = [1.5*truthTable.modisR17(1:n), 0.5*truthTable.modisR17(1:n), truthTable.modisT17(1:n)];
%         bayes_inputs.model.initialGuess = [truthTable.modisR17(1:n), truthTable.modisR17(1:n), truthTable.modisT17(1:n)];
%
%
%         %----------------------------------------------------------
%         % ----------- Define the Covariance Matrix ----------------
%         %----------------------------------------------------------
%
%         % For now lets claim the desired variables are independent
%         for ii = 1:n
%             bayes_inputs.model.covariance(:,:,ii) = diag(bayes_inputs.model.variance(ii,:));
%         end
%
%
%         %------------------------------------------------------
%         % ----------- Define cloud models inputs --------------
%         %------------------------------------------------------
%
%         % Define a custom cloud depth
%         %bayes_inputs.model.cloudDepth = 0.5;            % km
%
%         % Define cloud depth using Vocals Rex
%         bayes_inputs.model.cloudDepth = (vocalsRex.altitude(end) - vocalsRex.altitude(1))/1e3;      % km
%
%         % Define cloud top height using Vocals Rex
%         bayes_inputs.model.cloudTop_height = vocalsRex.altitude(end)/1e3;           % km
%
%         % Define cloud top height using MODIS data
%         %bayes_inputs.model.cloudTop_height = modis.cloud.topHeight(pixel_row, pixel_col)/1e3;        % km
%
%
%         % Define number of layers to use in libRadTran when defining
%         % vertically inhomogenous clouds
%         bayes_inputs.model.cloud_layers = 10;
%
%
%
%     end
%
%
%
%
%
% else
%
%     % We will only use what is in the truth table!
%     bayes_inputs.numPixels2Calculate = size(truthTable,1);
%     n = size(truthTable,1);
%     % the order of the values below: (r_top, r_bottom, tau_c)
%
%     %----------------------------------------------------
%     % ----------- Set the a priori value ----------------
%     %----------------------------------------------------
%
%     % The model mean is the a priori; our first guess
%     bayes_inputs.model.apriori = [truthTable.estR17, truthTable.estR17, truthTable.estT17]; % expected values for the effective radius (microns) and the optical depth
%
%
%     % lets create the variance and mean for each model parameter
%     % Using the same values defined by King and Vaughn (2012)
%     r_top_var = 3;
%     r_bot_var = 10;
%     percent_tau = 0.05;
%
%     bayes_inputs.model.variance = [linspace(r_top_var,r_top_var,n)',linspace(r_bot_var,r_bot_var,n)',percent_tau*truthTable.estT17]; % variance for the effective radius (microns squared) and optical thickness respectively
%
%
%     % For now lets claim the desired variables are independent
%     for ii = 1:n
%
%         bayes_inputs.model.covariance(:,:,ii) = diag(bayes_inputs.model.variance(ii,:));
%     end
%
%
% end
%
%
% end

