%% This function will write a .DAT ice cloud file for LibRadTran



% INPUTS:
%   (1) re - effective ice particle radius (microns) - this is either a single
%   value, a vector, or a matrix. A single value for re tells the function
%   to create a cloud with a single layer containing a constant particle
%   radius value. A vector tells the function to create a single ic file
%   with a particle profile. The length of the vector is equal to the number
%   of layers modeled. A matrix tells the function to create multiple ic
%   files, where the number of columns is equal to the number of ic files
%   created. The number of rows is equal to the number of layers modeled.
%   To create many ice cloud files at once that model a homogenous cloud,
%   simply set the column vectors of re to be identical values.
%   ***IMPORTANT*** the re values must start at the cloud bottom with the
%   first value (or row). The last value (or row) is the particle size at
%   cloud top.

%   (2) tau_c - cloud optical depth (unitless) - this is the cloud optical
%   depth, which is a monochromatic calculation. There is a single value
%   that defines the optical depth of the cloud. LibRadTran defines cloud
%   files in terms of two values that do not depend on wavelength: re and
%    ice content (IWC). But usually people talk about clouds as
%   having a certain particle size and a certain optical thickness. Enter a
%   single value for a single ic file, or a vector if you're creating
%   multiple ic files. If re is a matrix, and tau_c is a single value, each
%   ic file will have the entered tau_c. If each value in the re matrix
%   needs a unique tau value

%   (3) z_topBottom - altitude above sea level (kilometers) - this is a
%   vector with two values: [z_cloudTop, z_cloudBottom]. LibRadTran
%   constructs a cloud by creating layers, where each layer is homogenous
%   untill the next layer is defined. z_cloudTop defines where the cloud
%   ends; this is where the IWC should go to zero. This function will
%   compute a z vector equal in length to that of re using z_topBottom and
%   the geometric thickness H. If re is a matrix, the function expects
%   z_topBottom to be a matrix, where each column is a new ic file. If re
%   is a matrix, and z_topBottom is a vector, then this will be used for
%   every ic file.


%   (5) lambda - wavelength that defines the cloud optical depth
%   (nanometers) - If creating a single ic file, lambda is a single value. If
%   creating multiple ic files, lambda is a vector equal in length to the
%   number of columns of re. If re is a matrix and lambda is a single
%   value, this value will be used for each ic file created.

%   (6) distribution_str - a string telling the code which particle size
%   distribution to use  - One can chose from two options:
%       (a) 'mono' - monodispersed distribution
%       (b) 'gamma' - gamma particle distribution.
%       *** IMPORTANT *** For now, this function will NOT
%       use precomputed mie calculations using a gamma particle
%       distribution. The values returned by LibRadTran appear erroneously
%       high. Instead, if one wishes to use a gamma particle distribution,
%       the homogenous pre-computed mie table will be used to retrieve mie
%       properties, and then this code will integrate those values over the
%       size distribution.

%   (6) distribution_var - the variance of the size distribution, if
%   applicable. If one is modelling a homogenous cloud, this input will be
%   ignored, and the value can be anything. *** To match the optical
%   properties mie table precomputed by libRadtran, use a gamma
%   distribution alpha parameter of 7 ***

%   (7) vert_homogeneity_str - a string telling the code if the cloud is to be
%   modeled as vertically homogeneous. If vertically homogenous, the code will
%   assume every single r_e value represents a single cloud with a constant
%   particle radius. If vertically non-homogenous, each column of re is assumed
%   to be a single cloud with a particle profile.
%       (a) 'vert-homogeneous' - homogenous cloud where the entire cloud can be
%       modeled as a single layer with constant properties. If this option
%       is chosen, the code will expect a vector for re, where each entry
%       represents a different cloud.
%       (b) 'vert-non-homogeneous' - a non-homogeneous cloud implies a cloud
%       with multiple layers where the properties vary. If this option is
%       chosen, the code expects a single column vector for re, or a
%       matrix, where each column vector represents a cloud


%   (8) parameterization_str - a string telling the code which
%   parameterization to use for computing the mie scattering properties and
%   optical depth from ice ice content and effective particle radius:
%       (a) 'mie'
%       (b)

%   (9) index - this is the unique identifier that ensures files are not
%   written over one another. If one file is created, the index will be 1.
%   If many fiels are written in a for loop, each file will be tagged with
%   the number in the loop.

% OUTPUTS:
%   (1) .Dat file saved in the libRadTran folder:
%   /.../libRadtran-2.0.4/data/ic

% All look up tables were computed using the opensoure radiative transfer
% code, LibRadTran

% By Andrew John Buggee

%%

function [fileName] = write_ic_file(re, tau_c, z_topBottom, lambda, distribution_type, distribution_var,...
    vert_homogeneous_str, parameterization_str, index)

% ------------------------------------------------------------
% ---------------------- CHECK INPUTS ------------------------
% ------------------------------------------------------------

% Check to make sure there are 8 inputs, particle radius, cloud optical
% depth, and the altitude vector associated with this cloud


if nargin~=9
    error([newline,'Not enough inputs. Need 8: particle effective radius, optical depth, altitude,',...
        [' wavelength, particle distribution type, variance of the particle distribution,' ...
        ' homogeneity type, the parameterization used to compute IWC, and the unique file index.'], newline])
end

% Check to make sure re is the same length as the altitude vector

% first check to see if z_topBottom is a vector or a matrix
if size(z_topBottom,1)==1 || size(z_topBottom,2)==1
    % If true, then there must be two entries
    if length(z_topBottom)~=2
        error([newline,'Need two values for z_topBottom: altitude at cloud bottom top and cloud bottom', newline])
    end

    % make sure its a  column vector
    z_topBottom = reshape(z_topBottom,[],1);
    H = z_topBottom(1) - z_topBottom(2);            % geometric thickness

elseif size(z_topBottom,1)>1 && size(z_topBottom,2)>1
    % if true, then there can only be two rows, and it must be equal in
    % size to the r matrix
    if size(z_topBottom,1)~=2
        error([newline,'Need two values for z_topBottom: altitude at cloud bottom top and cloud bottom', newline])

    elseif size(z_topBottom,2)~=size(re,2) || size(z_topBottom,2)==1
        error([newline,'z_topBottom must have the same number of columns as re, or a single column that is used for all ic files.', newline])

    end

    H = z_topBottom(:,1) - z_topBottom(:,2);        % geometric cloud thickness

end



if length(lambda)>1 && length(lambda)~=size(re,2)

    error([newline,'Lambda must be either a single value or a vector equal in legnth to the number of columns in re.', newline])
end

if length(tau_c)>1 && length(tau_c)~=size(re,2)

    error([newline,'The optical depth must be either a single value or a vector equal in legnth to the number of columns in re.', newline])
end


% Check to make sure the distribution string is one of two possible values

if strcmp(distribution_type, 'mono')==false && strcmp(distribution_type, 'gamma')==false

    error([newline,'I dont recognize the particle distribution. Must be either "mono" or "gamma"', newline])
end


% Check to make sure the homogeneity string is one of two possible values

if strcmp(vert_homogeneous_str, 'vert-homogeneous')==false && strcmp(vert_homogeneous_str, 'vert-non-homogeneous')==false

    error([newline,'I dont recognize the homogeneity string. Must be either "vert-homogenous" or "vert-non-homogeneous"', newline])
end


% ----- Check to see if there are any NaNs in the radius vector -----

if any(isnan(re))==true

    error([newline, 'The effective radius has atleast one NaN value.', newline])
end







% Determine which computer you're using
computer_name = whatComputer;

% Find the folder where the mie calculations are stored
% find the folder where the ice cloud files are stored.
if strcmp(computer_name,'anbu8374')==true

    mie_calc_folder_path = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/Mie_Calculations/';
    ice_cloud_folder_path = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/ic/';

elseif strcmp(computer_name,'andrewbuggee')==true

    mie_calc_folder_path = '/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-particle-Retrieval/LibRadTran/libRadtran-2.0.4/Mie_Calculations/';
    ice_cloud_folder_path = '/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-particle-Retrieval/LibRadTran/libRadtran-2.0.4/data/ic/';

elseif strcmp(computer_name,'curc')==true

    mie_calc_folder_path = '/scratch/alpine/anbu8374/Mie_Calculations/';
    ice_cloud_folder_path = '/projects/anbu8374/software/libRadtran-2.0.5/data/ic/';


end

%%

% ------------------------------------------------------------
% ---------------------- COMPUTE IWC -------------------------
% ------------------------------------------------------------


rho_ice = 0.91675;                     % grams/cm^3 - density of ice at 0 C - wikipedia
% --------------------------------------------------


% --- STEP THROUGH COLUMNS OF re -----

% if re is a matrix, the code assumes that each column vector is a single
% cloud with a vertically-inhomogenous particle profile defined by the
% values in the column
if size(re,1)>1 && size(re,2)>1

    % We want to get all the mie properties we need before we loop through
    % the number of files

    num_files_2write = size(re,2);          % number of ic files to create



    % -------------------------------------------------------------------
    % ------ open the precomputed mie table and interpolate! ------------
    % -------------------------------------------------------------------

    % for writing ice cloud files, we only need the extinction efficiency
    % Since this function is used often, we've created a file with just Q_ext

    justQ = true;                       % Load the precomputed mie table of only Q_ext values



    if strcmp(parameterization_str, 'mie')==true

        % lambda has to be the same size as re
        if length(lambda)==num_files_2write

            % remember that vectorizing a matrix (:) will stack column vectors
            % on top of one another

            lambda = reshape(lambda,1,[]); % each column represents a new file

            % create identical column vecotrs
            lambda = repmat(lambda,size(re,1),1);


            % **** ONLY INTERPOLATING HOMOGENOUS MIE COMPUTATIONS ****
            % ********************************************************
            % IF GAMMA DISTRIBUTION DESIRED, CODE WILL MANUALLY
            % INTEGRATE OVER THE SIZE DISTRIBUTION DEFINED

            if strcmp(distribution_type,'gamma')==true

                % If we wish to estimate the mie properties of ice
                % for a distribution of particles, then we can skip the
                % pre-computed mie tables and estimate the values using the
                % average_mie_over_size_distribution directly

                % This function only deals with  ice clouds
                % define the index of refraction
                index_of_refraction = 'ice';

                % integrate over a size distribution to get an average
                [~, Qe_avg, ~] = average_mie_over_size_distribution(re, distribution_var, lambda,...
                    index_of_refraction, distribution_type, index);

            elseif strcmp(distribution_type,'mono')==true
                yq = interp_mie_computed_tables([repmat(lambda,numel(re),1), re], 'mono', justQ);

            else

                error([newline,'Invaled distribution type',newline])

            end



        elseif length(lambda)==1

            % **** ONLY INTERPOLATING HOMOGENOUS MIE COMPUTATIONS ****
            % ********************************************************
            % IF GAMMA DISTRIBUTION DESIRED, CODE WILL MANUALLY
            % INTEGRATE OVER THE SIZE DISTRIBUTION DEFINED

            if strcmp(distribution_type,'gamma')==true

                % If we wish to estimate the mie properties of  ice
                % for a distribution of particles, then we can skip the
                % pre-computed mie tables and estimate the values using the
                % average_mie_over_size_distribution directly

                % This function only deals with  ice clouds
                % define the index of refraction
                index_of_refraction = 'ice';

                % integrate over a size distribution to get an average
                [~, Qe_avg, ~] = average_mie_over_size_distribution(re, distribution_var, lambda,...
                    index_of_refraction, distribution_type, index);

            elseif strcmp(distribution_type,'mono')==true
                yq = interp_mie_computed_tables([repmat(lambda,numel(re),1), re], 'mono', justQ);

            else

                error([newline,'Invaled distribution type',newline])

            end

        end

    elseif strcmp(parameterization_str,'2limit')==true
        % set the value to be 2 for all calculations
        yq = 2*ones(length(re),5);

    end






    % if re is a vector and the vertically homogeneity is defined as
    % 'vert-non-homogeneous' then the the code assumes the vector defines a
    % single cloud with a particle profile
elseif ((size(re,1)==1 && size(re,2)>1) || (size(re,1)>1 && size(re,2)==1)) &&...
        strcmp(vert_homogeneous_str, 'vert-non-homogeneous')==true

    num_files_2write = 1;

    % re must be a column vector
    re = reshape(re,[],1);

    % the distribution variance should be a column vector
    distribution_var = reshape(distribution_var, [], 1);

    % -------------------------------------------------------------------
    % ------ open the precomputed mie table and interpolate! ------------
    % -------------------------------------------------------------------

    % for writing ice cloud files, we only need the extinction efficiency
    % Since this function is used often, we've created a file with just Q_ext

    justQ = true;                       % Load the precomputed mie table of only Q_ext values

    if strcmp(parameterization_str,'mie')==true


        if strcmp(distribution_type,'gamma')==true

            % -------------------------------------------------------------
            % --- MANUALLY INTEGRATE OVER THE SIZE DISTRIBUTION DEFINED ---
            % -------------------------------------------------------------
            % If we wish to estimate the mie properties of  ice
            % for a distribution of particles, then we can skip the
            % pre-computed mie tables and estimate the values using the
            % average_mie_over_size_distribution directly

            % This function only deals with  ice clouds
            % define the index of refraction
            %             index_of_refraction = 'ice';

            % this loop applies to a vertical particle profile. For now we
            % will apply the same distribution variance to each level in
            % the cloud.

            % integrate over a size distribution to get an average
            %             [~, Qe_avg, ~] = average_mie_over_size_distribution(re, distribution_var,...
            %                 lambda,index_of_refraction, distribution_type, index);


            % -------------------------------------------------------
            % ----------- USE LIBRADTRAN MIE CALCULATIONS -----------
            % -------------------------------------------------------
            % Libradtran doesn't compute the efficieny when a distribution
            % is specified. It computes the bulk coefficient per unit
            % concentration. For ice, since the density is 1 g/m^3, we
            % can simply multiply the output with the  ice content
            % and integrate over the path to get the optical depth.


            % What mie code should we use to compute the scattering properties?
            mie_program = 'MIEV0';               % type of mie algorithm to run

            % This function only deals with  ice clouds
            % define the index of refraction
            index_of_refraction = 'ice';

            size_distribution = {'gamma', distribution_var(1)};           % particle distribution

            % Do you want a long or short error file?
            err_msg_str = 'verbose';



            % The radius input is defined as [r_start, r_end, r_step].
            % where r_step is the interval between radii values (used only for
            % vectors of radii). A 0 tells the code there is no step. Finally, the
            % radius values have to be in increasing order.
            ext_bulk_coeff_per_IWC = zeros(length(re), 1);

            for rr = 1:length(re)

                mie_radius = [re(rr), re(rr), 0];    % microns


                % Create a mie file
                [input_filename, output_filename, mie_folder] = write_mie_file(mie_program, index_of_refraction,...
                    mie_radius, lambda, size_distribution, err_msg_str, index);

                % run the mie file
                [~] = runMIE(mie_folder,input_filename,output_filename);

                % Read the output of the mie file
                [ds,~,~] = readMIE(mie_folder,output_filename);

                ext_bulk_coeff_per_IWC(rr) = ds.Qext;       % km^-1 / (cm^3 / m^3)

            end
            % --------------------------------------------------------------




        elseif strcmp(distribution_type,'mono')==true
            yq = interp_mie_computed_tables([repmat(lambda,numel(re),1), re], 'mono', justQ);

        else

            error([newline,'Invaled distribution type',newline])

        end

    elseif strcmp(parameterization_str, '2limit')==true

        yq = 2*ones(length(re),5);

    end


    % if the re input is a vector and the homogenous string is defined as
    % vertically homogeneous, then the code assumes each value in the vector is
    % a single cloud, and the each value defines the homogenous particle size
    % for that cloud.
elseif (size(re,1)==1 || size(re,2)==1) && strcmp(vert_homogeneous_str, 'vert-homogeneous')==true

    num_files_2write = length(re);

    % re must be a column vector
    re = reshape(re,[],1);

    % the distribution variance should be a column vector
    distribution_var = reshape(distribution_var, [], 1);

    % -------------------------------------------------------------------
    % ------ open the precomputed mie table and interpolate! ------------
    % -------------------------------------------------------------------

    % for writing ice cloud files, we only need the extinction efficiency
    % Since this function is used often, we've created a file with just Q_ext

    justQ = true;                       % Load the precomputed mie table of only Q_ext values

    if strcmp(parameterization_str,'mie')==true




        if strcmp(distribution_type,'gamma')==true


            % -------------------------------------------------------
            % ----------- USE LIBRADTRAN MIE CALCULATIONS -----------
            % -------------------------------------------------------
            % Libradtran doesn't compute the efficieny when a distribution
            % is specified. It computes the bulk coefficient per unit
            % concentration. For ice, since the density is 1 g/m^3, we
            % can somply multiply the output with the  ice content
            % and integrate over the path to get the optical depth.


            % What mie code should we use to compute the scattering properties?
            mie_program = 'MIEV0';               % type of mie algorithm to run

            % This function only deals with  ice clouds
            % define the index of refraction
            index_of_refraction = 'ice';

            size_distribution = {'gamma', distribution_var(1)};           % particle distribution

            % Do you want a long or short error file?
            err_msg_str = 'verbose';


            % The radius input is defined as [r_start, r_end, r_step].
            % where r_step is the interval between radii values (used only for
            % vectors of radii). A 0 tells the code there is no step. Finally, the
            % radius values have to be in increasing order.

            % Libradtran doesn't compute the efficieny when a distribution
            % is specified. It computes the bulk coefficient per unit
            % concentration. For ice, since the density is not 1 g/m^3, we
            % have to multiply the output with the density of ice and
            % the ice water content, then integrate over the path to get the optical depth.
            ext_bulk_coeff_per_IWC = zeros(length(re), 1);

            for rr = 1:length(re)

                mie_radius = [re(rr), re(rr), 0];    % microns


                % Create a mie file
                [input_filename, output_filename, mie_folder] = write_mie_file(mie_program, index_of_refraction,...
                    mie_radius, lambda, size_distribution, err_msg_str, index);

                % run the mie file
                [~] = runMIE(mie_folder,input_filename,output_filename);

                % Read the output of the mie file
                [ds,~,~] = readMIE(mie_folder,output_filename);
                
                % bulk extinction coefficient per unit concentration of ice
                % particles
                ext_bulk_coeff_per_IWC(rr) = ds.Qext;       % km^-1 / (cm^3 / m^3)

            end
            % --------------------------------------------------------------





        
        else

            error([newline,'Invaled distribution type',newline])

        end


    else

        error([newline,'Invaled parameterization type',newline])


    end





else

    error([newline,'re is not a vector or a matrix! Check your inputs!', newline])

end



% grab the extinction efficiency values

if strcmp(distribution_type,'gamma')==true
%     Qext = Qe_avg';         % Extinction efficiency
    %Qext = linspace(2.0816, 2.0816, length(re))';        % value to match libRadTran

elseif strcmp(distribution_type,'mono')==true
    Qext = reshape(yq(:,3),[],num_files_2write);         % convert this back into a matrix corresponging to re

end




% now we will step through each ic file that needs to be created

% if H is a single value and num_files_2write is greater than 1, we will
% repeat it to create a vector with the same length
if length(H)==1 && num_files_2write>1
    H = repmat(H,num_files_2write,1);           % km - geometric thickness
end


% if tau_c is a single value and num_files_2write is greater than 1, we will
% repeat it to create a vector with the same length
if length(tau_c)==1 && num_files_2write>1
    tau_c = repmat(tau_c,num_files_2write,1);
end

if size(z_topBottom,2)==1 && num_files_2write>1
    z_topBottom = repmat(z_topBottom,1,num_files_2write);
end

% How many files are being created?
fileName = cell(1,num_files_2write);


% How many layers to model in the cloud?
if size(re,1)>1 && size(re,2)>1 && strcmp(vert_homogeneous_str, 'vert-non-homogeneous')==true
    nLayers = size(re,1)+1;             % Number of altitude levels we need to define a cloud
elseif (size(re,1)==1 || size(re,2)==1) && strcmp(vert_homogeneous_str, 'vert-non-homogeneous')==true
    nLayers = length(re)+1;             % Number of altitude levels we need to define a cloud
elseif (size(re,1)==1 || size(re,2)==1) && strcmp(vert_homogeneous_str, 'vert-homogeneous')==true
    nLayers = 1;

end


for nn = 1:num_files_2write


    % -------------------------------------------
    % ------ Create altitude vector! ------------
    % -------------------------------------------

    % the length of the altitude vector should be 1 unit longer than the length
    % of the effective radius. Thats because the last value in the altitude
    % vector is the altitude at cloud top, where the IWC has gone to zero


    % z must be a column vector
    if nLayers==1
        z = flipud(z_topBottom(:,nn));

    else
        z = linspace(z_topBottom(2,nn), z_topBottom(1,nn), nLayers)';                 % km - altitude vector
    end


    % -------------------------------------------------------------------
    % ------ compute number concentration and  ice content ------
    % -------------------------------------------------------------------

    % we need a number concentration for each file that is created

    if nLayers>1 && strcmp(distribution_type, 'gamma')==true

        % We could just integrate the size distributiion to get the total
        % number concentration, but we've chosen an arbitrary value for the
        % effective variance because it doesn't have much effect on
        % reflectance measurements over the solar spectral region. More
        % importantly, we want to connect two user defined variables, cloud
        % optical depth and effective radius, to the number concentration,
        % and thus the  ice content.
        %         z_meters = (z(1:end-1)-z(1))*1e3;       % meters - geometric depth, normalized
        %         re_meters = (re(:,nn)*1e-6);            % meters - effective radius converted to meters
        %
        %         Nc = tau_c(nn)./(pi*trapz(z_meters, Qext(:,nn).* re_meters.^2));                % m^(-3) - number concentration
        %
        %         % ------------------------------------------------------------------
        %         % --- Solve for the total  ice Content over the entire cloud ---
        %         % number concentration is constant with height. We make the
        %         % assumption that all particles can be modeled as the effective
        %         % radius. So the IWC simple changes with effective radius
        %         % ** LibRadTran requires IWC in units of grams/m^3 **
        %         lic = 4/3 * pi * rho_ice * re_meters.^3 .* Nc;                    % g/m^3 - grams of ice per meter cubed of air
        % -----------------------------------------------------------------



        % ----------------------------------------------------------------
        % ******** Integrating over monodispersed mie caluclation ********
        % ----------------------------------------------------------------
        % -- Assuming  ice content increases linearly with depth -
        %         re_meters = (re(:,nn)*1e-6);            % meters - effective radius converted to meters
        %         z_meters_midpoint = ((z(1:end-1)-z(1)) + (z(2)-z(1))/2)*1e3;       % meters - geometric depth, normalized
        %         dz = z_meters_midpoint(2)-z_meters_midpoint(1);           % meters
        %
        %
        %         slope = (4*rho_ice * tau_c) /(3*dz * sum(Qext .* z_meters_midpoint ./re_meters));     % g/m^3/m - slope of the lic profile
        %
        %         % solve for the linear  ice content profile
        %         lic = slope * z_meters_midpoint;                     % g/m^3 - grams of ice per meter cubed of air
        % ----------------------------------------------------------------



        % -----------------------------------------------------------------
        % ** Using libRadTran mie calculations with a size distribution ***
        % -----------------------------------------------------------------
        % ** Assuming  ice content increases linearly with depth **
        
        %z_kilometers_midpoint = ((z(1:end-1)-z(1)) + (z(2)-z(1))/2);       % kilometers - geometric depth at midpoint of each layer
        z_kilometers_upper_boundary = z(2:end) - z(1);                     % kilometers - geometric depth at upper boundary of each cloud layer
        dz_km = z(2) - z(1);           % kilometers

        %slope = tau_c /(dz_km * sum(ext_bluk_coeff_per_IWC .* z_kilometers_midpoint ));     % g/m^3/m - slope of the lic profile
        slope = tau_c(nn) /(dz_km * sum(ext_bulk_coeff_per_IWC .* z_kilometers_upper_boundary ));     % g/m^3/m - slope of the lic profile

        % solve for the linear  ice content profile
        %lic = slope * z_kilometers_midpoint;                     % g/m^3 - grams of ice per meter cubed of air
        iwc = slope * z_kilometers_upper_boundary;                     % g/m^3 - grams of ice per meter cubed of air
        % ----------------------------------------------------------------



        % -----------------------------------------------------------------
        % ******** compute IWC by integrated the size distribution ********
        % -----------------------------------------------------------------
        % *** There is another way to solve for the IWC ***
        % We've made an assumption about the particle size distribution and
        % we've computed the total number concentration. We can solve for
        % the  ice content by integrating the size distribution
        % *** IMPORTANT *** We have to play with the distribution width to
        % get the correct optical depth
        %         if strcmp(distribution_str,'gamma')==true
        %
        %             %distribution_var = 27;
        %             lic = zeros(size(re));
        %
        %             for zz = 1:length(re)
        %
        %                 [nr,r] = gamma_size_distribution_libRadTran2(re(zz), distribution_var(zz), Nc);       % [#/micron/m^3 , microns] - gamma particle size distribution
        %                 lic(zz) = trapz( r , 4/3 * pi * rho_ice * (r*1e-6).^3 .* nr);                % g/m^3 - grams of ice per meter cubed of air
        %
        %             end
        %
        %         end
        % ------------------------------------------------------------------
        % ------------------------------------------------------------------


        % create the ice cloud file name
        if index==0
            fileName{nn} = ['ic_rtop',num2str(round(re(end,nn))),'_rbot',num2str(round(re(1,nn))),'_T',num2str(round(tau_c(nn))),...
                '_', distribution_type,'_nn',num2str(nn), '.DAT'];
        elseif index>0
            fileName{nn} = ['ic_rtop',num2str(round(re(end,nn))),'_rbot',num2str(round(re(1,nn))),'_T',num2str(round(tau_c(nn))),...
                '_', distribution_type,'_nn',num2str(index), '.DAT'];
        end



    elseif nLayers==1 && strcmp(distribution_type, 'gamma')==true

        % ----------------------------------------------------------------
        % --- If there is one cloud layer, this is a homogensous cloud ---
        % ----------------------------------------------------------------

        %Nc = tau_c(nn)./(pi*(H(nn)*1e3)*Qext(nn).*(re(nn)*1e-6).^2);                 % m^(-3) - number concentration

        % Compute  ice Content
        %lic = 4/3 * pi * rho_ice * (re(nn)*1e-6).^3 .* Nc;                  % g/m^3 - grams of ice per meter cubed of air

        % Compute  ice Content
        iwc = (tau_c(nn) *rho_ice)./(ext_bulk_coeff_per_IWC(nn) .* H(nn));                           % g/m^3 - grams of ice per meter cubed of air


        % create the ice cloud file name
        if index==0
            fileName{nn} = ['ic_r',num2str(round(re(nn))),'_T',num2str(round(tau_c(nn))),'_', distribution_type,...
                '_nn',num2str(nn), '.DAT'];
        elseif index>0

            fileName{nn} = ['ic_r',num2str(round(re(nn))),'_T',num2str(round(tau_c(nn))),'_', distribution_type,...
                '_nn',num2str(index), '.DAT'];
        end




    elseif nLayers==1 && strcmp(distribution_type, 'mono')==true

        % ----------------------------------------------------------------
        % --- If there is one cloud layer, this is a homogensous cloud ---
        % ----------------------------------------------------------------

        %Nc = tau_c(nn)./(pi*(H(nn)*1e3)*Qext(nn).*(re(nn)*1e-6).^2);                 % m^(-3) - number concentration

        % Compute  ice Content
        %iwc = 4/3 * pi * rho_ice * (re(nn)*1e-6).^3 .* Nc;                  % g/m^3 - grams of ice per meter cubed of air

        % Compute  ice Content
        iwc = 4/3 * (re(nn)*1e-6) * rho_ice * tau_c(nn)./...
            (Qext(nn) * (H(nn)*1e3));                                                % g/m^3 - grams of ice per meter cubed of air

        % when assuming a di
        % create the ice cloud file name
        if index==0
            fileName{nn} = ['ic_r',num2str(round(re(nn))),'_T',num2str(round(tau_c(nn))),'_', distribution_type,...
                '_nn',num2str(nn), '.DAT'];
        elseif index>0

            fileName{nn} = ['ic_r',num2str(round(re(nn))),'_T',num2str(round(tau_c(nn))),'_', distribution_type,...
                '_nn',num2str(index), '.DAT'];
        end





    end



    % ------------------------------------------------------------
    % ----------------- WE NEED TO APPEND ZEROS ------------------
    % ------------------------------------------------------------

    % Wherever the cloud is, there needs to be zeros at the cloud top altitude,
    % and below the cloud bottom altitude. This information tells LibRadTran
    % where the boundaries of the cloud are

    % both the effective radius and the IWC need zeros on either boundary,
    % unless if the cloud is at the surface

    if (size(re,1)==1 || size(re,2)==1) && strcmp(vert_homogeneous_str, 'vert-non-homogeneous')==true

        if z(1)==0
            % If true, then the cloud starts at the surface and we only append
            % zeros above the cloud
            re_2write = [re; 0];
            iwc_2write = [iwc; 0];
            z_2write = z;

        else
            % In this case, we need zeros below the cloud bottom, and at cloud
            % top
            z_2write = [0; z];                 % create a value at the surface where the cloud parameters go to zero
            re_2write = [0; re; 0];
            iwc_2write = [0; iwc; 0];

        end

    elseif (size(re,1)==1 || size(re,2)==1) && strcmp(vert_homogeneous_str, 'vert-homogeneous')==true

        if z(1)==0
            % If true, then the cloud starts at the surface and we only append
            % zeros above the cloud
            re_2write = [re(nn); 0];
            iwc_2write = [iwc; 0];
            z_2write = z;

        else
            % In this case, we need zeros below the cloud bottom, and at cloud
            % top
            z_2write = [0; z];                 % create a value at the surface where the cloud parameters go to zero
            re_2write = [0; re(nn); 0];
            iwc_2write = [0; iwc; 0];

        end

    elseif (size(re,1)>1 && size(re,2)>1)

        % Cloud top height defines the altitude where there is no cloud.

        % if the minimum z value is 0 then the cloud is at the surface
        if z(1)==0
            % then we only append zeros above the cloud
            re_2write = [re(:,nn); 0];
            iwc_2write = [iwc; 0];
            z_2write = z;

        else
            % Then we need zeros on either end
            z_2write = [0; z];
            re_2write = [0; re(:,nn); 0];
            iwc_2write = [0; iwc; 0];

        end

    end


    % ------------------------------------------------------------
    % ---------------------- WRITE ic FILE -----------------------
    % ------------------------------------------------------------



    % Create the ice cloud file
    fileID = fopen([ice_cloud_folder_path,fileName{nn}], 'w');

    % fprintf writes lines in our text file from top to botom
    % ic.DAT files are written with the higher altitudes at the top, and the
    % surface at the bottom

    % to write column vectors in a text file, we have to store them as row
    % vectors

    toWrite = [flipud(z_2write)'; flipud(iwc_2write)'; flipud(re_2write)'];

    % Create the opening comment lines of the ic.DAT file

    fprintf(fileID, '%s %10s %7s %8s \n','#','z','IWC','R_eff');
    fprintf(fileID, '%s %10s %7s %8s \n','#','(km)','(g/m^3)','(micron)');

    % Write in the data
    fprintf(fileID,'%12.3f %7.4f %8.3f \n', toWrite);
    fclose(fileID);


end


end