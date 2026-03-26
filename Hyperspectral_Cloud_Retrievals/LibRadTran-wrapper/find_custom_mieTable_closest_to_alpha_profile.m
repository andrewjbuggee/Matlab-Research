


function out = find_custom_mieTable_closest_to_alpha_profile(alpha_profile, use_35_or_50, which_computer)




% load the set of VOCALS-REx in-situ observations
if strcmp(which_computer, 'anbu8374')==true


    % --- define the directory for where the custom pre-computed mie tables are ---
    custom_mie_tables_dir = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/wc/custom_mieTables/';


elseif strcmp(which_computer, 'andrewbuggee')==true



    % --- define the directory for where the custom pre-computed mie tables are ---
    custom_mie_tables_dir = '/Users/andrewbuggee/Documents/libRadtran-2.0.6/data/wc/custom_mieTables/';


elseif strcmp(which_computer, 'curc')==true



    % --- define the directory for where the custom pre-computed mie tables are ---
    custom_mie_tables_dir = '/projects/anbu8374/software/libRadtran-2.0.5/data/wc/custom_mieTables/';

end



% first, read all the filenames in the custom mie tables dir
% mie_table_filenames = dir(fullfile(custom_mie_tables_dir, '*.cdf'));

% % step through each filename and extract the alpha value, which is in
% % the filename
% % Extract alpha values from the filenames and store them
% mie_table_alpha_values = zeros( length(mie_table_filenames), 1);
%
% if use_35_or_50 == 35
%
%     for ff = 1:length(mie_table_filenames)
%
%         alpha_str = extractBetween(mie_table_filenames(ff).name, 'wc_mieTable_gamma_rEff_1-35microns_gammaDist_alpha_', '.cdf');
%
%         mie_table_alpha_values(ff) = str2double(alpha_str{1}); % Assuming alpha is in the filename
%
%     end
%
% elseif use_35_or_50 == 50
%
%     for ff = 1:length(mie_table_filenames)
%
%         alpha_str = extractBetween(mie_table_filenames(ff).name, 'wc_mieTable_gamma_rEff_1-50microns_gammaDist_alpha_', '.cdf');
%
%         mie_table_alpha_values(ff) = str2double(alpha_str{1}); % Assuming alpha is in the filename
%
%     end
%
%
% end



% The custom mie tables span the following alpha distribution width
% values
mie_table_alpha_values = [1  2	3	4	5	6	7	8	9	10	11	12	13	14	15,...
    16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34	35,...
    36	37	38	39	40	45	50	55	60	80	100	125	150];

% update the effective variance assumption for each cloud layer using
% the fit using log-normal fit of alpha values
% the distribution variance fits start at cloud base and move towards
% cloud top
out.closest_table_alpha_to_true = zeros(length(alpha_profile), 1);
out.mie_table_filename = cell(length(alpha_profile), 1);


% Lastly, take a mean of the vertical profile of effective variance
% This is the value that will be used in all calculations, since only a
% single mie table can be used in the libRadtran uvSpec calculations
mean_distribution_var = mean(alpha_profile);


if use_35_or_50 == 35

    for ll = 1:length(alpha_profile)

        % for each value, find the pre-computed mie table with the closest
        % alpha value. Save this filename
        % find the closest alpha value in the pre-computed mie tables
        [~, idx_min] = min( abs( alpha_profile(ll) - mie_table_alpha_values ) );

        % store the closest alpha value for use in the radiative transfer
        out.closest_table_alpha_to_true(ll) = mie_table_alpha_values(idx_min);

        out.mie_table_filename{ll} = ['wc_mieTable_gamma_rEff_1-35microns_gammaDist_alpha_',...
            num2str(mie_table_alpha_values(idx_min)), '.cdf'];


    end


    % find the pre-computed mie table with the closest
    % alpha value to the mean value and save it
    % find the closest alpha value in the pre-computed mie tables
    [~, idx_mean_min] = min( abs( mean_distribution_var - mie_table_alpha_values ) );

    % store the closest alpha value for use in the radiative transfer
    out.closest_table_alpha_to_mean_alpha = mie_table_alpha_values(idx_min);
    % store the filename of the closest alpha value
    out.mie_table_filename_closest_to_mean = [custom_mie_tables_dir, 'wc_mieTable_gamma_rEff_1-35microns_gammaDist_alpha_',...
        num2str(mie_table_alpha_values(idx_mean_min)), '.cdf'];







elseif use_35_or_50 == 50


    for ll = 1:length(alpha_profile)

        % for each value, find the pre-computed mie table with the closest
        % alpha value. Save this filename
        % find the closest alpha value in the pre-computed mie tables
        [~, idx_min] = min( abs( alpha_profile(ll) - mie_table_alpha_values ) );

        % store the closest alpha value for use in the radiative transfer
        out.closest_table_alpha_to_true(ll) = mie_table_alpha_values(idx_min);

        out.mie_table_filename{ll} = ['wc_mieTable_gamma_rEff_1-50microns_gammaDist_alpha_',...
            num2str(mie_table_alpha_values(idx_min)), '.cdf'];


    end




    % find the pre-computed mie table with the closest
    % alpha value to the mean value and save it
    % find the closest alpha value in the pre-computed mie tables
    [~, idx_mean_min] = min( abs( mean_distribution_var - mie_table_alpha_values ) );

    % store the closest alpha value for use in the radiative transfer
    out.closest_table_alpha_to_mean_alpha = mie_table_alpha_values(idx_min);
    % store the filename of the closest alpha value
    out.mie_table_filename_closest_to_mean = [custom_mie_tables_dir, 'wc_mieTable_gamma_rEff_1-50microns_gammaDist_alpha_',...
        num2str(mie_table_alpha_values(idx_mean_min)), '.cdf'];


end











end