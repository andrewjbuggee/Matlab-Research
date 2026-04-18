%% Compare two ensemble_profiles mat files
%
% Compares the two ensemble_profiles mat files produced by
% ensemble_vertical_statistics.m to:
%   (1) determine why the two files have a different number of profiles
%   (2) find the set of ensemble profiles common to both files
%   (3) find the set of ensemble profiles unique to each file
%
% Matching is done by fingerprinting each profile on its time vector
% (first time, last time, length). The time vector in each profile is
% seconds-since-midnight-UTC of the flight day, so combined with profile
% length and end time it uniquely identifies a profile across the 14
% source flight files.
%
% Andrew John Buggee

clear variables

%% Set paths to the two mat files

which_computer = whatComputer;

if strcmp(which_computer,'anbu8374')==true

    foldername = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/'];

elseif strcmp(which_computer,'andrewbuggee')==true

    foldername = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/'];

end

file_with_precip = 'ensemble_profiles_with_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25_04-Dec-2025.mat';
file_without_precip = 'ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25_drizzleLWP-threshold_5_10-Nov-2025.mat';


%% Load both mat files

A = load([foldername, file_with_precip]);       % 88 profiles (with precip)
B = load([foldername, file_without_precip]);    % 73 profiles (no precip)

ep_A = A.ensemble_profiles;
ep_B = B.ensemble_profiles;

n_A = length(ep_A);
n_B = length(ep_B);

fprintf('\nProfile counts:\n');
fprintf('  %s : %d profiles\n', file_with_precip, n_A);
fprintf('  %s : %d profiles\n', file_without_precip, n_B);
fprintf('  Difference (A - B) = %d\n\n', n_A - n_B);


%% Confirm that both files were built from the same set of source flights

fn_A = A.filename;
fn_B = B.filename;

if isequal(sort(fn_A), sort(fn_B))
    fprintf('Both files were generated from the same %d source flight files.\n\n', length(fn_A));
else
    fprintf('WARNING: the source filename lists differ between the two mat files!\n');
    only_in_A = setdiff(fn_A, fn_B);
    only_in_B = setdiff(fn_B, fn_A);
    fprintf('  Files only in A: %s\n', strjoin(only_in_A, ', '));
    fprintf('  Files only in B: %s\n\n', strjoin(only_in_B, ', '));
end


%% Build a fingerprint for each profile
%
% Each profile is uniquely identified by its time vector. We use:
%   [time(1), time(end), length(time)]
% rounded to the nearest 1e-3 seconds to avoid floating point mismatches.

fp_A = zeros(n_A, 3);
for nn = 1:n_A
    t = ep_A{nn}.time;
    fp_A(nn, :) = [round(t(1), 3), round(t(end), 3), numel(t)];
end

fp_B = zeros(n_B, 3);
for nn = 1:n_B
    t = ep_B{nn}.time;
    fp_B(nn, :) = [round(t(1), 3), round(t(end), 3), numel(t)];
end


%% Find overlapping profiles and profiles unique to each file

[common_fp, iA_common, iB_common] = intersect(fp_A, fp_B, 'rows', 'stable');

iA_unique = setdiff(1:n_A, iA_common);      % indices in A not found in B
iB_unique = setdiff(1:n_B, iB_common);      % indices in B not found in A

fprintf('Overlap (profiles in both files): %d\n', size(common_fp, 1));
fprintf('Unique to A (with_precip):        %d\n', length(iA_unique));
fprintf('Unique to B (without_precip):     %d\n\n', length(iB_unique));


%% Report the expected vs. actual difference
%
% If B is a strict subset of A (as the script logic suggests it should be),
% then length(iB_unique) should be 0 and length(iA_unique) should equal
% n_A - n_B = the number of drizzling profiles excluded from B.

if isempty(iB_unique)
    fprintf('B is a strict subset of A.\n');
    fprintf('The %d profiles unique to A are the ones excluded from B by the\n', length(iA_unique));
    fprintf('drizzle LWP threshold (> 5 g/m^2).\n\n');
else
    fprintf('B is NOT a strict subset of A. %d profiles in B are not in A.\n', length(iB_unique));
    fprintf('This means the upstream code or thresholds changed between the two runs\n');
    fprintf('(A was generated on 04-Dec-2025, B on 10-Nov-2025).\n\n');
end


%% Compute the drizzle/precipitation LWP for the profiles unique to A
%
% A profile in A but not in B is expected to have a 2DC LWP > 5 g/m^2
% (the drizzle threshold used when generating B). Report that LWP so we
% can confirm.

fprintf('Profiles unique to file A (with_precip) -- index in A, time(1), 2DC LWP:\n');
fprintf('%5s   %12s   %12s   %12s\n', 'idx_A', 'time(1) [s]', 'n_points', '2DC LWP [g/m^2]');

lwp_2DC_uniqueA = zeros(length(iA_unique), 1);

for kk = 1:length(iA_unique)

    nn = iA_unique(kk);
    p = ep_A{nn};

    % compute the 2DC liquid water path if the field exists
    if isfield(p, 'lwc_2DC') && ~isempty(p.lwc_2DC)
        lwp_2DC_uniqueA(kk) = abs(trapz(p.altitude, p.lwc_2DC));
    else
        lwp_2DC_uniqueA(kk) = NaN;
    end

    fprintf('%5d   %12.3f   %12d   %12.3f\n', ...
        nn, p.time(1), numel(p.time), lwp_2DC_uniqueA(kk));

end

fprintf('\n');

if all(~isnan(lwp_2DC_uniqueA))
    n_above = sum(lwp_2DC_uniqueA > 5);
    fprintf('Of the %d profiles unique to A, %d have 2DC LWP > 5 g/m^2 (the drizzle threshold).\n', ...
        length(iA_unique), n_above);
end


%% Save results for downstream use

profiles_unique_to_A = ep_A(iA_unique);
profiles_unique_to_B = ep_B(iB_unique);
profiles_common_fromA = ep_A(iA_common);
profiles_common_fromB = ep_B(iB_common);

save([foldername, 'ensemble_profile_comparison_', char(datetime("today")), '.mat'], ...
    'iA_common', 'iB_common', 'iA_unique', 'iB_unique', ...
    'profiles_unique_to_A', 'profiles_unique_to_B', ...
    'profiles_common_fromA', 'profiles_common_fromB', ...
    'file_with_precip', 'file_without_precip', 'lwp_2DC_uniqueA')

fprintf('\nResults saved to ensemble_profile_comparison_%s.mat\n', char(datetime("today")));
