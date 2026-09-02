%% Convert a cell array of per-profile structs into a struct of cell arrays
%
% DESCRIPTION:
%   Two layouts for ensemble profiles are used across this code base:
%
%     per-profile layout   ensemble_profiles{nn}.re      (what the collection
%                          loops build, and what gets saved to disk)
%
%     struct-of-cells      ensemble_profiles.re{nn}      (what the older
%                          analysis and plotting code expects)
%
%   This converts the first into the second. It is idempotent: handing it a
%   struct of cell arrays returns it unchanged, so it is safe to call on either.
%
% INPUTS:
%   profiles - cell array of per-profile structs, a struct array of profiles, or
%              a struct of cell arrays (returned unchanged)
%
%   inputs   - (optional) the thresholds structure used to build the ensemble.
%              Attached as ep.inputs, which the plotting and titling code reads
%              for the LWC and Nc thresholds.
%
% OUTPUT:
%   ep - struct whose every field is a 1 x n_profiles cell array
%
% NOTES:
%   - Fields absent from a given profile (for example 're' on VOCALS-REx
%     profiles whose 2DC data did not conform) are filled with [] for that
%     profile, so every field spans all profiles and indexing stays aligned.
%
% SEE ALSO:
%   ensemble_vertical_statistics, ensemble_vertical_statistics_ORACLES,
%   plot_multiple_median_ensemble_LWC_re_NC_vs_norm_optical_depth
%
% By Andrew John Buggee
%%

function ep = ensembleCell2StructOfCells(profiles, inputs)

%% Normalize the container to a cell array of profile structs

if iscell(profiles)
    % already a cell array of per-profile structs

elseif isstruct(profiles) && isscalar(profiles) && ~isempty(fieldnames(profiles))

    % Could already be a struct of cell arrays. Check the first field.
    fn = fieldnames(profiles);

    if iscell(profiles.(fn{1}))
        % already in the target layout -- return unchanged, but attach inputs
        ep = profiles;
        if nargin > 1 && ~isfield(ep, 'inputs')
            ep.inputs = inputs;
        end
        return
    end

    % a single profile struct
    profiles = {profiles};

elseif isstruct(profiles)
    % struct array of profiles
    profiles = num2cell(profiles);

else
    error([newline, 'Unrecognized ensemble layout: expected a cell array of ',...
        'profile structs, a struct array, or a struct of cell arrays.', newline])
end

if isempty(profiles)
    error([newline, 'No profiles to convert.', newline])
end

n_prof = numel(profiles);

%% Collect the union of field names across every profile

all_fields = {};

for nn = 1:n_prof
    all_fields = [all_fields; fieldnames(profiles{nn})];    %#ok<AGROW>
end

all_fields = unique(all_fields, 'stable');

%% Build one cell array per field

for ff = 1:length(all_fields)

    fname = all_fields{ff};
    ep.(fname) = cell(1, n_prof);

    for nn = 1:n_prof

        if isfield(profiles{nn}, fname)
            ep.(fname){nn} = profiles{nn}.(fname);
        else
            % Field missing for this profile; keep the slot so indexing by
            % profile number stays aligned across every field.
            ep.(fname){nn} = [];
        end

    end

end

%% Attach the threshold structure when it was provided

if nargin > 1
    ep.inputs = inputs;
end

end
