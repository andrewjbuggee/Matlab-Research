function check_mdcs_slurm()
    fprintf(' Checking for MATLAB Distributed Computing Server (MDCS) with SLURM integration...\n\n');
    
    % Get all available cluster profiles
    profiles = parallel.clusterProfiles;
    
    if isempty(profiles)
        fprintf(' No cluster profiles found.\n');
        return;
    else
        fprintf(' Found %d cluster profile(s):\n', numel(profiles));
        disp(profiles');
    end

    % Check each profile
    for i = 1:numel(profiles)
        profileName = profiles{i};
        fprintf('\n--- Checking profile: %s ---\n', profileName);
        try
            c = parcluster(profileName);

            % Check if it's not 'local'
            if strcmp(profileName, 'local')
                fprintf('  This is the default local profile. Limited to one node.\n');
                continue;
            end

            % Try to infer SLURM/MDCS involvement
            schedulerInfo = '';
            if isprop(c, 'AdditionalProperties')
                props = c.AdditionalProperties;
                if isstruct(props)
                    propNames = fieldnames(props);
                    for j = 1:numel(propNames)
                        if contains(lower(propNames{j}), 'slurm')
                            schedulerInfo = ' Likely SLURM';
                        end
                    end
                end
            end

            if isa(c, 'parallel.cluster.Generic') || isa(c, 'parallel.Cluster')
                schedulerHint = ' Can submit jobs via batch or parpool.';
            else
                schedulerHint = ' Not a generic cluster type.';
            end

            fprintf(' Profile loaded successfully.\n');
            fprintf('   Cluster class: %s\n', class(c));
            fprintf('   SLURM Hint: %s\n', ternary(~isempty(schedulerInfo), schedulerInfo, 'No clear SLURM markers'));
            fprintf('   %s\n', schedulerHint);

        catch ME
            fprintf('Could not load profile "%s": %s\n', profileName, ME.message);
        end
    end

    fprintf('\n Done.\n');
end

function out = ternary(condition, valTrue, valFalse)
    if condition
        out = valTrue;
    else
        out = valFalse;
    end
end
