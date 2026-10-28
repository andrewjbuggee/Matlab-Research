%% Read the ~1-s Mt. Soledad visibility record used to define cloud at the GCVI inlet
%
% DESCRIPTION:
%   Mt. Soledad (250 m MSL) had a visibility sensor at the aerosol inlets.
%   Its ~1-s record decided when the inlet switched to the ground-based
%   counterflow virtual impactor (GCVI; visibility < 5000 m and no detectable
%   precipitation) and, with the FM-120 liquid water content, defines the
%   5-min "cloud" intervals of Dedrick et al. (2026; supplement Table S4:
%   visibility < 1 km, and no 1-s value above 1 km in the interval). The
%   sensor's visibility is the distance at which the atmosphere's extinction
%   reduces contrast to ~5 %; in cloud it falls to tens to hundreds of
%   metres because the droplets scatter strongly (inside the cloud, not
%   below it, so this measures whether Mt. Soledad is IN cloud).
%
%   FILES (data/processed of the Python project, from the UCSD group):
%     msd_visibility_20230224-20230719.csv and msd_visibility_20230719-20240214.csv
%   each line: "yyyy-MM-dd HH:mm:ss",visibility_m (UTC; header line 1).
%   The two files overlap by about 10 h on 19 Jul 2023; repeated time
%   stamps are kept once (first file first).
%
%   CLEANING (as visibility_on_intervals in the Python quantities.py):
%   samples <= 0 m and > 48,270 m (the sensor's maximum, 30 statute miles)
%   are spikes or fill values and are set to NaN.
%
% INPUTS:
%   files - string array of the visibility CSV paths
%
% OUTPUTS:
%   vis - struct of column vectors sorted by time
%       time   datetime (UTC)
%       m      visibility (m), NaN where invalid
%
% REFERENCES:
%   Dedrick, J. L., and coauthors, 2026: J. Atmos. Sci., 83, 1061-1082,
%     doi:10.1175/JAS-D-26-0004.1 (p. 1064) and supplement Table S4.
%   Reader based on table_for_BAMS/readers/read_visibility_file.m.
%
% By Andrew John Buggee
%%

function vis = read_soledad_visibility(files)


VIS_SENSOR_MAX_M = 48270;    % 30 statute miles, the sensor's largest reported value

files = string(files);
t_all = cell(numel(files), 1);
v_all = cell(numel(files), 1);

for f = 1:numel(files)

    fid = fopen(files(f), 'r');
    if fid < 0
        error('read_soledad_visibility: cannot open %s', files(f));
    end
    cleaner = onCleanup(@() fclose(fid));   % closes the file even if an error occurs

    % "%f-%f-%f %f:%f:%f" reads the six parts of the quoted date and time
    % directly as numbers; readtable is far slower on ~10-15 million lines
    C = textscan(fid, '"%f-%f-%f %f:%f:%f",%f', 'HeaderLines', 1);
    clear cleaner

    t_all{f} = datetime(C{1}, C{2}, C{3}, C{4}, C{5}, C{6});
    v_all{f} = C{7};
end

t = vertcat(t_all{:});
v = vertcat(v_all{:});

% sort, and keep each time stamp once (unique keeps the first occurrence,
% i.e. the earlier file where the two overlap)
[vis.time, first] = unique(t, 'first');
vis.m = v(first);

% invalid samples -> NaN
vis.m(vis.m <= 0 | vis.m > VIS_SENSOR_MAX_M) = NaN;


end
