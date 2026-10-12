function vis = read_visibility_file(path)
% READ_VISIBILITY_FILE  Read one Mt. Soledad visibility file (~1-s samples).
%
%   vis = read_visibility_file(path)
%
%   File format (msd_visibility_*.csv), one line per sample:
%       "datetime (GMT)","visibility (m)"
%       "2023-02-24 00:00:01",40556
%   The two files hold ~10 and ~15 million lines. readtable is slow on files
%   this large, so textscan reads the date and time directly as numbers
%   (year, month, day, hour, minute, second); this takes a few seconds.
%
%   Output (struct of column vectors, sorted by time, NOTHING removed)
%     vis.time  datetime (UTC)
%     vis.m     visibility (m)
%   The record holds a few zeros and a few spikes far above the sensor's
%   maximum (48,270 m = 30 statute miles). Whether they are removed depends
%   on the calculation (see calc_rows31_32_FM120_hours.m).

fid = fopen(path, 'r');
if fid < 0
    error('Cannot open %s', path);
end
cleaner = onCleanup(@() fclose(fid));   % closes the file even if an error occurs

% "%f-%f-%f %f:%f:%f" reads the six parts of the quoted date/time;
% the quotes and the comma are matched as literal text.
C = textscan(fid, '"%f-%f-%f %f:%f:%f",%f', 'HeaderLines', 1);

vis.time = datetime(C{1}, C{2}, C{3}, C{4}, C{5}, C{6});
vis.m    = C{7};

% Make sure the samples are in time order ('sort' keeps the file order of equal times)
[vis.time, order] = sort(vis.time);
vis.m = vis.m(order);

end
