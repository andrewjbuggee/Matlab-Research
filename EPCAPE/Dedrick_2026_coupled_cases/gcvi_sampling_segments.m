%% Continuous GCVI sampling periods at Mt. Soledad, from the 1-s enhancement-factor file
%
% DESCRIPTION:
%   The ground-based counterflow virtual impactor (GCVI; Brechtel
%   Manufacturing) at Mt. Soledad sampled cloud-droplet residuals when the
%   visibility was < 5000 m and no precipitation was detected (Dedrick et al.
%   2026, p. 1064). Its 1-s record in the UC San Diego Library item
%   "Aerosol Quantification Products" (file EPCAPE_CVI_EF.csv, local name
%   bb2743661p_36_1_CVI_EF.csv) holds the inlet's flows, airspeed, cut size
%   and the enhancement factor (EF) that corrects residual concentrations
%   for the inlet's concentration of droplets (Shingler et al. 2012, Eq. 2).
%   The record is written only while the GCVI runs, and EF <= 0 (or inf)
%   marks start-up and idle states (EF ~ 5 when sampling).
%
%   A GCVI "segment" is a run of 1-s samples with a finite EF > EfMin and no
%   gap longer than MaxGap_s. This is the Python gcvi_segments
%   (comparisons/seasonal_averages/quantities.py) in MATLAB. Its segments
%   are the raw material of the Mt. Soledad cloud events: Dedrick et al.
%   start from the GCVI periods longer than 1 h (154 of them; supplement
%   Table S4).
%
%   FILE FORMAT: 9 comment lines ("# ..."), a header line, then
%   "yyyy-MM-dd HH:mm:ss,lfeflow,lfetemp,airspeed,cntflow,cnttemp,EF,cutsize"
%   with Windows (CRLF) line ends; times are UTC (GMT).
%
% INPUTS:
%   ef_file - path to bb2743661p_36_1_CVI_EF.csv
%   Optional name-value pairs:
%   'EfMin'     a sample counts as sampling if EF > EfMin (default 1)
%   'MaxGap_s'  a gap longer than this starts a new segment (default 60 s)
%
% OUTPUTS:
%   segs - table, one row per segment:
%       start_utc, end_utc (last sample + 1 s), duration_h, mean_EF,
%       median_cut_um (median GCVI cut size, um)
%   gcvi - struct with the 1-s samples that count as sampling: time, EF,
%          cut_um, airspeed_ms
%
% REFERENCES:
%   Shingler, T., and coauthors, 2012: Characterisation and airborne
%     deployment of a new counterflow virtual impactor inlet. Atmos. Meas.
%     Tech., 5, 1259-1269, doi:10.5194/amt-5-1259-2012.
%   Data: Berta, V., S. Han, A. Osawa, J. Dedrick, N. Maneenoi, and L. M.
%     Russell, 2025: Aerosol quantification products at Mt. Soledad and
%     Scripps Pier during EPCAPE 2023-24. UC San Diego Library Digital
%     Collections, doi:10.6075/J0WM1DSC (as cited by Dedrick et al. 2026;
%     readme entry EPCAPE_CVI_EF.csv).
%   Dedrick, J. L., and coauthors, 2026: J. Atmos. Sci., 83, 1061-1082,
%     doi:10.1175/JAS-D-26-0004.1; supplement Text S6, S7 and Table S4.
%
% By Andrew John Buggee
%%

function [segs, gcvi] = gcvi_sampling_segments(ef_file, varargin)


p = inputParser;
addParameter(p, 'EfMin', 1);
addParameter(p, 'MaxGap_s', 60);
parse(p, varargin{:});
opt = p.Results;



% -------------------------------------------------------------------
% ------------------- Read the 1-s GCVI record ----------------------
% -------------------------------------------------------------------

fid = fopen(ef_file, 'r');
if fid < 0
    error('gcvi_sampling_segments: cannot open %s', ef_file);
end
cleaner = onCleanup(@() fclose(fid));

% 9 comment lines + 1 header line. Some fields are empty (",,", e.g. the EF
% while the counterflow is off), so the line is split at commas and an
% empty field becomes NaN. EndOfLine is detected automatically, so the CRLF
% line ends are fine.
C = textscan(fid, '%{yyyy-MM-dd HH:mm:ss}D %f %f %f %f %f %f %f', ...
             'Delimiter', ',', 'EmptyValue', NaN, 'HeaderLines', 10);
clear cleaner

t        = C{1};
airspeed = C{4};       % wind-tunnel airspeed (m s-1)
EF       = C{7};       % enhancement factor (corrected)
cut_um   = C{8};       % GCVI cut size (um)
if numel(EF) ~= numel(t)
    error('gcvi_sampling_segments: the file could not be read to the end (stopped near %s)', string(t(end)));
end

% sampling = finite EF above EfMin
on = isfinite(EF) & EF > opt.EfMin;
[gcvi.time, order] = sort(t(on));
EF_on = EF(on);        gcvi.EF          = EF_on(order);
cut_on = cut_um(on);   gcvi.cut_um      = cut_on(order);
as_on = airspeed(on);  gcvi.airspeed_ms = as_on(order);



% -------------------------------------------------------------------
% -------------- Split the sampling samples into segments -----------
% -------------------------------------------------------------------

gap_s = seconds(diff(gcvi.time));
new_segment = [true; gap_s > opt.MaxGap_s];   % a long gap starts a new segment
seg_id = cumsum(new_segment);
n_seg = seg_id(end);

first = find(new_segment);
last  = [first(2:end) - 1; numel(gcvi.time)];

start_utc     = gcvi.time(first);
end_utc       = gcvi.time(last) + seconds(1);          % each sample covers 1 s
duration_h    = hours(end_utc - start_utc);
mean_EF       = accumarray(seg_id, gcvi.EF, [n_seg 1], @mean);
median_cut_um = accumarray(seg_id, gcvi.cut_um, [n_seg 1], @median);

segs = table(start_utc, end_utc, duration_h, mean_EF, median_cut_um);


end
