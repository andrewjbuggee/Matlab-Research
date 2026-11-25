function P = interpsonde_30min_profiles(folder)
% INTERPSONDE_30MIN_PROFILES  30-min mean temperature and pressure profiles from ARM INTERPSONDE daily files.
%
%   P = interpsonde_30min_profiles(folder)
%
%   Reads every epcinterpolatedsondeM1.c1.*.nc file in folder (the variable subset
%   temp, bar_pres, rh downloaded with the Python tool EPCAPE/download_data,
%   product thermo_interpsonde_M1) and averages each height over 30-min windows
%   [t, t + 30 min), the resolution of s_ARG00 in Dedrick et al. (2026, supplement
%   Table S7) and of the posted w* file.
%
%   INTERPSONDE is ARM's interpolated-sonde value-added product: radiosonde
%   profiles interpolated in time onto a fixed 1-min x 332-height grid, scaled to the
%   microwave-radiometer water vapour and tied to the surface met station (input
%   datastreams in the file's global attributes: griddedsonde, mwrlos, met; ARM
%   DOI 10.5439/1095316; Fairless et al. 2016, as cited by Dedrick et al. 2026, p. 1065).
%
%   QC: samples with any bit set in qc_temp / qc_bar_pres (below valid_min, above
%   valid_max, or missing in the input) are set to NaN before averaging.
%
%   Heights: the file gives `height` in km without saying above what. The first
%   level, 0.007 km, equals the site altitude alt = 7 m, which suggests heights
%   above mean sea level (my inference, not stated in the file). The 7-m question
%   is below the 20-m grid spacing.
%
%   Output struct
%     P.time      start of each 30-min window (datetime column, UTC), only windows with data
%     P.height_m  height of each level, m (as in the file, see above), 1 x 332
%     P.alt_m     site altitude from the file, m
%     P.T_K       mean temperature, K (windows x levels)
%     P.p_Pa      mean pressure, Pa (windows x levels)

files = dir(fullfile(folder, 'epcinterpolatedsondeM1.c1.*.nc'));
assert(~isempty(files), 'no epcinterpolatedsondeM1.c1 files in %s', folder);
files = sort({files.name});

t_all = NaT(0, 1); T_all = zeros(0, 0); p_all = zeros(0, 0);
for f = 1:numel(files)
    fn = fullfile(folder, files{f});
    % time is "seconds since <day>" (the file's units attribute); base_time is not in the subset,
    % so the day comes from the file name: epcinterpolatedsondeM1.c1.YYYYMMDD.hhmmss.nc
    day = datetime(regexp(files{f}, '\.(\d{8})\.', 'tokens', 'once'), 'InputFormat', 'yyyyMMdd');
    t = day + seconds(double(ncread(fn, 'time')));
    if f == 1
        P.height_m = 1000 * double(ncread(fn, 'height'))';   % km -> m
        P.alt_m = double(ncread(fn, 'alt'));
    end
    T_C = double(ncread(fn, 'temp'))';      % (time x height) after the transpose; degC
    p_kPa = double(ncread(fn, 'bar_pres'))';  % kPa
    T_C(double(ncread(fn, 'qc_temp'))' ~= 0) = NaN;
    p_kPa(double(ncread(fn, 'qc_bar_pres'))' ~= 0) = NaN;

    % 30-min means of each level: window index of every 1-min sample
    w = floor(minutes(t - day) / 30);       % 0..47
    n_w = 48;
    T_w = NaN(n_w, numel(P.height_m)); p_w = T_w;
    for k = 0:n_w - 1
        in = w == k;
        if any(in)
            T_w(k + 1, :) = mean(T_C(in, :), 1, 'omitnan');
            p_w(k + 1, :) = mean(p_kPa(in, :), 1, 'omitnan');
        end
    end
    t_all = [t_all; day + minutes(30 * (0:n_w - 1)')];   %#ok<AGROW>
    T_all = [T_all; T_w + 273.15];                         %#ok<AGROW> degC -> K
    p_all = [p_all; p_w * 1000];                           %#ok<AGROW> kPa -> Pa
end

has = any(isfinite(T_all), 2) & any(isfinite(p_all), 2);
P.time = t_all(has);
P.T_K = T_all(has, :);
P.p_Pa = p_all(has, :);

end
