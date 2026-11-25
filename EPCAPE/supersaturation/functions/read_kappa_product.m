function K = read_kappa_product(file_path)
% READ_KAPPA_PRODUCT  Read the posted EPCAPE hygroscopicity file (bb2743661p_37_1_kappa.nc) into a timetable.
%
%   K = read_kappa_product(file_path)
%
%   The file is one item of the UC San Diego Library collection "Aerosol
%   Quantification Products at Mt. Soledad and Scripps Pier during EPCAPE 2023-24"
%   (https://doi.org/10.6075/J0WM1DSC). Its title attribute is "EPCAPE kappa_calc
%   retrieved from Mt. Soledad and Scripps Pier". It holds four series of the
%   composition-based hygroscopicity parameter kappa (unitless), one value per
%   composition sample:
%     kappa_pier_01  Scripps Pier, organic kappa = 0.10 (ACSM, ~30-min samples)
%     kappa_sol_01   Mt. Soledad, organic kappa = 0.10 (AMS, 5-min samples)
%     kappa_sol_f44  Mt. Soledad, organic kappa from f44
%     kappa_sol_OC   Mt. Soledad, organic kappa from the O/C ratio
%   (descriptions from the file's long_name attributes; sampling intervals are my
%   measurement from the time stamps). The method is the volume-weighted mixing
%   rule kappa = sum_i eps_i kappa_i (Petters and Kreidenweis 2007, p. 1962, eq. 7)
%   with the component densities and kappas of Dedrick et al. (2026), supplement
%   Text S4, eq. S1 and Table S2 (organics: 1.50 g cm-3, kappa = 0.10).
%   Whether this is the "kcalc.nc" that Jeramy Dedrick used for s_kappa is not
%   stated anywhere I have found.
%
%   The file has one time axis for all four series; each series is NaN at the
%   other instrument's sample times.
%
%   Output: timetable with the four series as variables (names as in the file),
%   row times in UTC (the file's time is "seconds since 1970-01-01", i.e. POSIX
%   time, which is UTC by definition).

t = datetime(ncread(file_path, 'time'), 'ConvertFrom', 'posixtime');   % UTC
names = {'kappa_pier_01', 'kappa_sol_01', 'kappa_sol_f44', 'kappa_sol_OC'};
K = timetable('RowTimes', t);
for k = 1:numel(names)
    K.(names{k}) = double(ncread(file_path, names{k}));
    K.Properties.VariableDescriptions{names{k}} = ncreadatt(file_path, names{k}, 'long_name');
end
K.Properties.DimensionNames{1} = 'time';

end
