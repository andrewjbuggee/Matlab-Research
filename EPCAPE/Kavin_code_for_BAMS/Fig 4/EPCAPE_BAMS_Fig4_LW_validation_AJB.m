%% EPCAPE_BAMS_Fig4_LW_validation_AJB  -- does RADFLUX reproduce Dan's LW files?
% BAMS Fig. 4 (EPCAPE_BAMS_Fig4_AJB_v1.m) draws the monthly downwelling LW
% line from the 1-min ARM RADFLUX file instead of Dan's <MON>_LW_Stepplot.csv
% files, because Dan's FEB file pools Feb 2023 and Feb 2024 and cannot give
% two February points. This script shows that the RADFLUX recipe reproduces
% Dan's numbers, so the swap changes only the February points.
%
% RECIPE (same as Fig. 4 v1, extended to all five statistics in Dan's files)
%   1. RADFLUX downwelling_longwave, 1-min, qc_downwelling_longwave == 0 only.
%   2. Hourly means, hour-start stamps, PST (UTC-8, no daylight saving).
%   3. For each calendar month and each hour of the day (0-23), the Min,
%      25th percentile, Median, 75th percentile and Max of those hourly means.
%      Percentiles: linear interpolation between order statistics (= Excel
%      PERCENTILE.INC / numpy 'linear'; local function prctileLinear, because
%      R2024b prctile has no 'inclusive' method). MATLAB's default prctile is
%      also computed and printed for comparison.
%   4. Monthly value plotted in Fig. 4 = mean over the 24 hours of the hourly
%      medians (the statistic Kavin averaged from Dan's files).
%   Dan's FEB file is compared with RADFLUX Feb 2023 + Feb 2024 POOLED, which
%   is what his file contains; Fig. 4 then splits them.
%
% DAN'S FILES: <MON>_LW_Stepplot.csv, columns Hour, Min, 25%ile, Median,
%   75%ile, Max (24 rows, W m^-2). How Dan made them is not documented; the
%   match below is the evidence for the recipe.
%
% FIGURE (EPCAPE_BAMS_Fig4_LW_validation_AJB.png/.pdf)
%   a  1:1 plot, all 12 months x 24 hours x 5 statistics (1440 pairs)
%   b  difference (RADFLUX - Dan) of the hourly medians vs hour, per month
%   c  monthly LW as plotted in Fig. 4: Dan vs RADFLUX, plus RADFLUX Feb 2023
%      and Feb 2024 separately
%   d  diurnal cycle of the hourly median: Dan (lines) vs RADFLUX (dots)
%
% Andrew J. Buggee, Scripps Institution of Oceanography, 2026-10-14

clear; close all

%% 1. Settings

danDir_andrew     = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/EPCAPE/Kavin_code_for_BAMS/Annual Echo Plot';
radfluxDir_andrew = '/Users/andrewbuggee/Documents/VS_CODE/Python-Research/EPCAPE/data/processed';
if isfolder(danDir_andrew)
    danDir = danDir_andrew;  radfluxDir = radfluxDir_andrew;   % Andrew's machine
else
    danDir = pwd;            radfluxDir = pwd;                 % anyone else
end
% Outputs (PNG/PDF) go next to this script, in the 'Fig 4' folder; the input
% data stay in 'Annual Echo Plot' (Kavin's Fig. 5 code reads them too).
outDir = fileparts(mfilename('fullpath'));
if isempty(outDir), outDir = pwd; end   % e.g. code pasted into the Command Window
exportFigs = true;

tz = '-08:00';                                     % PST, as in Fig. 4 / Fig. 5 v2
monthNames = ["JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"];
statNames  = ["Min","25th pct","Median","75th pct","Max"];

%% 2. Dan's monthly hour-of-day statistics: dan(hour, stat, calendar month)

dan = nan(24, 5, 12);                              % [W m^-2]
for mo = 1:12
    M = readmatrix(fullfile(danDir, monthNames(mo) + "_LW_Stepplot.csv"), 'NumHeaderLines',1);
    % column 1 = hour 0..23; columns 2-6 = Min, 25%ile, Median, 75%ile, Max
    dan(M(:,1)+1, :, mo) = M(:,2:6);
end

%% 3. RADFLUX: 1-min -> hourly means (PST)

radFile = fullfile(radfluxDir,'radflux_M1_20230215_20240214.nc');
tRaw  = double(ncread(radFile,'time'));            % seconds since 1970-01-01 UTC
tRad  = datetime(tRaw,'ConvertFrom','posixtime','TimeZone','UTC');
tRad.TimeZone = tz;
lw    = double(ncread(radFile,'downwelling_longwave'));    % [W m^-2]
qc    = ncread(radFile,'qc_downwelling_longwave');
fprintf('RADFLUX LW: %d of %d 1-min samples fail QC (qc ~= 0) and are removed\n', ...
        nnz(qc ~= 0 & ~isnan(lw)), nnz(~isnan(lw)));
lw(qc ~= 0) = NaN;
hourly = retime(timetable(tRad, lw), 'hourly', 'mean');    % hour-start stamps
moH = month(hourly.tRad);
hrH = hour(hourly.tRad);

%% 4. Same statistics from RADFLUX: rad(hour, stat, calendar month)
% Calendar month = Dan's grouping (Feb 2023 and Feb 2024 pooled).

rad      = nan(24, 5, 12);                         % linear (numpy/Excel) percentiles
radExact = nan(24, 5, 12);                         % prctile default ('exact')
for mo = 1:12
    for h = 0:23
        v = hourly.lw(moH == mo & hrH == h);
        v = v(~isnan(v));
        if isempty(v), continue, end
        rad(h+1,:,mo)      = [min(v), prctileLinear(v,[25 50 75]), max(v)];
        radExact(h+1,:,mo) = [min(v), prctile(v,[25 50 75]),                      max(v)];
    end
end

% Monthly Fig. 4 value: mean over hours of the hourly medians.
danMonthly = squeeze(mean(dan(:,3,:),1,'omitmissing'));    % [12 x 1]
radMonthly = squeeze(mean(rad(:,3,:),1,'omitmissing'));

% Feb 2023 and Feb 2024 separately (what Fig. 4 plots).
febYear = nan(1,2);
for y = 1:2
    yr = 2022 + y;
    med = nan(24,1);
    for h = 0:23
        v = hourly.lw(moH == 2 & year(hourly.tRad) == yr & hrH == h);
        med(h+1) = median(v,'omitmissing');
    end
    febYear(y) = mean(med,'omitmissing');
end

%% 5. Agreement statistics

d      = rad - dan;                                % [24 x 5 x 12]
dExact = radExact - dan;
fprintf('\nRADFLUX - Dan, all 12 months x 24 hours (W m^-2):\n');
fprintf('%-10s %8s %8s %10s %10s\n','statistic','bias','RMSE','max|diff|','RMSE exact');
for s = 1:5
    x  = d(:,s,:);       x = x(~isnan(x));
    xe = dExact(:,s,:);  xe = xe(~isnan(xe));
    fprintf('%-10s %8.3f %8.3f %10.3f %10.3f\n', statNames(s), mean(x), ...
            sqrt(mean(x.^2)), max(abs(x)), sqrt(mean(xe.^2)));
end
monthlyTbl = table(monthNames(:), round(danMonthly,2), round(radMonthly,2), ...
                   round(radMonthly-danMonthly,2), 'VariableNames', ...
                   {'Month','Dan_Wm2','RADFLUX_Wm2','diff_Wm2'});
disp(monthlyTbl)
fprintf('RADFLUX Feb 2023 = %.1f, Feb 2024 = %.1f W m^-2 (Dan FEB, pooled = %.1f)\n', ...
        febYear(1), febYear(2), danMonthly(2));

%% 6. Figure

fig = figure('Visible','on','Units','pixels','Position',[50 50 1000 780]);
resetToFactoryDefaults(fig);                       % undo startup.m dark theme
fig.Color = 'w';
set(fig,'DefaultAxesFontSize',11,'DefaultLegendFontSize',9,'DefaultTextFontSize',11);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

statColors = [0.35 0.35 0.35; 0 0.447 0.741; 0.85 0.325 0.098; 0.494 0.184 0.556; 0.466 0.674 0.188];
monthColors = turbo(12);

% ---- a: 1:1 ---------------------------------------------------------------
axA = nexttile(tl,1); hold(axA,'on'); box(axA,'on')
lims = [240 430];
plot(axA, lims, lims, 'k-', 'HandleVisibility','off')
for s = 1:5
    x = dan(:,s,:);  y = rad(:,s,:);
    ds = d(:,s,:);   ds = ds(~isnan(ds));
    scatter(axA, x(:), y(:), 14, statColors(s,:), 'filled', 'MarkerFaceAlpha',0.6, ...
            'DisplayName', sprintf('%s  (RMSE %.2f, max %.2f)', statNames(s), ...
            sqrt(mean(ds.^2)), max(abs(ds))));
end
axis(axA,'equal'); xlim(axA,lims); ylim(axA,lims)
xlabel(axA,'Dan''s LW\_Stepplot value [W m^{-2}]')
ylabel(axA,'RADFLUX recipe [W m^{-2}]')
legend(axA,'Location','southeast')
title(axA,'a   12 months \times 24 h \times 5 statistics','FontWeight','normal', ...
      'HorizontalAlignment','left','Units','normalized','Position',[0 1.01])

% ---- b: median differences vs hour ---------------------------------------
axB = nexttile(tl,2); hold(axB,'on'); box(axB,'on')
for mo = 1:12
    plot(axB, 0:23, d(:,3,mo), '.-', 'Color',monthColors(mo,:), 'DisplayName',monthNames(mo))
end
yline(axB,0,'k-','HandleVisibility','off')
xlim(axB,[-0.5 23.5]); set(axB,'XTick',0:3:23)
dMax = max(abs(d(:,3,:)),[],'all');
ylim(axB, 1.15*dMax*[-1 1])
xlabel(axB,'Hour [PST]'); ylabel(axB,'RADFLUX - Dan, hourly median [W m^{-2}]')
legend(axB,'Location','eastoutside')
title(axB,'b   difference of the hourly medians','FontWeight','normal', ...
      'HorizontalAlignment','left','Units','normalized','Position',[0 1.01])

% ---- c: monthly values used in Fig. 4 -------------------------------------
axC = nexttile(tl,3); hold(axC,'on'); box(axC,'on')
plot(axC, 1:12, danMonthly, 'ko', 'MarkerSize',9, 'LineWidth',1.2, 'DisplayName','Dan (calendar month)')
plot(axC, 1:12, radMonthly, '.', 'Color',statColors(3,:), 'MarkerSize',16, ...
     'DisplayName','RADFLUX, same pooling')
plot(axC, [1.8 2.2], febYear, 's', 'Color',statColors(2,:), 'MarkerFaceColor',statColors(2,:), ...
     'MarkerSize',6, 'DisplayName','RADFLUX Feb 2023 / Feb 2024 (Fig. 4)')
text(axC, 1.8, febYear(1), ' 2023', 'FontSize',9, 'VerticalAlignment','top')
text(axC, 2.2, febYear(2), ' 2024', 'FontSize',9, 'VerticalAlignment','bottom')
xlim(axC,[0.5 12.5]); set(axC,'XTick',1:12,'XTickLabel',extractBefore(monthNames,2) + lower(extractAfter(monthNames,1)))
ylabel(axC,'Mean of hourly medians [W m^{-2}]')
legend(axC,'Location','south')   % northwest hid the May point
title(axC,sprintf('c   monthly LW (max |diff| = %.2f W m^{-2})', max(abs(radMonthly-danMonthly))), ...
      'FontWeight','normal','HorizontalAlignment','left','Units','normalized','Position',[0 1.01])

% ---- d: diurnal medians ----------------------------------------------------
axD = nexttile(tl,4); hold(axD,'on'); box(axD,'on')
for mo = 1:12
    plot(axD, 0:23, dan(:,3,mo), '-', 'Color',monthColors(mo,:), 'LineWidth',1.2, ...
         'HandleVisibility','off')
    plot(axD, 0:23, rad(:,3,mo), '.', 'Color',monthColors(mo,:), 'MarkerSize',9, ...
         'HandleVisibility','off')
end
plot(axD, NaN, NaN, 'k-', 'LineWidth',1.2, 'DisplayName','Dan')
plot(axD, NaN, NaN, 'k.', 'MarkerSize',12, 'DisplayName','RADFLUX')
xlim(axD,[-0.5 23.5]); set(axD,'XTick',0:3:23)
xlabel(axD,'Hour [PST]'); ylabel(axD,'Hourly median LW [W m^{-2}]')
legend(axD,'Location','southwest')
title(axD,'d   diurnal cycle of the median (colours as in b)','FontWeight','normal', ...
      'HorizontalAlignment','left','Units','normalized','Position',[0 1.01])

if exportFigs
    exportgraphics(fig, fullfile(outDir,'EPCAPE_BAMS_Fig4_LW_validation_AJB.png'),'Resolution',200);
    exportgraphics(fig, fullfile(outDir,'EPCAPE_BAMS_Fig4_LW_validation_AJB.pdf'),'ContentType','vector');
end

%% Local functions

function p = prctileLinear(v, pct)
% prctileLinear  Percentiles by linear interpolation between sorted values at
% rank 1 + (n-1)*pct/100 (Hyndman & Fan 1996 type 7; numpy 'linear', Excel
% PERCENTILE.INC). MATLAB's prctile default uses rank 0.5 + n*pct/100 instead.
v = sort(v(:));
n = numel(v);
r = 1 + (n-1)*pct(:).'/100;          % fractional rank
lo = floor(r);  hi = min(lo+1, n);
p = v(lo).' + (r - lo).*(v(hi).' - v(lo).');
end

function resetToFactoryDefaults(fig)
% resetToFactoryDefaults  Give one figure MATLAB's factory graphics defaults
% (Andrew's startup.m sets a dark theme at the root level).
rootDefaults = fieldnames(get(groot,'default'));
for k = 1:numel(rootDefaults)
    name = rootDefaults{k};
    if startsWith(name,'defaultFigure','IgnoreCase',true), continue, end
    set(fig, name, get(groot, regexprep(name,'^default','factory','ignorecase')));
end
fig.Color = get(groot,'factoryFigureColor');
end
