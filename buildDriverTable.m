function drivers = buildDriverTable(epochs, workDir, opts)
%BUILDDRIVERTABLE Hourly solar/geomagnetic driver table (Sect. 2.2, 3.4).
%   drivers = BUILDDRIVERTABLE(epochs, workDir, ...)
%     epochs  : K x 1 datetime (UTC, hourly), aligned 1:1 to the GIM record
%     workDir : download/cache folder (default 'index_cache')
%
%   Output table (K rows) with
%     F107     daily 10.7-cm flux, interpolated to 1 h [sfu]        (GFZ)
%     F107bar  81-day running mean of F107 [sfu]  -- TRAILING mean by
%              default (revised manuscript, Sect. 2.2; reply to the PDF
%              review, comment 2: the centred mean needs 40 days of future
%              flux and is not available to a forecaster)
%     Kp       3-hourly planetary index, interpolated to 1 h        (GFZ)
%     Dst      hourly storm-time index [nT]                    (WDC Kyoto)
%   and, when 'mode' is 'operational' (Sect. 3.4), the latency-constrained
%   quantities that a forecaster has when issuing a 24-h forecast for
%   epoch t+24 at time t:
%     KpFc     NOAA SWPC day-1 Kp forecast valid for the 3-h interval
%              containing the epoch (issued the day before)
%     KpFcErr  KpFc - Kp (for reporting the forecast error, Sect. 2.2)
%   The persisted F10.7, the trailing mean and the last observed Dst
%   needed for the target-epoch block are taken from the same table by
%   buildTECDataset (it uses the row of the CURRENT epoch t).
%
%   Options
%     'mode'        'hindcast' (default) | 'operational'
%     'meanType'    'trailing' (default) | 'centered'   (81-day F10.7 mean)
%     'kpForecastDir'  folder with archived NOAA SWPC 3-day Kp forecast
%                   text files named yyyymmdd_3-day-forecast.txt (the
%                   "NOAA Kp index breakdown" table). Files not present
%                   are downloaded from the NOAA archive URL template in
%                   fetchKpForecast when possible. Epochs without an
%                   archived forecast are filled by an empirical error
%                   model (observed Kp + N(0, 0.9), clipped), flagged in
%                   drivers.KpFcSource (1 = archived, 0 = error model).
%
%   Sources: GFZ Kp_ap_Ap_SN_F107_since_1932.txt; WDC Kyoto monthly Dst
%   (final -> provisional -> realtime); NOAA SWPC 3-day forecast archive.

arguments
    epochs (:,1) datetime
    workDir (1,:) char = 'index_cache'
    opts.mode (1,:) char {mustBeMember(opts.mode, {'hindcast','operational'})} = 'hindcast'
    opts.meanType (1,:) char {mustBeMember(opts.meanType, {'trailing','centered'})} = 'trailing'
    opts.kpForecastDir (1,:) char = ''
    opts.kpErrSigma (1,1) double = 0.9
end

if ~isfolder(workDir), mkdir(workDir); end
epochs.TimeZone = 'UTC';
K = numel(epochs);

% =============== 1) GFZ: daily F10.7 + 3-hourly Kp =======================
gfzFile = fullfile(workDir, 'Kp_ap_Ap_SN_F107_since_1932.txt');
if ~isfile(gfzFile)
    websave(gfzFile, 'https://kp.gfz-potsdam.de/app/files/Kp_ap_Ap_SN_F107_since_1932.txt', weboptions('Timeout', 120));
end
[dayDates, F107daily, Kp3h, kpTimes] = parseGFZ(gfzFile);

switch opts.meanType
    case 'trailing'
        F107bar_daily = movmean(F107daily, [80 0], 'omitnan');   % past 81 days only
    case 'centered'
        F107bar_daily = movmean(F107daily, 81, 'omitnan');
end

dayNoon = dayDates + hours(12);
F107    = interp1(dayNoon, F107daily,    epochs, 'linear', 'extrap');
F107bar = interp1(dayNoon, F107bar_daily, epochs, 'linear', 'extrap');
Kp = interp1(kpTimes, Kp3h, epochs, 'linear', 'extrap');
Kp = min(max(Kp, 0), 9);

% =============== 2) Kyoto: hourly Dst ====================================
months = unique(dateshift(epochs, 'start', 'month'));
dstT = datetime.empty(0,1); dstT.TimeZone = 'UTC'; dstV = [];
for m = months'
    [tm, vm] = fetchDstMonth(year(m), month(m), workDir);
    dstT = [dstT; tm]; dstV = [dstV; vm]; %#ok<AGROW>
end
[dstT, iu] = unique(dstT); dstV = dstV(iu);
[tf, loc] = ismember(dateshift(epochs, 'start', 'hour'), dstT);
Dst = nan(K, 1); Dst(tf) = dstV(loc(tf));
if any(~tf)
    Dst = fillmissing(Dst, 'linear', 'EndValues', 'nearest');
    warning('%d Dst epochs missing at source; filled by interpolation.', nnz(~tf));
end

drivers = table(F107, F107bar, Kp, Dst);

% =============== 3) Operational: NOAA day-1 Kp forecast ==================
if strcmp(opts.mode, 'operational')
    [KpFc, src] = fetchKpForecast(epochs, Kp, opts.kpForecastDir, workDir, opts.kpErrSigma);
    drivers.KpFc = KpFc;
    drivers.KpFcErr = KpFc - Kp;
    drivers.KpFcSource = src;
    nA = nnz(src == 1);
    fprintf('Operational Kp forecast: %d/%d epochs archived (mean |err| %.2f, >2 units in %.1f %% of 3-h intervals)\n', ...
        nA, K, mean(abs(drivers.KpFcErr(src == 1))), 100 * mean(abs(drivers.KpFcErr(src == 1)) > 2));
end

fprintf(['Driver table (%s, %s F10.7 mean): %d hourly rows | F10.7 %.0f-%.0f sfu | ', ...
    'Kp max %.1f | Dst min %d nT\n'], opts.mode, opts.meanType, K, ...
    min(F107), max(F107), max(Kp), round(min(Dst)));
end

% =========================================================================
function [dayDates, F107, Kp3h, kpTimes] = parseGFZ(fileName)
raw = readmatrix(fileName, 'FileType', 'text', 'CommentStyle', '#');
dayDates = datetime(raw(:,1), raw(:,2), raw(:,3), 'TimeZone', 'UTC');
Kp8   = raw(:, 8:15);
F107  = raw(:, 26);                % F10.7 observed [sfu]
F107(F107 < 0) = NaN;              % -1 = missing
% burst-contaminated daily values (e.g. 939 sfu on 7 Mar 2011): outside
% [60, 350] sfu or > 30 % away from the 5-day running median -> NaN
med5 = movmedian(F107, 5, 'omitnan');
burst = F107 < 60 | F107 > 350 | abs(F107 - med5) > 0.30 * med5;
nB = nnz(burst & ~isnan(F107));
F107(burst) = NaN;
F107  = fillmissing(F107, 'linear', 'EndValues', 'nearest');
fprintf('GFZ F10.7: %d burst-contaminated daily values replaced by interpolation (%.2f %% of days)\n', nB, 100 * nB / numel(F107));
Kp8(Kp8 < 0) = NaN;
kpTimes = reshape((dayDates + hours(1.5) + hours(3)*(0:7))', [], 1);
Kp3h = reshape(Kp8', [], 1);
good = ~isnan(Kp3h);
kpTimes = kpTimes(good); Kp3h = Kp3h(good);
fprintf('GFZ file parsed: %d days (%s .. %s)\n', numel(dayDates), ...
    datestr(dayDates(1)), datestr(dayDates(end)));
end

% =========================================================================
function [t, v] = fetchDstMonth(yr, mo, workDir)
% One month of hourly Dst. Source 1: WDC Kyoto (final -> provisional ->
% realtime, 3 attempts each with a pause). Source 2 (fallback): NASA
% OMNI2 hourly yearly file (spdf.gsfc.nasa.gov, no login), field 41 = Dst.
persistent kyotoFails
if isempty(kyotoFails), kyotoFails = 0; end
tag = sprintf('%04d%02d', yr, mo);
localFile = fullfile(workDir, ['dst_' tag '.txt']);
if ~isfile(localFile) && kyotoFails >= 2          % Kyoto unreachable this session: go straight to OMNI2
    [t, v] = fetchDstOMNI(yr, mo, workDir); return
end
if ~isfile(localFile)
    stems = {'dst_final', 'dst_provisional', 'dst_realtime'};
    ok = false;
    for s = 1:numel(stems)
        url = sprintf('https://wdc.kugi.kyoto-u.ac.jp/%s/%s/dst%02d%02d.for.request', stems{s}, tag, mod(yr, 100), mo);
        for attempt = 1:3
            try
                websave(localFile, url, weboptions('Timeout', 60));
                txt = fileread(localFile);
                if contains(txt, 'DST') && ~contains(lower(txt), '<html'), ok = true; break; end
                delete(localFile);
            catch
                if isfile(localFile), delete(localFile); end
                pause(2 * attempt);
            end
        end
        if ok, break; end
    end
    if ok
        kyotoFails = 0; [t, v] = parseWDCDst(localFile); return
    end
    kyotoFails = kyotoFails + 1;
    warning('buildDriverTable:dst', 'Kyoto Dst not retrievable for %s; using NASA OMNI2%s.', tag, ...
        tern(kyotoFails >= 2, ' (Kyoto skipped for the rest of this session)', ''));
end
if isfile(localFile)
    [t, v] = parseWDCDst(localFile);
else
    [t, v] = fetchDstOMNI(yr, mo, workDir);
end
end

function [t, v] = fetchDstOMNI(yr, mo, workDir)
% NASA OMNI2 hourly yearly file: year, day-of-year, hour, ... , Dst (field 41)
f = fullfile(workDir, sprintf('omni2_%04d.dat', yr));
if ~isfile(f)
    websave(f, sprintf('https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2_%04d.dat', yr), weboptions('Timeout', 120));
end
raw = readmatrix(f, 'FileType', 'text');
t = datetime(raw(:, 1), 1, 1, 'TimeZone', 'UTC') + days(raw(:, 2) - 1) + hours(raw(:, 3));
v = raw(:, 41); v(v >= 9999) = NaN;
keep = month(t) == mo; t = t(keep); v = v(keep);
v = fillmissing(v, 'linear', 'EndValues', 'nearest');
end

% =========================================================================
function [t, v] = parseWDCDst(fileName)
lines = readlines(fileName);
t = datetime.empty(0,1); t.TimeZone = 'UTC'; v = [];
for L = lines'
    s = char(L);
    if numel(s) < 24 || ~strcmp(s(1:3), 'DST'), continue; end
    s = [s, repmat(' ', 1, max(0, 120 - numel(s)))];
    yy = str2double(s(4:5)); mo = str2double(s(6:7)); dd = str2double(s(9:10));
    cc = str2double(s(13:14)); if isnan(cc) || cc < 19, cc = 19 + (yy < 50); end
    best = []; bestErr = inf;
    for startCol = [21, 17]
        seg = s(startCol : startCol + 95);
        vals = sscanf(seg, '%4d');
        dm   = str2double(s(startCol + 96 : startCol + 99));
        if numel(vals) ~= 24, continue; end
        vals = double(vals); vals(vals == 9999) = NaN;
        if isnan(dm) || dm == 9999
            err = 0; off = 0;
        else
            off = round((dm - mean(vals, 'omitnan')) / 100) * 100;
            err = abs(dm - (mean(vals, 'omitnan') + off));
        end
        if err < bestErr, bestErr = err; best = vals + off; end
    end
    if isempty(best), continue; end
    d0 = datetime(100*cc + yy, mo, dd, 'TimeZone', 'UTC');
    t = [t; d0 + hours(0:23)']; v = [v; best]; %#ok<AGROW>
end
v = fillmissing(v, 'linear', 'EndValues', 'nearest');
end

% =========================================================================
function [KpFc, src] = fetchKpForecast(epochs, KpObs, fcDir, workDir, sigmaErr)
% Day-1 Kp forecast valid at each epoch, from the NOAA SWPC 3-day forecast
% issued on the previous day. Archived files are expected as
%   <fcDir>/yyyymmdd_3-day-forecast.txt      (yyyymmdd = issue date)
% containing the "NOAA Kp index breakdown" table:
%      NOAA Kp index breakdown May 09-May 11 2024
%                   May 09       May 10       May 11
%      00-03UT       2.67         3.00         4.33
%      ...
% If a file is missing it is requested from the NOAA archive URL template
% below; if that fails, the empirical error model is used (src = 0).
if isempty(fcDir), fcDir = fullfile(workDir, 'kp_forecast'); end
if ~isfolder(fcDir), mkdir(fcDir); end
urlTemplate = 'https://services.swpc.noaa.gov/text/archive/3-day-forecast/%s_3-day-forecast.txt';

K = numel(epochs);
KpFc = nan(K, 1); src = zeros(K, 1);
dayList = unique(dateshift(epochs, 'start', 'day'));
for d = dayList'
    issue = d - days(1);                              % forecast issued the day before
    fName = fullfile(fcDir, [datestr(issue, 'yyyymmdd') '_3-day-forecast.txt']);
    if ~isfile(fName)
        try, websave(fName, sprintf(urlTemplate, datestr(issue, 'yyyymmdd'))); catch, end
    end
    kp8 = [];
    if isfile(fName), kp8 = parseKpBreakdown(fName, d); end
    idx = find(dateshift(epochs, 'start', 'day') == d);
    if numel(kp8) == 8
        slot = floor(hour(epochs(idx)) / 3) + 1;
        KpFc(idx) = kp8(slot); src(idx) = 1;
    end
end
miss = isnan(KpFc);
if any(miss)
    rng(20240510);                                    % reproducible error model
    KpFc(miss) = min(max(KpObs(miss) + sigmaErr * randn(nnz(miss), 1), 0), 9);
    warning('buildDriverTable:kpfc', ...
        '%d epochs without archived Kp forecast: empirical error model used (sigma = %.1f).', nnz(miss), sigmaErr);
end
end

function s = tern(c, a, b), if c, s = a; else, s = b; end, end
