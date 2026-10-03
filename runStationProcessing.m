%RUNSTATIONPROCESSING Receiver-level VTEC of all stations (Sect. 2.3), standalone and resumable.
%   Runs in its own (light) MATLAB session, independently of runExperiments:
%   it needs only the epoch vector and the Dst series to select the days
%   (every cfg.stationDayStride-th day of 2023-2024 plus all days within
%   +/-48 h of a storm), then downloads and processes the RINEX files in
%   parallel (one station-day per worker; downloads are I/O-bound, so more
%   workers than cores help) and writes
%       stations/<code>_<yyyy>_<mm>.mat   one file per station-month (resume)
%       station_<code>.mat                all observations of a station
%       stations/positions.csv            receiver positions from the headers
%   runExperiments (do.stations / do.anatolia = true) then finds everything
%   in place and only evaluates. Stop with Ctrl-C at any time; re-running
%   continues with the unfinished months.
if ~usejava('desktop'), diary('log_stations.txt'); end       % unattended (matlab -batch): keep a log
cfg.stationYears = 2023:2024;
cfg.stationDayStride = 3;
cfg.rinexReader = 'rinexread';
cfg.workers = 16;                                     % pool size
cfg.loopWorkers = 16;                                 % concurrent station-days (benchmarkDownloads decides)
cfg.alternateSources = true;                          % alternate CDDIS / BKG to spread the load
cfg.netrc = fullfile(pwd, 'earthdata_netrc');
if ~isfile(cfg.netrc), error('runStationProcessing:netrc', 'Missing %s', cfg.netrc); end
% ---- days: systematic subsample + all storm days ----
load('TECdatabase.mat', 'epochs'); load('drivers_hindcast.mat', 'drivers');
sel = ismember(year(epochs), cfg.stationYears);
[~, stormMask] = stormAnalysis({}, {}, epochs(sel), drivers.Dst(sel), -100, 48);
allDays = (datetime(cfg.stationYears(1),1,1,'TimeZone','UTC') : datetime(cfg.stationYears(end),12,31,'TimeZone','UTC'))';
ep = epochs(sel); stormDays = unique(dateshift(ep(stormMask), 'start', 'day'));
stationDays = union(allDays(1:cfg.stationDayStride:end), stormDays);
clear epochs drivers ep
fprintf('Station validation: %d days per station (%d storm-window days)\n', numel(stationDays), nnz(ismember(stationDays, stormDays)));
save('stationDays.mat', 'stationDays');
% ---- worker pool ----
if isempty(gcp('nocreate'))
    try, parpool('local', cfg.workers); catch, parpool('local'); end
end
% ---- stations ----
Sall = stationList('all'); dcbCache = containers.Map(); t0 = tic;
months = unique(dateshift(stationDays, 'start', 'month')); missing = cell(0, 2);
for s = 1:height(Sall)
    code = char(Sall.code(s));
    fprintf('=== %s (%d/%d), elapsed %.1f h ===\n', code, s, height(Sall), toc(t0) / 3600);
    try
        obsAll = processStationYears(code, stationDays, dcbCache, cfg.rinexReader, cfg.loopWorkers, cfg.alternateSources); %#ok<NASGU>
    catch ME
        fprintf(2, '!!! %s failed: %s\n', code, ME.message);         % continue with the next station
    end
    for mo = months'
        if ~isfile(fullfile('stations', sprintf('%s_%04d_%02d.mat', code, year(mo), month(mo))))
            missing(end+1, :) = {code, datestr(mo, 'yyyy-mm')}; %#ok<AGROW>
        end
    end
end
T = cell2table(missing, 'VariableNames', {'station', 'month'}); writetable(T, 'stations_missing_months.csv');
fprintf('Station processing pass finished in %.1f h; %d station-months missing.\n', toc(t0) / 3600, size(missing, 1));
if isempty(missing)
    fid = fopen('STATIONS_DONE.txt', 'w'); fprintf(fid, 'complete %s\n', datestr(now)); fclose(fid);
end
