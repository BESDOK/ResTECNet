function T = benchmarkDownloads(levels)
%BENCHMARKDOWNLOADS Throughput of the RINEX download chain versus concurrency.
%   T = BENCHMARKDOWNLOADS([1 2 4 8]) downloads (and unpacks, without any
%   VTEC processing) daily files of the station ONSA for n distinct,
%   uncached days with n parallel workers, separately from CDDIS and from
%   BKG, and reports the wall time and the throughput in station-days per
%   minute. The best level (highest days/min) and the faster server are the
%   settings for runStationProcessing (cfg.loopWorkers / alternateSources).
%   Takes about 10 minutes; downloaded files are deleted again.
arguments
    levels (1,:) double = [1 2 4 8]
end
netrc = fullfile(pwd, 'earthdata_netrc');
if isempty(gcp('nocreate')), parpool('Processes', max(levels)); end
base = datetime(2022, 3, 1, 'TimeZone', 'UTC'); used = 0; rows = {};
for src = {'cddis', 'bkg'}
    for n = levels
        days = base + caldays(used + (1:n)'); used = used + n;
        t0 = tic;
        ok = false(n, 1);
        parfor (k = 1:n, n)
            o = downloadRINEX('ONSA', days(k), 'rinex_bench', 'cddisNetrc', netrc, 'country', 'SWE', 'sources', src{1});
            ok(k) = ~isempty(o);
        end
        el = toc(t0);
        rows(end+1, :) = {src{1}, n, sum(ok), el, sum(ok) / (el / 60)}; %#ok<AGROW>
        fprintf('%-6s %2d parallel: %d/%d days in %.0f s -> %.1f days/min\n', src{1}, n, sum(ok), n, el, sum(ok) / (el / 60));
        if isfolder('rinex_bench'), delete(fullfile('rinex_bench', '*ONSA*')); end
    end
end
T = cell2table(rows, 'VariableNames', {'server', 'parallel', 'daysOK', 'seconds', 'daysPerMin'});
disp(T);
end
