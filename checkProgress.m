function checkProgress()
%CHECKPROGRESS Progress of the station processing (run from any MATLAB session in the code folder).
%   Counts the finished station-months (stations/<code>_<yyyy>_<mm>.mat), estimates the
%   rate from the file time stamps of the last hour and prints an ETA.
S = stationList('all'); nMonths = 24; total = height(S) * nMonths;
f = dir(fullfile('stations', '*_20??_??.mat'));
done = numel(f);
perStation = zeros(height(S), 1);
for i = 1:height(S), perStation(i) = nnz(startsWith({f.name}, [char(S.code(i)) '_'])); end
fprintf('Station-months finished: %d of %d (%.0f %%)\n', done, total, 100 * done / total);
for i = 1:height(S), fprintf('  %-5s %2d/%d\n', S.code(i), perStation(i), nMonths); end
if ~isempty(f)
    t = [f.datenum]; recent = t > now - 1/24;
    if nnz(recent) >= 3
        rate = nnz(recent);                                   % station-months per hour
        fprintf('Last hour: %d station-months -> remaining about %.1f h\n', rate, (total - done) / rate);
    else
        fprintf('Fewer than 3 station-months in the last hour (finished blocks are written in chunks).\n');
    end
    fprintf('Newest file: %s\n', datestr(max(t)));
end
for n = {'STATIONS_DONE.txt', 'run_unattended_status.txt'}
    if isfile(n{1}), fprintf('--- %s ---\n%s\n', n{1}, fileread(n{1})); end
end
end
