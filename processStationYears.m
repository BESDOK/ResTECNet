function obsAll = processStationYears(code, dayList, dcbCache, reader, workers, alternate)
%PROCESSSTATIONYEARS Receiver-level VTEC of one station for an explicit day list.
%   Resumable: one file per station-month (stations/<code>_yyyy_mm.mat) is
%   written after each year-block; finished months are skipped. All pending
%   days of a year are processed in ONE parfor (download + VTEC processing
%   per station-day), so that the workers stay busy instead of waiting for
%   the slowest day of every month. A month is saved only if at least half
%   of its days could be processed (otherwise it is retried on the next run).
%     workers   concurrent station-days (default 8)
%     alternate true: odd/even days use the BKG / CDDIS mirror first
if nargin < 5, workers = 8; end
if nargin < 6, alternate = false; end
St = stationList('all'); iso = char(St.iso(St.code == string(code)));
if ~isfolder('stations'), mkdir('stations'); end
netrc = fullfile(pwd, 'earthdata_netrc');
months = unique(dateshift(dayList, 'start', 'month'));
pieces = {};
for yr = unique(year(months))'
    mo = months(year(months) == yr);
    pend = false(numel(mo), 1);
    for i = 1:numel(mo)                                          % finished months
        f = monthFile(code, mo(i));
        if isfile(f)
            M = load(f);
            if height(M.obsM) > 0, pieces{end+1} = M.obsM; end %#ok<AGROW>
        else
            pend(i) = true;
        end
    end
    mo = mo(pend);
    if isempty(mo), continue; end
    % DCBs of the pending months (small, broadcast to the workers)
    dcbCell = cell(numel(mo), 1); dcbC1Cell = cell(numel(mo), 1); okM = true(numel(mo), 1);
    for i = 1:numel(mo)
        key = sprintf('%d-%02d', year(mo(i)), month(mo(i)));
        if ~isKey(dcbCache, key)
            try, dcbCache(key) = readCODEDCB(year(mo(i)), month(mo(i)), 'dcb_cache', 'P1P2');
            catch ME, warning('%s: %s', key, ME.message); end
        end
        if ~isKey(dcbCache, key), okM(i) = false; continue; end
        keyC = [key '-c1'];                                      % P1-C1 DCB (C/A-code receivers)
        if ~isKey(dcbCache, keyC)
            try, dcbCache(keyC) = readCODEDCB(year(mo(i)), month(mo(i)), 'dcb_cache', 'P1C1'); catch, dcbCache(keyC) = []; end
        end
        dcbCell{i} = dcbCache(key); dcbC1Cell{i} = dcbCache(keyC);
    end
    allD = dayList(ismember(dateshift(dayList, 'start', 'month'), mo));
    [~, mIdx] = ismember(dateshift(allD, 'start', 'month'), mo);
    keep = okM(mIdx); allD = allD(keep); mIdx = mIdx(keep);
    nD = numel(allD); if nD == 0, continue; end
    parts = cell(nD, 1); metas = cell(nD, 1); tdl = nan(nD, 1); tpr = nan(nD, 1);
    tYear = tic;
    parfor (k = 1:nD, workers)
        try
            t1 = tic;
            [oFk, nFk] = downloadRINEX(code, allD(k), 'rinex_cache', 'cddisNetrc', netrc, 'country', iso, 'preferBKG', alternate && mod(k, 2) == 0);
            tdl(k) = toc(t1);
            if isempty(oFk) || isempty(nFk), continue; end
            t2 = tic;
            [parts{k}, metas{k}] = processStationVTEC(oFk, nFk, dcbCell{mIdx(k)}, 'dcbP1C1', dcbC1Cell{mIdx(k)}, 'reader', reader);
            tpr(k) = toc(t2);
            cleanupObs(oFk);
        catch ME
            warning('%s %s: %s', code, datestr(allD(k)), ME.message);
        end
    end
    km = find(~cellfun(@isempty, metas), 1);
    if ~isempty(km), recordPosition(code, metas{km}); end        % exact coordinates for the tables
    for i = 1:numel(mo)
        sel = mIdx == i; if ~any(sel), continue; end
        okParts = parts(sel & ~cellfun(@isempty, parts));
        if numel(okParts) < max(1, ceil(0.5 * nnz(sel)))
            warning('%s %s: only %d of %d days processed -- month not saved, will be retried on the next run.', code, datestr(mo(i), 'yyyy-mm'), numel(okParts), nnz(sel));
            continue
        end
        obsM = vertcat(okParts{:});
        save(monthFile(code, mo(i)), 'obsM', '-v7.3');
        pieces{end+1} = obsM; %#ok<AGROW>
        fprintf('%s %s: %d observations from %d of %d days saved\n', code, datestr(mo(i), 'yyyy-mm'), height(obsM), numel(okParts), nnz(sel));
    end
    el = toc(tYear) / 60;
    fprintf('%s %d: %d station-days in %.1f min = %.1f days/min (per day: download %.0f s, processing %.0f s)\n', ...
        code, yr, nD, el, nD / el, mean(tdl, 'omitnan'), mean(tpr, 'omitnan'));
end
if isempty(pieces), obsAll = timetable(); else, obsAll = sortrows(vertcat(pieces{:})); end
save(sprintf('station_%s.mat', code), 'obsAll', '-v7.3');
end

function f = monthFile(code, d)
f = fullfile('stations', sprintf('%s_%04d_%02d.mat', code, year(d), month(d)));
end

function recordPosition(code, meta)
% append the receiver position read from the RINEX header (once per station)
f = fullfile('stations', 'positions.csv');
if isfile(f) && any(startsWith(readlines(f), [code ',']))
    return
end
fid = fopen(f, 'a'); fprintf(fid, '%s,%.4f,%.4f\n', code, meta.rxLat, meta.rxLon); fclose(fid);
end

function cleanupObs(oF)
% delete the observation file and its compressed / Hatanaka intermediates
[p, n] = fileparts(oF);
try, delete(fullfile(p, [n '.*'])); catch, end
end
