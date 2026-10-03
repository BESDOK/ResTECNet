function [Tst, Treg, hourly, raw] = evaluateStationVTEC(statsList, names, ds, lat, lon, stations, opts)
%EVALUATESTATIONVTEC Forecast vs dual-frequency station VTEC (Sect. 4.6/4.7).
%   [Tst, Treg, hourly] = EVALUATESTATIONVTEC(statsList, names, ds, lat, lon, stations, ...)
%     statsList : cell of stats structs (each with .Yhat H x W x 1 x N, TECU,
%                 aligned to ds.epochsY). The CODE final map (ds.Y) is
%                 evaluated automatically as the accuracy floor "CODE GIM".
%     names     : cell of model names
%     stations  : struct array with fields code, regime and obs (timetable
%                 from processStationVTEC, all days concatenated)
%   Options
%     'elevMin'   elevation cut-off for the VTEC comparison (default 30)
%     'window'    half-window around each hourly epoch [min] (default 30)
%     'stormMask' N x 1 logical (storm windows) -> extra rows "…(storm)"
%   Method (Sect. 2.3 step 5): every map is bilinearly interpolated to the
%   IPP coordinates of each accepted observation (longitude wrapped
%   circularly on the 72-column grid); the hourly station value is the mean
%   over the accepted IPPs within +/- window; RMSE is computed over the hours
%   with observations.
%   Tst  : table, one row per station, one column per model (+ CODE GIM), with the
%          station's mean offset to the CODE GIM removed (absorbs the receiver-DCB
%          level, uncertain by a few TECU with thin-shell processing)
%   raw  : struct with the same tables WITHOUT the offset removal (raw.Tst, raw.Treg)
%   Treg : same aggregated by regime (pooled squared errors) and over all
%   hourly : struct array with the hourly series per station (for
%            regionalGradients and figures)

arguments
    statsList cell
    names cell
    ds struct
    lat (:,1) double
    lon (:,1) double
    stations struct
    opts.elevMin (1,1) double = 30
    opts.window (1,1) double = 30
    opts.stormMask logical = false(0,1)
end

M = numel(statsList);
allNames = [names(:)', {'CODE GIM'}];
maps = [statsList, {struct('Yhat', single(double(ds.Y) * ds.Tmax))}];
N = numel(ds.epochsY);
nS = numel(stations);
hourly = struct('code', {}, 'regime', {}, 'time', {}, 'obs', {}, 'model', {}, 'names', {});

for s = 1:nS
    O = stations(s).obs;
    O = O(O.ok & O.el >= opts.elevMin, :);
    obsH = nan(N, 1); modH = nan(N, M + 1);
    if height(O) > 0
        tt = O.Properties.RowTimes; tt.TimeZone = 'UTC';
        h0 = dateshift(tt, 'start', 'hour');
        hr = h0 + hours(minutes(tt - h0) >= 30);                 % nearest hourly epoch (+/- 30 min)
        [tf, loc] = ismember(hr, ds.epochsY);
        idx = find(tf & abs(minutes(tt - hr)) <= opts.window);
        [ue, ~, gi] = unique(loc(idx));                           % epoch index of every group
        cntG = accumarray(gi, 1);
        [~, ord] = sort(gi); idxS = idx(ord); ends = cumsum(cntG); starts = [1; ends(1:end-1) + 1];
        for g = 1:numel(ue)
            if cntG(g) < 3, continue; end
            ii = idxS(starts(g):ends(g)); n = ue(g);
            obsH(n) = mean(O.VTEC(ii));
            for m = 1:M + 1
                v = interpMapAtPoints(double(maps{m}.Yhat(:, :, 1, n)), lat, lon, O.ippLat(ii), O.ippLon(ii));
                modH(n, m) = mean(v);
            end
        end
    end
    hourly(s).code = stations(s).code; hourly(s).regime = stations(s).regime;
    hourly(s).time = ds.epochsY; hourly(s).obs = obsH; hourly(s).model = modH; hourly(s).names = allNames;
    valid = ~isnan(obsH);
    if any(valid), hourly(s).biasVsCODE = mean(obsH(valid) - modH(valid, end)); else, hourly(s).biasVsCODE = 0; end
    hourly(s).corrVsCODE = NaN;
    if nnz(valid) > 10, c = corrcoef(obsH(valid), modH(valid, end)); hourly(s).corrVsCODE = c(1, 2); end
    if nnz(valid) > 10 && hourly(s).corrVsCODE < 0.9
        warning('evaluateStationVTEC:corr', 'Station %s: correlation with the CODE GIM is only %.2f -- check the station processing.', stations(s).code, hourly(s).corrVsCODE);
    end
end

[Tst, Treg] = summarize(hourly, allNames, opts.stormMask, true);
[Tr, Rr] = summarize(hourly, allNames, opts.stormMask, false);
raw.Tst = Tr; raw.Treg = Rr;
fprintf('RMSE (TECU) vs station VTEC by regime, station-mean offset removed:\n'); disp(Treg);
fprintf('RMSE (TECU) vs station VTEC by regime, raw (receiver-DCB level included):\n'); disp(Rr);
fprintf('Station VTEC - CODE GIM offsets (TECU): %s\n', strjoin(compose('%s %+.1f', string({hourly.code})', [hourly.biasVsCODE]'), ', '));
end

% -------------------------------------------------------------------------
function [Tst, Treg] = summarize(hourly, allNames, stormMask, debias)
% per-station and per-regime RMSE; debias = remove the station's mean
% (station VTEC - CODE GIM) offset, which absorbs the receiver-DCB level.
% Epochs for which a model has no forecast (NaN, e.g. missing C1PG maps)
% are excluded for that model only.
nS = numel(hourly); nM = numel(allNames);
sumSq = zeros(nS, nM); cnt = zeros(nS, nM); sumSqS = zeros(nS, nM); cntS = zeros(nS, nM);
useStorm = ~isempty(stormMask);
for s = 1:nS
    obsH = hourly(s).obs; modH = hourly(s).model;
    b = 0; if debias, b = hourly(s).biasVsCODE; end
    e = modH - obsH + b;                                   % model minus (offset-corrected) station
    ok = ~isnan(e); e(~ok) = 0;
    sumSq(s, :) = sum(e.^2, 1); cnt(s, :) = sum(ok, 1);
    if useStorm
        oks = ok & stormMask(:);
        sumSqS(s, :) = sum((e .* oks).^2, 1); cntS(s, :) = sum(oks, 1);
    end
end
vn = matlab.lang.makeValidName(allNames);
R = sqrt(sumSq ./ max(cnt, 1)); R(cnt == 0) = NaN;
Tst = array2table(R, 'VariableNames', vn, 'RowNames', cellstr(string({hourly.code})));
Tst.nHours = cnt(:, end);                                  % hours with station data (CODE GIM column)
if useStorm
    Rs = sqrt(sumSqS ./ max(cntS, 1)); Rs(cntS == 0) = NaN;
    Tst = [Tst; array2table([Rs, cntS(:, end)], 'VariableNames', [vn, {'nHours'}], 'RowNames', cellstr(string({hourly.code}) + " (storm)"))]; %#ok<AGROW>
end
regimes = unique(string({hourly.regime}), 'stable');
Rr = zeros(numel(regimes) + 1, nM);
for r = 1:numel(regimes)
    i = string({hourly.regime}) == regimes(r);
    Rr(r, :) = sqrt(sum(sumSq(i, :), 1) ./ max(sum(cnt(i, :), 1), 1));
end
Rr(end, :) = sqrt(sum(sumSq, 1) ./ max(sum(cnt, 1), 1));
rn = [cellstr(regimes(:)); {'All stations'}];
if useStorm
    Rr(end+1, :) = sqrt(sum(sumSqS, 1) ./ max(sum(cntS, 1), 1)); rn{end+1} = 'All stations (storm)';
end
Treg = array2table(Rr, 'VariableNames', vn, 'RowNames', rn);
end

% -------------------------------------------------------------------------
function v = interpMapAtPoints(map, lat, lon, pLat, pLon)
% bilinear interpolation with circular longitude on the W unique meridians
[H, W] = size(map);
dlon = lon(2) - lon(1);
fj = (pLon - lon(1)) / dlon;                         % fractional column (0-based)
j1 = floor(fj); b = fj - j1; j1 = mod(j1, W) + 1; j2 = mod(j1, W) + 1;
fi = (pLat - lat(1)) / (lat(2) - lat(1));
i1 = floor(fi); a = fi - i1; i1 = min(max(i1 + 1, 1), H - 1); i2 = i1 + 1;
a = min(max(a, 0), 1);
v = (1 - a) .* (1 - b) .* map(sub2ind([H W], i1, j1)) + (1 - a) .* b .* map(sub2ind([H W], i1, j2)) ...
  + a .* (1 - b) .* map(sub2ind([H W], i2, j1)) + a .* b .* map(sub2ind([H W], i2, j2));
end
