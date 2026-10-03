function T = evaluateSlantDelay(statsList, names, ds, lat, lon, obs, opts)
%EVALUATESLANTDELAY Slant-delay error of map-based corrections over real geometry.
%   T = EVALUATESLANTDELAY(statsList, names, ds, lat, lon, obs, ...)
%     statsList/names : as in evaluateStationVTEC (CODE final map added
%                       automatically); a Klobuchar entry is added when
%                       'klobuchar' (struct with alpha, beta from the nav
%                       header, and rxLat, rxLon) is given
%     obs             : timetable from processStationVTEC (all days), one
%                       station (e.g. MERS, year 2024)
%   For every observation above 'elevMin' (default 10 deg) the map VTEC is
%   interpolated to the IPP, mapped to slant TEC with the MSLM factor
%   stored in obs.mf and compared with the levelled dual-frequency STEC.
%   Returns RMS slant-delay error at L1 [m] (1 TECU = 0.162 m) for all
%   elevations and for el >= 30 deg (Sect. 4.7.1), with the station's mean
%   offset to the CODE GIM (receiver-DCB level) removed, plus the raw values.

arguments
    statsList cell
    names cell
    ds struct
    lat (:,1) double
    lon (:,1) double
    obs timetable
    opts.elevMin (1,1) double = 10
    opts.klobuchar = []
end

mPerTECU = 0.162;
allNames = [names(:)', {'CODE GIM'}];
maps = [statsList, {struct('Yhat', single(double(ds.Y) * ds.Tmax))}];
O = obs(obs.ok & obs.el >= opts.elevMin, :);
tt = O.Properties.RowTimes; tt.TimeZone = 'UTC';
% map valid at each observation: linear interpolation in time between the
% enclosing hourly target epochs (nearest epoch within 30 min if the
% record has a gap); observations farther than 30 min from any epoch are dropped
eY = ds.epochsY; eY.TimeZone = 'UTC';
kA = interp1(datenum(eY), 1:numel(eY), datenum(tt), 'previous');      % index of epoch <= tt
kB = interp1(datenum(eY), 1:numel(eY), datenum(tt), 'next');          % index of epoch >= tt
kN = interp1(datenum(eY), 1:numel(eY), datenum(tt), 'nearest', 'extrap');
kN = min(max(round(kN), 1), numel(eY));
okT = abs(minutes(tt - eY(kN))) <= 30;
O = O(okT, :); tt = tt(okT); kA = kA(okT); kB = kB(okT); kN = kN(okT);
encl = ~isnan(kA) & ~isnan(kB) & (kB - kA <= 1);                     % enclosed by consecutive epochs
kA(~encl) = kN(~encl); kB(~encl) = kN(~encl);
wB = zeros(size(kA));
i2 = kB > kA;
wB(i2) = seconds(tt(i2) - eY(kA(i2))) ./ seconds(eY(kB(i2)) - eY(kA(i2)));
nO = height(O); M = numel(maps);
errT = nan(nO, M);                                          % model STEC minus station STEC [TECU]
for m = 1:M
    Yall = double(maps{m}.Yhat);
    vA = zeros(nO, 1); vB = zeros(nO, 1);
    for n = unique([kA; kB])'
        iA = kA == n; iB = kB == n;
        if any(iA), vA(iA) = interpMapAtPoints(Yall(:, :, 1, n), lat, lon, O.ippLat(iA), O.ippLon(iA)); end
        if any(iB), vB(iB) = interpMapAtPoints(Yall(:, :, 1, n), lat, lon, O.ippLat(iB), O.ippLon(iB)); end
    end
    v = (1 - wB) .* vA + wB .* vB;
    errT(:, m) = v .* O.mf - O.STEC;
end
off = mean(errT(:, M), 'omitnan');                          % CODE GIM minus station: receiver-DCB level [TECU]
errRaw = errT * mPerTECU; errDeb = (errT - off) * mPerTECU; % offset removed (station mean vs CODE GIM)
if ~isempty(opts.klobuchar)
    k = opts.klobuchar;
    dI = klobucharModel(k.alpha, k.beta, k.rxLat, k.rxLon, O.el, O.az, tt);            % m at L1
    errRaw(:, end+1) = dI - O.STEC * mPerTECU;
    errDeb(:, end+1) = dI - (O.STEC + off) * mPerTECU;       % station STEC shifted by the offset: obs' = obs + off
    allNames{end+1} = 'Klobuchar';
end
hi = O.el >= 30;
rms = @(E, m) sqrt(mean(E(m, :).^2, 1, 'omitnan'))';
T = table(allNames(:), rms(errDeb, true(nO, 1)), rms(errDeb, hi), rms(errRaw, true(nO, 1)), rms(errRaw, hi), ...
    'VariableNames', {'Correction', 'RMS_slant_m_all', 'RMS_slant_m_el30', 'RMS_raw_m_all', 'RMS_raw_m_el30'});
fprintf('Slant-delay RMS error at L1 (m), %d observations, station offset %+.2f TECU removed (raw columns keep it):\n', nO, off); disp(T);
end

function v = interpMapAtPoints(map, lat, lon, pLat, pLon)
[H, W] = size(map);
dlon = lon(2) - lon(1);
fj = (pLon - lon(1)) / dlon; j1 = floor(fj); b = fj - j1; j1 = mod(j1, W) + 1; j2 = mod(j1, W) + 1;
fi = (pLat - lat(1)) / (lat(2) - lat(1)); i1 = floor(fi); a = min(max(fi - i1, 0), 1);
i1 = min(max(i1 + 1, 1), H - 1); i2 = i1 + 1;
v = (1-a).*(1-b).*map(sub2ind([H W], i1, j1)) + (1-a).*b.*map(sub2ind([H W], i1, j2)) ...
  + a.*(1-b).*map(sub2ind([H W], i2, j1)) + a.*b.*map(sub2ind([H W], i2, j2));
end
