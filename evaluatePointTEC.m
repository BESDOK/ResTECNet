function T = evaluatePointTEC(statsList, names, ds, lat, lon, siteLat, siteLon, Dst, RefY)
%EVALUATEPOINTTEC Point-forecast statistics over a site (Kayseri case study).
%   T = EVALUATEPOINTTEC(statsList, names, ds, lat, lon, siteLat, siteLon, Dst)
%     statsList : cell array of stats structs from evaluateTEC
%                 (each must contain .Yhat, H x W x 1 x N, in TECU)
%     names     : cell array of model names
%     ds        : dataset struct from buildTECDataset (uses .Y, .Tmax,
%                 .epochsY)
%     lat, lon  : grid vectors [deg]
%     siteLat, siteLon : site coordinates, e.g. Kayseri 38.72, 35.49
%     Dst       : N x 1 hourly Dst [nT] aligned with ds.epochsY
%
%   Returns a table with RMSE (TECU) per model, stratified as in the
%   manuscript (Table: Kayseri): All, Winter (DJF), Equinoxes (MAM+SON),
%   Summer (JJA), Quiet (Dst > -30 nT), Storm windows (Dst <= -100 nT,
%   +/- 48 h). Bilinear interpolation of the TEC maps to the site
%   coordinates (interpolation of TEC at the nodes). Optional RefY
%   (H x W x 1 x N, TECU) replaces the CODE reference, e.g. the UQRG maps
%   interpolated to Kayseri (Sect. 4.7: 1.83 vs CODE, 1.97 vs UQRG).

arguments
    statsList cell
    names cell
    ds struct
    lat (:,1) double
    lon (:,1) double
    siteLat (1,1) double = 38.72     % Kayseri
    siteLon (1,1) double = 35.49
    Dst (:,1) double = []
    RefY = []
end

N = size(ds.Y, 4);
epochsY = ds.epochsY;

% ---- bilinear interpolation weights for the site ----
[wIdx, wVal] = bilinearWeights(lat, lon, siteLat, siteLon);

if isempty(RefY), RefY = double(ds.Y) * ds.Tmax; end
Tref = interpSite(double(RefY), wIdx, wVal);              % N x 1, TECU

% ---- strata masks ----
mo = month(epochsY);
mask.All      = true(N, 1);
mask.Winter   = ismember(mo, [12 1 2]);
mask.Equinox  = ismember(mo, [3 4 5 9 10 11]);
mask.Summer   = ismember(mo, [6 7 8]);
if ~isempty(Dst)
    stormMask = false(N, 1);
    below = find(Dst <= -100);
    while ~isempty(below)
        [~, k] = min(Dst(below));
        c = below(k);
        stormMask = stormMask | abs(hours(epochsY - epochsY(c))) <= 48;
        below(abs(hours(epochsY(below) - epochsY(c))) < 72) = [];
    end
    mask.Quiet = (Dst > -30) & ~stormMask;
    mask.Storm = stormMask;
end

strata = fieldnames(mask);
M = numel(statsList);
R = zeros(numel(strata), M);

for m = 1:M
    Tp = interpSite(double(statsList{m}.Yhat), wIdx, wVal);  % N x 1
    e2 = (Tp - Tref).^2;
    for s = 1:numel(strata)
        R(s, m) = sqrt(mean(e2(mask.(strata{s})), 'omitnan'));        % epochs without a forecast (e.g. CODE C1PG gaps) do not poison the stratum
    end
end

T = array2table(R, 'VariableNames', matlab.lang.makeValidName(names), ...
    'RowNames', strata);
fprintf('Point-forecast RMSE (TECU) at (%.2f N, %.2f E):\n', siteLat, siteLon);
disp(T);
end

% -------------------------------------------------------------------------
function [wIdx, wVal] = bilinearWeights(lat, lon, sLat, sLon)
% Indices and weights of the four surrounding nodes.
[~, iLat] = sort(abs(lat - sLat)); i1 = min(iLat(1:2)); i2 = max(iLat(1:2));
[~, iLon] = sort(abs(lon - sLon)); j1 = min(iLon(1:2)); j2 = max(iLon(1:2));
a = (sLat - lat(i1)) / (lat(i2) - lat(i1));
b = (sLon - lon(j1)) / (lon(j2) - lon(j1));
wIdx = [i1 j1; i1 j2; i2 j1; i2 j2];
wVal = [(1-a)*(1-b); (1-a)*b; a*(1-b); a*b];
end

function v = interpSite(A, wIdx, wVal)
% A: H x W x 1 x N -> v: N x 1 bilinear value at the site.
N = size(A, 4);
v = zeros(N, 1);
for k = 1:4
    v = v + wVal(k) * squeeze(A(wIdx(k,1), wIdx(k,2), 1, :));
end
end
