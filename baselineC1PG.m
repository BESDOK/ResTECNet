function stats = baselineC1PG(c1pgFile, ds, lat)
%BASELINEC1PG CODE one-day predicted GIM as an operational baseline (Sect. 3.6, ii).
%   stats = BASELINEC1PG(c1pgFile, ds, lat)
%     c1pgFile : MAT-file produced by
%                buildTECDatabase(..., 'product', 'C1PG', 'fillGaps', false)
%                for the test years (TEC, epochs, lat, lon)
%     ds       : test dataset from buildTECDataset (uses .epochsY, .Tmax)
%     lat      : H x 1 latitudes
%   The predicted map valid at each target epoch ds.epochsY is compared
%   with the CODE final map (ds.Y). Epochs without a C1PG map are NaN and
%   are ignored by the statistics (their count is reported).

S = load(c1pgFile, 'TEC', 'epochs');
[tf, loc] = ismember(dateshift(ds.epochsY, 'start', 'hour'), S.epochs);
N = numel(ds.epochsY);
Yhat = nan([size(ds.Y, 1:2), 1, N], 'single');
Yhat(:, :, 1, tf) = single(S.TEC(:, :, loc(tf))) / ds.Tmax;
if any(~tf)
    warning('baselineC1PG:missing', '%d of %d target epochs have no C1PG map.', nnz(~tf), N);
end
dsP = ds; dsP.Yhat = Yhat;
stats = evaluateTEC([], dsP, lat);
stats.available = tf;
end
