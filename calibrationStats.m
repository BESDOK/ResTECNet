function R = calibrationStats(statsU, ds, lat, stormMask, opts)
%CALIBRATIONSTATS Calibration of the Res-TECNet-U uncertainty map (Sect. 4.9).
%   R = CALIBRATIONSTATS(statsU, ds, lat, stormMask)
%     statsU    : from evaluateTEC on the 'uncertainty' network (.Yhat and
%                 .Sigma in TECU)
%     ds        : test dataset (reference .Y)
%     stormMask : N x 1 logical of storm-window epochs (from stormAnalysis
%                 logic; [] to skip)
%   Options: 'nBins' (default 10) equally populated sigma classes.
%   R fields: cov1 (fraction of residuals within +/-1 sigma), cov2,
%   cov1_storm, binSigma / binRMSE (nBins x 1), ratioRange [min max] of
%   RMSE/sigma across classes, corrBins, sigmaInflationStorm (mean sigma in
%   storm windows / mean sigma elsewhere). Latitude weights are applied
%   through weighted sampling of the coverage statistics.

arguments
    statsU struct
    ds struct
    lat (:,1) double
    stormMask logical = false(0,1)
    opts.nBins (1,1) double = 10
end
E = double(statsU.Yhat) - double(ds.Y) * ds.Tmax;
S = double(statsU.Sigma);
w = repmat(makeLatWeights(lat), [1, size(E, 2), 1, size(E, 4)]);
z = abs(E) ./ S;
R.cov1 = sum(w(:) .* (z(:) <= 1)) / sum(w(:));
R.cov2 = sum(w(:) .* (z(:) <= 2)) / sum(w(:));
if ~isempty(stormMask)
    zs = z(:, :, :, stormMask); ws = w(:, :, :, stormMask);
    R.cov1_storm = sum(ws(:) .* (zs(:) <= 1)) / sum(ws(:));
    R.sigmaInflationStorm = mean(S(:, :, :, stormMask), 'all') / mean(S(:, :, :, ~stormMask), 'all');
end
% equally populated sigma classes
[~, order] = sort(S(:));
nb = opts.nBins; edges = round(linspace(0, numel(order), nb + 1));
R.binSigma = zeros(nb, 1); R.binRMSE = zeros(nb, 1);
Ev = E(:); Sv = S(:);
for b = 1:nb
    i = order(edges(b) + 1 : edges(b + 1));
    R.binSigma(b) = mean(Sv(i)); R.binRMSE(b) = sqrt(mean(Ev(i).^2));
end
cc = corrcoef(R.binSigma, R.binRMSE); R.corrBins = cc(1, 2);
ratio = R.binRMSE ./ R.binSigma; R.ratioRange = [min(ratio), max(ratio)];
fprintf('Calibration: 68%% -> %.1f %%, 95%% -> %.1f %%, corr(sigma,RMSE) over %d classes = %.2f, ratio %.2f-%.2f\n', ...
    100 * R.cov1, 100 * R.cov2, nb, R.corrBins, R.ratioRange(1), R.ratioRange(2));
if isfield(R, 'sigmaInflationStorm')
    fprintf('Calibration in storm windows: sigma inflation x%.2f, 68%% -> %.1f %%\n', R.sigmaInflationStorm, 100 * R.cov1_storm);
end
end
