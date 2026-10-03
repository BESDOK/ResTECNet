function R = interpolationInfluence(TEC, epochs, lat, lon, drivers, splits, cfg, opts)
%INTERPOLATIONINFLUENCE Effect of the pre-2015 2 h -> 1 h resampling (Sect. 4.8.1).
%   R = INTERPOLATIONINFLUENCE(TEC, epochs, lat, lon, drivers, splits, cfg, ...)
%   Test 1 (decimation): the native 1-h maps of 'decimYears' (default
%     2015-2018) are decimated to 2 h and re-interpolated with the same
%     sun-fixed scheme; R.decimRMSE_global and R.decimRMSE_EIA give the RMS
%     difference between interpolated and original odd-hour maps [TECU]
%     (EIA = |dipole magnetic latitude| <= 15 deg).
%   Test 2 (native-only training, optional 'trainNative' = true): a model
%     trained only on native 1-h data from 19 Oct 2014 to the end of the
%     training interval is evaluated on the test set; R.nativeOnly holds
%     the test statistics (compare with the full-record model).

arguments
    TEC (:,:,:) {mustBeNumeric}
    epochs (:,1) datetime
    lat (:,1) double
    lon (:,1) double
    drivers table
    splits struct
    cfg struct
    opts.decimYears double = 2015:2018
    opts.trainNative (1,1) logical = false
    opts.trainOpts cell = {}
    opts.DstTest double = []
    opts.stride (1,1) double = 1
    opts.scaled (1,1) logical = true   % TEC already divided by Tmax (runExperiments)
end
epochs.TimeZone = 'UTC';
dlon = lon(2) - lon(1);
sc = 1; if opts.scaled, sc = cfg.Tmax; end          % test 1 reports TECU

% ---------- Test 1 ----------
sel = find(ismember(year(epochs), opts.decimYears));
even = sel(mod(hour(epochs(sel)), 2) == 0);
odd  = sel(mod(hour(epochs(sel)), 2) == 1);
[LON, LAT] = meshgrid(lon, lat);
magLat = asind(sind(LAT) .* sind(80.65) + cosd(LAT) .* cosd(80.65) .* cosd(LON + 72.68));
eia = abs(magLat) <= 15;
w = repmat(makeLatWeights(lat), 1, numel(lon));
sq = 0; sqE = 0; n = 0;
for k = odd'
    a = k - 1; b = k + 1;
    if ~ismember(a, even) || ~ismember(b, even), continue; end
    A = circshift(TEC(:, :, a), -round(15 / dlon), 2);
    B = circshift(TEC(:, :, b),  round(15 / dlon), 2);
    interp = 0.5 * A + 0.5 * B;
    e = double(interp - TEC(:, :, k)) * sc;
    sq = sq + mean(w .* e.^2, 'all'); sqE = sqE + mean(e(eia).^2); n = n + 1;
end
R.decimRMSE_global = sqrt(sq / n); R.decimRMSE_EIA = sqrt(sqE / n); R.nMaps = n;
fprintf('2h->1h re-interpolation error over %d maps: %.2f TECU global, %.2f TECU EIA\n', ...
    n, R.decimRMSE_global, R.decimRMSE_EIA);

% ---------- Test 2 ----------
if opts.trainNative
    nativeStart = datetime(2014, 10, 19, 'TimeZone', 'UTC');
    trN = splits.tr & (epochs >= nativeStart);
    mk = @(m, varargin) buildTECDataset(TEC, epochs, drivers, 'L', cfg.L, 'Delta', cfg.Delta, ...
        'Tmax', cfg.Tmax, 'mask', m, 'scaled', opts.scaled, varargin{:});
    dsTr = mk(trN, 'stride', opts.stride, 'targets', false); dsVa = mk(splits.va, 'zmu', dsTr.zmu, 'zsig', dsTr.zsig);
    dsTe = mk(splits.te, 'zmu', dsTr.zmu, 'zsig', dsTr.zsig);
    H = dsTr.H; W = dsTr.W; C = dsTr.C;
    net = createResTECNet(H, W, C, 'L', cfg.L);
    net = trainResTECNet(net, dsTr, dsVa, lat, 'lambda', cfg.lambda, opts.trainOpts{:});
    R.nativeOnly = evaluateTEC(net, dsTe, lat);
    if ~isempty(opts.DstTest)
        R.nativeOnlyStorm = stormAnalysis({R.nativeOnly}, {'nativeOnly'}, dsTe.epochsY, opts.DstTest, -100, 48);
    end
    fprintf('Native-only training: test RMSE %.3f TECU\n', R.nativeOnly.RMSE);
end
end
