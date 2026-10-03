function stats = evaluateTEC(net, ds, lat, opts)
%EVALUATETEC Latitude-weighted RMSE / MAE / bias / correlation (Sect. 3.7).
%   stats = EVALUATETEC(net, ds, lat, ...)
%     net : trained dlnetwork; or a function handle f(X) returning
%           H x W x 1(2) x N numeric (scaled) predictions for the
%           model-function baselines (ConvLSTM, Transformer); or [] to
%           evaluate a precomputed ds.Yhat (H x W x 1 x N, scaled like ds.Y)
%     ds  : dataset struct from buildTECDataset
%     lat : H x 1 latitudes [deg]
%   Options
%     'RefY'   alternative reference H x W x 1 x N in TECU (e.g. the UQRG
%              maps at ds.epochsY, Sect. 4.6); default: ds.Y * ds.Tmax
%     'batch'  prediction batch size (default 64)
%
%   stats: RMSE, MAE, Bias [TECU], Corr, perEpochRMSE (N x 1),
%          Yhat (H x W x 1 x N, TECU) and, for a 2-channel network,
%          Sigma (H x W x 1 x N, predicted std in TECU, Sect. 4.9).
%   All statistics are computed on the W unique meridians of the grid.

arguments
    net
    ds struct
    lat (:,1) double
    opts.RefY = []
    opts.batch (1,1) double = 64
end

wlat = makeLatWeights(lat);
Tmax = ds.Tmax;
N = size(ds.Y, 4);

if isempty(net)
    Yhat = ds.Yhat;
else
    Yhat = [];
    useGPU = ~isa(net, 'function_handle') && canUseGPU;        % feed the GPU explicitly: predict() on a plain array runs on the CPU
    t0 = tic; nextPrint = 0;
    for i = 1:opts.batch:N
        j = min(i + opts.batch - 1, N);
        X = tecGetBatch(ds, i:j);
        if isa(net, 'function_handle')
            P = net(X);
        else
            if useGPU
                try, P = predict(net, gpuArray(X)); catch, useGPU = false; P = predict(net, X); end
            else
                P = predict(net, X);
            end
            if isa(P, 'dlarray'), P = extractdata(P); end
            P = gather(P);
        end
        if isempty(Yhat), Yhat = zeros([size(P, 1:3), N], 'single'); end
        Yhat(:, :, :, i:j) = single(P);
        if i >= nextPrint && ~isa(net, 'function_handle')
            fprintf('    evaluateTEC: %d / %d (%s, %.0f s)\n', i, N, ternary(useGPU, 'GPU', 'CPU'), toc(t0)); nextPrint = i + 5000;
        end
    end
end

Sigma = [];
if size(Yhat, 3) == 2
    Sigma = exp(0.5 * Yhat(:, :, 2, :)) * single(Tmax);
    Yhat = Yhat(:, :, 1, :);
end
Yhat = Yhat * single(Tmax);                                   % TECU, single

% ---- chunked statistics (memory-lean: no full double copies) ----
H = size(Yhat, 1); W = size(Yhat, 2);
w = repmat(wlat(:), 1, W);
sumSq = zeros(H, W); sumAbs = zeros(H, W); sumE = zeros(H, W); cnt = zeros(H, W);
perEpochRMSE = zeros(N, 1); perEpochCorr = zeros(N, 1);
chunk = 512;
for i = 1:chunk:N
    j = min(i + chunk - 1, N);
    P = double(Yhat(:, :, 1, i:j));
    if isempty(opts.RefY), T = double(ds.Y(:, :, 1, i:j)) * Tmax; else, T = double(opts.RefY(:, :, 1, i:j)); end
    E = P - T; valid = ~isnan(E);
    E(~valid) = 0;
    sumSq = sumSq + sum(E.^2, 4); sumAbs = sumAbs + sum(abs(E), 4); sumE = sumE + sum(E, 4); cnt = cnt + sum(valid, 4);
    for n = 1:size(E, 4)
        e = E(:, :, 1, n); v = valid(:, :, 1, n);
        if ~any(v, 'all'), perEpochRMSE(i+n-1) = NaN; perEpochCorr(i+n-1) = NaN; continue; end
        perEpochRMSE(i+n-1) = sqrt(sum(w(v) .* e(v).^2) / sum(w(v)));
        p = P(:, :, 1, n); t = T(:, :, 1, n);
        c = corrcoef(p(v), t(v)); perEpochCorr(i+n-1) = c(1, 2);
    end
end
cnt(cnt == 0) = NaN;
stats.RMSE = sqrt(mean(w .* (sumSq ./ cnt), 'all', 'omitnan'));
stats.MAE  = mean(w .* (sumAbs ./ cnt), 'all', 'omitnan');
stats.Bias = mean(w .* (sumE ./ cnt), 'all', 'omitnan');
stats.Corr = mean(perEpochCorr, 'omitnan');
stats.perEpochRMSE = perEpochRMSE;
stats.Yhat = Yhat;
if ~isempty(Sigma), stats.Sigma = Sigma; end
end

function r = ternary(c, a, b)
if c, r = a; else, r = b; end
end
