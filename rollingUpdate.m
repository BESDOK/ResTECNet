function [stats, log] = rollingUpdate(base, kind, mk, epochs, drivers, dsTe, lat, cfg, testYears, opts)
%ROLLINGUPDATE Annually updated evaluation (Sect. 3.3, "operational retraining").
%   [stats, log] = ROLLINGUPDATE(base, kind, mk, epochs, drivers, dsTe, lat, cfg, testYears, ...)
%   For each test year Y the model is fine-tuned, starting from the static
%   network trained on 1999-2018, on the data of the years max(2019, Y-2)..Y-1
%   (never on year Y or later), and evaluated on year Y. The predictions of
%   all years are assembled into one stats struct aligned with dsTe.epochsY,
%   so that the result is directly comparable with the static model. This
%   mirrors what an operational system does (periodic retraining) and what
%   the adaptive baselines (SH-AR, C1PG) do implicitly.
%     base : dlnetwork ('dlnet') or struct params ('fcn')
%     kind : 'dlnet' | 'fcn'
%     mk   : @(mask, drv, varargin) buildTECDataset(...) with 'scaled'/'badMask' set
%   Options: 'epochs' (4), 'lr' (3e-4), 'stride' (3), 'loss' ('mse'|'nll'),
%            'modelFcn' (for 'fcn'), 'drvTrain' (drivers table for training,
%            default = drivers), 'tag' (checkpoint tag), 'mode'.

arguments
    base
    kind (1,:) char
    mk function_handle
    epochs (:,1) datetime
    drivers table
    dsTe struct
    lat (:,1) double
    cfg struct
    testYears double
    opts.epochs (1,1) double = 4
    opts.lr (1,1) double = 3e-4
    opts.stride (1,1) double = 3
    opts.loss (1,:) char = 'mse'
    opts.modelFcn = []
    opts.tag (1,:) char = 'roll'
    opts.mode (1,:) char = 'hindcast'
end
N = dsTe.N;
Yhat = nan([dsTe.H, dsTe.W, 1 + strcmp(opts.loss, 'nll'), N], 'single');
yTe = year(dsTe.epochsY);
log = table('Size', [0 3], 'VariableTypes', {'double','double','double'}, 'VariableNames', {'Year','TrainSamples','EpochsRun'});
for Y = testYears
    upYears = max(2019, Y - 2):(Y - 1);
    trMask = ismember(year(epochs), upYears);
    vaMask = ismember(year(epochs), Y - 1) & month(epochs) >= 11;       % last 2 months of Y-1 as monitor
    dsUp = mk(trMask, drivers, 'zmu', cfg.zmu, 'zsig', cfg.zsig, 'stride', opts.stride, 'targets', false, 'mode', opts.mode);
    dsMo = mk(vaMask, drivers, 'zmu', cfg.zmu, 'zsig', cfg.zsig, 'stride', 6, 'mode', opts.mode);
    rowFile = sprintf('roll_%s_%d.mat', opts.tag, Y);
    if isfile(rowFile)
        R = load(rowFile); P = R.P;
    else
        switch kind
            case 'dlnet'
                net = trainResTECNet(base, dsUp, dsMo, lat, 'lambda', cfg.lambda, 'loss', opts.loss, 'maxEpochs', opts.epochs, ...
                    'patience', opts.epochs, 'initialLR', opts.lr, 'dropPeriod', opts.epochs, 'plots', 'none', 'augment', false, 'background', true);
                pred = @(X) predict(net, X);
            case 'fcn'
                p = trainModelFunction(base, opts.modelFcn, dsUp, dsMo, lat, 'lambda', cfg.lambda, 'maxEpochs', opts.epochs, ...
                    'patience', opts.epochs, 'initialLR', opts.lr, 'miniBatch', 32, 'augment', false);
                pred = @(X) gather(extractdata(opts.modelFcn(p, dlarray(single(X), 'SSCB'))));
        end
        idxY = find(yTe == Y);
        P = zeros([dsTe.H, dsTe.W, size(Yhat, 3), numel(idxY)], 'single');
        for i = 1:64:numel(idxY)
            j = min(i + 63, numel(idxY));
            P(:, :, :, i:j) = single(pred(tecGetBatch(dsTe, idxY(i:j))));
        end
        save(rowFile, 'P', '-v7.3');
    end
    idxY = find(yTe == Y);
    Yhat(:, :, :, idxY) = P;
    log(end+1, :) = {Y, dsUp.N, opts.epochs}; %#ok<AGROW>
    fprintf('rollingUpdate %s: year %d done (fine-tuned on %s, %d samples)\n', opts.tag, Y, mat2str(upYears), dsUp.N);
end
dsP = dsTe; dsP.Yhat = Yhat;
if size(Yhat, 3) == 2
    dsP.Yhat = Yhat(:, :, 1, :);
    stats = evaluateTEC([], dsP, lat);
    stats.Sigma = exp(0.5 * double(Yhat(:, :, 2, :))) * dsTe.Tmax;
else
    stats = evaluateTEC([], dsP, lat);
end
end
