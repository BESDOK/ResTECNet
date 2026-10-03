function [params, info] = trainModelFunction(params, modelFcn, dsTrain, dsVal, lat, opts)
%TRAINMODELFUNCTION Custom Adam loop for the model-function baselines.
%   [params, info] = TRAINMODELFUNCTION(params, modelFcn, dsTrain, dsVal, lat, ...)
%     params   : learnable struct from initConvLSTM / initTransformer
%     modelFcn : @(params, X) -> dlarray H x W x 1 x B (SSCB, scaled TEC)
%   Options identical to trainResTECNet ('lambda', 'maxEpochs', 'miniBatch',
%   'initialLR', 'patience', 'augment', 'seed') so that the ConvLSTM and
%   Transformer baselines use the same loss (Eq. 6), optimizer, schedule
%   and early-stopping rule as Res-TECNet (Sect. 3.6, vi-vii).
%   info.valLoss (per epoch), info.bestEpoch, info.numParams.

arguments
    params struct
    modelFcn function_handle
    dsTrain struct
    dsVal struct
    lat (:,1) double
    opts.lambda (1,1) double = 0.1
    opts.maxEpochs (1,1) double = 200
    opts.miniBatch (1,1) double = 32
    opts.initialLR (1,1) double = 1e-3
    opts.patience (1,1) double = 15
    opts.augment (1,1) logical = true
    opts.seed (1,1) double = 1
    opts.valSubset (1,1) double = 2000     % validation samples used per epoch
    opts.checkpointFile (1,:) char = ''    % .mat written every epoch; resumed if present
end

rng(opts.seed);
wlat = makeLatWeights(lat);
Ntr = dsTrain.N; Nva = dsVal.N;
useGPU = canUseGPU();
info.numParams = countParams(params);
fprintf('trainModelFunction: %d learnable parameters, %d training samples\n', info.numParams, Ntr);

avgG = []; avgSqG = []; iter = 0;
best = inf; bestParams = params; bestEpoch = 0; wait = 0;
info.valLoss = nan(opts.maxEpochs, 1);
lr = opts.initialLR;
valIdx = sort(randperm(Nva, min(Nva, opts.valSubset)));
startEpoch = 1;
if ~isempty(opts.checkpointFile) && isfile(opts.checkpointFile)
    C = load(opts.checkpointFile);
    params = C.params; avgG = C.avgG; avgSqG = C.avgSqG; iter = C.iter; best = C.best;
    bestParams = C.bestParams; bestEpoch = C.bestEpoch; wait = C.wait; lr = C.lr;
    info.valLoss(1:numel(C.valLoss)) = C.valLoss; startEpoch = C.epoch + 1;
    fprintf('trainModelFunction: resuming after epoch %d from %s\n', C.epoch, opts.checkpointFile);
end

for epoch = startEpoch:opts.maxEpochs
    if epoch > 1 && mod(epoch - 1, 20) == 0, lr = lr * 0.5; end
    order = randperm(Ntr);
    for i = 1:opts.miniBatch:Ntr
        idx = order(i : min(i + opts.miniBatch - 1, Ntr));
        [X, Y] = tecGetBatch(dsTrain, idx);
        if opts.augment
            s = randi(size(X, 2)) - 1;
            X = circshift(X, s, 2); Y = circshift(Y, s, 2);
        end
        X = dlarray(X, 'SSCB'); Y = dlarray(Y, 'SSCB');
        if useGPU, X = gpuArray(X); Y = gpuArray(Y); end
        [loss, grads] = dlfeval(@modelGradients, params, modelFcn, X, Y, wlat, opts.lambda);
        iter = iter + 1;
        [params, avgG, avgSqG] = adamupdate(params, grads, avgG, avgSqG, iter, lr);
    end
    % ---- validation ----
    vl = 0; nb = 0;
    for i = 1:64:numel(valIdx)
        idx = valIdx(i : min(i + 63, numel(valIdx)));
        [Xv, Yv] = tecGetBatch(dsVal, idx); X = dlarray(Xv, 'SSCB'); Y = dlarray(Yv, 'SSCB');
        if useGPU, X = gpuArray(X); Y = gpuArray(Y); end
        vl = vl + double(gather(extractdata(lossTEC(modelFcn(params, X), Y, wlat, opts.lambda)))) * numel(idx);
        nb = nb + numel(idx);
    end
    vl = vl / nb; info.valLoss(epoch) = vl;
    fprintf('epoch %3d | lr %.2e | train loss %.4e | val loss %.4e\n', epoch, lr, ...
        double(gather(extractdata(loss))), vl);
    if vl < best
        best = vl; bestParams = params; bestEpoch = epoch; wait = 0;
    else
        wait = wait + 1;
    end
    if ~isempty(opts.checkpointFile)
        valLoss = info.valLoss(1:epoch); %#ok<NASGU>
        save(opts.checkpointFile, 'params', 'avgG', 'avgSqG', 'iter', 'best', 'bestParams', ...
            'bestEpoch', 'wait', 'lr', 'valLoss', 'epoch', '-v7.3');
    end
    if wait >= opts.patience, fprintf('Early stopping at epoch %d\n', epoch); break; end
end
params = bestParams; info.bestEpoch = bestEpoch; info.bestValLoss = best;
end

function [loss, grads] = modelGradients(params, modelFcn, X, Y, wlat, lambda)
Yhat = modelFcn(params, X);
loss = lossTEC(Yhat, Y, wlat, lambda);
grads = dlgradient(loss, params);
end

function n = countParams(p)
n = 0;
f = fieldnames(p);
for i = 1:numel(f)
    v = p.(f{i});
    if isstruct(v)
        if isscalar(v), n = n + countParams(v);
        else, for j = 1:numel(v), n = n + countParams(v(j)); end
        end
    elseif iscell(v)
        for j = 1:numel(v), n = n + countParams(struct('x', v{j})); end
    else
        n = n + numel(v);
    end
end
end
