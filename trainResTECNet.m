function [net, info] = trainResTECNet(net, dsTrain, dsVal, lat, opts)
%TRAINRESTECNET Train Res-TECNet / plain CNN / Res-TECNet-U with trainnet.
%   [net, info] = TRAINRESTECNET(net, dsTrain, dsVal, lat, ...)
%     'lambda'     gradient-loss weight                (default 0.1)
%     'loss'       'mse' (Eq. 6) | 'nll' (Eq. 7, uncertainty head)
%     'maxEpochs'  (default 200)   'miniBatch' (default 64)
%     'initialLR'  (default 1e-3)  'patience'  (default 15)
%     'dropPeriod' LR halving period in epochs (default 20)
%     'background' assemble mini-batches in the thread-based background
%                  pool (default false; try true for speed once memory is
%                  known to be sufficient)
%     'augment'    random circular longitude shifts    (default true)
%     'seed'       rng seed for reproducibility        (default 1)
%     'plots'      'training-progress' | 'none'
%     'checkpointDir'  folder for per-epoch checkpoints; if it already
%                  holds checkpoints, training RESUMES from the latest one
%                  (power-failure safe). '' = no checkpoints.
%   Sect. 3.3 of the revised manuscript: Adam, LR halved every 20 epochs,
%   He initialization, early stopping on the validation loss (patience
%   15), best-validation network retained.

arguments
    net dlnetwork
    dsTrain struct
    dsVal struct
    lat (:,1) double
    opts.lambda (1,1) double = 0.1
    opts.loss (1,:) char {mustBeMember(opts.loss, {'mse','nll'})} = 'mse'
    opts.maxEpochs (1,1) double = 200
    opts.miniBatch (1,1) double = 64
    opts.initialLR (1,1) double = 1e-3
    opts.patience (1,1) double = 15
    opts.dropPeriod (1,1) double = 20
    opts.background (1,1) logical = false
    opts.augment (1,1) logical = true
    opts.seed (1,1) double = 1
    opts.plots (1,:) char = 'training-progress'
    opts.checkpointDir (1,:) char = ''
end

rng(opts.seed);
wlat = makeLatWeights(lat);

% on-demand sample assembly (tecGetBatch) through transformed index datastores
cdsTr = transform(arrayDatastore((1:dsTrain.N)'), @(c) readSample(dsTrain, c{1}, opts.augment));
cdsVa = transform(arrayDatastore((1:dsVal.N)'),   @(c) readSample(dsVal,   c{1}, false));

switch opts.loss
    case 'mse', lossFcn = @(Y, T) lossTEC(Y, T, wlat, opts.lambda);
    case 'nll', lossFcn = @(Y, T) lossTECNLL(Y, T, wlat, opts.lambda);
end

nIter = ceil(dsTrain.N / opts.miniBatch);
preEnv = 'serial';
if opts.background, preEnv = 'background'; end   % thread pool; the dataset (ds.Ts) is shared, not copied per worker
options = trainingOptions('adam', ...
    'InitialLearnRate', opts.initialLR, ...
    'LearnRateSchedule', 'piecewise', ...
    'LearnRateDropFactor', 0.5, ...
    'LearnRateDropPeriod', opts.dropPeriod, ...
    'MaxEpochs', opts.maxEpochs, ...
    'MiniBatchSize', opts.miniBatch, ...
    'Shuffle', 'every-epoch', ...
    'ValidationData', cdsVa, ...
    'ValidationFrequency', max(1, nIter), ...     % once per epoch
    'ValidationPatience', opts.patience, ...
    'OutputNetwork', 'best-validation', ...
    'ExecutionEnvironment', 'auto', ...
    'PreprocessingEnvironment', preEnv, ...
    'Plots', opts.plots, ...
    'Verbose', true);

% ---- resume from the latest checkpoint, if any ----
epochsDone = 0;
if ~isempty(opts.checkpointDir)
    if ~isfolder(opts.checkpointDir), mkdir(opts.checkpointDir); end
    ck = dir(fullfile(opts.checkpointDir, 'net_checkpoint__*.mat'));
    stateFile = fullfile(opts.checkpointDir, 'resume_state.mat');
    base = 0; if isfile(stateFile), B = load(stateFile); base = B.epochsBase; end
    if ~isempty(ck)
        [~, k] = max([ck.datenum]);                          % newest checkpoint
        it = str2double(regexp(ck(k).name, '__(\d+)__', 'tokens', 'once'));
        C = load(fullfile(opts.checkpointDir, ck(k).name));
        net = C.net; epochsDone = base + round(it / nIter);
        fprintf('trainResTECNet: resuming from %s (%d epochs done)\n', ck(k).name, epochsDone);
        epochsBase = epochsDone; save(stateFile, 'epochsBase');  % iterations restart at 1 below
        save(fullfile(opts.checkpointDir, 'resume_copy.mat'), 'net');     % keep a copy before clearing
        delete(fullfile(opts.checkpointDir, 'net_checkpoint__*.mat'));
        movefile(fullfile(opts.checkpointDir, 'resume_copy.mat'), fullfile(opts.checkpointDir, 'net_checkpoint__0__resume.mat'));
    end
    options.CheckpointPath = opts.checkpointDir;
    options.CheckpointFrequency = 1;
    options.CheckpointFrequencyUnit = 'epoch';
end
remaining = opts.maxEpochs - epochsDone;
if remaining <= 0
    info = struct('resumed', true, 'epochsDone', epochsDone); return
end
options.MaxEpochs = remaining;
options.InitialLearnRate = opts.initialLR * 0.5^floor(epochsDone / opts.dropPeriod);
[net, tinfo] = trainnet(cdsTr, net, lossFcn, options);
info = struct('trainnet', tinfo, 'epochsDone', epochsDone);
try, info.epochsDone = epochsDone + max(tinfo.TrainingHistory.Epoch); catch, end
try, info.valLoss = tinfo.ValidationHistory.Loss; catch, end
if ~isempty(opts.checkpointDir)
    epochsBase = info.epochsDone; save(fullfile(opts.checkpointDir, 'resume_state.mat'), 'epochsBase');
end
end

function out = readSample(ds, idx, augment)
% one observation {X, Y} for trainnet; random circular longitude shift
[X, Y] = tecGetBatch(ds, idx);
if augment
    s = randi(size(X, 2)) - 1;
    if s > 0, X = circshift(X, s, 2); Y = circshift(Y, s, 2); end
end
out = {X, Y};
end
