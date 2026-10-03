%RUNEXPERIMENTS Driver script for all experiments of the REVISED manuscript.
%   Reproduces (numbers refer to the revised manuscript AISR-D-26-01824 R1):
%     Table 1  hyperparameter selection            (hyperparameterSearch)
%     Table 2  overall accuracy, all baselines      (Sect. 4.2)
%     Table 3  storm-time accuracy                  (Sect. 4.3)
%     Table 4  16 IGS stations by regime + UQRG     (Sect. 4.6)
%     Table 5  six Anatolian stations, gradients    (Sect. 4.7)
%     Table 6  slant-delay and SPP positioning      (Sect. 4.7.1)
%     Table 7  architectural + per-driver ablation  (Sect. 4.8)
%     Sect. 4.8.1 interpolation influence; Sect. 4.9 calibration
%     Figs. 1-6 (fig4 forecast maps, fig1-3 May 2024, fig0 error map)
%   All numeric results are written to results/*.csv and results/all.mat.
%   MATLAB R2024a; Deep Learning Toolbox; Navigation Toolbox for the
%   RINEX-based sections (9-10). All GNSS data come from open archives
%   (IGS/EPN via CDDIS and BKG); no licensed network data are needed.
%   Set the flags below to run subsets.

%% ---------------- 0. Configuration ----------------
try, g_ = gpuDevice; fprintf('GPU: %s, %.1f GB free\n', g_.Name, g_.AvailableMemory / 1e9); catch, fprintf(2, '!!! NO GPU available: network inference and training would fall back to the CPU and be extremely slow\n'); end
fprintf('[%s] %s\n', datestr(now, 'HH:MM:SS'), '0. Configuration');
if ~usejava('desktop'), diary('log_experiments.txt'); end          % unattended (matlab -batch): keep a log
cfg.L = 12; cfg.Delta = 24; cfg.Tmax = 200; cfg.lambda = 0.1; cfg.B = 8; cfg.filters = 64;
cfg.posStation = 'MERS';
cfg.stationYears = 2023:2024;      % receiver-level validation window (every 3rd day + all storm-window days)
cfg.stationDayStride = 3;
cfg.rinexReader = 'rinexread';    % 'rinexread' (default; 2.5 s per daily file) | 'fast' (validated identical but slower: 4.4 s)
cfg.netrc = fullfile(pwd, 'earthdata_netrc');   % Earthdata credentials for CDDIS (GIMs + RINEX)
if ~isfile(cfg.netrc), error('runExperiments:netrc', 'Create %s (machine urs.earthdata.nasa.gov login U password P)', cfg.netrc); end   % open IGS station nearest to Kayseri (~260 km) for the slant / SPP experiments
trainYears = 1999:2018; valYears = 2019:2020; testYears = 2021:2024;
do.download = true;  do.hyper = true;   do.trainMain = true; do.baselines = true;
stationsReady = isfile('STATIONS_DONE.txt');                 % written by runStationProcessing when all station-months exist
do.stations = stationsReady; do.anatolia = stationsReady; do.ablation = true; do.uncertainty = true;
do.rolling = false;                % annual-update evaluation: tested, no gain (kept for the record)
do.adaptive = true;                % Res-TECNet-A: SH-AR adaptive base + learned correction (final model, Sect. 3.2)
resDir = 'results'; if ~isfolder(resDir), mkdir(resDir); end
% ---- compute budget (Sect. 3.3 / Table 1 must quote the values actually used) ----
budget.preset = 'overnight';     % 'overnight' (~9-10 GPU-h) | 'fast' (~22 GPU-h) | 'full' (days)
switch budget.preset
    case 'overnight'
        % second pass (26 Sep): the first 15 epochs left the network at persistence level
        % (validation loss still falling); continue from the checkpoints with a slower
        % LR decay and WITHOUT the random longitude-shift augmentation (it breaks the
        % consistency between the map and the UT/DOY channels).
        budget.trainStride = 6;   budget.maxEpochs = 40; budget.patience = 8; budget.dropPeriod = 12;
        budget.augment     = false;
        budget.warmStart   = true;  budget.fineEpochs = 8;
        budget.baselineYears  = 2011:2018;
        budget.baselineStride = 6;  budget.baselineEpochs = 8;
        budget.ablStride   = 6;   budget.ablEpochs = 8;  budget.ablTrainYears = 2011:2018;
        budget.ablRows     = {'full','persistBase','noGlobalSkip','zeroPad','noDrivers','noTargetBlock','noKp','noDst','noisyTarget','operational'};   % plainCNN: see Table 2 (from scratch)
        budget.hyper       = struct('L', [], 'lambda', [0 0.05 0.1 0.2 0.5], 'B', [4 12 16], 'filters', []);   % rows that admit a warm start; lambda = 0.1 is the retained configuration (B = 8); finished rows are cached
    case 'fast'
        budget.trainStride = 4;   budget.maxEpochs = 25; budget.patience = 6; budget.dropPeriod = 8;
        budget.augment     = false;
        budget.warmStart   = true;  budget.fineEpochs = 10;
        budget.baselineYears  = 1999:2018; budget.baselineStride = 4; budget.baselineEpochs = 10;
        budget.ablStride   = 6;   budget.ablEpochs = 10; budget.ablTrainYears = 2011:2018;
        budget.ablRows     = {'full','noGlobalSkip','plainCNN','zeroPad','noDrivers','noTargetBlock','noKp','noDst','noisyTarget','operational'};
        budget.hyper       = struct('L', [6 24], 'lambda', [0 0.5], 'B', [4 16], 'filters', [32 128]);
    case 'full'
        budget.trainStride = 3;   budget.maxEpochs = 40; budget.patience = 8; budget.dropPeriod = 10;
        budget.augment     = false;
        budget.warmStart   = false; budget.fineEpochs = 40;
        budget.baselineYears  = 1999:2018; budget.baselineStride = 3; budget.baselineEpochs = 40;
        budget.ablStride   = 6;   budget.ablEpochs = 15; budget.ablTrainYears = 1999:2018;
        budget.ablRows     = {};
        budget.hyper       = struct('L', [6 12 24 48], 'lambda', [0 0.05 0.1 0.2 0.5], 'B', [4 8 12 16], 'filters', [32 64 128]);
end
trainOpts    = {'plots', 'none', 'maxEpochs', budget.maxEpochs, 'patience', budget.patience, 'dropPeriod', budget.dropPeriod, 'background', true, 'augment', budget.augment};   % 32 GB RAM: background batch assembly
trainOptsFT  = {'plots', 'none', 'maxEpochs', budget.fineEpochs, 'patience', 4, 'dropPeriod', 4, 'initialLR', 3e-4, 'background', true, 'augment', budget.augment};
trainOptsAbl = {'plots', 'none', 'maxEpochs', budget.ablEpochs, 'patience', 4, 'dropPeriod', 4, 'initialLR', 3e-4, 'background', true, 'augment', budget.augment};   % fine-tuning budget (warm start)

%% ---------------- 1. Databases and drivers ----------------
fprintf('[%s] %s\n', datestr(now, 'HH:MM:SS'), '1. Databases and drivers');
if do.download                                    % each database is built only once
    if ~isfile('TECdatabase.mat'),  buildTECDatabase('1999-01-01', '2024-12-31', 'TECdatabase.mat',  'ionex_cache', 'product', 'CODG', 'cddisNetrc', cfg.netrc); end
    if ~isfile('C1PGdatabase.mat'), buildTECDatabase('2021-01-01', '2024-12-31', 'C1PGdatabase.mat', 'ionex_cache', 'product', 'C1PG', 'fillGaps', false, 'cddisNetrc', cfg.netrc); end
    if ~isfile('UQRGdatabase.mat'), buildTECDatabase('2021-01-01', '2024-12-31', 'UQRGdatabase.mat', 'ionex_cache', 'product', 'UQRG', 'fillGaps', false, 'cddisNetrc', cfg.netrc); end
end
load('TECdatabase.mat', 'TEC', 'epochs', 'lat', 'lon', 'badMask');
TEC = TEC / single(cfg.Tmax);                   % scaled once, in place (~4.6 GB); every dataset references it without copying
[H, W, ~] = size(TEC); C = cfg.L + 16; cfg.H = H; cfg.W = W; cfg.C = C;
assert(W == 72, 'Grid must contain the 72 unique meridians (readIONEXTEC dropDuplicateLon).');
if ~isfile('drivers_hindcast.mat')
    drivers = buildDriverTable(epochs, 'index_cache', 'mode', 'hindcast', 'meanType', 'trailing');
    save('drivers_hindcast.mat', 'drivers');
    driversOp = buildDriverTable(epochs, 'index_cache', 'mode', 'operational', 'meanType', 'trailing');
    save('drivers_operational.mat', 'driversOp');
else
    load('drivers_hindcast.mat', 'drivers'); load('drivers_operational.mat', 'driversOp');
end

%% ---------------- 2. Datasets ----------------
fprintf('[%s] %s\n', datestr(now, 'HH:MM:SS'), '2. Datasets');
splits.tr = ismember(year(epochs), trainYears);
splits.va = ismember(year(epochs), valYears);
splits.te = ismember(year(epochs), testYears);
mk = @(mask, drv, varargin) buildTECDataset(TEC, epochs, drv, 'L', cfg.L, 'Delta', cfg.Delta, 'Tmax', cfg.Tmax, ...
    'badMask', badMask, 'mask', mask, 'scaled', true, varargin{:});
dsTr = mk(splits.tr, drivers, 'stride', budget.trainStride, 'targets', false);
cfg.zmu = dsTr.zmu; cfg.zsig = dsTr.zsig;
dsVa = mk(splits.va, drivers, 'zmu', cfg.zmu, 'zsig', cfg.zsig);
dsTe = mk(splits.te, drivers, 'zmu', cfg.zmu, 'zsig', cfg.zsig);
% operational datasets (Sect. 3.4): archived Kp forecasts exist from 2011
dsTrOp = mk(splits.tr, driversOp, 'zmu', cfg.zmu, 'zsig', cfg.zsig, 'mode', 'operational', 'stride', budget.trainStride, 'targets', false);
dsVaOp = mk(splits.va, driversOp, 'zmu', cfg.zmu, 'zsig', cfg.zsig, 'mode', 'operational');
dsTeOp = mk(splits.te, driversOp, 'zmu', cfg.zmu, 'zsig', cfg.zsig, 'mode', 'operational');
DstAll = drivers.Dst;
DstTest = DstAll(dsTe.idxT + dsTe.Delta);                                  % aligned with epochsY (indices into the full record)

%% ---------------- 4. Main models ----------------
fprintf('[%s] %s\n', datestr(now, 'HH:MM:SS'), '4. Main models');
% Every trained network is saved to its own .mat and skipped on restart;
% trainResTECNet resumes from per-epoch checkpoints in checkpoints/<tag>.
% A saved network is reused; if budget.maxEpochs exceeds the epochs already
% trained (info.epochsDone), training CONTINUES from the checkpoints.
[net, info, cont] = trainOrContinue('resTECNet_hindcast.mat', 'net', @() createResTECNet(H, W, C, 'L', cfg.L), ...
    dsTr, dsVa, lat, cfg, trainOpts, 'checkpoints/hindcast', budget.maxEpochs);
if cont            % hindcast changed -> redo the fine-tuned networks
    if isfile('resTECNet_operational.mat'), delete('resTECNet_operational.mat'); end
    if isfile('resTECNet_uncertainty.mat'), delete('resTECNet_uncertainty.mat'); end
    if isfolder('checkpoints/operational'), rmdir('checkpoints/operational', 's'); end
    if isfolder('checkpoints/uncertainty'), rmdir('checkpoints/uncertainty', 's'); end
end
if isfile('resTECNet_operational.mat')
    load('resTECNet_operational.mat', 'netOp');
else
    netOp = createResTECNet(H, W, C, 'L', cfg.L);
    if budget.warmStart, netOp = transferWeights(netOp, net); to = trainOptsFT; else, to = trainOpts; end
    netOp = trainResTECNet(netOp, dsTrOp, dsVaOp, lat, 'lambda', cfg.lambda, to{:}, 'checkpointDir', 'checkpoints/operational');
    save('resTECNet_operational.mat', 'netOp');
end
if do.uncertainty
    if isfile('resTECNet_uncertainty.mat')
        load('resTECNet_uncertainty.mat', 'netU');
    else
        netU = createResTECNet(H, W, C, 'L', cfg.L, 'uncertainty', true);
        if budget.warmStart, netU = transferWeights(netU, net); to = trainOptsFT; else, to = trainOpts; end
        netU = trainResTECNet(netU, dsTr, dsVa, lat, 'lambda', cfg.lambda, 'loss', 'nll', to{:}, 'checkpointDir', 'checkpoints/uncertainty');
        save('resTECNet_uncertainty.mat', 'netU');
    end
end
statsRes = evaluateTEC(net, dsTe, lat);
statsOp  = evaluateTEC(netOp, dsTeOp, lat);

%% ---------------- 4A. Res-TECNet-A: adaptive base + learned correction ----------------
fprintf('[%s] %s\n', datestr(now, 'HH:MM:SS'), '4A. Res-TECNet-A: adaptive base + learned correction');
% The SH-AR 24-h forecast (computed once for the whole record, from past data
% only) is appended as input channel L+17 and replaces T_t in the global skip.
if do.adaptive
    if isfile('SHARforecasts.mat')
        load('SHARforecasts.mat', 'Fsh');
    else
        Fsh = computeSHARForecasts(TEC, epochs, lat, lon, 'Delta', cfg.Delta, 'order', 48);
        save('SHARforecasts.mat', 'Fsh', '-v7.3');
    end
    mkA = @(mask, drv, varargin) mk(mask, drv, 'base', Fsh, varargin{:});
    dsTrA = mkA(splits.tr, drivers, 'zmu', cfg.zmu, 'zsig', cfg.zsig, 'stride', budget.trainStride, 'targets', false);
    dsVaA = mkA(splits.va, drivers, 'zmu', cfg.zmu, 'zsig', cfg.zsig);
    dsTeA = mkA(splits.te, drivers, 'zmu', cfg.zmu, 'zsig', cfg.zsig);
    dsTrOpA = mkA(splits.tr, driversOp, 'zmu', cfg.zmu, 'zsig', cfg.zsig, 'mode', 'operational', 'stride', budget.trainStride, 'targets', false);
    dsVaOpA = mkA(splits.va, driversOp, 'zmu', cfg.zmu, 'zsig', cfg.zsig, 'mode', 'operational');
    dsTeOpA = mkA(splits.te, driversOp, 'zmu', cfg.zmu, 'zsig', cfg.zsig, 'mode', 'operational');
    assert(isequal(dsTeA.epochsY, dsTe.epochsY) && isequal(dsTeOpA.epochsY, dsTe.epochsY) && isequal(dsTeOp.epochsY, dsTe.epochsY), ...
        'Test epochs of the adaptive-base datasets differ from the baseline datasets: results would be misaligned.');
    CA = dsTrA.C; cfg.CA = CA; cfg.baseChannel = CA;
    trainOptsA = {'plots', 'none', 'maxEpochs', 15, 'patience', 5, 'dropPeriod', 6, 'initialLR', 5e-4, 'background', true, 'augment', false};
    [netA, infoA, ~] = trainOrContinue('resTECNetA_hindcast.mat', 'netA', ...
        @() transferWeights(createResTECNet(H, W, CA, 'L', cfg.L, 'baseChannel', CA), net), ...   % warm start from the static net
        dsTrA, dsVaA, lat, cfg, trainOptsA, 'checkpoints/A_hindcast', 15);
    if isfile('resTECNetA_operational.mat'), load('resTECNetA_operational.mat', 'netOpA');
    else
        netOpA = transferWeights(createResTECNet(H, W, CA, 'L', cfg.L, 'baseChannel', CA), netA);
        netOpA = trainResTECNet(netOpA, dsTrOpA, dsVaOpA, lat, 'lambda', cfg.lambda, trainOptsFT{:}, 'checkpointDir', 'checkpoints/A_operational');
        save('resTECNetA_operational.mat', 'netOpA');
    end
    if do.uncertainty
        if isfile('resTECNetA_uncertainty.mat'), load('resTECNetA_uncertainty.mat', 'netUA');
        else
            netUA = transferWeights(createResTECNet(H, W, CA, 'L', cfg.L, 'baseChannel', CA, 'uncertainty', true), netA);
            netUA = trainResTECNet(netUA, dsTrA, dsVaA, lat, 'lambda', cfg.lambda, 'loss', 'nll', trainOptsFT{:}, 'checkpointDir', 'checkpoints/A_uncertainty');
            save('resTECNetA_uncertainty.mat', 'netUA');
        end
    end
    statsResA = evaluateTEC(netA, dsTeA, lat);
    statsOpA  = evaluateTEC(netOpA, dsTeOpA, lat);
    % SH-AR baseline directly from the same forecast field (identical method to baselineSHAR)
    dsP = dsTe; dsP.Yhat = reshape(Fsh(:, :, dsTe.idxT + dsTe.Delta), H, W, 1, dsTe.N);
    statsSHARfield = evaluateTEC([], dsP, lat);
end

%% ---------------- 5. Baselines (Table 2) ----------------
fprintf('[%s] %s\n', datestr(now, 'HH:MM:SS'), '5. Baselines (Table 2)');
statsPers = persistenceBaseline(dsTe, lat);
statsC1PG = baselineC1PG('C1PGdatabase.mat', dsTe, lat);
leadIn = splits.te | ismember(year(epochs), valYears(end));            % 30-day lead-in for SH-AR
statsSHAR = baselineSHAR(TEC(:,:,leadIn) * cfg.Tmax, epochs(leadIn), dsTe, lat, lon, 'order', 48);
if do.baselines
    if false                                   % (no combined cache: each baseline has its own file)
    else
        netPlain = trainOrContinue('baseline_plain.mat', 'netPlain', ...
            @() createResTECNet(H, W, C, 'L', cfg.L, 'residual', false, 'globalSkip', false), ...
            dsTr, dsVa, lat, cfg, trainOpts, 'checkpoints/plain', budget.maxEpochs);
        splitsB = splits; splitsB.tr = ismember(year(epochs), budget.baselineYears);
        dsTrB = mk(splitsB.tr, drivers, 'zmu', cfg.zmu, 'zsig', cfg.zsig, 'stride', budget.baselineStride, 'targets', false);
        [pCL, cfgCL] = initConvLSTM(C, cfg.L, 64, 2);
        fCL = @(p, X) modelConvLSTM(p, X, cfgCL);
        if isfile('baseline_convlstm.mat'), load('baseline_convlstm.mat', 'pCL', 'cfgCL');
        else
            [pCL, infoCL] = trainModelFunction(pCL, fCL, dsTrB, dsVa, lat, 'lambda', cfg.lambda, 'maxEpochs', budget.baselineEpochs, ...
                'patience', budget.patience, 'miniBatch', 32, 'checkpointFile', 'checkpoints/convlstm_ckpt.mat');
            save('baseline_convlstm.mat', 'pCL', 'cfgCL', 'infoCL');
        end
        [pTR, cfgTR] = initTransformer(H, W, C, cfg.L, 'd', 256, 'numLayers', 6, 'numHeads', 8);
        fTR = @(p, X) modelTransformer(p, X, cfgTR);
        if isfile('baseline_transformer.mat'), load('baseline_transformer.mat', 'pTR', 'cfgTR');
        else
            [pTR, infoTR] = trainModelFunction(pTR, fTR, dsTrB, dsVa, lat, 'lambda', cfg.lambda, 'maxEpochs', budget.baselineEpochs, ...
                'patience', budget.patience, 'miniBatch', 32, 'checkpointFile', 'checkpoints/transformer_ckpt.mat');
            save('baseline_transformer.mat', 'pTR', 'cfgTR', 'infoTR');
        end
    end
    fCL = @(p, X) modelConvLSTM(p, X, cfgCL); fTR = @(p, X) modelTransformer(p, X, cfgTR);
    statsPlain = evaluateTEC(netPlain, dsTe, lat);
    predCL = @(X) gather(extractdata(fCL(pCL, dlarray(single(X), 'SSCB'))));
    statsCL = evaluateTEC(predCL, dsTe, lat);
    predTR = @(X) gather(extractdata(fTR(pTR, dlarray(single(X), 'SSCB'))));
    statsTR = evaluateTEC(predTR, dsTe, lat);
end
% IRI-2020 (climatological reference): sample externally onto the grid at
% dsTe.epochsY -> YhatIRI (H x W x 1 x N, TECU); dsI = dsTe; dsI.Yhat = single(YhatIRI/cfg.Tmax);
% statsIRI = evaluateTEC([], dsI, lat);
%% ---------------- 5b. Annually updated models (Sect. 3.3) ----------------
fprintf('[%s] %s\n', datestr(now, 'HH:MM:SS'), '5b. Annually updated models (Sect. 3.3)');
% Every learned model is fine-tuned once per test year on the two preceding
% years (>= 2019, never on the test year), mirroring operational retraining
% and the adaptivity of SH-AR / C1PG. Static and updated results are both reported.
if do.rolling
    statsResR = rollingUpdate(net,   'dlnet', mk, epochs, drivers,   dsTe,   lat, cfg, testYears, 'tag', 'hindcast');
    statsOpR  = rollingUpdate(netOp, 'dlnet', mk, epochs, driversOp, dsTeOp, lat, cfg, testYears, 'tag', 'operational', 'mode', 'operational');
    statsPlainR = rollingUpdate(netPlain, 'dlnet', mk, epochs, drivers, dsTe, lat, cfg, testYears, 'tag', 'plain');
    statsCLR = rollingUpdate(pCL, 'fcn', mk, epochs, drivers, dsTe, lat, cfg, testYears, 'tag', 'convlstm', 'modelFcn', fCL);
    statsTRR = rollingUpdate(pTR, 'fcn', mk, epochs, drivers, dsTe, lat, cfg, testYears, 'tag', 'transformer', 'modelFcn', fTR);
    if do.uncertainty
        statsUR = rollingUpdate(netU, 'dlnet', mk, epochs, drivers, dsTe, lat, cfg, testYears, 'tag', 'uncertainty', 'loss', 'nll');
    end
end

models = {'IRI-2020', 'Persistence-24h', 'CODE C1PG', 'SH-AR', 'Plain CNN', 'ConvLSTM', ...
          'Transformer', 'Res-TECNet operational', 'Res-TECNet hindcast'};
S = {[], statsPers, statsC1PG, statsSHAR, statsPlain, statsCL, statsTR, statsOp, statsRes};
if do.adaptive
    models = [models, {'Res-TECNet-A operational', 'Res-TECNet-A hindcast'}];
    S = [S, {statsOpA, statsResA}];
    iSH = find(strcmp(models, 'SH-AR'), 1);
    relDiff = abs(statsSHARfield.RMSE - S{iSH}.RMSE) / S{iSH}.RMSE;
    fprintf('SH-AR cross-check: forecast field %.3f vs independent baselineSHAR %.3f TECU (%.1f %% apart)\n', statsSHARfield.RMSE, S{iSH}.RMSE, 100 * relDiff);
    if relDiff > 0.03, warning('SH-AR implementations disagree by more than 3 %% -- check computeSHARForecasts / baselineSHAR.'); end
    S{iSH} = statsSHARfield;                          % same forecast field as the base channel
end
if do.rolling
    models = [models, {'Plain CNN (annual update)', 'ConvLSTM (annual update)', 'Transformer (annual update)', ...
                       'Res-TECNet operational (annual update)', 'Res-TECNet hindcast (annual update)'}];
    S = [S, {statsPlainR, statsCLR, statsTRR, statsOpR, statsResR}];
end
if exist('statsIRI', 'var'), S{1} = statsIRI; end
keep = ~cellfun(@isempty, S); models = models(keep); S = S(keep);
T2 = table(models(:), cellfun(@(s) s.RMSE, S)', cellfun(@(s) s.MAE, S)', ...
    cellfun(@(s) s.Bias, S)', cellfun(@(s) s.Corr, S)', 'VariableNames', {'Model','RMSE','MAE','Bias','Corr'});
bud = repmat({'-'}, numel(models), 1);
bud(ismember(models, {'Plain CNN','Res-TECNet operational','Res-TECNet hindcast'})) = {sprintf('%d-%d, %d-h stride, %d ep', ...
    trainYears(1), trainYears(end), budget.trainStride, budget.maxEpochs)};
bud(ismember(models, {'ConvLSTM','Transformer'})) = {sprintf('%d-%d, %d-h stride, %d ep', ...
    budget.baselineYears(1), budget.baselineYears(end), budget.baselineStride, budget.baselineEpochs)};
bud(contains(models, 'annual update')) = {'static + yearly fine-tune on the 2 preceding years (4 ep)'};
bud(contains(models, 'Res-TECNet-A')) = {'static net + SH-AR base channel, 15 ep fine-tune'};
T2.TrainingBudget = bud;
% Res-TECNet under the baselines' budget (for a like-for-like comparison) is
% the 'full' row of Table 7 (driverAblation, same window/stride/epochs);
% report it in the Table 2 footnote.
disp(T2); writetable(T2, fullfile(resDir, 'table2_overall.csv'));

%% ---------------- 6. Storm-time (Table 3) ----------------
fprintf('[%s] %s\n', datestr(now, 'HH:MM:SS'), '6. Storm-time (Table 3)');
[T3, stormMask, quietMask, events] = stormAnalysis(S, models, dsTe.epochsY, DstTest, -100, 48);
may = abs(hours(dsTe.epochsY - datetime(2024,5,11,2,0,0,'TimeZone','UTC'))) <= 48;
T3.May2024_TECU = cellfun(@(s) sqrt(mean(s.perEpochRMSE(may).^2)), S)';
disp(T3); writetable(T3, fullfile(resDir, 'table3_storm.csv'));
% per-year RMSE and relative error (Sect. 4.2): global-mean TEC of the year as denominator
yrs = unique(year(dsTe.epochsY))'; Ty = table(yrs', 'VariableNames', {'Year'});
meanTEC = arrayfun(@(y) mean(dsTe.Y(:,:,:,year(dsTe.epochsY) == y), 'all') * cfg.Tmax, yrs)';
Ty.MeanTEC_TECU = meanTEC;
for k = 1:numel(models)
    Ty.(matlab.lang.makeValidName(models{k})) = arrayfun(@(y) sqrt(mean(S{k}.perEpochRMSE(year(dsTe.epochsY) == y).^2, 'omitnan')), yrs)';
end
disp(Ty); writetable(Ty, fullfile(resDir, 'table2b_per_year.csv'));

%% ---------------- 6b. Reference results for the remaining sections ----------------
fprintf('[%s] %s\n', datestr(now, 'HH:MM:SS'), '6b. Reference results for the remaining sections');
% Figures, station comparison and positioning use the annually updated networks when available
if do.rolling, statsRes = statsResR; statsOp = statsOpR; if do.uncertainty, statsU_forCal = statsUR; end, end
if do.adaptive                      % the final model of the revised manuscript
    statsRes = statsResA; statsOp = statsOpA; net = netA; dsTe = dsTeA; dsTeOp = dsTeOpA;
    if do.uncertainty, statsU_forCal = evaluateTEC(netUA, dsTeA, lat); end
end

%% ---------------- 7. Figures: forecast maps (Fig. 2) and May 2024 (Figs. 3-5) ----------------
fprintf('[%s] %s\n', datestr(now, 'HH:MM:SS'), '7. Figures: forecast maps (Fig. 2) and May 2024 (Figs. 3-5)');
plotForecastMaps(statsRes, dsTe, lat, lon, 'outFile', 'fig4.jpg');
mayWin = epochs >= datetime(2024,5,6,'TimeZone','UTC') & epochs < datetime(2024,5,15,'TimeZone','UTC');
fc = struct('name', {'Res-TECNet 24-h forecast (hindcast)', 'Res-TECNet 24-h forecast (operational)'}, ...
            'epochs', {dsTe.epochsY, dsTeOp.epochsY}, 'Yhat', {statsRes.Yhat, statsOp.Yhat}, ...
            'color', {[0.85 0.1 0.1], [1.0 0.55 0.0]});
stationObs = [];  if isfile(sprintf('station_%s.mat', cfg.posStation)), Lst = load(sprintf('station_%s.mat', cfg.posStation)); stationObs = Lst.obsAll; end
Sst = stationList('anatolia'); iS = find(Sst.code == cfg.posStation, 1);
may24 = visualizeMay2024(TEC(:,:,mayWin) * cfg.Tmax, epochs(mayWin), lat, lon, drivers(mayWin,:), ...
    'forecasts', fc, 'station', stationObs, 'stationName', cfg.posStation, 'stationLat', Sst.lat(iS), 'stationLon', Sst.lon(iS));

%% ---------------- 8. Spatial error map (Fig. 6) ----------------
fprintf('[%s] %s\n', datestr(now, 'HH:MM:SS'), '8. Spatial error map (Fig. 6)');
plotErrorMap(statsRes, dsTe, lat, lon, 'stations', stationList('all'), 'outFile', 'fig0.jpg');

%% ---------------- 9. Independent observations (Table 4, Sect. 4.6) ----------------
fprintf('[%s] %s\n', datestr(now, 'HH:MM:SS'), '9. Independent observations (Table 4, Sect. 4.6)');
% ---- station-day list: every stationDayStride-th day of cfg.stationYears plus every day inside a storm window ----
allDays = (datetime(cfg.stationYears(1),1,1,'TimeZone','UTC') : datetime(cfg.stationYears(end),12,31,'TimeZone','UTC'))';
stormDays = unique(dateshift(dsTe.epochsY(stormMask), 'start', 'day'));
stationDays = union(allDays(1:cfg.stationDayStride:end), stormDays(ismember(year(stormDays), cfg.stationYears)));
fprintf('Station validation: %d days per station (%d storm-window days)\n', numel(stationDays), nnz(ismember(stationDays, stormDays)));
if do.stations
  try
    Sigs = stationList('igs'); stations = struct('code', {}, 'regime', {}, 'obs', {});
    dcbCache = containers.Map();
    for s = 1:height(Sigs)
        obsAll = processStationYears(char(Sigs.code(s)), stationDays, dcbCache, cfg.rinexReader);
        if height(obsAll) == 0, warning('%s: no observations processed -- station skipped', Sigs.code(s)); continue; end
        [obsAll, dInfo] = harmonizeReceiverDCB(obsAll);
        fprintf('%s: %d days (%d without DCB estimate), receiver-DCB day-to-day scatter %.2f TECU RMS (max %.1f), median %.1f ns\n', Sigs.code(s), dInfo.nDays, dInfo.nNoEstimate, dInfo.dailyScatterTECU, dInfo.maxAbsShiftTECU, dInfo.medianDCB_ns);
        stations(end+1) = struct('code', char(Sigs.code(s)), 'regime', char(Sigs.regime(s)), 'obs', obsAll); %#ok<SAGROW>
    end
    idx4 = ismember(models, {'Persistence-24h','CODE C1PG','SH-AR','Transformer','Res-TECNet-A operational','Res-TECNet-A hindcast'});
    if nnz(idx4) < 4, idx4 = ismember(models, {'Persistence-24h','CODE C1PG','Transformer','Res-TECNet operational','Res-TECNet hindcast'}); end
    [T4st, T4, ~, raw4] = evaluateStationVTEC(S(idx4), models(idx4), dsTe, lat, lon, stations, 'stormMask', stormMask);
    nm = T4st.Properties.RowNames; isS = endsWith(nm, ' (storm)');
    low = ~isS & T4st.nHours < 200;
    if any(low), warning('Stations with fewer than 200 valid hours (check data / processing): %s', strjoin(nm(low), ', ')); end
    if T4{'All stations', 'CODEGIM'} > 6, warning('CODE GIM vs station VTEC (offset removed) is %.1f TECU overall: unusually large, inspect the station processing.', T4{'All stations', 'CODEGIM'}); end
    writetable(T4, fullfile(resDir, 'table4_regimes.csv'), 'WriteRowNames', true);
    writetable(T4st, fullfile(resDir, 'table4_stations.csv'), 'WriteRowNames', true);
    writetable(raw4.Treg, fullfile(resDir, 'table4_regimes_raw.csv'), 'WriteRowNames', true);
    writetable(raw4.Tst, fullfile(resDir, 'table4_stations_raw.csv'), 'WriteRowNames', true);
    % global comparison against the independent UQRG product
    U = load('UQRGdatabase.mat', 'TEC', 'epochs');
    [tf, loc] = ismember(dateshift(dsTe.epochsY, 'start', 'hour'), U.epochs);
    Yu = nan([H W 1 numel(tf)]); Yu(:,:,1,tf) = U.TEC(:,:,loc(tf));
    rmseOf = @(st) st.RMSE;
    uq.ResTECNet   = rmseOf(evaluateTEC(net, dsTe, lat, 'RefY', Yu));
    uq.Persistence = rmseOf(evaluateTEC([], setfield(dsTe, 'Yhat', dsTe.Tlast), lat, 'RefY', Yu)); %#ok<SFLD>
    uq.C1PG        = rmseOf(evaluateTEC([], setfield(dsTe, 'Yhat', single(statsC1PG.Yhat / cfg.Tmax)), lat, 'RefY', Yu)); %#ok<SFLD>
    uq.CODEvsUQRG  = rmseOf(evaluateTEC([], setfield(dsTe, 'Yhat', dsTe.Y), lat, 'RefY', Yu)); %#ok<SFLD>
    if exist('statsResA', 'var')                                  % the final model of the manuscript and SH-AR, same reference
        uq.SHAR          = rmseOf(evaluateTEC([], setfield(dsTe, 'Yhat', single(statsSHARfield.Yhat / cfg.Tmax)), lat, 'RefY', Yu)); %#ok<SFLD>
        uq.ResTECNetA_hind = rmseOf(evaluateTEC([], setfield(dsTe, 'Yhat', single(statsResA.Yhat / cfg.Tmax)), lat, 'RefY', Yu)); %#ok<SFLD>
        uq.ResTECNetA_oper = rmseOf(evaluateTEC([], setfield(dsTe, 'Yhat', single(statsOpA.Yhat / cfg.Tmax)), lat, 'RefY', Yu)); %#ok<SFLD>
    end
    disp(uq); save(fullfile(resDir, 'uqrg_global.mat'), 'uq');
  catch ME
    fprintf(2, '!!! Section 9 (stations) failed: %s\n%s\n', ME.message, getReport(ME, 'basic'));
  end
end

%% ---------------- 10. Anatolia (Table 5) + gradients + slant / SPP (Table 6) ----------------
fprintf('[%s] %s\n', datestr(now, 'HH:MM:SS'), '10. Anatolia (Table 5) + gradients + slant / SPP (Table 6)');
if do.anatolia
  try
    % six open-access IGS/EPN stations (ISTA, TUBI, ANKR, MERS, NICO, NSSP); processed like the IGS set
    Sana = stationList('anatolia'); stAna = struct('code', {}, 'regime', {}, 'obs', {});
    if ~exist('dcbCache', 'var'), dcbCache = containers.Map(); end
    for s = 1:height(Sana)
        obsAll = processStationYears(char(Sana.code(s)), stationDays, dcbCache, cfg.rinexReader);
        if height(obsAll) == 0, warning('%s: no observations processed -- station skipped', Sana.code(s)); continue; end
        [obsAll, dInfo] = harmonizeReceiverDCB(obsAll);
        fprintf('%s: %d days (%d without DCB estimate), receiver-DCB day-to-day scatter %.2f TECU RMS (max %.1f), median %.1f ns\n', Sana.code(s), dInfo.nDays, dInfo.nNoEstimate, dInfo.dailyScatterTECU, dInfo.maxAbsShiftTECU, dInfo.medianDCB_ns);
        stAna(end+1) = struct('code', char(Sana.code(s)), 'regime', 'ANATOLIA', 'obs', obsAll); %#ok<SAGROW>
    end
    idx5 = ismember(models, {'IRI-2020','Persistence-24h','CODE C1PG','SH-AR','Res-TECNet-A operational','Res-TECNet-A hindcast'});
    if nnz(idx5) < 3, idx5 = ismember(models, {'IRI-2020','Persistence-24h','CODE C1PG','Res-TECNet operational','Res-TECNet hindcast'}); end
    [T5st, T5, hourly, raw5] = evaluateStationVTEC(S(idx5), models(idx5), dsTe, lat, lon, stAna, 'stormMask', stormMask);
    writetable(T5st, fullfile(resDir, 'table5_anatolia.csv'), 'WriteRowNames', true);
    writetable(raw5.Tst, fullfile(resDir, 'table5_anatolia_raw.csv'), 'WriteRowNames', true);
    writetable(raw5.Treg, fullfile(resDir, 'table5_regimes_raw.csv'), 'WriteRowNames', true);
    Tgrad = regionalGradients(hourly, {'ARUC','AUT1','East-West (ARUC - AUT1)'; 'ISTA','BSHM','North-South (ISTA - BSHM)'});
    writetable(Tgrad, fullfile(resDir, 'table5_gradients.csv'));
    % GIM-based point statistics at Kayseri for comparison (1.83 vs CODE / vs UQRG)
    Tkay = evaluatePointTEC(S(idx5), models(idx5), dsTe, lat, lon, 38.72, 35.49, DstTest);
    writetable(Tkay, fullfile(resDir, 'kayseri_point_vsCODE.csv'), 'WriteRowNames', true);
    if exist('Yu', 'var')
        TkayU = evaluatePointTEC(S(idx5), models(idx5), dsTe, lat, lon, 38.72, 35.49, DstTest, Yu);
        writetable(TkayU, fullfile(resDir, 'kayseri_point_vsUQRG.csv'), 'WriteRowNames', true);
    end
    % slant-delay errors and positioning at the open station nearest to Kayseri (cfg.posStation), 2024
    k = find(string({stAna.code}) == string(cfg.posStation), 1);
    if ~isempty(k)
        obsK = stAna(k).obs; obsK = obsK(year(obsK.Properties.RowTimes) == 2024, :);
        [~, nF] = downloadRINEX(cfg.posStation, datetime(2024,5,10,'TimeZone','UTC'), 'rinex_cache', 'cddisNetrc', cfg.netrc, 'country', 'TUR');
        [a_, b_] = readNavIono(nF); ip_ = find(Sana.code == string(cfg.posStation), 1);
        kl = struct('alpha', a_, 'beta', b_, 'rxLat', Sana.lat(ip_), 'rxLon', Sana.lon(ip_));
        T6a = evaluateSlantDelay(S(idx5), models(idx5), dsTe, lat, lon, obsK, 'klobuchar', kl);
        writetable(T6a, fullfile(resDir, 'table6_slant.csv'));
        mkMap = @(Yhat, ep) struct('Yhat', Yhat, 'epochs', ep, 'lat', lat, 'lon', lon);
        corr = struct('name', {'None','Klobuchar (broadcast)','Persistence-24h','CODE C1PG','SH-AR','Res-TECNet-A, operational','Res-TECNet-A, hindcast','CODE final GIM (post-processed)'}, ...
            'type', {'none','klobuchar','map','map','map','map','map','map'}, ...
            'maps', {[],[], mkMap(statsPers.Yhat, dsTe.epochsY), mkMap(statsC1PG.Yhat, dsTe.epochsY), ...
                     mkMap(statsSHARfield.Yhat, dsTe.epochsY), mkMap(statsOpA.Yhat, dsTeOpA.epochsY), mkMap(statsResA.Yhat, dsTeA.epochsY), mkMap(single(double(dsTe.Y)*cfg.Tmax), dsTe.epochsY)}, ...
            'hShell', 506.7, 'alpha', 0.9782);
        errAll = []; 
        for d = datetime(2024,1,1,'TimeZone','UTC') : days(7) : datetime(2024,12,31,'TimeZone','UTC')   % weekly sampling
            [oF, nF] = downloadRINEX(cfg.posStation, d, 'rinex_cache', 'cddisNetrc', cfg.netrc, 'country', 'TUR');
            if isempty(oF), continue; end
            [~, sol] = sppPositioning(oF, nF, corr);
            errAll = cat(1, errAll, sol.errENU);
        end
        h = squeeze(hypot(errAll(:,1,:), errAll(:,2,:))); v = squeeze(abs(errAll(:,3,:)));
        bad = any(h > 100 | v > 100, 2);                          % gross blunders (same for every correction, e.g. bad code data) are not ionospheric errors
        fprintf('Positioning: %d of %d epochs removed as gross blunders (>100 m in any solution)\n', nnz(bad), numel(bad));
        h(bad, :) = NaN; v(bad, :) = NaN;
        T6b = table({corr.name}', sqrt(mean(h.^2,'omitnan'))', prctile(h,95)', sqrt(mean(v.^2,'omitnan'))', prctile(v,95)', ...
            'VariableNames', {'Correction','Horizontal_RMS_m','Horizontal_95_m','Vertical_RMS_m','Vertical_95_m'});
        disp(T6b); writetable(T6b, fullfile(resDir, 'table6_positioning.csv'));
    end
  catch ME
    fprintf(2, '!!! Section 10 (Anatolia / slant / SPP) failed: %s\n%s\n', ME.message, getReport(ME, 'basic'));
  end
end

%% ---------------- 11. Ablation (Table 7) and interpolation influence (Sect. 4.8.1) ----------------
fprintf('[%s] %s\n', datestr(now, 'HH:MM:SS'), '11. Ablation (Table 7) and interpolation influence (Sect. 4.8.1)');
if do.ablation
    stationEval = [];
    if exist('stAna', 'var') && ~isempty(stAna)
        stationEval = @(st) anatoliaStorm(st, dsTe, lat, lon, stAna, stormMask);
    end
    drvAbl = drivers; drvAbl.KpFc = driversOp.KpFc;      % operational row needs KpFc
    splitsAbl = splits; splitsAbl.tr = ismember(year(epochs), budget.ablTrainYears);
    % every row (including 'full') is trained under the same reduced budget so that the differences are comparable
    ablBase = []; ablTag = ''; init = struct();
    if do.adaptive
        ablBase = Fsh; ablTag = '_A2';                                    % '_A' rows (from scratch, not converged) are superseded
        init.default = netA; init.persistBase = net; init.operational = netOpA;
    else
        init.default = net; init.operational = netOp;
    end
  try
    T7 = driverAblation(TEC, epochs, drvAbl, splitsAbl, lat, cfg, DstTest, 'trainOpts', trainOptsAbl, 'stationEval', stationEval, ...
        'stride', budget.ablStride, 'rows', budget.ablRows, 'base', ablBase, 'tag', ablTag, 'init', init);
    writetable(T7, fullfile(resDir, 'table7_ablation.csv'));
  catch ME
    fprintf(2, '!!! Ablation (Table 7) failed: %s\n%s\n', ME.message, getReport(ME, 'basic'));
  end
end

%% ---------------- 12. Predictive uncertainty (Sect. 4.9) ----------------
fprintf('[%s] %s\n', datestr(now, 'HH:MM:SS'), '12. Predictive uncertainty (Sect. 4.9)');
if do.uncertainty
    if exist('statsU_forCal', 'var'), statsU = statsU_forCal; else, statsU = evaluateTEC(netU, dsTe, lat); end
    Rcal = calibrationStats(statsU, dsTe, lat, stormMask);
    fprintf('Res-TECNet-U RMSE %.3f TECU\n', statsU.RMSE);
    save(fullfile(resDir, 'calibration.mat'), 'Rcal');
end
save(fullfile(resDir, 'all.mat'), 'T2', 'T3', 'cfg', 'budget', '-v7.3');

%% ---------------- 13. Heavy runs last: interpolation influence (Sect. 4.8.1) and hyperparameters (Table 1) ----------------
fprintf('[%s] %s\n', datestr(now, 'HH:MM:SS'), '13. Interpolation influence (Sect. 4.8.1) and hyperparameters (Table 1)');
if do.ablation
  try
    intFile = fullfile(resDir, 'interpolation_influence.mat');
    if isfile(intFile)
        load(intFile, 'Rint');
    else
        Rint = interpolationInfluence(TEC, epochs, lat, lon, drivers, splits, cfg, 'trainNative', true, 'scaled', true, ...
            'trainOpts', {'plots', 'none', 'maxEpochs', 40, 'patience', 8, 'dropPeriod', 12, 'background', true, 'augment', false}, ...
            'DstTest', DstTest, 'stride', 3);   % from scratch, same schedule as the main network
        save(intFile, 'Rint');
    end
    summarizeStruct(Rint, fullfile(resDir, 'interpolation_influence_summary.txt'));   % small text file (the .mat is too large to share)
  catch ME
    fprintf(2, '!!! Interpolation influence failed: %s\n%s\n', ME.message, getReport(ME, 'basic'));
  end
end
try
if do.hyper && do.adaptive && ~exist('Fsh', 'var')
    if isfile('SHARforecasts.mat'), load('SHARforecasts.mat', 'Fsh'); else, Fsh = computeSHARForecasts(TEC, epochs, lat, lon, 'Delta', cfg.Delta, 'order', 48); save('SHARforecasts.mat', 'Fsh', '-v7.3'); end
end
if do.hyper
    splitsAbl = splits; splitsAbl.tr = ismember(year(epochs), budget.ablTrainYears);
    hBase = []; hTag = '_P2'; hInit = net; if do.adaptive, hBase = Fsh; hTag = '_A2'; hInit = netA; end
    T1 = hyperparameterSearch(TEC, epochs, drivers, splitsAbl, lat, cfg, DstAll, 'trainOpts', trainOptsAbl, 'stride', budget.ablStride, ...
        'L', budget.hyper.L, 'lambda', budget.hyper.lambda, 'B', budget.hyper.B, 'filters', budget.hyper.filters, 'base', hBase, 'tag', hTag, 'init', hInit);
    writetable(T1, fullfile(resDir, 'table1_hyperparameters.csv'));
end
catch ME
  fprintf(2, '!!! Hyperparameter search failed: %s\n%s\n', ME.message, getReport(ME, 'basic'));
end

%% ---------------- local helper ----------------
fprintf('[%s] %s\n', datestr(now, 'HH:MM:SS'), 'local helper');
function r = anatoliaStorm(st, dsTe, lat, lon, stAna, stormMask)
[~, Treg] = evaluateStationVTEC({st}, {'m'}, dsTe, lat, lon, stAna, 'stormMask', stormMask);
r = Treg{'All stations (storm)', 'm'};
end

function [net, info, continued] = trainOrContinue(matFile, varName, makeNet, dsTr, dsVa, lat, cfg, trainOpts, ckDir, maxEpochs)
% load a saved network; continue training from its checkpoints when the
% requested maxEpochs exceeds the epochs already done; otherwise train new
continued = false; info = struct('epochsDone', 0);
if isfile(matFile)
    S = load(matFile); net = S.(varName);
    done = 0; if isfield(S, 'info') && isfield(S.info, 'epochsDone'), done = S.info.epochsDone;
    elseif isfile(fullfile(ckDir, 'resume_state.mat')), R = load(fullfile(ckDir, 'resume_state.mat')); done = R.epochsBase; end
    if done >= maxEpochs, info = struct('epochsDone', done); return; end
    fprintf('%s: %d epochs done, continuing to %d\n', matFile, done, maxEpochs);
    continued = true;
else
    net = makeNet();
end
[net, info] = trainResTECNet(net, dsTr, dsVa, lat, 'lambda', cfg.lambda, trainOpts{:}, 'checkpointDir', ckDir);
S = struct(varName, net, 'info', info, 'cfg', cfg); save(matFile, '-struct', 'S');
end
