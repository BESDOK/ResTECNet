%RUNSMOKETEST Reduced end-to-end rehearsal of runExperiments (step 1).
%   Downloads the CODE final GIMs for 2017-2021 (short-name .Z files for
%   2017-2022 -> checks the 7-Zip/uncompress path), the C1PG predicted
%   GIMs and UQRG maps for 2021, builds hindcast + operational drivers,
%   trains Res-TECNet for a few epochs on 2017-2018 (val 2019, test 2021),
%   evaluates persistence / C1PG / SH-AR, produces Table 2/3 style output,
%   fig4.jpg (forecast maps) and fig0.jpg (error map), and re-runs the
%   May-2024 figures with the (weak) network overlaid.
%   Expected wall time: downloads 20-40 min, training 15-60 min (GPU).
%   Results in results_smoke/. Nothing here is publication-grade: it only
%   proves that every stage runs on real data before the full run.

cfg.L = 12; cfg.Delta = 24; cfg.Tmax = 200; cfg.lambda = 0.1; cfg.B = 8; cfg.filters = 64;
trainYears = 2017:2018; valYears = 2019; testYears = 2021;
maxEpochs = 3;                                   % rehearsal only
trainStride = 3;                                 % same thinning as the full run (budget.trainStride)
resDir = 'results_smoke'; if ~isfolder(resDir), mkdir(resDir); end
t0 = tic;

%% 1. databases
if ~isfile('TEC_smoke.mat')
    buildTECDatabase('2017-01-01', '2021-12-31', 'TEC_smoke.mat', 'ionex_cache', 'product', 'CODG');
end
if ~isfile('C1PG_smoke.mat')
    buildTECDatabase('2021-01-01', '2021-12-31', 'C1PG_smoke.mat', 'ionex_cache', 'product', 'C1PG', 'fillGaps', false);
end
if ~isfile('UQRG_smoke.mat')
    try
        buildTECDatabase('2021-01-01', '2021-12-31', 'UQRG_smoke.mat', 'ionex_cache', 'product', 'UQRG', 'fillGaps', false);
    catch ME, warning('UQRG skipped: %s', ME.message);
    end
end
load('TEC_smoke.mat', 'TEC', 'epochs', 'lat', 'lon', 'badMask');
TEC = TEC / single(cfg.Tmax);                   % scaled once; datasets reference it
[H, W, K] = size(TEC); C = cfg.L + 16; cfg.H = H; cfg.W = W; cfg.C = C;
fprintf('Database: %d maps %s..%s, %d unfilled gaps (%.1f min)\n', K, datestr(epochs(1)), datestr(epochs(end)), nnz(badMask), toc(t0)/60);

%% 2. drivers
if ~isfile('drivers_smoke.mat')
    drivers = buildDriverTable(epochs, 'index_cache', 'mode', 'hindcast');
    driversOp = buildDriverTable(epochs, 'index_cache', 'mode', 'operational');
    save('drivers_smoke.mat', 'drivers', 'driversOp');
else
    load('drivers_smoke.mat', 'drivers', 'driversOp');
end
fprintf('Kp forecast source: %d archived / %d error-model epochs\n', nnz(driversOp.KpFcSource == 1), nnz(driversOp.KpFcSource == 0));

%% 3. datasets
splits.tr = ismember(year(epochs), trainYears);
splits.va = ismember(year(epochs), valYears);
splits.te = ismember(year(epochs), testYears);
mk = @(mask, drv, varargin) buildTECDataset(TEC, epochs, drv, 'L', cfg.L, 'Delta', cfg.Delta, 'Tmax', cfg.Tmax, ...
    'badMask', badMask, 'mask', mask, 'scaled', true, varargin{:});
dsTr = mk(splits.tr, drivers, 'stride', trainStride, 'targets', false); cfg.zmu = dsTr.zmu; cfg.zsig = dsTr.zsig;
dsVa = mk(splits.va, drivers, 'zmu', cfg.zmu, 'zsig', cfg.zsig);
dsTe = mk(splits.te, drivers, 'zmu', cfg.zmu, 'zsig', cfg.zsig);
dsTeOp = mk(splits.te, driversOp, 'zmu', cfg.zmu, 'zsig', cfg.zsig, 'mode', 'operational');
DstTest = drivers.Dst(dsTe.idxT + dsTe.Delta);
fprintf('Datasets: train %d, val %d, test %d samples (%.1f GB scaled maps)\n', ...
    dsTr.N, dsVa.N, dsTe.N, numel(dsTr.Ts)*4/1e9);

%% 4. short training
net = createResTECNet(H, W, C, 'L', cfg.L);
[net, info] = trainResTECNet(net, dsTr, dsVa, lat, 'lambda', cfg.lambda, 'maxEpochs', maxEpochs, 'plots', 'none');
% -> note the seconds per epoch here; the full run has ~10x more training samples
save(fullfile(resDir, 'net_smoke.mat'), 'net', 'info', 'cfg');
fprintf('Training done (%.1f min)\n', toc(t0)/60);
statsRes = evaluateTEC(net, dsTe, lat);
statsOp  = evaluateTEC(net, dsTeOp, lat);          % same net, operational inputs (rehearsal only)

%% 5. baselines and tables
statsPers = persistenceBaseline(dsTe, lat);
statsC1PG = baselineC1PG('C1PG_smoke.mat', dsTe, lat);
leadIn = splits.te | ismember(year(epochs), valYears(end));
statsSHAR = baselineSHAR(TEC(:,:,leadIn) * cfg.Tmax, epochs(leadIn), dsTe, lat, lon, 'order', 48);
models = {'Persistence-24h','CODE C1PG','SH-AR','Res-TECNet (3 epochs, op. inputs)','Res-TECNet (3 epochs)'};
S = {statsPers, statsC1PG, statsSHAR, statsOp, statsRes};
T2 = table(models(:), cellfun(@(s) s.RMSE, S)', cellfun(@(s) s.MAE, S)', cellfun(@(s) s.Bias, S)', cellfun(@(s) s.Corr, S)', ...
    'VariableNames', {'Model','RMSE','MAE','Bias','Corr'});
disp(T2); writetable(T2, fullfile(resDir, 'table2_smoke.csv'));
[T3, stormMask] = stormAnalysis(S, models, dsTe.epochsY, DstTest, -100, 48);
disp(T3); writetable(T3, fullfile(resDir, 'table3_smoke.csv'));

%% 6. figures
ex = [datetime(2021,3,13,18,0,0,'TimeZone','UTC'); datetime(2021,11,4,18,0,0,'TimeZone','UTC')]; % quiet / Nov-2021 storm
plotForecastMaps(statsRes, dsTe, lat, lon, ex, 'outFile', fullfile(resDir, 'fig4_smoke.jpg'), ...
    'labels', {'Quiet (13 Mar 2021, 18:00 UT)', 'Storm (4 Nov 2021, 18:00 UT)'});
plotErrorMap(statsRes, dsTe, lat, lon, 'stations', stationList('all'), 'outFile', fullfile(resDir, 'fig0_smoke.jpg'));

%% 7. UQRG global check (if available)
if isfile('UQRG_smoke.mat')
    U = load('UQRG_smoke.mat', 'TEC', 'epochs');
    [tf, loc] = ismember(dateshift(dsTe.epochsY, 'start', 'hour'), U.epochs);
    Yu = nan([H W 1 numel(tf)]); Yu(:,:,1,tf) = U.TEC(:,:,loc(tf));
    rmseOf = @(st) st.RMSE;
    fprintf('vs UQRG: Res-TECNet %.2f | persistence %.2f | CODE-UQRG %.2f TECU (%d/%d epochs)\n', ...
        rmseOf(evaluateTEC(net, dsTe, lat, 'RefY', Yu)), ...
        rmseOf(evaluateTEC([], setfield(dsTe, 'Yhat', dsTe.Tlast), lat, 'RefY', Yu)), ...
        rmseOf(evaluateTEC([], setfield(dsTe, 'Yhat', dsTe.Y), lat, 'RefY', Yu)), nnz(tf), numel(tf)); %#ok<SFLD>
end

%% 8. May 2024 figures with the rehearsal network
copyfile(fullfile(resDir, 'net_smoke.mat'), 'resTECNet_hindcast.mat');   % runFullPaperSims looks for this name
runFullPaperSims
disp(out.dailyDaytime);
fprintf('Smoke test finished in %.1f min. Results in %s/\n', toc(t0)/60, resDir);
