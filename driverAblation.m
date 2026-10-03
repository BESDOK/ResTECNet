function [Tabl, statsAbl] = driverAblation(TEC, epochs, drivers, splits, lat, cfg, DstTest, opts)
%DRIVERABLATION Architectural and per-driver ablation (Sect. 4.8, Table 7).
%   [Tabl, statsAbl] = DRIVERABLATION(TEC, epochs, drivers, splits, lat, cfg, DstTest, ...)
%     splits  : struct with logical masks .tr, .va, .te over epochs
%     cfg     : struct with L, Delta, Tmax, lambda, zmu, zsig (training
%               statistics), H, W, C
%     DstTest : Dst aligned with the test targets (for storm windows)
%   Options
%     'rows'  cellstr subset of the configurations below (default all)
%     'trainOpts' cell of name/value pairs passed to trainResTECNet
%     'stride'  training-sample stride (reduced budget)
%     'init'    warm-start networks: every variant is initialised from the
%               final network (layers with matching names/sizes; zeroed or
%               removed parts start from their initialisation) and fine-tuned
%               under the reduced budget, so that all rows are converged and
%               comparable (a from-scratch short training is not)
%     'reuse'   struct row -> .mat file of an already trained network
%               (e.g. plainCNN -> baseline_plain.mat); evaluated, not retrained
%     'stationEval' function handle @(stats) -> RMSE at the Anatolian
%                   stations in the storm windows (optional last column)
%   Configurations (each retrained from scratch, hindcast inputs):
%     'full', 'persistBase' (global skip from T_t instead of the SH-AR base,
%     no base channel: the original Res-TECNet), 'noGlobalSkip', 'plainCNN', 'zeroPad', 'noF107', 'noKp',
%     'noDst', 'noTime', 'noDrivers', 'noTargetBlock', 'noisyTarget'
%     (inference-only), 'operational' (needs drivers.KpFc)

arguments
    TEC (:,:,:) {mustBeNumeric}
    epochs (:,1) datetime
    drivers table
    splits struct
    lat (:,1) double
    cfg struct
    DstTest (:,1) double
    opts.rows cell = {}
    opts.trainOpts cell = {}
    opts.stationEval = []
    opts.stride (1,1) double = 1
    opts.scaled (1,1) logical = true   % TEC already divided by Tmax (runExperiments)
    opts.reuse struct = struct()
    opts.base = []                     % SH-AR base field (Res-TECNet-A); [] = persistence base
    opts.tag (1,:) char = ''           % file-name tag (e.g. '_A') so that rows of different model variants do not collide
    opts.init struct = struct()        % warm start: .default (final net), .persistBase (static net), .operational (op net)
end
allRows = {'full','persistBase','noGlobalSkip','plainCNN','zeroPad','noF107','noKp','noDst', ...
           'noTime','noDrivers','noTargetBlock','noisyTarget','operational'};
if isempty(opts.rows), opts.rows = allRows; end
L = cfg.L; D = cfg.Delta;
mk = @(mask, varargin) buildTECDataset(TEC, epochs, drivers, 'L', L, 'Delta', D, 'Tmax', cfg.Tmax, ...
    'zmu', cfg.zmu, 'zsig', cfg.zsig, 'mask', mask, 'scaled', opts.scaled, 'base', opts.base, varargin{:});
mkP = @(mask, varargin) buildTECDataset(TEC, epochs, drivers, 'L', L, 'Delta', D, 'Tmax', cfg.Tmax, ...
    'zmu', cfg.zmu, 'zsig', cfg.zsig, 'mask', mask, 'scaled', opts.scaled, varargin{:});   % no base (persistence)
useBase = ~isempty(opts.base);
netOptsBase = {}; if useBase, netOptsBase = {'baseChannel', cfg.L + 17}; end

statsAbl = struct(); rows = {};
netFull = [];
for r = 1:numel(opts.rows)
    name = opts.rows{r};
    netOpts = [{'L', L}, netOptsBase]; dsOpts = {}; mkRow = mk;
    switch name
        case 'full'
        case 'persistBase',   netOpts = {'L', L}; mkRow = mkP;
        case 'noGlobalSkip',  netOpts = [netOpts, {'globalSkip', false}];
        case 'plainCNN',      netOpts = [netOpts, {'residual', false, 'globalSkip', false}];
        case 'zeroPad',       netOpts = [netOpts, {'circular', false}];
        case 'noF107',        dsOpts = {'zeroChannels', {'F107'}};
        case 'noKp',          dsOpts = {'zeroChannels', {'Kp'}};
        case 'noDst',         dsOpts = {'zeroChannels', {'Dst'}};
        case 'noTime',        dsOpts = {'zeroChannels', {'time'}};
        case 'noDrivers',     dsOpts = {'zeroChannels', {'all'}};
        case 'noTargetBlock', dsOpts = {'zeroChannels', {'target'}};
        case 'operational',   dsOpts = {'mode', 'operational'};
        case 'noisyTarget'
            if isempty(netFull), error('driverAblation:order', 'Run ''full'' before ''noisyTarget''.'); end
            dsTe = mk(splits.te);
            n = dsTe.N;
            dsTe.Ztd(:, 1) = dsTe.Ztd(:, 1) + single(10 / cfg.zsig(1)) * randn(n, 1, 'single');  % F10.7(t+D) +/-10 sfu
            dsTe.Ztd(:, 3) = dsTe.Ztd(:, 3) + single(1  / cfg.zsig(3)) * randn(n, 1, 'single');  % Kp(t+D) +/-1
            st = evaluateTEC(netFull, dsTe, lat);
            [rows, statsAbl] = addRow(rows, statsAbl, name, st, dsTe, DstTest, opts.stationEval);
            continue
    end
    rowFile = sprintf('ablation%s_%s.mat', opts.tag, name);
    dsTe = mkRow(splits.te, dsOpts{:});
    if isfile(rowFile)
        R = load(rowFile); net = R.net; if isfield(R, 'st'), st = R.st; else, st = R.st0; end   % finished on a previous run
        if ~isempty(opts.stationEval) && (~isfield(st, 'Yhat') || isempty(st.Yhat)), st = evaluateTEC(net, dsTe, lat); end   % cached rows are stored without the 0.7 GB map array
    elseif isfield(opts.reuse, name) && isfile(opts.reuse.(name))
        R = load(opts.reuse.(name)); f = fieldnames(R); net = R.(f{find(cellfun(@(x) isa(R.(x), 'dlnetwork'), f), 1)});
        st = evaluateTEC(net, dsTe, lat);                              % reuse the main-run network
    else
        dsTr = mkRow(splits.tr, 'stride', opts.stride, 'targets', false, dsOpts{:}); dsVa = mkRow(splits.va, dsOpts{:});
        net = createResTECNet(cfg.H, cfg.W, dsTr.C, netOpts{:});
        initNet = [];
        if isfield(opts.init, name), initNet = opts.init.(name);
        elseif isfield(opts.init, 'default'), initNet = opts.init.default; end
        if ~isempty(initNet), net = transferWeights(net, initNet); end
        net = trainResTECNet(net, dsTr, dsVa, lat, 'lambda', cfg.lambda, opts.trainOpts{:}, ...
            'checkpointDir', fullfile('checkpoints', ['abl' opts.tag '_' name]));
        st = evaluateTEC(net, dsTe, lat);
    end
    if strcmp(name, 'full'), netFull = net; end
    [rows, statsAbl] = addRow(rows, statsAbl, name, st, dsTe, DstTest, opts.stationEval);
    if ~isfile(rowFile), stSave = st; stSave.Yhat = []; st0 = stSave; save(rowFile, 'net', 'st0', '-v7.3'); end
end
Tabl = cell2table(rows, 'VariableNames', {'Configuration', 'Overall_TECU', 'Storm_TECU', 'AnatoliaStorm_TECU'});
disp(Tabl);
end

function [rows, S] = addRow(rows, S, name, st, dsTe, DstTest, stationEval)
Ts = stormAnalysis({st}, {name}, dsTe.epochsY, DstTest, -100, 48);
anat = NaN; if ~isempty(stationEval), anat = stationEval(st); end
rows(end+1, :) = {name, st.RMSE, Ts.StormRMSE_TECU(1), anat};
st.Yhat = [];                                   % do not accumulate 0.7 GB per row in memory
S.(matlab.lang.makeValidName(name)) = st;
end
