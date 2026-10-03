function Thyp = hyperparameterSearch(TEC, epochs, drivers, splits, lat, cfg, DstVal, opts)
%HYPERPARAMETERSEARCH One-factor-at-a-time search on the validation set (Table 1).
%   Thyp = HYPERPARAMETERSEARCH(TEC, epochs, drivers, splits, lat, cfg, DstVal, ...)
%     Factors (defaults = grids of the revised manuscript, Sect. 3.3):
%       'L'       {6, 12, 24, 48}
%       'lambda'  {0, 0.05, 0.1, 0.2, 0.5}
%       'B'       {4, 8, 12, 16}
%       'filters' {32, 64, 128}
%     Each row changes one factor relative to the retained configuration
%     cfg (L = 12, lambda = 0.1, B = 8, 64 filters). Criteria: latitude-
%     weighted RMSE over all validation epochs and over the validation-set
%     storm windows (Dst <= -100 nT, +/- 48 h). DstVal must be aligned to
%     the validation epochs (full hourly vector, masked by splits.va).
%   Options: 'trainOpts' passed to trainResTECNet (e.g. {'maxEpochs', 60}).

arguments
    TEC (:,:,:) {mustBeNumeric}
    epochs (:,1) datetime
    drivers table
    splits struct
    lat (:,1) double
    cfg struct
    DstVal (:,1) double
    opts.L double = [6 12 24 48]        % pass [] to skip a factor
    opts.lambda double = [0 0.05 0.1 0.2 0.5]
    opts.B double = [4 8 12 16]
    opts.filters double = [32 64 128]
    opts.trainOpts cell = {}
    opts.stride (1,1) double = 1
    opts.scaled (1,1) logical = true   % TEC already divided by Tmax (runExperiments)
    opts.base = []
    opts.tag (1,:) char = ''
    opts.init = []                     % warm start from the final network (rows that keep the input layout)
end
rows = {};
grid = [num2cell(opts.L), num2cell(opts.lambda), num2cell(opts.B), num2cell(opts.filters)];
fac  = [repmat({'L'}, 1, numel(opts.L)), repmat({'lambda'}, 1, numel(opts.lambda)), ...
        repmat({'B'}, 1, numel(opts.B)), repmat({'filters'}, 1, numel(opts.filters))];
for r = 1:numel(grid)
    L = cfg.L; lam = cfg.lambda; B = cfg.B; nf = cfg.filters;
    switch fac{r}
        case 'L', L = grid{r};
        case 'lambda', lam = grid{r};
        case 'B', B = grid{r};
        case 'filters', nf = grid{r};
    end
    dsTr = buildTECDataset(TEC, epochs, drivers, 'L', L, 'Delta', cfg.Delta, 'Tmax', cfg.Tmax, ...
        'mask', splits.tr, 'scaled', opts.scaled, 'stride', opts.stride, 'targets', false, 'base', opts.base);
    dsVa = buildTECDataset(TEC, epochs, drivers, 'L', L, 'Delta', cfg.Delta, 'Tmax', cfg.Tmax, ...
        'mask', splits.va, 'scaled', opts.scaled, 'zmu', dsTr.zmu, 'zsig', dsTr.zsig, 'base', opts.base);
    H = dsTr.H; W = dsTr.W;
    rowFile = sprintf('hyper%s_%s_%g.mat', opts.tag, fac{r}, grid{r});
    if isfile(rowFile)
        R = load(rowFile); net = R.net; st = R.st;                     % finished on a previous run
    else
        bc = {}; if ~isempty(opts.base), bc = {'baseChannel', L + 17}; end
        net = createResTECNet(H, W, dsTr.C, 'L', L, 'numBlocks', B, 'numFilters', nf, bc{:});
        if ~isempty(opts.init) && L == cfg.L && nf == cfg.filters, net = transferWeights(net, opts.init); end
        net = trainResTECNet(net, dsTr, dsVa, lat, 'lambda', lam, 'plots', 'none', opts.trainOpts{:}, ...
            'checkpointDir', fullfile('checkpoints', sprintf('hyper%s_%s_%g', opts.tag, fac{r}, grid{r})));
        st = evaluateTEC(net, dsVa, lat); st.Yhat = [];
        save(rowFile, 'net', 'st', '-v7.3');
    end
    dv = DstVal(dsVa.idxT + dsVa.Delta);
    Ts = stormAnalysis({st}, {'net'}, dsVa.epochsY, dv, -100, 48);
    rows(end+1, :) = {fac{r}, grid{r}, st.RMSE, Ts.StormRMSE_TECU(1), countParams(net)}; %#ok<AGROW>
    fprintf('%s = %g : val RMSE %.3f, storm %.3f\n', fac{r}, grid{r}, st.RMSE, Ts.StormRMSE_TECU(1));
end
Thyp = cell2table(rows, 'VariableNames', {'Factor', 'Value', 'ValRMSE_TECU', 'ValStormRMSE_TECU', 'NumParams'});
disp(Thyp);
end

function n = countParams(net)
n = 0;
for i = 1:height(net.Learnables), n = n + numel(net.Learnables.Value{i}); end
end
