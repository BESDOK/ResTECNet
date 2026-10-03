function ds = buildTECDataset(TEC, epochs, drivers, opts)
%BUILDTECDATASET Sample index set for Res-TECNet (Eqs. 1-3, Sect. 3.1/3.4).
%   ds = BUILDTECDATASET(TEC, epochs, drivers, ...)
%
%   Memory-lean design: the input tensor X (H x W x (L+16) per sample) is
%   NOT materialised. The dataset stores the scaled maps once (ds.Ts,
%   single), the per-sample driver rows (ds.Zt, ds.Ztd, N x 8) and the
%   sample indices; X is assembled on demand by TECGETBATCH(ds, idx),
%   which trainResTECNet / trainModelFunction / evaluateTEC call per
%   mini-batch. This keeps a 20-year training set at ~4 GB instead of
%   ~100 GB.
%
%   TEC     : H x W x K hourly maps [TECU] (W = 72 unique meridians; NaN
%             marks unfilled gaps from fillTECGaps)
%   epochs  : K x 1 datetime (UTC), hourly, regular
%   drivers : table from buildDriverTable (F107, F107bar, Kp, Dst [, KpFc])
%   Options
%     'L', 'Delta', 'Tmax'      (12, 24, 200)
%     'zmu','zsig'  1x4 training statistics of [F107 F107bar Kp Dst]
%                   (computed from THIS data when empty)
%     'mode'   'hindcast'    : target block z_{t+D} = observed indices at t+D
%              'operational' : F107(t), F107bar(t), drivers.KpFc(t+D), Dst(t)
%                              + harmonic time of t+D (Sect. 3.4)
%     'zeroChannels'  cellstr among {'F107','Kp','Dst','time','target','all'}
%                   (per-driver ablation, Sect. 4.8) -- zeroes the columns
%                   of ds.Zt / ds.Ztd
%     'badMask' K x 1 logical (unfilled gaps)
%     'stride'  keep every stride-th current epoch t (default 1 = hourly);
%               used to thin the TRAINING set (consecutive hourly samples
%               are highly redundant) -- validation/test stay at 1 h
%     'mask'    K x 1 logical selecting the split (e.g. training years) inside
%               the FULL record; samples whose window [t-L+1, t+D] leaves the
%               mask are dropped. Pass the full TEC array (not TEC(:,:,mask))
%               so that no copy is made.
%     'scaled'  true if TEC is already single and divided by Tmax; then
%               ds.Ts references it without copying (memory-lean full run)
%     'base'    H x W x K single, SAME scaling as TEC: the adaptive base map
%               (SH-AR forecast for epoch k, from computeSHARForecasts). When
%               given, tecGetBatch appends it as channel C = L+17 (base for
%               epoch t+Delta) and Res-TECNet-A uses it for the global skip
%     'targets' materialise ds.Y / ds.Tlast (default true; set false for
%               training sets -- tecGetBatch reads them from ds.Ts)
%
%   Output ds (struct):
%     .Ts      H x W x K single   scaled maps (shared by all samples)
%     .idxT    N x 1              index of the current epoch t in Ts
%     .Zt, .Ztd N x 8 single      driver blocks z_t, z_{t+D}
%                                 order [F107 F107bar Kp Dst sinUT cosUT sinDOY cosDOY]
%     .Y       H x W x 1 x N single  targets (scaled)   -- materialised
%     .Tlast   H x W x 1 x N single  T_t (scaled)       -- materialised
%     .epochsY N x 1 datetime, .N, .H, .W, .C (= L+16), .L, .Delta,
%     .zmu, .zsig, .Tmax, .mode
%   X for samples idx: [X, Y] = tecGetBatch(ds, idx)  (H x W x C x numel(idx))

arguments
    TEC (:,:,:) {mustBeNumeric}
    epochs (:,1) datetime
    drivers table
    opts.L (1,1) double = 12
    opts.Delta (1,1) double = 24
    opts.Tmax (1,1) double = 200
    opts.zmu double = []
    opts.zsig double = []
    opts.mode (1,:) char {mustBeMember(opts.mode, {'hindcast','operational'})} = 'hindcast'
    opts.zeroChannels cell = {}
    opts.badMask logical = false(0,1)
    opts.stride (1,1) double = 1
    opts.mask logical = false(0,1)
    opts.scaled (1,1) logical = false
    opts.targets (1,1) logical = true
    opts.base = []
end

[H, W, K] = size(TEC);
L = opts.L; D = opts.Delta;
if height(drivers) ~= K
    error('buildTECDataset:align', 'drivers must have one row per epoch.');
end

Z = [drivers.F107, drivers.F107bar, drivers.Kp, drivers.Dst];
if isempty(opts.zmu)
    opts.zmu = mean(Z, 1); opts.zsig = std(Z, 0, 1);
end
zeroVar = opts.zsig == 0 | ~isfinite(opts.zsig);
if any(zeroVar)
    warning('buildTECDataset:zeroVariance', ...
        'Driver column(s) %s have zero variance; scale set to 1.', mat2str(find(zeroVar)));
    opts.zsig(zeroVar) = 1;
end
Zn = (Z - opts.zmu) ./ opts.zsig;
ut  = hour(epochs) + minute(epochs)/60;
doy = day(epochs, 'dayofyear');
Henc = [sin(2*pi*ut/24), cos(2*pi*ut/24), sin(2*pi*doy/365.25), cos(2*pi*doy/365.25)];
Zcur = [Zn, Henc];

% ---- valid samples ----
if ~isempty(opts.badMask), mapBad = opts.badMask(:);
else, mapBad = squeeze(isnan(TEC(1, 1, :))) | squeeze(isnan(TEC(end, end, :)));   % cheap NaN test (whole-map NaN), K x 1
end
mapBad = mapBad(:);
if ~isempty(opts.mask), mapBad = mapBad | ~opts.mask(:); end
baseBad = false(K, 1);
if ~isempty(opts.base), baseBad = squeeze(isnan(opts.base(1, 1, :))); baseBad = baseBad(:); end
idx0 = (L : opts.stride : (K - D))';
cumBad = cumsum([0; double(mapBad(:))]);
histBad = cumBad(idx0 + 1) - cumBad(idx0 - L + 1);
ok = histBad == 0 & ~mapBad(idx0 + D) & ~baseBad(idx0 + D);
idx0 = idx0(ok);
N = numel(idx0);

% ---- per-sample driver rows ----
Zt = Zcur(idx0, :);
switch opts.mode
    case 'hindcast'
        Ztd = Zcur(idx0 + D, :);
    case 'operational'
        if ~ismember('KpFc', drivers.Properties.VariableNames)
            error('buildTECDataset:kpfc', ...
                'Operational mode needs drivers.KpFc (buildDriverTable with ''mode'',''operational'').');
        end
        KpFcN = (drivers.KpFc - opts.zmu(3)) ./ opts.zsig(3);
        Ztd = [Zn(idx0, 1:2), KpFcN(idx0 + D), Zn(idx0, 4), Henc(idx0 + D, :)];
end
if ~isempty(opts.zeroChannels)
    ch = [];
    for c = 1:numel(opts.zeroChannels)
        switch lower(opts.zeroChannels{c})
            case 'f107',   ch = [ch, 1 2]; %#ok<AGROW>
            case 'kp',     ch = [ch, 3];   %#ok<AGROW>
            case 'dst',    ch = [ch, 4];   %#ok<AGROW>
            case 'time',   ch = [ch, 5:8]; %#ok<AGROW>
            case 'all',    ch = [ch, 1:8]; %#ok<AGROW>
            case 'target', Ztd(:) = 0;
            otherwise, error('buildTECDataset:zero', 'Unknown channel group %s', opts.zeroChannels{c});
        end
    end
    ch = unique(ch); Zt(:, ch) = 0; Ztd(:, ch) = 0;
end

if opts.scaled, ds.Ts = TEC; else, ds.Ts = single(TEC) / single(opts.Tmax); end
ds.idxT = idx0;
ds.Zt = single(Zt); ds.Ztd = single(Ztd);
if opts.targets
    ds.Y = reshape(ds.Ts(:, :, idx0 + D), H, W, 1, N);
    ds.Tlast = reshape(ds.Ts(:, :, idx0), H, W, 1, N);
else
    ds.Y = []; ds.Tlast = [];
end
ds.epochsY = epochs(idx0 + D);
ds.Base = opts.base;                              % [] or shared reference
ds.N = N; ds.H = H; ds.W = W; ds.C = L + 16 + double(~isempty(opts.base)); ds.L = L; ds.Delta = D;
ds.zmu = opts.zmu; ds.zsig = opts.zsig; ds.Tmax = opts.Tmax; ds.mode = opts.mode;
ds.zeroChannels = opts.zeroChannels;
fprintf('buildTECDataset (%s): %d samples, %d x %d grid, C = %d, %d dropped for gaps/mask (%.2f GB own)\n', ...
    opts.mode, N, H, W, ds.C, nnz(~ok), (numel(ds.Ts) * double(~opts.scaled) + 2*numel(ds.Y)) * 4 / 1e9);
end
