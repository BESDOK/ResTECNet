function report = verifyRevisionCoverage()
%VERIFYREVISIONCOVERAGE Check that the code base implements every revision item.
%   report = VERIFYREVISIONCOVERAGE()
%   Part A maps each reviewer/editor requirement of AISR-D-26-01824 to the
%   implementing function(s) and checks that they exist on the path.
%   Part B runs synthetic, download-free smoke tests of the pipeline
%   (small grids, few epochs): IONEX reader incl. duplicate-column drop,
%   gap filling, hindcast/operational dataset assembly, channel zeroing,
%   network construction incl. uncertainty head, both losses, ConvLSTM /
%   Transformer forward passes, SH-AR, persistence, evaluation, storm
%   masks, calibration, station comparison, gradients, Klobuchar, geodesy.
%   Prints PASS/FAIL per item and returns a table. Run in MATLAB R2024a
%   from the code folder: >> report = verifyRevisionCoverage
%   Items that need toolboxes not installed are reported as SKIP.

rows = {};
% ---------------------------------------------------------------- Part A
req = {
 'R1.2 / A2.3'   'Literature 2024-2026 + metric comparison (manuscript)'      {}
 'R2.M1 / A2.4'  'Independent GNSS station VTEC, 4 regimes + other GIM (UQRG)' {'stationList','downloadRINEX','processStationVTEC','evaluateStationVTEC','buildTECDatabase'}
 'R2.M2 / A2.11' 'Operational (latency-constrained) drivers vs hindcast'        {'buildDriverTable','buildTECDataset'}
 'R2.M3 / A2.19' 'Baselines: persistence, C1PG, SH-AR, ConvLSTM, Transformer; adaptive base (Res-TECNet-A)'   {'persistenceBaseline','baselineC1PG','baselineSHAR','computeSHARForecasts','modelConvLSTM','initConvLSTM','modelTransformer','initTransformer','trainModelFunction'}
 'R2.M4 / A2.28' 'Full dual-frequency CORS processing, multiple stations, gradients' {'processStationVTEC','readCODEDCB','gpsSatPos','geoTools','regionalGradients'}
 'R2.M5 / A2.32' 'Slant-domain error + positioning experiment'                 {'evaluateSlantDelay','sppPositioning','klobucharModel'}
 'R2.M6 / A2.36' 'Predicted uncertainty map (heteroscedastic head) + calibration' {'createResTECNet','lossTECNLL','calibrationStats'}
 'R2.m1 / A2.9'  '73rd longitude column handled (W = 72)'                      {'readIONEXTEC'}
 'R2.m2 / A2.8'  '2 h -> 1 h interpolation influence quantified'               {'interpolationInfluence','fillTECGaps'}
 'R2.m3 / A2.17' 'Hyperparameter selection documented'                         {'hyperparameterSearch'}
 'R2.m4 / A2.38' 'Per-driver ablation + zero padding + target block'           {'driverAblation','buildTECDataset'}
 'A1.1 / A2.26'  'Res-TECNet forecast added to the Kayseri May-2024 figure'    {'visualizeMay2024'}
 'A1.2 / A2.10'  '81-day F10.7 mean trailing (no future flux)'                 {'buildDriverTable'}
 'A2.41 / A2.42' 'Colour-bar labels; new forecast-map figure; stations on error map' {'plotForecastMaps','plotErrorMap','visualizeMay2024'}
 'A2.45 / E-code' 'Public code: download chain + all experiments in one driver' {'runExperiments','runFullPaperSims','buildTECDatabase'}
};
for i = 1:size(req, 1)
    fns = req{i, 3}; missing = fns(cellfun(@(f) exist(f, 'file') ~= 2, fns));
    ok = isempty(missing);
    rows(end+1, :) = {req{i,1}, req{i,2}, tern(ok, 'PASS', ['FAIL: missing ' strjoin(missing, ', ')])}; %#ok<AGROW>
end

% ---------------------------------------------------------------- Part B
H = 71; W = 72; K = 400; L = 12; D = 24; Tmax = 200;
lat = (87.5:-2.5:-87.5)'; lon = (-180:5:175)';
rng(7);
[LON, LAT] = meshgrid(lon, lat);
ep = (datetime(2024,5,1,'TimeZone','UTC') + hours(0:K-1))';
TEC = zeros(H, W, K);
for k = 1:K
    ut = hour(ep(k));
    TEC(:,:,k) = 20 + 30 * cosd(LAT).^4 .* max(0, cosd(LON + 15*ut - 180)) + 2*randn(H, W);
end
TEC = max(TEC, 0.5);
F107s = 150 + 5*randn(K,1);
drv = table(F107s, movmean(F107s, [80 0]), 2 + abs(randn(K,1)), -20 + 10*randn(K,1), ...
    'VariableNames', {'F107','F107bar','Kp','Dst'});
drv.Dst(200:230) = -150; drv.KpFc = min(max(drv.Kp + 0.9*randn(K,1), 0), 9);

hasDL = license('test', 'Neural_Network_Toolbox') && exist('dlnetwork', 'class') == 8;

% B1 IONEX reader with duplicate column
rows(end+1, :) = runTest('B1 readIONEXTEC drops duplicate +180 column', @() testIonex(lat, lon));
% B2 gap filling
rows(end+1, :) = runTest('B2 fillTECGaps hourly grid / short gap fill / long gap flagged', @() testGaps(TEC, ep, lon));
% B3 datasets
rows(end+1, :) = runTest('B3 buildTECDataset + tecGetBatch hindcast/operational/zeroChannels', @() testDataset(TEC, ep, drv, L, D, Tmax));
% B4 network + losses
if hasDL
    rows(end+1, :) = runTest('B4 createResTECNet (+uncertainty, zero-pad) forward pass', @() testNet(H, W, L, TEC, ep, drv, D, Tmax));
    rows(end+1, :) = runTest('B5 lossTEC / lossTECNLL finite on dlarrays', @() testLoss(H, W, lat));
    rows(end+1, :) = runTest('B6 ConvLSTM and Transformer forward shapes', @() testModelFcns(H, W, L, D, TEC, ep, drv, Tmax));
else
    rows(end+1, :) = {'B4-B6 network tests', 'need Deep Learning Toolbox', 'SKIP'};
end
% B7 SH-AR + persistence + evaluateTEC + storm masks
rows(end+1, :) = runTest('B7 baselineSHAR / persistence / evaluateTEC / stormAnalysis', @() testBaselines(TEC, ep, drv, lat, lon, L, D, Tmax));
% B8 calibration
rows(end+1, :) = runTest('B8 calibrationStats coverage in [0,1]', @() testCalib(TEC, ep, drv, lat, L, D, Tmax));
% B9 station comparison + gradients
rows(end+1, :) = runTest('B9 evaluateStationVTEC / regionalGradients with synthetic IPPs', @() testStations(TEC, ep, drv, lat, lon, L, D, Tmax));
% B10 geodesy + Klobuchar
rows(end+1, :) = runTest('B10 geoTools ecef2lla/azel/ipp/mslm + klobucharModel', @() testGeo());
% B11 driver-table parsers (Kp breakdown)
rows(end+1, :) = runTest('B11 NOAA Kp-breakdown parser (synthetic file)', @() testKpParser());
% B15/B16 station-chain building blocks
rows(end+1, :) = runTest('B15 fast RINEX-3 reader (synthetic file: columns, blanks, events, GNSS filter)', @() testRinexReader());
rows(end+1, :) = runTest('B16 gpsSatPos circular orbit + Sagnac direction', @() testOrbitSagnac());
rows(end+1, :) = runTest('B17 harmonizeReceiverDCB removes day-to-day receiver-DCB noise', @() testHarmonize());
% B13/B14 numerical cross-checks against brute-force computations
rows(end+1, :) = runTest('B13 evaluateTEC equals an independent brute-force computation', @() testEvalBrute(TEC, ep, drv, lat, L, D, Tmax));
rows(end+1, :) = runTest('B14 evaluateStationVTEC: offset removal, NaN models, hourly binning', @() testStationOffsets(TEC, ep, drv, lat, lon, L, D, Tmax));
% B12 adaptive base: SH-AR forecast field + base channel + Res-TECNet-A skip
if hasDL
    rows(end+1, :) = runTest('B12 computeSHARForecasts + base channel + Res-TECNet-A forward', @() testAdaptive(TEC, ep, drv, lat, lon, L, D, Tmax));
else
    rows(end+1, :) = {'B12 adaptive base', 'need Deep Learning Toolbox', 'SKIP'};
end

report = cell2table(rows, 'VariableNames', {'Item', 'Description', 'Result'});
disp(report);
nF = nnz(startsWith(report.Result, 'FAIL'));
fprintf('\n%d items, %d FAIL, %d SKIP\n', height(report), nF, nnz(strcmp(report.Result, 'SKIP')));
end

% ======================================================================
function row = runTest(name, f)
try
    msg = f(); if isempty(msg), msg = 'PASS'; end
    row = {name, '', msg};
catch ME
    row = {name, '', ['FAIL: ' ME.message]};
end
end
function s = tern(c, a, b), if c, s = a; else, s = b; end, end

function msg = testIonex(lat, lon)
f = [tempname '.INX']; fid = fopen(f, 'w');
hdr = @(txt, lab) fprintf(fid, '%-60s%s\n', txt, lab);
hdr('     1.0            IONOSPHERE MAPS     GPS', 'IONEX VERSION / TYPE');
hdr('    -1', 'EXPONENT');
hdr('    87.5 -87.5  -2.5', 'LAT1 / LAT2 / DLAT');
hdr('  -180.0 180.0   5.0', 'LON1 / LON2 / DLON');
hdr('', 'END OF HEADER');
hdr('     1', 'START OF TEC MAP');
hdr('  2024     5    10     0     0     0', 'EPOCH OF CURRENT MAP');
for i = 1:numel(lat)
    hdr(sprintf('  %6.1f-180.0 180.0   5.0 450.0', lat(i)), 'LAT/LON1/LON2/DLON/H');
    vals = [round(100*(1:72)/72) round(100/72)];        % 73 values, last == first
    for j = 1:16:73, fprintf(fid, '%s\n', sprintf('%5d', vals(j:min(j+15,73)))); end
end
hdr('     1', 'END OF TEC MAP');
fclose(fid);
[T, e, la, lo] = readIONEXTEC(f);
[T73, ~, ~, lo73] = readIONEXTEC(f, 'dropDuplicateLon', false);
delete(f);
assert(size(T, 2) == 72 && numel(lo) == 72 && lo(end) == 175, 'W must be 72');
assert(size(T73, 2) == 73 && lo73(end) == 180, 'full read must give 73');
assert(abs(T(1,1) - 0.1) < 1e-9 && numel(e) == 1 && numel(la) == 71, 'values/epoch');
msg = '';
end

function msg = testGaps(TEC, ep, lon)
keep = true(size(ep)); keep(50:52) = false; keep(100:110) = false;   % 3-h and 11-h gaps
[Th, eh, bad] = fillTECGaps(TEC(:,:,keep), ep(keep), lon);
assert(numel(eh) == numel(ep), 'hourly grid restored');
assert(all(~isnan(Th(:,:,50:52)), 'all') && ~any(bad(50:52)), 'short gap filled');
assert(all(bad(100:110)) && all(isnan(Th(:,:,105)), 'all'), 'long gap flagged');
msg = '';
end

function msg = testDataset(TEC, ep, drv, L, D, Tmax)
ds = buildTECDataset(TEC, ep, drv, 'L', L, 'Delta', D, 'Tmax', Tmax);
[Xh, Yh] = tecGetBatch(ds, 1:3);
assert(size(Xh, 3) == L + 16 && ds.N == numel(ep) - L - D + 1 && isequal(size(Yh), [size(Xh,1) size(Xh,2) 1 3]), 'hindcast size');
dso = buildTECDataset(TEC, ep, drv, 'L', L, 'Delta', D, 'Tmax', Tmax, 'zmu', ds.zmu, 'zsig', ds.zsig, 'mode', 'operational');
t = 30; n = find(ds.idxT == t);
% operational target block: F107 and Dst from row t, Kp from KpFc(t+D)
Xo = tecGetBatch(dso, n); Xh = tecGetBatch(ds, n);
assert(abs(Xo(1,1,L+9) - Xh(1,1,L+1)) < 1e-6, 'op F107(t+D) == F107(t)');
kpfc = single((drv.KpFc(t+D) - ds.zmu(3)) / ds.zsig(3));
assert(abs(Xo(1,1,L+11) - kpfc) < 1e-5, 'op Kp(t+D) == KpFc');
dsz = buildTECDataset(TEC, ep, drv, 'L', L, 'Delta', D, 'Tmax', Tmax, 'zeroChannels', {'Kp','time'});
Xz = tecGetBatch(dsz, 1:5);
assert(all(Xz(:,:,[L+3 L+11 L+5:L+8 L+13:L+16],:) == 0, 'all') && any(Xz(:,:,L+1,:) ~= 0, 'all'), 'zeroChannels');
dst = buildTECDataset(TEC, ep, drv, 'L', L, 'Delta', D, 'Tmax', Tmax, 'zeroChannels', {'target'});
Xt = tecGetBatch(dst, 1:5);
assert(all(Xt(:,:,L+9:L+16,:) == 0, 'all') && any(Xt(:,:,L+1,:) ~= 0, 'all'), 'target block');
T2 = TEC; T2(:,:,60) = NaN;
dsn = buildTECDataset(T2, ep, drv, 'L', L, 'Delta', D, 'Tmax', Tmax);
assert(dsn.N == ds.N - (L + 1) && ~any(dsn.idxT >= 60 & dsn.idxT <= 60 + L - 1) && ~any(dsn.idxT + D == 60), 'NaN samples dropped (exactly L+1)');
msg = '';
end

function msg = testNet(H, W, L, TEC, ep, drv, D, Tmax)
ds = buildTECDataset(TEC(:,:,1:80), ep(1:80), drv(1:80,:), 'L', L, 'Delta', D, 'Tmax', Tmax);
C = ds.C; X2 = tecGetBatch(ds, 1:2);
net = createResTECNet(H, W, C, 'L', L, 'numBlocks', 2, 'numFilters', 8);
Y = predict(net, X2); assert(isequal(size(Y), [H W 1 2]), 'output size');
netU = createResTECNet(H, W, C, 'L', L, 'numBlocks', 2, 'numFilters', 8, 'uncertainty', true);
Y = predict(netU, X2); assert(isequal(size(Y), [H W 2 2]), 'uncertainty output size');
netZ = createResTECNet(H, W, C, 'L', L, 'numBlocks', 2, 'numFilters', 8, 'circular', false);
Y = predict(netZ, X2); assert(isequal(size(Y), [H W 1 2]), 'zero-pad output size');
netP = createResTECNet(H, W, C, 'L', L, 'numBlocks', 2, 'numFilters', 8, 'residual', false, 'globalSkip', false);
Y = predict(netP, X2); assert(isequal(size(Y), [H W 1 2]), 'plain CNN output size');
st = evaluateTEC(netU, ds, (87.5:-2.5:-87.5)');
assert(isfield(st, 'Sigma') && all(st.Sigma(:) > 0), 'Sigma from uncertainty head');
msg = '';
end

function msg = testLoss(H, W, lat)
wl = makeLatWeights(lat);
Y = dlarray(rand(H, W, 1, 3, 'single'), 'SSCB'); T = dlarray(rand(H, W, 1, 3, 'single'), 'SSCB');
l1 = lossTEC(Y, T, wl, 0.1); assert(isfinite(extractdata(l1)) && extractdata(l1) > 0, 'lossTEC');
Y2 = dlarray(cat(3, rand(H, W, 1, 3, 'single'), 0.1*randn(H, W, 1, 3, 'single')), 'SSCB');
l2 = lossTECNLL(Y2, T, wl, 0.1); assert(isfinite(extractdata(l2)), 'lossTECNLL');
assert(abs(mean(wl) - 1) < 1e-12, 'weights unit mean');
msg = '';
end

function msg = testModelFcns(H, W, L, D, TEC, ep, drv, Tmax)
ds = buildTECDataset(TEC(:,:,1:60), ep(1:60), drv(1:60,:), 'L', L, 'Delta', D, 'Tmax', Tmax);
C = ds.C; [X2, Y2] = tecGetBatch(ds, 1:2);
X = dlarray(X2, 'SSCB');
[p, cfg] = initConvLSTM(C, L, 4, 2);
Y = modelConvLSTM(p, X, cfg); Yc = extractdata(Y);
assert(isequal(size(Yc), [H W 1 2]) && all(isfinite(Yc(:))), sprintf('ConvLSTM shape: size %s, %d non-finite', mat2str(size(Yc)), nnz(~isfinite(Yc(:)))));
assert(all(isfinite(X2(:))), 'test input contains non-finite values');
[pt, cft] = initTransformer(H, W, C, L, 'd', 16, 'numLayers', 1, 'numHeads', 2);
Y = modelTransformer(pt, X, cft); Yn = extractdata(Y);
assert(isequal(size(Yn), [H W 1 2]) && all(isfinite(Yn(:))), ...
    sprintf('Transformer shape: size %s, dims %s, %d non-finite', mat2str(size(Yn)), dims(Y), nnz(~isfinite(Yn(:)))));
% gradient flows through both
wl = makeLatWeights((87.5:-2.5:-87.5)');
Yt = dlarray(Y2, 'SSCB');
g = dlfeval(@(pp) gradOf(pp, @(q, XX) modelTransformer(q, XX, cft), X, Yt, wl), pt);
assert(isstruct(g) && ~isempty(extractdata(g.Wout)), 'Transformer gradient');
msg = '';
end
function g = gradOf(p, f, X, Y, wl)
loss = lossTEC(f(p, X), Y, wl, 0.1); g = dlgradient(loss, p);
end

function msg = testBaselines(TEC, ep, drv, lat, lon, L, D, Tmax)
ds = buildTECDataset(TEC(:,:,301:end), ep(301:end), drv(301:end,:), 'L', L, 'Delta', D, 'Tmax', Tmax);
[st, info] = baselineSHAR(TEC, ep, ds, lat, lon, 'degree', 3, 'order', 6, 'window', 120, 'refit', 48);
assert(info.nCoef == 16 && isfinite(st.RMSE), 'SH-AR');
sp = persistenceBaseline(ds, lat); assert(isfinite(sp.RMSE) && numel(sp.perEpochRMSE) == ds.N, 'persistence');
Dst = drv.Dst(301:end); Dst = Dst(ds.idxT + D);
[T, sm, qm] = stormAnalysis({sp, st}, {'p','s'}, ds.epochsY, Dst, -100, 48);
assert(height(T) == 2 && ~any(sm & qm), 'stormAnalysis');
% RefY path
st2 = evaluateTEC([], setfield(ds, 'Yhat', ds.Tlast), lat, 'RefY', double(ds.Y) * Tmax + 1); %#ok<SFLD>
assert(abs(st2.Bias - (sp.Bias - 1)) < 1e-3, 'RefY offsets bias by -1');
msg = '';
end

function msg = testCalib(TEC, ep, drv, lat, L, D, Tmax)
ds = buildTECDataset(TEC(:,:,1:100), ep(1:100), drv(1:100,:), 'L', L, 'Delta', D, 'Tmax', Tmax);
N = ds.N;
stU.Yhat = single(double(ds.Y) * Tmax + 2 * randn(size(ds.Y)));
stU.Sigma = single(2 * ones(size(ds.Y)));
sm = false(N, 1); sm(1:10) = true;
R = calibrationStats(stU, ds, lat, sm);
assert(R.cov1 > 0.5 && R.cov1 < 0.85 && R.cov2 > 0.9 && abs(R.sigmaInflationStorm - 1) < 1e-6, 'coverage');
msg = '';
end

function msg = testStations(TEC, ep, drv, lat, lon, L, D, Tmax)
ds = buildTECDataset(TEC(:,:,1:120), ep(1:120), drv(1:120,:), 'L', L, 'Delta', D, 'Tmax', Tmax);
st.Yhat = single(double(ds.Y) * Tmax + 1);                    % model = reference + 1 TECU
codes = {'AAAA', 'BBBB', 'CCCC'}; lat0 = [38.7 40 36]; lon0 = [35.5 33 30];
stations = struct('code', {}, 'regime', {}, 'obs', {});
for s = 1:3
    tt = []; V = []; la = []; lo = []; el = [];
    for n = 1:numel(ds.epochsY)
        for j = 1:5
            tq = ds.epochsY(n) + minutes(5*j - 15);
            k = find(ep(1:120) == ds.epochsY(n), 1);         % obs value from the map valid at the hour
            ipla = lat0(s) + randn; iplo = lon0(s) + randn;
            v = interp2(lon, flipud(lat), flipud(TEC(:,:,k)), iplo, ipla);
            tt = [tt; tq]; V = [V; v]; la = [la; ipla]; lo = [lo; iplo]; el = [el; 40 + 10*rand]; %#ok<AGROW>
        end
    end
    obs = timetable(tt, V, la, lo, el, true(size(V)), ones(size(V)), 'VariableNames', {'VTEC','ippLat','ippLon','el','ok','mf'});
    obs.STEC = obs.VTEC; obs.prn = ones(height(obs), 1); obs.arc = ones(height(obs), 1); obs.az = zeros(height(obs), 1);
    stations(end+1) = struct('code', codes{s}, 'regime', tern(s < 3, 'MID', 'EIA'), 'obs', obs); %#ok<AGROW>
end
[Tst, Treg, hourly] = evaluateStationVTEC({st}, {'M'}, ds, lat, lon, stations, 'elevMin', 30);
assert(abs(Tst.M(1) - 1) < 0.3 && Tst.CODEGIM(1) < 0.3, 'station RMSE ~1 for +1 model, ~0 for reference');
assert(height(Treg) == 3 && all(isfinite(Treg.M)), 'regime table');
Tg = regionalGradients(hourly, {'AAAA','BBBB','A-B'});
assert(height(Tg) == 2 && Tg.RMSE_TECU(1) < 0.3, 'gradient of constant offset cancels');
Ts = evaluateSlantDelay({st}, {'M'}, ds, lat, lon, stations(1).obs, 'elevMin', 10);
assert(height(Ts) == 2 && abs(Ts.RMS_slant_m_all(1) - 0.162) < 0.05, 'slant error ~0.162 m for +1 TECU');
msg = '';
end

function msg = testGeo()
lla = geoTools('ecef2lla', [4.1e6 3.2e6 3.9e6]);
assert(abs(lla(1) - 37.98) < 1 && abs(lla(2) - 37.98) < 1, 'ecef2lla');
rx = [4.1e6 3.2e6 3.9e6]; sx = rx * 1.2 + [1e6 0 5e6];
ae = geoTools('azel', rx, sx); assert(ae{2} > 0 && ae{2} <= 90 && ae{1} >= 0 && ae{1} < 360, 'azel');
ip = geoTools('ipp', 38.7, 35.5, 180, 90, 450); assert(abs(ip{1} - 38.7) < 1e-6 && abs(ip{2} - 35.5) < 1e-6, 'zenith IPP');
ip = geoTools('ipp', 38.7, 35.5, 0, 30, 450); assert(ip{1} > 38.7, 'northward IPP');
M = geoTools('mslm', [90 30 10], 450, 0.9782); assert(abs(M(1) - 1) < 1e-9 && M(2) > 1 && M(3) > M(2), 'mslm');
a = [1e-8 0 -6e-8 0]; b = [9e4 0 -2e5 0];
dI = klobucharModel(a, b, 38.7, 35.5, [30; 80], [90; 90], repmat(datetime(2024,5,10,12,0,0,'TimeZone','UTC'), 2, 1));
assert(all(dI > 0) && dI(1) > dI(2), 'Klobuchar positive, larger at low elevation');
msg = '';
end

function msg = testKpParser()
f = [tempname '.txt']; fid = fopen(f, 'w');
fprintf(fid, ':Product: 3-Day Forecast\nNOAA Kp index breakdown May 09-May 11 2024\n\n');
fprintf(fid, '             May 09       May 10       May 11\n');
vals = [2.67 3 4.33; 2.33 3.33 5; 2 4 6.67; 2 4.67 7.33; 3 5.33 8; 3.67 6 8.33; 4 6.67 7; 4.33 7 6];
slots = {'00-03UT','03-06UT','06-09UT','09-12UT','12-15UT','15-18UT','18-21UT','21-00UT'};
for r = 1:8, fprintf(fid, '%s     %5.2f        %5.2f (G1)   %5.2f (G4)\n', slots{r}, vals(r,:)); end
fclose(fid);
kp8 = parseKpBreakdown(f, datetime(2024,5,10,'TimeZone','UTC'));
kpNone = parseKpBreakdown(f, datetime(2024,6,1,'TimeZone','UTC'));
delete(f);
assert(numel(kp8) == 8 && abs(kp8(6) - 6) < 1e-9 && abs(kp8(1) - 3) < 1e-9, 'day-2 column, storm tags removed');
assert(isempty(kpNone), 'missing day -> []');
msg = '';
end

function msg = testAdaptive(TEC, ep, drv, lat, lon, L, D, Tmax)
Ts = single(TEC) / Tmax;
F = computeSHARForecasts(Ts, ep, lat, lon, 'Delta', D, 'degree', 3, 'order', 6, 'window', 120, 'refit', 48);
assert(isequal(size(F), size(Ts)) && all(isnan(F(1,1,1:120))) && ~isnan(F(1,1,end)), 'forecast field shape / lead-in');
e = double(F(:,:,200:end) - Ts(:,:,200:end)) * Tmax;
assert(sqrt(mean(e(:).^2)) < 15, 'SH-AR field error unreasonable');
ds = buildTECDataset(Ts, ep, drv, 'L', L, 'Delta', D, 'Tmax', Tmax, 'scaled', true, 'base', F);
assert(ds.C == L + 17 && ds.N > 0 && all(~isnan(F(1,1,ds.idxT + D))), 'base channel dataset');
X = tecGetBatch(ds, 1:2);
assert(size(X, 3) == L + 17 && isequal(X(:,:,end,1), F(:,:,ds.idxT(1) + D)), 'base channel content');
net = createResTECNet(ds.H, ds.W, ds.C, 'L', L, 'numBlocks', 2, 'numFilters', 8, 'baseChannel', ds.C);
Y = predict(net, X); assert(isequal(size(Y), [ds.H ds.W 1 2]) && all(isfinite(Y(:))), 'Res-TECNet-A forward');
net0 = createResTECNet(ds.H, ds.W, L + 16, 'L', L, 'numBlocks', 2, 'numFilters', 8);
netT = transferWeights(net, net0);                       % partial copy of the stem conv (28 -> 29 channels)
Lt = netT.Learnables; L0 = net0.Learnables;
i = find(strcmp(Lt.Layer, 'stem_conv') & strcmp(Lt.Parameter, 'Weights'), 1);
j = find(strcmp(L0.Layer, 'stem_conv') & strcmp(L0.Parameter, 'Weights'), 1);
assert(isequal(extractdata(Lt.Value{i}(:,:,1:L+16,:)), extractdata(L0.Value{j})), 'transferWeights partial copy');
msg = '';
end

function msg = testEvalBrute(TEC, ep, drv, lat, L, D, Tmax)
ds = buildTECDataset(TEC(:,:,1:140), ep(1:140), drv(1:140,:), 'L', L, 'Delta', D, 'Tmax', Tmax);
rng(3); ds.Yhat = ds.Y + 0.02 * randn(size(ds.Y), 'single') + 0.01;     % scaled prediction
st = evaluateTEC([], ds, lat);
[H, W, ~, N] = size(ds.Y);
w = cosd(lat); w = w / mean(w);
P = double(ds.Yhat) * Tmax; T = double(ds.Y) * Tmax;
acc = 0; accA = 0; accB = 0; per = zeros(N, 1);
for n = 1:N
    a2 = 0; sw = 0;
    for i = 1:H
        for j = 1:W
            e = P(i, j, 1, n) - T(i, j, 1, n);
            acc = acc + w(i) * e^2; accA = accA + w(i) * abs(e); accB = accB + w(i) * e;
            a2 = a2 + w(i) * e^2; sw = sw + w(i);
        end
    end
    per(n) = sqrt(a2 / sw);
end
n3 = N * H * W;
assert(abs(st.RMSE - sqrt(acc / n3)) < 1e-6 * sqrt(acc / n3), 'RMSE');
assert(abs(st.MAE - accA / n3) < 1e-6 * (accA / n3), 'MAE');
assert(abs(st.Bias - accB / n3) < 1e-6 + 1e-6 * abs(accB / n3), 'Bias');
assert(max(abs(st.perEpochRMSE - per)) < 1e-6, 'perEpochRMSE');
msg = '';
end

function msg = testStationOffsets(TEC, ep, drv, lat, lon, L, D, Tmax)
ds = buildTECDataset(TEC(:,:,1:130), ep(1:130), drv(1:130,:), 'L', L, 'Delta', D, 'Tmax', Tmax);
N = ds.N; ref = double(ds.Y) * Tmax;
stA.Yhat = single(ref + 1);                                  % model = GIM + 1 TECU
stB.Yhat = single(ref + 2); stB.Yhat(:, :, :, 1:20) = NaN;   % model with missing forecasts
stations = struct('code', {}, 'regime', {}, 'obs', {}); offs = [-3 -5];
for s = 1:2
    tt = []; V = []; la = []; lo = []; el = [];
    for n = 1:N
        for j = 1:5
            tq = ds.epochsY(n) + minutes(5*j - 15);
            ipla = 38 + 0.5*randn; iplo = 35 + 0.5*randn;
            v = interp2(lon, flipud(lat), flipud(ref(:,:,1,n)), iplo, ipla) + offs(s);
            tt = [tt; tq]; V = [V; v]; la = [la; ipla]; lo = [lo; iplo]; el = [el; 45]; %#ok<AGROW>
        end
    end
    obs = timetable(tt, V, la, lo, el, true(size(V)), ones(size(V)), 'VariableNames', {'VTEC','ippLat','ippLon','el','ok','mf'});
    obs.STEC = obs.VTEC; obs.prn = ones(height(obs), 1); obs.arc = ones(height(obs), 1); obs.az = zeros(height(obs), 1);
    stations(end+1) = struct('code', sprintf('S%d', s), 'regime', 'MID', 'obs', obs); %#ok<AGROW>
end
[Tst, Treg, hourly, raw] = evaluateStationVTEC({stA, stB}, {'A', 'B'}, ds, lat, lon, stations, 'elevMin', 30);
assert(all(isfinite(Treg{:, :}), 'all'), 'NaN forecasts must not contaminate the pooled statistics');
assert(abs(hourly(1).biasVsCODE - offs(1)) < 0.3 && abs(hourly(2).biasVsCODE - offs(2)) < 0.3, 'station offsets');
% raw: model - station = 1 + 3 = 4 (A), 2 + 3 = 5 (B at S1); GIM: 3
assert(abs(raw.Tst{'S1', 'A'} - 4) < 0.4 && abs(raw.Tst{'S1', 'CODEGIM'} - 3) < 0.4 && abs(raw.Tst{'S2', 'A'} - 6) < 0.4, 'raw');
% offset removed: A -> 1, B -> 2, GIM -> 0
assert(abs(Tst{'S1', 'A'} - 1) < 0.4 && abs(Tst{'S2', 'B'} - 2) < 0.4 && Tst{'S1', 'CODEGIM'} < 0.4, 'offset removed');
assert(Tst.nHours(1) >= 0.9 * N && Tst.nHours(1) <= N, 'hourly binning found the observations (one group per hour)');
msg = '';
end

function msg = testRinexReader()
f = [tempname '.rnx']; fid = fopen(f, 'w');
hd = @(s, l) fprintf(fid, '%-60s%s\n', s, l);
hd('     3.04           OBSERVATION DATA    M', 'RINEX VERSION / TYPE');
hd(sprintf('%c  %3d%s', 'G', 6, sprintf(' %s', 'C1C', 'L1C', 'C2W', 'L2W', 'S1C', 'S2W')), 'SYS / # / OBS TYPES');
hd(sprintf('%c  %3d%s', 'E', 2, sprintf(' %s', 'C1X', 'L1X')), 'SYS / # / OBS TYPES');
hd('', 'END OF HEADER');
ep = @(sec, fl, n) fprintf(fid, '> %4d %02d %02d %02d %02d%11.7f  %d%3d\n', 2024, 5, 10, 0, floor(sec/60), mod(sec, 60), fl, n);
ob = @(id, x) fprintf(fid, '%3s%s\n', id, strjoin(arrayfun(@fld, x, 'UniformOutput', false), ''));
g01 = {[22000000.123 115600000.456 22000010.789 90100000.321 45.5 40.2], ...
       [22000100.123 115601000.456 22000110.789 90101000.321 45.6 40.3], ...
       [22000200.123 115602000.456 22000210.789 90102000.321 45.7 40.4]};
g05 = {[23000000.111 120000000.222 23000010.333 93500000.444 44.0 39.0], ...
       [23000100.111 120001000.222 NaN 93501000.444 44.1 NaN], ...
       [23000200.111 120002000.222 23000210.333 93502000.444 44.2 39.2]};
ep(0, 0, 3);  ob('G01', g01{1}); ob('E11', [25000000.5 131000000.5]); ob('G05', g05{1});
ep(30, 0, 3); ob('G01', g01{2}); ob('E11', [25000100.5 131001000.5]); ob('G05', g05{2});
ep(45, 4, 1); fprintf(fid, '%-60s%s\n', 'GNSS EVENT TEXT THAT STARTS WITH G', 'COMMENT');
ep(60, 0, 2); ob('G05', g05{3}); ob('G01', g01{3});
fclose(fid);
G = readRinexObsGPS(f); delete(f);
assert(height(G) == 6, 'GPS records only (Galileo and event line excluded)');
assert(isequal(G.SatelliteID, [1; 5; 1; 5; 5; 1]), 'satellite ids / order');
assert(all(ismember({'C1C', 'L1C', 'C2W', 'L2W'}, G.Properties.VariableNames)) && ~ismember('S1C', G.Properties.VariableNames), 'variables');
assert(abs(G.C1C(1) - 22000000.123) < 1e-3 && abs(G.L1C(2) - 120000000.222) < 1e-3 && abs(G.L2W(6) - 90102000.321) < 1e-3, 'values');
assert(isnan(G.C2W(4)) && ~isnan(G.C2W(3)) && isnan(G.C2W(3)) == false, 'blank field -> NaN');
assert(abs(seconds(G.Time(3) - G.Time(1)) - 30) < 1e-6 && abs(seconds(G.Time(5) - G.Time(1)) - 60) < 1e-6, 'epoch times');
msg = '';
end
function s = fld(x)
if isnan(x), s = repmat(' ', 1, 16); else, s = sprintf('%14.3f  ', x); end
end

function msg = testOrbitSagnac()
a = 5153.6^2; OmE = 7.2921151467e-5; mu = 3.986005e14;
eph = struct('Eccentricity', 0, 'sqrtA', 5153.6, 'Toe', 0, 'M0', 0, 'Delta_n', 0, 'omega', 0, 'OMEGA0', 0, 'OMEGA_DOT', 0, ...
    'i0', deg2rad(55), 'IDOT', 0, 'Cuc', 0, 'Cus', 0, 'Crc', 0, 'Crs', 0, 'Cic', 0, 'Cis', 0, 'SVClockBias', 1e-4, ...
    'SVClockDrift', 0, 'SVClockDriftRate', 0, 'TGD', 2e-9, 'Time', datetime(2024, 5, 5, 0, 0, 0));
P = 2 * pi / sqrt(mu / a^3);
[xyz, dts] = gpsSatPos(eph, [0; P/4; P/2; 3*P/4]);
assert(max(abs(sqrt(sum(xyz.^2, 2)) - a)) < 1e-3, 'circular orbit radius');
assert(abs(xyz(1, 1) - a) < 1e-3 && abs(xyz(1, 2)) < 1e-3 && abs(xyz(1, 3)) < 1e-3, 'start at the ascending node');
assert(abs(xyz(2, 3) - a * sin(deg2rad(55))) < 1e-2, 'maximum latitude after a quarter period');
assert(abs(xyz(3, 3)) < 1e-2, 'back at the equator after half a period');
assert(abs(dts(1) - (1e-4 - 2e-9)) < 1e-12, 'clock: bias minus TGD (no relativistic term for e = 0)');
% Earth rotation: a satellite on +x at travel time 0.07 s must appear rotated towards -y in the reception frame
r = geoTools('sagnac', [2.6e7 0 0], 0.07);
assert(r(2) < 0 && abs(r(2) + 2.6e7 * sin(OmE * 0.07)) < 1e-6 && abs(r(1) - 2.6e7 * cos(OmE * 0.07)) < 1e-6, 'Sagnac direction');
assert(abs(norm(r) - 2.6e7) < 1e-6, 'rotation preserves the range');
msg = '';
end

function msg = testHarmonize()
rng(11); nd = 40; per = 12; D = 5; noise = 1.5 * randn(nd, 1);          % true DCB 5 TECU; daily estimates D + noise
mf = 1.2 + 0.3 * rand(nd * per, 1); day = repelem((1:nd)', per);
t = datetime(2023, 1, 1, 'TimeZone', 'UTC') + days(day - 1) + hours(rand(nd * per, 1) * 20);
[t, o] = sort(t); mf = mf(o); day = day(o);
vTrue = 20 + 5 * rand(nd * per, 1);
est = D + noise(day);                                                   % estimate used on each day
stec = vTrue .* mf + D - est;                                           % processing subtracted the estimate: error D - est
% day 7: no estimate possible (DCB applied = 0 exactly), its VTEC is wrong by the full receiver DCB
d7 = day == 7; est(d7) = 0; stec(d7) = vTrue(d7) .* mf(d7) + D - 0;
obs = timetable(t, stec, stec ./ mf, mf, est, true(size(mf)), 'VariableNames', {'STEC', 'VTEC', 'mf', 'dRec', 'ok'});
[h, info] = harmonizeReceiverDCB(obs, 'halfWindowDays', 15, 'minDays', 5);
assert(info.nNoEstimate == 1, 'the day without estimate is counted');
errH7 = h.VTEC(d7) - vTrue(d7); errOthers = h.VTEC(~d7) - vTrue(~d7);
assert(abs(mean(errH7) - mean(errOthers)) < 1.0, 'the day without estimate is corrected with the neighbours'' median');
errRaw = obs.VTEC - vTrue; errH = h.VTEC - vTrue;
dayMean = @(e) accumarray(day, e, [], @mean);
dm = dayMean(errH); dr = dayMean(errRaw); keepD = (1:nd)' ~= 7;
assert(std(dm(keepD)) < 0.35 * std(dr(keepD)), 'day-to-day scatter of the VTEC error must shrink');
assert(abs(info.dailyScatterTECU - 1.5) < 0.45, 'reported daily scatter ~ noise level (day without estimate excluded)');
assert(max(abs(h.STEC - h.VTEC .* h.mf)) < 1e-9, 'STEC and VTEC stay consistent');
msg = '';
end
