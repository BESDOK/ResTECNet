function [obs, meta] = processStationVTEC(obsFile, navFile, dcbSat, opts)
%PROCESSSTATIONVTEC Dual-frequency GPS slant/vertical TEC at one station, one day.
%   [obs, meta] = PROCESSSTATIONVTEC(obsFile, navFile, dcbSat, ...)
%     obsFile : daily RINEX (2 or 3) observation file, 30-s sampling
%     navFile : broadcast navigation file (brdc)
%     dcbSat  : table from readCODEDCB (P1P2, ns); optionally with a
%               second table for P1C1 in opts.dcbP1C1 when the receiver
%               tracks C/A code instead of P1
%   Options (defaults implement Sect. 2.3 of the revised manuscript)
%     'elevCut'    elevation cut-off kept in the output (default 10 deg;
%                  the VTEC comparison applies 30 deg, the slant one 10 deg)
%     'hShell'     single-layer height [km]              (default 506.7,
%                  the shell height of the CODE MSLM, for consistency with the GIMs)
%     'alpha'      MSLM scaling factor                    (default 0.9782)
%     'minArcMin'  minimum arc length [min]               (default 30)
%     'dcbRec'     receiver DCB [ns]; [] -> estimated daily by minimizing
%                  the VTEC scatter across satellites in the 02-05 LT window
%     'rxXYZ'      receiver ECEF [m]; [] -> RINEX header APPROX POSITION
%
%   Processing steps (numbered as in Sect. 2.3):
%    (1) geometry-free phase L4 per continuous arc (cycle slips by the
%        Melbourne-Wubbena and L4-rate tests, arcs < 30 min discarded),
%        levelled to the code combination P4 by the sin^2(el)-weighted arc mean;
%    (2) satellite DCBs from CODE monthly solutions, receiver DCB estimated
%        daily (returned in meta.dcbRec_ns and meta.dcbRecScatter);
%    (3) STEC -> VTEC with the MSLM at hShell (alpha), IPP at the same height;
%    (4) residual screening at 3 x daily MAD (per satellite, on VTEC);
%    (5) output per observation, ready for evaluateStationVTEC.
%
%   Output timetable obs: Time, prn, el, az, ippLat, ippLon, STEC, VTEC,
%   arc (id), mf (mapping function), ok (screening flag), dRec (the daily
%   receiver DCB subtracted from STEC, in TECU; harmonizeReceiverDCB replaces
%   it by a running median).
%   Requires Navigation Toolbox (rinexread); header fields are read by
%   readRinexApproxXYZ.

arguments
    obsFile (1,:) char
    navFile (1,:) char
    dcbSat table
    opts.elevCut (1,1) double = 10
    opts.hShell (1,1) double = 506.7
    opts.alpha (1,1) double = 0.9782
    opts.minArcMin (1,1) double = 30
    opts.dcbRec double = []
    opts.rxXYZ double = []
    opts.dcbP1C1 = []
    opts.sampling (1,1) double = 60      % observations kept every 'sampling' s (RINEX is 30 s; 60 s is ample for hourly means)
    opts.reader (1,:) char {mustBeMember(opts.reader, {'fast','rinexread'})} = 'rinexread'   % measured: rinexread 2.5 s vs fast parser 4.4 s per daily file
    opts.dcbElev (1,1) double = 20       % lowest elevation used in the receiver-DCB estimation [deg]
    opts.dcbSign (1,1) double = 1        % P4_true = P4_obs + dcbSign*c*DCB_sat(P1-P2); validate with validateStationVTEC
end

c  = 299792458; f1 = 1575.42e6; f2 = 1227.60e6;
lam1 = c/f1; lam2 = c/f2;
Kion = 40.3e16 * (f1^2 - f2^2) / (f1^2 * f2^2);           % m per TECU on P4 (~0.105)

% ---------------- read RINEX ----------------
G = readObsGPS(obsFile, opts.reader);
tG = G.Time; G = G(mod(round(seconds(tG - dateshift(tG, 'start', 'day'))), opts.sampling) == 0, :);   % decimate
if isempty(opts.rxXYZ)
    rx = readRinexApproxXYZ(obsFile);
else
    rx = opts.rxXYZ(:)';
end
lla = geoTools('ecef2lla', rx); rxLat = lla(1); rxLon = lla(2);

Eph = loadGpsNav(navFile);

% observable selection (RINEX-3 codes; RINEX-2 mapped by rinexread)
[P1, nameP1] = pickCol(G, {'C1W','C1P','C1C'});
[P2, nameP2] = pickCol(G, {'C2W','C2P','C2L','C2X','C2S'});
[L1, nameL1] = pickCol(G, {'L1C','L1W','L1P'}); [L2, nameL2] = pickCol(G, {'L2W','L2P','L2L','L2X','L2S'});
if isempty(nameP1) || isempty(nameP2) || isempty(nameL1) || isempty(nameL2)
    error('processStationVTEC:obs', 'Dual-frequency GPS observables not found in %s.', obsFile);
end
usedC1 = strcmp(nameP1, 'C1C');                               % C/A code on L1: P1-C1 DCB needed
if usedC1 && isempty(opts.dcbP1C1)
    warning('processStationVTEC:c1c', 'C1C used on L1 but no P1-C1 DCB given (''dcbP1C1''): satellite-dependent errors up to ~6 TECU.');
end
if ~ismember(nameP2, {'C2W','C2P'})
    warning('processStationVTEC:c2', 'L2 code %s is not P(Y): the P2-C2 bias is not corrected.', nameP2);
end
ok = ~isnan(P1) & ~isnan(P2) & ~isnan(L1) & ~isnan(L2);
G = G(ok, :); P1 = P1(ok); P2 = P2(ok); L1 = L1(ok) * lam1; L2 = L2(ok) * lam2;
tg = G.Time; tg.TimeZone = '';                              % GPS time without zone (matches the navigation records)
t = tg; t.TimeZone = 'UTC';
prn = double(G.SatelliteID);

% ---------------- satellite geometry ----------------
N = height(G);
el = nan(N, 1); az = nan(N, 1);
% rinexread returns epochs in the file's time system (GPS time); weeks
% start on Sunday, so TOW follows directly (leap seconds are irrelevant here)
tow = seconds(t - dateshift(t, 'start', 'week'));
for s = unique(prn)'                                          % per satellite and ephemeris record (vectorized)
    iobs = find(prn == s);
    E = Eph(double(Eph.SatelliteID) == s, :); if isempty(E), continue; end
    [~, k] = min(abs(seconds(E.Time(:)' - tg(iobs))), [], 2);   % nearest broadcast record
    for ke = unique(k)'
        ii = iobs(k == ke);
        sxyz = gpsSatPos(E(ke, :), tow(ii) - P1(ii) / c);
        sxyz = geoTools('sagnac', sxyz, P1(ii) / c);
        ae = geoTools('azel', rx, sxyz);
        az(ii) = ae{1}; el(ii) = ae{2};
    end
end
keep = el >= opts.elevCut & ~isnan(el);
t = t(keep); prn = prn(keep); P1 = P1(keep); P2 = P2(keep); L1 = L1(keep); L2 = L2(keep);
az = az(keep); el = el(keep); N = numel(t);

% ---------------- combinations ----------------
P4 = P2 - P1;                                              % m (ionosphere delays P2 more)
L4 = L1 - L2;                                              % m (= P4 + arc ambiguity, phase advance)
MW = (f1*L1 - f2*L2)/(f1 - f2) - (f1*P1 + f2*P2)/(f1 + f2); % Melbourne-Wubbena, m

% ---------------- arcs and cycle slips (step 1) ----------------
arc = zeros(N, 1); STECphase = nan(N, 1); nextArc = 0;
l4max = 0.5 * opts.sampling / 30;                          % L4 rate threshold scaled with the sampling interval [m]
for s = unique(prn)'
    idx = find(prn == s); [~, o] = sort(t(idx)); idx = idx(o);
    ti = t(idx); brk = [true; seconds(diff(ti)) > 300];   % gap > 5 min
    mwRun = []; l4prev = NaN;
    for j = 1:numel(idx)
        slip = false;
        if ~brk(j)
            if numel(mwRun) >= 10 && abs(MW(idx(j)) - mean(mwRun)) > max(3 * std(mwRun), 1.0), slip = true; end
            if ~isnan(l4prev) && abs(L4(idx(j)) - l4prev) > l4max, slip = true; end    % L4 rate test
        end
        if brk(j) || slip
            nextArc = nextArc + 1; mwRun = [];
        end
        arc(idx(j)) = nextArc; mwRun(end+1) = MW(idx(j)); %#ok<AGROW>
        l4prev = L4(idx(j));
    end
end
% levelling: code-levelled phase, weight sin^2(el); discard short arcs
for a = unique(arc)'
    ia = find(arc == a);
    if minutes(max(t(ia)) - min(t(ia))) < opts.minArcMin, arc(ia) = 0; continue; end
    w = sind(el(ia)).^2;
    bias = sum(w .* (L4(ia) - P4(ia))) / sum(w);           % <L4 - P4> = arc ambiguity
    STECphase(ia) = L4(ia) - bias;                        % code-levelled phase = smooth P4, m
end
keep = arc > 0; t = t(keep); prn = prn(keep); az = az(keep); el = el(keep);
arc = arc(keep); P4lev = STECphase(keep); N = numel(t);

% ---------------- DCBs (step 2) ----------------
dsat = zeros(N, 1); prnStr = compose("G%02d", prn);
[tf, loc] = ismember(prnStr, string(dcbSat.prn)); dsat(tf) = dcbSat.dcb_ns(loc(tf));
if usedC1 && ~isempty(opts.dcbP1C1)                       % C/A receivers: P1 = C1 + DCB(P1C1)
    [tf, loc] = ismember(prnStr, string(opts.dcbP1C1.prn)); dsat(tf) = dsat(tf) - opts.dcbP1C1.dcb_ns(loc(tf));
end
mf = geoTools('mslm', el, opts.hShell, opts.alpha);
stecNoRec = (P4lev + opts.dcbSign * c * 1e-9 * dsat) / Kion;              % TECU, receiver DCB still inside
if isempty(opts.dcbRec)
    lt = mod(hour(t) + minute(t)/60 + rxLon/15, 24);
    night = lt >= 2 & lt <= 5 & el >= opts.dcbElev;
    if nnz(night) > 50
        % minimum within-epoch scatter across satellites (thin-shell model):
        % aa = VTEC + d*bb, bb = 1/mf; d minimizes the variance of aa - d*bb at each epoch
        aa = stecNoRec(night) ./ mf(night); bb = 1 ./ mf(night);
        g = findgroups(t(night)); cnt = accumarray(g, 1);
        ga = accumarray(g, aa, [], @mean); gb = accumarray(g, bb, [], @mean);
        ad = aa - ga(g); bd = bb - gb(g);
        use = cnt(g) >= 3;
        if nnz(use) > 30 && sum(bd(use).^2) > 0
            drecTECU = sum(ad(use) .* bd(use)) / sum(bd(use).^2);
        else
            warning('processStationVTEC:dcb', 'Too few multi-satellite night epochs; receiver DCB set to 0.'); drecTECU = 0;
        end
    else
        warning('processStationVTEC:dcb', 'Too few night-time observations; receiver DCB set to 0.');
        drecTECU = 0;
    end
    dcbRec_ns = drecTECU * Kion / (c * 1e-9);
else
    dcbRec_ns = opts.dcbRec; drecTECU = dcbRec_ns * c * 1e-9 / Kion;
end
STEC = stecNoRec - drecTECU;

% ---------------- VTEC, IPP (step 3) ----------------
VTEC = STEC ./ mf;
ip = geoTools('ipp', rxLat, rxLon, az, el, opts.hShell);
ippLat = ip{1}; ippLon = ip{2};

% ---------------- screening (step 4) ----------------
okFlag = true(N, 1);
for s = unique(prn)'
    i = prn == s; v = VTEC(i);
    madv = 1.4826 * median(abs(v - median(v)));
    okFlag(i) = abs(v - median(v)) <= 3 * max(madv, 1) | (numel(v) < 20);
end
okFlag = okFlag & VTEC > -5 & VTEC < 250;

dRec = repmat(drecTECU, N, 1);                              % receiver DCB used for this day [TECU] (see harmonizeReceiverDCB)
obs = timetable(t, prn, el, az, ippLat, ippLon, STEC, VTEC, arc, mf, okFlag, dRec, ...
    'VariableNames', {'prn','el','az','ippLat','ippLon','STEC','VTEC','arc','mf','ok','dRec'});
meta.rxXYZ = rx; meta.rxLat = rxLat; meta.rxLon = rxLon;
meta.dcbRec_ns = dcbRec_ns; meta.Kion = Kion; meta.hShell = opts.hShell; meta.alpha = opts.alpha;
meta.nArcs = numel(unique(arc)); meta.nObs = N; meta.usedC1 = usedC1;
meta.codes = {nameP1, nameP2, nameL1, nameL2};
fprintf('%s: %d obs, %d arcs, receiver DCB %.2f ns, %.1f %% screened\n', ...
    obsFile, N, meta.nArcs, dcbRec_ns, 100 * mean(~okFlag));
end

% -------------------------------------------------------------------------
function [v, name] = pickCol(G, names)
% first observable among 'names' whose column has data (> 50 % non-NaN)
v = nan(height(G), 1); name = '';
vn = G.Properties.VariableNames;
for k = 1:numel(names)
    if ismember(names{k}, vn) && mean(~isnan(G.(names{k}))) > 0.5
        v = G.(names{k}); name = names{k}; return
    end
end
end
