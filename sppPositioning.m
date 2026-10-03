function [T, sol] = sppPositioning(obsFile, navFile, corrections, opts)
%SPPPOSITIONING Single-frequency GPS L1 code positioning with ionospheric corrections.
%   [T, sol] = SPPPOSITIONING(obsFile, navFile, corrections, ...)
%     corrections : struct array with fields
%        name : label, e.g. 'None', 'Klobuchar', 'Res-TECNet operational'
%        type : 'none' | 'klobuchar' | 'map'
%        maps : for 'map' -> struct with .Yhat (H x W x 1 x N, TECU),
%               .epochs (N x 1 datetime), .lat, .lon  (hourly maps of the
%               forecast/reference product; interpolated in space to the
%               IPP and in time between the enclosing hours)
%        hShell, alpha : for 'map' -> MSLM parameters (Sect. 2.3)
%     Options: 'elevCut' (default 10), 'truthXYZ' (ECEF [m]; default RINEX
%              header position), 'interval' (s, default 30), 'weighting'
%              elevation-dependent 1/sin^2(el) (default true)
%   Epoch-wise weighted least squares on C/A pseudoranges with broadcast
%   orbits/clocks, relativistic and TGD corrections, Saastamoinen
%   troposphere (standard atmosphere) and the selected ionospheric
%   correction. Returns horizontal/vertical RMS and 95th-percentile
%   errors [m] per correction (Sect. 4.7.1, Table 6) and the solutions.
%   Requires Navigation Toolbox (rinexread / rinexinfo).

arguments
    obsFile (1,:) char
    navFile (1,:) char
    corrections struct
    opts.elevCut (1,1) double = 10
    opts.truthXYZ double = []
    opts.interval (1,1) double = 60
    opts.weighting (1,1) logical = true
end
c = 299792458;
G = readObsGPS(obsFile);
Eph = loadGpsNav(navFile);
truth = opts.truthXYZ; if isempty(truth), truth = readRinexApproxXYZ(obsFile); end
lla = geoTools('ecef2lla', truth);
try, [alphaK, betaK] = readNavIono(navFile); catch, alphaK = []; betaK = []; end
if ismember('C1C', G.Properties.VariableNames), C1 = G.C1C; else, C1 = G.C1W; end
ok = ~isnan(C1); G = G(ok, :); C1 = C1(ok);
tg = G.Time; tg.TimeZone = ''; t = tg; t.TimeZone = 'UTC'; prn = double(G.SatelliteID);
epochs = unique(t); epochs = epochs(mod(seconds(epochs - epochs(1)), opts.interval) == 0);
nC = numel(corrections);
err = nan(numel(epochs), 3, nC);
tow = seconds(t - dateshift(t, 'start', 'week'));

for e = 1:numel(epochs)
    i = find(t == epochs(e));
    if numel(i) < 5, continue; end
    % satellite states
    sx = nan(numel(i), 3); dts = nan(numel(i), 1);
    for k = 1:numel(i)
        E = Eph(Eph.SatelliteID == prn(i(k)), :);
        if isempty(E), continue; end
        [~, j] = min(abs(seconds(E.Time - tg(i(k))))); E = E(j, :);
        [p, dt] = gpsSatPos(E, tow(i(k)) - C1(i(k)) / c);
        p = geoTools('sagnac', p, C1(i(k)) / c);                % Sagnac (frame of reception)
        sx(k, :) = p; dts(k) = dt;
    end
    good = ~isnan(dts);
    i = i(good); sx = sx(good, :); dts = dts(good);
    for ci = 1:nC
        x = [truth, 0]; % initial guess (position + clock)
        for it = 1:8
            d = sx - x(1:3); rho = sqrt(sum(d.^2, 2));
            ae = geoTools('azel', x(1:3), sx); az = ae{1}; el = ae{2};
            use = el >= opts.elevCut;
            if nnz(use) < 4, break; end
            trop = saastamoinen(el, lla(3));
            switch corrections(ci).type
                case 'none',      iono = zeros(size(el));
                case 'klobuchar', iono = klobucharModel(alphaK, betaK, lla(1), lla(2), el, az, repmat(epochs(e), numel(el), 1));
                case 'map'
                    Mp = corrections(ci).maps;
                    ip = geoTools('ipp', lla(1), lla(2), az, el, corrections(ci).hShell);
                    v = mapValueAtTime(Mp, epochs(e), ip{1}, ip{2});
                    iono = v .* geoTools('mslm', el, corrections(ci).hShell, corrections(ci).alpha) * 0.162;
            end
            pr = C1(i) + c * dts - trop - iono;
            r = pr - (rho + x(4));
            Hm = [-d ./ rho, ones(numel(rho), 1)];
            w = ones(size(el)); if opts.weighting, w = sind(el).^2; end
            w(~use) = 0;
            dx = (Hm' * (w .* Hm)) \ (Hm' * (w .* r));
            x = x + dx';
            if norm(dx(1:3)) < 1e-3, break; end
        end
        if nnz(use) >= 4
            enu = ecef2enu(x(1:3) - truth, lla(1), lla(2));
            err(e, :, ci) = enu;
        end
    end
end

names = {corrections.name}';
hRMS = zeros(nC, 1); h95 = hRMS; vRMS = hRMS; v95 = hRMS;
for ci = 1:nC
    h = hypot(err(:, 1, ci), err(:, 2, ci)); v = abs(err(:, 3, ci));
    hRMS(ci) = sqrt(mean(h.^2, 'omitnan')); h95(ci) = prctile(h, 95);
    vRMS(ci) = sqrt(mean(v.^2, 'omitnan')); v95(ci) = prctile(v, 95);
end
T = table(names, hRMS, h95, vRMS, v95, 'VariableNames', ...
    {'Correction', 'Horizontal_RMS_m', 'Horizontal_95_m', 'Vertical_RMS_m', 'Vertical_95_m'});
sol.epochs = epochs; sol.errENU = err; sol.truthXYZ = truth;
disp(T);
end

% -------------------------------------------------------------------------
function trop = saastamoinen(el, h)
h = max(h, 0);
P = 1013.25 * (1 - 2.2557e-5 * h)^5.2568; Tk = 15 - 6.5e-3 * h + 273.15;
e = 6.108 * 0.5 * exp((17.15 * Tk - 4684) / (Tk - 38.45));
z = deg2rad(90 - el);
trop = 0.002277 ./ cos(z) .* (P + (1255 ./ Tk + 0.05) * e - tan(z).^2);
end

function enu = ecef2enu(d, lat, lon)
lat = deg2rad(lat); lon = deg2rad(lon);
Rm = [-sin(lon) cos(lon) 0; -sin(lat)*cos(lon) -sin(lat)*sin(lon) cos(lat); ...
      cos(lat)*cos(lon) cos(lat)*sin(lon) sin(lat)];
enu = (Rm * d(:))';
end

function v = mapValueAtTime(Mp, tq, pLat, pLon)
% linear interpolation in time between the two enclosing hourly maps
k = find(Mp.epochs <= tq, 1, 'last');
if isempty(k), k = 1; end
k2 = min(k + 1, numel(Mp.epochs));
a = 0; if k2 > k, a = seconds(tq - Mp.epochs(k)) / seconds(Mp.epochs(k2) - Mp.epochs(k)); end
v1 = interpMapAtPoints(double(Mp.Yhat(:, :, 1, k)), Mp.lat, Mp.lon, pLat, pLon);
v2 = interpMapAtPoints(double(Mp.Yhat(:, :, 1, k2)), Mp.lat, Mp.lon, pLat, pLon);
v = (1 - a) * v1 + a * v2;
end

function v = interpMapAtPoints(map, lat, lon, pLat, pLon)
[H, W] = size(map);
dlon = lon(2) - lon(1);
fj = (pLon - lon(1)) / dlon; j1 = floor(fj); b = fj - j1; j1 = mod(j1, W) + 1; j2 = mod(j1, W) + 1;
fi = (pLat - lat(1)) / (lat(2) - lat(1)); i1 = floor(fi); a = min(max(fi - i1, 0), 1);
i1 = min(max(i1 + 1, 1), H - 1); i2 = i1 + 1;
v = (1-a).*(1-b).*map(sub2ind([H W], i1, j1)) + (1-a).*b.*map(sub2ind([H W], i1, j2)) ...
  + a.*(1-b).*map(sub2ind([H W], i2, j1)) + a.*b.*map(sub2ind([H W], i2, j2));
end
