function R = validateStationVTEC(obs, TEC, epochs, lat, lon, opts)
%VALIDATESTATIONVTEC Sanity check of a processed station day against the CODE GIM.
%   R = VALIDATESTATIONVTEC(obs, TEC, epochs, lat, lon)
%     obs    : timetable from processStationVTEC
%     TEC    : H x W x K CODE maps [TECU] (e.g. TECmay2024.mat), epochs: K x 1 UTC
%   Every accepted observation (el >= 'elevMin') is compared with the GIM
%   interpolated to its ionospheric pierce point (bilinear in space with
%   circular longitude, linear in time). Expected for a healthy chain at a
%   mid-latitude station: |bias| below ~2 TECU, RMS of 2-4 TECU. A wrong
%   sign of the satellite DCBs or an orbit/geometry error gives RMS > 8 TECU.
%   R: n, bias, rms, corr (TECU), and R.bySat = [prn, bias, rms] per satellite.
arguments
    obs timetable
    TEC (:,:,:) {mustBeNumeric}
    epochs (:,1) datetime
    lat (:,1) double
    lon (:,1) double
    opts.elevMin (1,1) double = 30
end
epochs.TimeZone = 'UTC';
O = obs(obs.ok & obs.el >= opts.elevMin, :);
tt = O.Properties.RowTimes; tt.TimeZone = 'UTC';
kf = interp1(datenum(epochs), 1:numel(epochs), datenum(tt), 'linear');
valid = ~isnan(kf); O = O(valid, :); kf = kf(valid);
k0 = floor(kf); w = kf - k0; k1 = min(k0 + 1, numel(epochs));
gim = nan(height(O), 1);
for kk = unique(k0)'
    i = k0 == kk;
    v0 = interpPoints(double(TEC(:, :, kk)), lat, lon, O.ippLat(i), O.ippLon(i));
    v1 = interpPoints(double(TEC(:, :, kk + (kk < numel(epochs)))), lat, lon, O.ippLat(i), O.ippLon(i));
    gim(i) = (1 - w(i)) .* v0 + w(i) .* v1;
end
d = O.VTEC - gim;
R.n = numel(d); R.bias = mean(d, 'omitnan'); R.rms = sqrt(mean(d.^2, 'omitnan'));
c = corrcoef(O.VTEC, gim, 'Rows', 'complete'); R.corr = c(1, 2);
sats = unique(O.prn)'; bySat = zeros(numel(sats), 3);          % columns: prn, bias, rms
for j = 1:numel(sats)
    m = O.prn == sats(j);
    bySat(j, :) = [sats(j), mean(d(m), 'omitnan'), sqrt(mean(d(m).^2, 'omitnan'))];
end
R.bySat = bySat;
fprintf('Station VTEC - CODE GIM (%d IPPs, el >= %d): bias %.2f, RMS %.2f TECU, corr %.3f; per-satellite bias range %.1f .. %.1f TECU\n', ...
    R.n, opts.elevMin, R.bias, R.rms, R.corr, min(bySat(:, 2)), max(bySat(:, 2)));
end

function v = interpPoints(map, lat, lon, pLat, pLon)
[H, W] = size(map); dlon = lon(2) - lon(1);
fj = (pLon - lon(1)) / dlon; j1 = floor(fj); b = fj - j1; j1 = mod(j1, W) + 1; j2 = mod(j1, W) + 1;
fi = (pLat - lat(1)) / (lat(2) - lat(1)); i1 = floor(fi); a = min(max(fi - i1, 0), 1);
i1 = min(max(i1 + 1, 1), H - 1); i2 = i1 + 1;
v = (1-a).*(1-b).*map(sub2ind([H W], i1, j1)) + (1-a).*b.*map(sub2ind([H W], i1, j2)) ...
  + a.*(1-b).*map(sub2ind([H W], i2, j1)) + a.*b.*map(sub2ind([H W], i2, j2));
end
