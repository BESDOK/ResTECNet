function F = computeSHARForecasts(Ts, epochs, lat, lon, opts)
%COMPUTESHARFORECASTS SH-AR 24-h forecast for every epoch of the record.
%   F = COMPUTESHARFORECASTS(Ts, epochs, lat, lon, ...) returns F, H x W x K
%   single, where F(:,:,k) is the spherical-harmonic autoregressive forecast
%   for epoch k issued at epoch k - Delta from data up to k - Delta only
%   (NaN where the lead-in is insufficient). Same method as baselineSHAR
%   (degree 15, sun-fixed frame, direct AR of order p over a 30-day window,
%   refitted every 'refit' hours), computed once for the whole record so that
%   it can serve as the adaptive base of Res-TECNet-A (Sect. 3.2) and as the
%   SH-AR baseline. Ts may be scaled (single); F is returned in the same
%   unit. Memory: two nCoef x K coefficient arrays (~1 GB for 26 years).
arguments
    Ts (:,:,:) {mustBeNumeric}
    epochs (:,1) datetime
    lat (:,1) double
    lon (:,1) double
    opts.Delta (1,1) double = 24
    opts.degree (1,1) double = 15
    opts.order (1,1) double = 48
    opts.window (1,1) double = 720
    opts.refit (1,1) double = 168
end
[H, W, K] = size(Ts); D = opts.Delta; p = opts.order; win = opts.window;
dlon = lon(2) - lon(1);
t0 = tic;
% ---- basis and weighted least-squares solver ----
[LON, LAT] = meshgrid(lon, lat);
x = sind(LAT(:)); lam = deg2rad(LON(:)); nCoef = (opts.degree + 1)^2;
A = zeros(H * W, nCoef); c = 0;
for n = 0:opts.degree
    P = legendre(n, x', 'norm')';
    for m = 0:n
        c = c + 1; A(:, c) = P(:, m + 1) .* cos(m * lam);
        if m > 0, c = c + 1; A(:, c) = P(:, m + 1) .* sin(m * lam); end
    end
end
w = repmat(makeLatWeights(lat), 1, W); w = w(:);
S = (A' * (w .* A)) \ (A' .* w');                       % nCoef x (H*W)
% ---- coefficients in the sun-fixed frame, grouped by UT hour ----
ut = hour(epochs); shift = round(15 * ut / dlon);
Cf = zeros(nCoef, K);
for h = 0:23
    idx = find(ut == h); if isempty(idx), continue; end
    for b = 1:5000:numel(idx)
        j = idx(b:min(b+4999, numel(idx)));
        M = circshift(double(Ts(:, :, j)), -shift(j(1)), 2);
        Cf(:, j) = S * reshape(M, H * W, []);
    end
end
bad = any(isnan(Cf), 1);
fprintf('SH-AR: %d coefficients per map computed for %d epochs (%.1f min)\n', nCoef, K, toc(t0)/60);
% ---- direct AR forecast, refitted every 'refit' hours ----
Chat = nan(nCoef, K);                                   % Chat(:,k) = forecast for epoch k
tStart = win + p + 1;
for tb = tStart:opts.refit:(K - D)
    tt = (tb - win + 1 : tb)'; tt = tt(tt - p >= 1 & tt + D <= tb);
    ok = true(size(tt)); for k = 1:p, ok = ok & ~bad(tt - k + 1)'; end; ok = ok & ~bad(tt + D)';
    tt = tt(ok); if numel(tt) < 2 * p, continue; end
    Xl = zeros(numel(tt), p); a = zeros(p, nCoef);
    for cc = 1:nCoef
        for k = 1:p, Xl(:, k) = Cf(cc, tt - k + 1)'; end
        a(:, cc) = lsqminnorm(Xl, Cf(cc, tt + D)');
    end
    for t = tb:min(tb + opts.refit - 1, K - D)
        if any(bad(t - p + 1 : t)), continue; end
        lags = Cf(:, t : -1 : t - p + 1);
        Chat(:, t + D) = sum(lags .* a', 2);
    end
end
fprintf('SH-AR: autoregression done (%.1f min)\n', toc(t0)/60);
% ---- reconstruct and rotate back ----
F = nan(H, W, K, 'single');
for h = 0:23
    idx = find(ut == h & ~isnan(Chat(1, :))'); if isempty(idx), continue; end
    for b = 1:5000:numel(idx)
        j = idx(b:min(b+4999, numel(idx)));
        M = reshape(A * Chat(:, j), H, W, []);
        F(:, :, j) = single(circshift(M, shift(j(1)), 2));
    end
end
fprintf('SH-AR: %d forecast maps reconstructed (%.1f min)\n', nnz(~isnan(F(1,1,:))), toc(t0)/60);
end
