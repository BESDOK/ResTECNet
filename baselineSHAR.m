function [stats, info] = baselineSHAR(TEC, epochs, ds, lat, lon, opts)
%BASELINESHAR Spherical-harmonic coefficient autoregression (Sect. 3.6, iii).
%   [stats, info] = BASELINESHAR(TEC, epochs, ds, lat, lon, ...)
%     TEC, epochs : hourly CODE record covering ds.epochsY AND a lead-in of
%                   at least opts.window + opts.order + Delta hours
%     ds          : test dataset (uses .epochsY, .Delta, .Tmax, .Y)
%     lat, lon    : grid vectors (W = 72 unique meridians)
%   Options
%     'degree'  SH degree/order (default 15 -> 256 coefficients)
%     'order'   AR order p in hours (default 48, selected on the
%               validation set among {24, 48, 96})
%     'window'  fitting window in hours (default 720 = 30 days)
%     'refit'   re-estimate the AR coefficients every 'refit' hours
%               (default 168)
%   Method: every map is rotated to a sun-fixed frame (15 deg/h, exact
%   column shifts on the 5-deg grid), projected onto real spherical
%   harmonics by weighted least squares, and each coefficient series is
%   extrapolated 24 h ahead by a direct AR predictor
%        c(t+D) = sum_{k=1..p} a_k c(t-k+1)
%   fitted by least squares over the sliding window; the predicted map
%   is reconstructed and rotated back to the Earth-fixed frame.

arguments
    TEC (:,:,:) {mustBeNumeric}
    epochs (:,1) datetime
    ds struct
    lat (:,1) double
    lon (:,1) double
    opts.degree (1,1) double = 15
    opts.order (1,1) double = 48
    opts.window (1,1) double = 720
    opts.refit (1,1) double = 168
end

[H, W, K] = size(TEC);
D = ds.Delta; p = opts.order; win = opts.window;
dlon = lon(2) - lon(1);

% ---- real SH basis on the grid ----
[A, nCoef] = shBasis(lat, lon, opts.degree);            % (H*W) x nCoef
w = repmat(makeLatWeights(lat), 1, W); w = w(:);
S = (A' * (w .* A)) \ (A' .* w');                        % nCoef x (H*W) solver

% ---- coefficient series in the sun-fixed frame ----
ut = hour(epochs) + minute(epochs) / 60;
shift = round(15 * ut / dlon);                           % cells
Cf = zeros(nCoef, K);
for k = 1:K
    m = circshift(TEC(:, :, k), -shift(k), 2);           % to sun-fixed
    Cf(:, k) = S * m(:);
end

% ---- direct AR forecast, refitted every 'refit' hours ----
[tf, kt] = ismember(dateshift(ds.epochsY, 'start', 'hour'), epochs);
if ~all(tf), error('baselineSHAR:cover', 'TEC record must cover all target epochs.'); end
N = numel(kt);
Yhat = nan(H, W, 1, N, 'single');
a = []; lastFit = -inf;
for n = 1:N
    kTar = kt(n); t = kTar - D;                          % current epoch index
    if t - p - win < 1, continue; end
    if isempty(a) || (t - lastFit) >= opts.refit
        % fit a_k for every coefficient: rows = issue epochs in window
        tt = (t - win + 1 : t)';                         % issue epochs
        tt = tt(tt - p >= 1 & tt + D <= t);              % targets inside history
        Xlag = zeros(numel(tt), p);
        a = zeros(p, nCoef);
        for c = 1:nCoef
            for k = 1:p, Xlag(:, k) = Cf(c, tt - k + 1)'; end
            a(:, c) = lsqminnorm(Xlag, Cf(c, tt + D)');
        end
        lastFit = t;
    end
    lags = Cf(:, t : -1 : t - p + 1);                    % nCoef x p
    cHat = sum(lags .* a', 2);
    m = reshape(A * cHat, H, W);
    Yhat(:, :, 1, n) = single(circshift(m, shift(kTar), 2)) / ds.Tmax;
end
if any(isnan(Yhat(1, 1, 1, :)))
    warning('baselineSHAR:leadin', '%d targets skipped (insufficient lead-in).', nnz(isnan(Yhat(1,1,1,:))));
end
dsP = ds; dsP.Yhat = Yhat;
stats = evaluateTEC([], dsP, lat);
info.nCoef = nCoef; info.order = p; info.window = win;
end

% -------------------------------------------------------------------------
function [A, nCoef] = shBasis(lat, lon, nmax)
[LON, LAT] = meshgrid(lon, lat);
x = sind(LAT(:)); lam = deg2rad(LON(:));
nCoef = (nmax + 1)^2;
A = zeros(numel(x), nCoef); c = 0;
for n = 0:nmax
    P = legendre(n, x', 'norm')';                        % numel x (n+1), m = 0..n
    for m = 0:n
        c = c + 1; A(:, c) = P(:, m + 1) .* cos(m * lam);
        if m > 0
            c = c + 1; A(:, c) = P(:, m + 1) .* sin(m * lam);
        end
    end
end
end
