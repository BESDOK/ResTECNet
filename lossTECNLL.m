function loss = lossTECNLL(Y, T, wlat, lambda)
%LOSSTECNLL Latitude-weighted heteroscedastic Gaussian NLL + gradient L1
%   (Eq. 7 of the revised manuscript, Sect. 3.5; Kendall & Gal, 2017).
%   loss = LOSSTECNLL(Y, T, wlat, lambda)
%     Y      : network output, dlarray H x W x 2 x B (SSCB):
%              channel 1 = predicted (scaled) TEC, channel 2 = log-variance
%     T      : target, dlarray H x W x 1 x B
%     wlat   : H x 1 latitude weights (unit mean)
%     lambda : gradient-term weight (paper: 0.1)
%
%   The gradient term acts on the mean channel only. Use with trainnet as
%   @(Y,T) lossTECNLL(Y, T, wlat, 0.1) together with
%   createResTECNet(..., 'uncertainty', true).

wl = dlarray(single(wlat(:)));
mu = Y(:, :, 1, :);
s  = Y(:, :, 2, :);                          % log sigma^2
s  = max(min(s, 10), -10);                   % numerical guard

nll = 0.5 * exp(-s) .* (mu - T).^2 + 0.5 * s;
nllTerm = mean(wl .* nll, 'all');

dYi = mu(2:end, :, :, :) - mu(1:end-1, :, :, :);
dTi = T(2:end, :, :, :)  - T(1:end-1, :, :, :);
dYj = cat(2, mu(:, 2:end, :, :), mu(:, 1, :, :)) - mu;
dTj = cat(2, T(:, 2:end, :, :),  T(:, 1, :, :))  - T;
gradTerm = mean(abs(dYi - dTi), 'all') + mean(abs(dYj - dTj), 'all');

loss = nllTerm + lambda * gradTerm;
end
