function loss = lossTEC(Y, T, wlat, lambda)
%LOSSTEC Latitude-weighted MSE + gradient-difference L1 loss (Eq. 6).
%   loss = LOSSTEC(Y, T, wlat, lambda)
%     Y      : network prediction, dlarray H x W x 1 x B (SSCB)
%     T      : target,             dlarray H x W x 1 x B
%     wlat   : H x 1 latitude weights, w_i = cosd(lat_i), normalized so
%              that mean(wlat) == 1 (see makeLatWeights)
%     lambda : weight of the gradient term (paper: 0.1)
%
%   Use with trainnet as:  @(Y,T) lossTEC(Y, T, wlat, 0.1)

wl = dlarray(single(wlat(:)));            % H x 1, broadcasts over W,C,B

% latitude-weighted MSE
e2 = (Y - T).^2;
mseTerm = mean(wl .* e2, 'all');

% forward-difference gradients (latitude and longitude, circular in lon)
dYi = Y(2:end, :, :, :) - Y(1:end-1, :, :, :);
dTi = T(2:end, :, :, :) - T(1:end-1, :, :, :);
dYj = cat(2, Y(:, 2:end, :, :), Y(:, 1, :, :)) - Y;
dTj = cat(2, T(:, 2:end, :, :), T(:, 1, :, :)) - T;

gradTerm = mean(abs(dYi - dTi), 'all') + mean(abs(dYj - dTj), 'all');

loss = mseTerm + lambda * gradTerm;
end
