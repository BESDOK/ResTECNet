function Y = modelConvLSTM(p, X, cfg)
%MODELCONVLSTM Two-layer ConvLSTM baseline (Sect. 3.6, vi).
%   Y = MODELCONVLSTM(params, X, cfg), X dlarray H x W x C x B ('SSCB') built by
%   buildTECDataset (C = L + 16). The L historical maps are consumed as a
%   sequence; at every step the 16 broadcast driver channels (z_t and
%   z_{t+D}) are concatenated to the current map so that the recurrent
%   baseline receives exactly the same information as Res-TECNet.
%   Circular longitude padding, a 3x3 linear head and the global skip
%   (T_t + R) are identical to Res-TECNet. Output H x W x 1 x B.
%   [params, cfg] from initConvLSTM(C, L, numFilters, numLayers); use as
%   modelFcn = @(p, X) modelConvLSTM(p, X, cfg) (params holds dlarrays only).

L = cfg.L; nl = numel(p.layers);
[H, W, ~, B] = size(X);
drv = X(:, :, L+1:end, :);                      % 16 driver channels
Tt  = X(:, :, L, :);
h = cell(nl, 1); c = cell(nl, 1);
for l = 1:nl
    z = zeros(H, W, cfg.nf, B, 'like', extractdata(X(1, 1, 1, 1)));
    h{l} = dlarray(z, 'SSCB'); c{l} = h{l};
end
for k = 1:L
    inp = cat(3, X(:, :, k, :), drv);
    for l = 1:nl
        [h{l}, c{l}] = cellStep(inp, h{l}, c{l}, p.layers(l), cfg.nf);
        inp = h{l};
    end
end
R = dlconv(circPad(h{nl}), p.Whead, p.bhead, 'Padding', [1 1 0 0]);
Y = Tt + R;
end

function [h, c] = cellStep(x, h, c, q, nf)
g = dlconv(circPad(cat(3, x, h)), q.W, q.b, 'Padding', [1 1 0 0]);   % 4*nf channels
i = sigmoid(g(:, :, 1:nf, :));
f = sigmoid(g(:, :, nf+1:2*nf, :));
o = sigmoid(g(:, :, 2*nf+1:3*nf, :));
u = tanh(g(:, :, 3*nf+1:4*nf, :));
c = f .* c + i .* u;
h = o .* tanh(c);
end

function Z = circPad(X)
Z = cat(2, X(:, end, :, :), X, X(:, 1, :, :));
end
