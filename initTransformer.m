function [p, cfg] = initTransformer(H, W, C, L, opts)
%INITTRANSFORMER Learnable parameters for modelTransformer.
%   [p, cfg] = INITTRANSFORMER(H, W, C, L, 'd', 256, 'numLayers', 6, 'numHeads', 8,
%                        'patch', [5 4], 'mlpRatio', 4)
%   Retained configuration of the revised manuscript: d = 256, 6 layers,
%   8 heads, 5 x 4 patches (selected on the validation set among
%   d in {128,256,512} and depth in {4,6,8}); about 4.3e6 parameters.
arguments
    H (1,1) double
    W (1,1) double
    C (1,1) double
    L (1,1) double
    opts.d (1,1) double = 256
    opts.numLayers (1,1) double = 6
    opts.numHeads (1,1) double = 8
    opts.patch (1,2) double = [5 4]
    opts.mlpRatio (1,1) double = 4
end
ph = opts.patch(1); pw = opts.patch(2); d = opts.d;
if mod(W, pw) ~= 0, error('initTransformer:patch', 'W must be a multiple of the patch width.'); end
Hp = ceil(H / ph) * ph; nT = (Hp / ph) * (W / pw);

cfg.L = L; cfg.patch = opts.patch; cfg.numLayers = opts.numLayers; cfg.numHeads = opts.numHeads;
g = @(varargin) dlarray(single(randn(varargin{:}) * 0.02));
z = @(n) dlarray(zeros(n, 1, 'single'));
o = @(n) dlarray(ones(n, 1, 'single'));

p.Wpatch = dlarray(single(randn(ph, pw, C, d) * sqrt(2 / (ph * pw * C))));
p.bpatch = z(d);
p.pos = dlarray(single(randn(d, nT) * 0.02));         % d x nT, unformatted (broadcast over B)
dm = opts.mlpRatio * d;
layers = struct([]);
for l = 1:opts.numLayers
    layers(l).ln1b = z(d); layers(l).ln1g = o(d);
    layers(l).Wq = g(d, d); layers(l).bq = z(d);
    layers(l).Wk = g(d, d); layers(l).bk = z(d);
    layers(l).Wv = g(d, d); layers(l).bv = z(d);
    layers(l).Wo = g(d, d); layers(l).bo = z(d);
    layers(l).ln2b = z(d); layers(l).ln2g = o(d);
    layers(l).W1 = g(dm, d); layers(l).b1 = z(dm);
    layers(l).W2 = g(d, dm); layers(l).b2 = z(d);
end
p.layers = layers;
p.lnfb = z(d); p.lnfg = o(d);
p.Wout = g(ph * pw, d); p.bout = z(ph * pw);
end
