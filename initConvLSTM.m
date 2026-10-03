function [p, cfg] = initConvLSTM(C, L, numFilters, numLayers)
%INITCONVLSTM Learnable parameters for modelConvLSTM (He initialization).
%   [p, cfg] = INITCONVLSTM(C, L, numFilters, numLayers)
%   p holds the learnable dlarrays only; cfg the fixed settings (L, nf).
%   Retained configuration of the revised manuscript: 64 filters, 2
%   layers (selected on the validation set among {32,64,128} x {1,2,3}).
arguments
    C (1,1) double
    L (1,1) double
    numFilters (1,1) double = 64
    numLayers (1,1) double = 2
end
cfg.L = L; cfg.nf = numFilters;
cin = 1 + (C - L);                               % one map + 16 drivers
layers = struct('W', {}, 'b', {});
for l = 1:numLayers
    ci = cin + numFilters;                       % [x, h] concatenated
    layers(l).W = dlarray(single(randn(3, 3, ci, 4 * numFilters) * sqrt(2 / (9 * ci))));
    b = zeros(4 * numFilters, 1, 'single');
    b(numFilters+1 : 2*numFilters) = 1;          % forget-gate bias = 1
    layers(l).b = dlarray(b);
    cin = numFilters;
end
p.layers = layers;
p.Whead = dlarray(single(randn(3, 3, numFilters, 1) * sqrt(2 / (9 * numFilters))));
p.bhead = dlarray(zeros(1, 1, 'single'));
end
