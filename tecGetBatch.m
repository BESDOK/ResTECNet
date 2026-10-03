function [X, Y, Tlast] = tecGetBatch(ds, idx)
%TECGETBATCH Assemble the input tensor of Eq. (3) for a set of samples.
%   [X, Y, Tlast] = TECGETBATCH(ds, idx)
%     ds  : struct from buildTECDataset
%     idx : sample indices (1..ds.N)
%   X : H x W x C x B single, channels 1..L = T_{t-L+1..t}, L+1..L+8 = z_t,
%       L+9..L+16 = z_{t+Delta}, [L+17 = SH-AR base for t+Delta when ds.Base
%       is set]; Y, Tlast : H x W x 1 x B (scaled).
idx = idx(:)';
B = numel(idx); H = ds.H; W = ds.W; L = ds.L;
X = zeros(H, W, ds.C, B, 'single');
for b = 1:B
    t = ds.idxT(idx(b));
    X(:, :, 1:L, b) = ds.Ts(:, :, t-L+1:t);
    X(:, :, L+1:L+8, b) = repmat(reshape(ds.Zt(idx(b), :), 1, 1, 8), H, W, 1);
    X(:, :, L+9:L+16, b) = repmat(reshape(ds.Ztd(idx(b), :), 1, 1, 8), H, W, 1);
    if ds.C > L + 16, X(:, :, L+17, b) = ds.Base(:, :, t + ds.Delta); end   % adaptive base (SH-AR) channel
end
if nargout > 1
    tt = ds.idxT(idx);
    Y = reshape(ds.Ts(:, :, tt + ds.Delta), H, W, 1, B);        % read from the shared map array
end
if nargout > 2, Tlast = reshape(ds.Ts(:, :, tt), H, W, 1, B); end
end
