function Y = modelTransformer(p, X, cfg)
%MODELTRANSFORMER Spatiotemporal vision-Transformer baseline (Sect. 3.6, vii).
%   Y = MODELTRANSFORMER(params, X, cfg), X dlarray H x W x C x B ('SSCB').
%   Every input channel of Eq. (3) is divided into patches of size
%   cfg.patch = [ph pw] (default [5 4]); the 71 x 72 grid is zero-padded in
%   latitude to a multiple of ph (75 rows), giving nT = (75/5)*(72/4) = 270
%   tokens of dimension ph*pw*C = 560, linearly embedded to d (256),
%   followed by cfg.numLayers pre-norm encoder layers (multi-head
%   self-attention with cfg.numHeads heads, GELU MLP), a linear
%   de-patching head to ph*pw pixels per token and the same global skip
%   T_t + R as Res-TECNet. Follows the design principle of Shih et al.
%   (2024) with the inputs, period and loss of Res-TECNet.
%   [params, cfg] from initTransformer; use
%   modelFcn = @(p, X) modelTransformer(p, X, cfg).
%
%   Implementation note: apart from the patch-embedding dlconv, all
%   operations are written explicitly on unformatted dlarrays (matrix
%   products, pagemtimes, manual layer normalization / softmax / GELU) so
%   that the forward pass and its gradient do not depend on the dimension
%   semantics of the attention / layernorm / fullyconnect functions.

L = cfg.L; ph = cfg.patch(1); pw = cfg.patch(2); nHeads = cfg.numHeads;
H = size(X, 1); W = size(X, 2); B = size(X, 4);
Tt = X(:, :, L, :);
Hp = ceil(H / ph) * ph;
if Hp > H
    X = cat(1, X, 0 * X(1:Hp - H, :, :, :));          % zero rows, same format/device
end
nh = Hp / ph; nw = W / pw; nT = nh * nw;

% ---- patch embedding: conv with kernel = stride = patch ----
E = dlconv(X, p.Wpatch, p.bpatch, 'Stride', [ph pw]);  % nh x nw x d x B ('SSCB')
E = stripdims(E);
d = size(E, 3);
Z = reshape(E, [nT, d, B]);                            % token = i + nh*(j-1)
Z = permute(Z, [2 1 3]);                               % d x nT x B (unformatted)
Z = Z + p.pos;                                         % d x nT (learnable)

% ---- encoder ----
for l = 1:cfg.numLayers
    q = p.layers(l);
    Zn = layerNorm(Z, q.ln1g, q.ln1b);
    Q = linear(Zn, q.Wq, q.bq); K = linear(Zn, q.Wk, q.bk); V = linear(Zn, q.Wv, q.bv);
    A = multiHeadAttention(Q, K, V, nHeads);
    Z = Z + linear(A, q.Wo, q.bo);
    Zn = layerNorm(Z, q.ln2g, q.ln2b);
    Z = Z + linear(gelu_(linear(Zn, q.W1, q.b1)), q.W2, q.b2);
end
Z = layerNorm(Z, p.lnfg, p.lnfb);

% ---- de-patching head ----
P = linear(Z, p.Wout, p.bout);                         % (ph*pw) x nT x B
R = reshape(P, [ph, pw, nh, nw, B]);                   % (pi,pj,i,j,b)
R = permute(R, [1 3 2 4 5]);                           % (pi,i,pj,j,b)
R = reshape(R, [Hp, W, 1, B]);
R = R(1:H, :, :, :);
Y = Tt + dlarray(R, 'SSCB');
end

% -------------------------------------------------------------------------
function Y = linear(Z, Wm, b)
% Z: din x nT x B  ->  dout x nT x B
din = size(Z, 1); nT = size(Z, 2); B = size(Z, 3);
Y = Wm * reshape(Z, [din, nT * B]) + b;
Y = reshape(Y, [size(Wm, 1), nT, B]);
end

function Zn = layerNorm(Z, g, b)
mu = mean(Z, 1);
v  = mean((Z - mu).^2, 1);
Zn = (Z - mu) ./ sqrt(v + 1e-5) .* g + b;
end

function A = multiHeadAttention(Q, K, V, nHeads)
% Q,K,V: d x nT x B ; scaled dot-product attention per head
d = size(Q, 1); nT = size(Q, 2); B = size(Q, 3);
dh = d / nHeads;
split = @(M) reshape(permute(reshape(M, [dh, nHeads, nT, B]), [1 3 2 4]), [dh, nT, nHeads * B]);
Qh = split(Q); Kh = split(K); Vh = split(V);
S = pagemtimes(Kh, 'transpose', Qh, 'none') / sqrt(dh);   % nT(keys) x nT(queries) x (heads*B)
S = S - max(S, [], 1);
Ex = exp(S);
Sm = Ex ./ sum(Ex, 1);                                     % softmax over keys
Ah = pagemtimes(Vh, Sm);                                   % dh x nT x (heads*B)
A = reshape(permute(reshape(Ah, [dh, nT, nHeads, B]), [1 3 2 4]), [d, nT, B]);
end

function Y = gelu_(X)
Y = 0.5 * X .* (1 + tanh(sqrt(2 / pi) * (X + 0.044715 * X.^3)));
end
