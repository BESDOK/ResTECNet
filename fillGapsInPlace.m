function [TECh, badMask] = fillGapsInPlace(TECh, have, lon, maxGapHours)
%FILLGAPSINPLACE Sun-fixed linear filling of short gaps in an hourly map array.
%   [TECh, badMask] = FILLGAPSINPLACE(TECh, have, lon, maxGapHours)
%     TECh : H x W x Kh single, NaN where 'have' is false
%     have : Kh x 1 logical (map present)
%   Gaps of at most maxGapHours consecutive hours enclosed by present maps
%   are filled; others are flagged in badMask. Operates in place (MATLAB
%   copy-on-write: call as TECh = fillGapsInPlace(TECh, ...)).
Kh = numel(have); dlon = lon(2) - lon(1);
badMask = false(Kh, 1);
missing = find(~have);
if isempty(missing), return; end
brk = [true; diff(missing) > 1];
starts = missing(brk); ends = missing([brk(2:end); true]);
for g = 1:numel(starts)
    a = starts(g) - 1; b = ends(g) + 1; len = ends(g) - starts(g) + 1;
    if a < 1 || b > Kh || len > maxGapHours
        badMask(starts(g):ends(g)) = true; continue
    end
    A0 = TECh(:, :, a); B0 = TECh(:, :, b);
    for n = starts(g):ends(g)
        fa = (n - a) / (b - a);
        A = circshift(A0, -round(15 * (n - a) / dlon), 2);
        B = circshift(B0,  round(15 * (b - n) / dlon), 2);
        TECh(:, :, n) = (1 - fa) * A + fa * B;
    end
end
end
