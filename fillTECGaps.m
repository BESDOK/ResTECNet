function [TECh, epochsH, badMask] = fillTECGaps(TEC, epochs, lon, opts)
%FILLTECGAPS Regularize to an hourly grid and fill short gaps (Sect. 2.1).
%   [TECh, epochsH, badMask] = FILLTECGAPS(TEC, epochs, lon, 'maxGapHours', 6, 'tecMax', 250)
%   Screening: maps with NaN or nonphysical values (< 0, > tecMax) are
%   discarded; the record is placed on a gap-free hourly grid; gaps
%   <= maxGapHours are filled by sun-fixed linear interpolation; longer gaps
%   stay NaN and are flagged in badMask. Memory-lean: the hourly array is
%   allocated once and filled in place (used for the 26-year record).
arguments
    TEC (:,:,:) {mustBeNumeric}
    epochs (:,1) datetime
    lon (:,1) double
    opts.maxGapHours (1,1) double = 6
    opts.tecMax (1,1) double = 250
end
epochs.TimeZone = 'UTC';
[H, W, K] = size(TEC);
e0 = dateshift(min(epochs), 'start', 'hour'); e1 = dateshift(max(epochs), 'start', 'hour');
epochsH = (e0:hours(1):e1)'; Kh = numel(epochsH);
TECh = nan(H, W, Kh, 'single');
have = false(Kh, 1); nRej = 0;
[tf, loc] = ismember(dateshift(epochs, 'start', 'hour'), epochsH);
for k = 1:K
    if ~tf(k) || have(loc(k)), continue; end
    m = TEC(:, :, k);
    if any(isnan(m), 'all') || any(m < 0, 'all') || any(m > opts.tecMax, 'all'), nRej = nRej + 1; continue; end
    TECh(:, :, loc(k)) = single(m); have(loc(k)) = true;
end
fprintf('fillTECGaps: %d of %d maps rejected by screening (%.3f %%)\n', nRej, K, 100 * nRej / max(K, 1));
[TECh, badMask] = fillGapsInPlace(TECh, have, lon, opts.maxGapHours);
fprintf('fillTECGaps: %d hourly epochs, %d filled, %d left as gaps\n', Kh, nnz(~have) - nnz(badMask), nnz(badMask));
end
