function [T, stormMask, quietMask, events] = stormAnalysis(statsList, names, epochsY, Dst, thresholdDst, winHours)
%STORMANALYSIS Storm-time RMSE over +/- windows around Dst minima (Table 2).
%   T = STORMANALYSIS(statsList, names, epochsY, Dst, thresholdDst, winHours)
%     statsList : cell array of stats structs from evaluateTEC
%                 (must contain .perEpochRMSE aligned with epochsY)
%     names     : cell array of model names (same length)
%     epochsY   : N x 1 datetime of target epochs
%     Dst       : N x 1 hourly Dst [nT] aligned with epochsY
%     thresholdDst : storm threshold, e.g. -100
%     winHours  : half-window, e.g. 48
%
%   Returns a table with one row per model: mean storm-time RMSE,
%   quiet-time (Dst > -30 nT, outside storm windows) RMSE, and storm count,
%   plus the storm/quiet masks and the event indices for reuse
%   (evaluateStationVTEC, calibrationStats, driverAblation).

arguments
    statsList cell
    names cell
    epochsY (:,1) datetime
    Dst (:,1) double
    thresholdDst (1,1) double = -100
    winHours (1,1) double = 48
end

% --- identify storm events: local Dst minima below threshold,
%     separated by at least 3 days ---
below = find(Dst <= thresholdDst);
events = [];
while ~isempty(below)
    [~, k] = min(Dst(below));
    center = below(k);
    events(end+1) = center; %#ok<AGROW>
    below(abs(hours(epochsY(below) - epochsY(center))) < 72) = [];
end
events = sort(events);

stormMask = false(numel(epochsY), 1);
for e = events
    stormMask = stormMask | ...
        abs(hours(epochsY - epochsY(e))) <= winHours;
end
quietMask = Dst > -30 & ~stormMask;

M = numel(statsList);
stormRMSE = zeros(M, 1); quietRMSE = zeros(M, 1);
for m = 1:M
    r2 = statsList{m}.perEpochRMSE.^2;
    stormRMSE(m) = sqrt(mean(r2(stormMask), 'omitnan'));
    quietRMSE(m) = sqrt(mean(r2(quietMask), 'omitnan'));   % NaN epochs (e.g. missing C1PG maps) ignored
end

T = table(names(:), stormRMSE, quietRMSE, ...
    repmat(numel(events), M, 1), ...
    'VariableNames', {'Model', 'StormRMSE_TECU', 'QuietRMSE_TECU', 'NumStorms'});

fprintf('Detected %d storm events (Dst <= %g nT):\n', numel(events), thresholdDst);
disp(epochsY(events));
end
