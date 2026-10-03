function T = regionalGradients(hourly, pairs)
%REGIONALGRADIENTS Observed vs predicted VTEC differences between station pairs (Sect. 4.7).
%   T = REGIONALGRADIENTS(hourly, pairs)
%     hourly : struct array from evaluateStationVTEC (fields code, obs,
%              model, names)
%     pairs  : cell array of {codeA, codeB, label}, e.g.
%              {'VAAN','ISTA','East-West (VAAN - ISTA)';
%               'ANKR','ADAN','North-South (ANKR - ADAN)'}
%   For each pair the difference series obs(A) - obs(B) is compared with
%   model(A) - model(B) for every model: standard deviation of the observed
%   difference, correlation and RMSE of the predicted difference.

arguments
    hourly struct
    pairs cell
end
codes = string({hourly.code});
rows = {};
for p = 1:size(pairs, 1)
    a = find(codes == pairs{p, 1}, 1); b = find(codes == pairs{p, 2}, 1);
    if isempty(a) || isempty(b), warning('regionalGradients:pair', 'Pair %s missing.', pairs{p, 3}); continue; end
    dObs = hourly(a).obs - hourly(b).obs;
    valid = ~isnan(dObs);
    for m = 1:numel(hourly(a).names)
        dMod = hourly(a).model(:, m) - hourly(b).model(:, m);
        v = valid & ~isnan(dMod);
        cc = corrcoef(dObs(v), dMod(v));
        rows(end+1, :) = {pairs{p, 3}, hourly(a).names{m}, std(dObs(v)), cc(1, 2), ...
            sqrt(mean((dMod(v) - dObs(v)).^2))}; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'Pair', 'Model', 'ObsDiffStd_TECU', 'Corr', 'RMSE_TECU'});
disp(T);
end
