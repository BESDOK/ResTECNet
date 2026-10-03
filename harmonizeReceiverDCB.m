function [obs, info] = harmonizeReceiverDCB(obs, opts)
%HARMONIZERECEIVERDCB Replace daily receiver-DCB estimates by a running median.
%   [obs, info] = HARMONIZERECEIVERDCB(obs, 'halfWindowDays', 15)
%   The receiver DCB is a hardware property that changes slowly (days to
%   weeks), whereas the day-by-day estimates from the minimum-scatter
%   thin-shell method fluctuate by about 1 ns (3 TECU) on quiet days and
%   fail on geomagnetically disturbed days (the night-time ionosphere is
%   then not homogeneous). Every station-day's estimate d_k (obs.dRec,
%   TECU) is therefore replaced by the median d_ref of the valid estimates
%   within +/- halfWindowDays (at least 'minDays' days):
%       STEC' = STEC + (d_k - d_ref),  VTEC' = VTEC + (d_k - d_ref) / mf.
%   A day for which no estimate was possible (too few night-time
%   observations: processStationVTEC then applies d_k = 0 exactly) is
%   excluded from the reference and corrected with its neighbours' median;
%   if no reference exists the observations of the day are flagged
%   ok = false. info: nDays, nNoEstimate, dailyScatterTECU (RMS of
%   d_k - d_ref over valid days), maxAbsShiftTECU and the median receiver
%   DCB in ns. The absolute level of d_ref remains uncertain; it is handled
%   by the station-offset removal of evaluateStationVTEC.
arguments
    obs timetable
    opts.halfWindowDays (1,1) double = 15
    opts.minDays (1,1) double = 5
end
info = struct('nDays', 0, 'nNoEstimate', 0, 'dailyScatterTECU', NaN, 'maxAbsShiftTECU', NaN, 'medianDCB_ns', NaN);
if ~ismember('dRec', obs.Properties.VariableNames) || height(obs) == 0
    warning('harmonizeReceiverDCB:none', 'No dRec variable: nothing to harmonize.'); return
end
t = obs.Properties.RowTimes; day = dateshift(t, 'start', 'day');
[ud, ~, gi] = unique(day);
dDay = accumarray(gi, obs.dRec, [], @(x) x(1));
okDay = dDay ~= 0 & ~isnan(dDay);                              % 0 = no estimate possible
ref = nan(size(dDay));
for k = 1:numel(ud)
    win = abs(days(ud - ud(k))) <= opts.halfWindowDays & okDay;
    if nnz(win) >= opts.minDays, ref(k) = median(dDay(win)); end
end
own = isnan(ref) & okDay; ref(own) = dDay(own);                % sparse neighbourhood: keep the day's own estimate
shift = dDay - ref;                                            % NaN where no reference exists
sh = shift(gi); bad = isnan(sh); sh(bad) = 0;
obs.STEC = obs.STEC + sh;
obs.VTEC = obs.VTEC + sh ./ obs.mf;
obs.dRec = obs.dRec - sh;
if any(bad) && ismember('ok', obs.Properties.VariableNames), obs.ok(bad) = false; end
Kion = 40.3e16 * (1575.42e6^2 - 1227.60e6^2) / (1575.42e6^2 * 1227.60e6^2);
v = okDay & ~isnan(shift);
info.nDays = numel(ud); info.nNoEstimate = nnz(~okDay);
info.dailyScatterTECU = sqrt(mean(shift(v).^2)); info.maxAbsShiftTECU = max(abs(shift(v)));
info.medianDCB_ns = median(ref, 'omitnan') * Kion / (299792458 * 1e-9);
end
