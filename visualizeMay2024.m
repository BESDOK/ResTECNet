function out = visualizeMay2024(TEC, epochs, lat, lon, drivers, opts)
%VISUALIZEMAY2024 Figures 3-5 of the revised manuscript for the May 2024 superstorm.
%   out = VISUALIZEMAY2024(TEC, epochs, lat, lon, drivers, ...)
%   Options
%     'forecasts'  struct array with fields name, epochs (target epochs),
%                  Yhat (H x W x 1 x N, TECU), color -- e.g. the Res-TECNet
%                  hindcast and operational forecasts (Fig. 5 overlays;
%                  reply to the PDF review, comment 1)
%     'station'    timetable from processStationVTEC of the open station
%                  nearest to Kayseri (MERS by default; grey dots in Fig. 5
%                  = hourly means of ok observations with el>=30) drawn
%                  together with the GIM interpolated to that station, so
%                  that the station-GIM agreement is location-consistent
%     'stationName','stationLat','stationLon'  (default MERS, 36.56, 34.26)
%     'persistence' logical: also draw the 24-h persistence series
%     'outPrefix'  file prefix (default '') -> fig1.jpg, fig2.jpg, fig3.jpg
%   Produces fig1.jpg (drivers), fig2.jpg (maps, colour bars labelled),
%   fig3.jpg (Kayseri series with forecasts and station) and prints the
%   key numbers quoted in Sect. 4.4, including the forecast RMSE over the
%   window and the station-GIM agreement. Returns them in 'out'.

arguments
    TEC (:,:,:) {mustBeNumeric}
    epochs (:,1) datetime
    lat (:,1) double
    lon (:,1) double
    drivers table
    opts.forecasts struct = struct('name', {}, 'epochs', {}, 'Yhat', {}, 'color', {})
    opts.station = []
    opts.stationName (1,:) char = 'MERS'
    opts.stationLat (1,1) double = 36.56
    opts.stationLon (1,1) double = 34.26
    opts.persistence (1,1) logical = true
    opts.outPrefix (1,:) char = ''
end
epochs.TimeZone = 'UTC';
kayLat = 38.72; kayLon = 35.49;
tQuiet = datetime(2024,5,8,18,0,0,'TimeZone','UTC');
tMain  = datetime(2024,5,10,18,0,0,'TimeZone','UTC');
tRecov = datetime(2024,5,12,18,0,0,'TimeZone','UTC');
kQ = find(epochs == tQuiet); kM = find(epochs == tMain); kR = find(epochs == tRecov);
stormWin = [datetime(2024,5,10,17,0,0,'TimeZone','UTC'), datetime(2024,5,12,6,0,0,'TimeZone','UTC')];

%% ---------------- Fig. 3 (file fig1): drivers ----------------
fA = figure('Color','w','Position',[80 80 860 420]);
tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile; area(epochs, drivers.Dst, 'FaceColor',[0.85 0.33 0.10], 'FaceAlpha',0.25, ...
    'EdgeColor',[0.85 0.33 0.10], 'LineWidth',1.2);
hold on; yline(-100,'--k','Dst = -100 nT','LabelHorizontalAlignment','left'); shadeStorm(stormWin);
ylabel('Dst (nT)'); grid on; set(gca,'FontSize',11); title('(a) Storm-time disturbance index');
nexttile; stairs(epochs, drivers.Kp, 'Color',[0 0.35 0.65], 'LineWidth',1.4); hold on; shadeStorm(stormWin);
ylabel('Kp'); ylim([0 9.5]); grid on; set(gca,'FontSize',11); title('(b) Planetary geomagnetic index');
xlabel(tl, 'Universal time, 6-14 May 2024', 'FontSize', 11);
exportgraphics(fA, [opts.outPrefix 'fig1.jpg'], 'Resolution', 300);

%% ---------------- Fig. 4 (file fig2): maps with labelled colour bars ----------------
fB = figure('Color','w','Position',[60 60 980 640]);
tlB = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
cmax = max(TEC(:,:,kM), [], 'all');
dQ = TEC(:,:,kM) - TEC(:,:,kQ); dR = TEC(:,:,kR) - TEC(:,:,kQ);
drawMap(TEC(:,:,kQ), [0 cmax], turbo, sprintf('(a) Quiet reference: %s UT', datestr(tQuiet,'dd mmm HH:MM')), 'Vertical TEC (TECU)');
drawMap(TEC(:,:,kM), [0 cmax], turbo, sprintf('(b) Main phase: %s UT', datestr(tMain,'dd mmm HH:MM')), 'Vertical TEC (TECU)');
drawMap(dQ, max(abs(dQ),[],'all')*[-1 1], parula, '(c) Main phase - quiet (same UT)', 'TEC difference (TECU)');
drawMap(dR, max(abs(dR),[],'all')*[-1 1], parula, '(d) Recovery - quiet (same UT)', 'TEC difference (TECU)');
title(tlB, 'CODE GIM vertical TEC, May 2024 superstorm', 'FontSize', 12);
exportgraphics(fB, [opts.outPrefix 'fig2.jpg'], 'Resolution', 300);

%% ---------------- Fig. 5 (file fig3): Kayseri series + forecasts + station ----------------
tecKay = interpKay(TEC, lat, lon, kayLat, kayLon);
refMask = epochs < datetime(2024,5,9,'TimeZone','UTC');
refCurve = accumarray(hour(epochs(refMask))+1, tecKay(refMask), [24 1], @mean);
tecRef = refCurve(hour(epochs)+1);

fC = figure('Color','w','Position',[80 80 900 420]);
yyaxis left; hold on
hG = plot(epochs, tecKay, 'k-', 'LineWidth', 1.8);
hRef = plot(epochs, tecRef, 'k--', 'LineWidth', 1.0);
hh = [hG, hRef]; leg = {'CODE GIM, Kayseri', 'Quiet-day diurnal reference (6-8 May)'};
out.kayseriRMSE = table('Size', [0 2], 'VariableTypes', {'string','double'}, 'VariableNames', {'Series','RMSE_TECU'});
if opts.persistence
    persist = nan(size(tecKay)); [tf, loc] = ismember(epochs - hours(24), epochs); persist(tf) = tecKay(loc(tf));
    hP = plot(epochs, persist, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2); hh(end+1) = hP; leg{end+1} = '24-h persistence';
    out.kayseriRMSE(end+1, :) = {"Persistence", sqrt(mean((persist - tecKay).^2, 'omitnan'))};
end
for f = 1:numel(opts.forecasts)
    F = opts.forecasts(f); F.epochs.TimeZone = 'UTC';
    fk = interpKay(reshape(F.Yhat, size(F.Yhat,1), size(F.Yhat,2), []), lat, lon, kayLat, kayLon);
    [tf, loc] = ismember(epochs, F.epochs);
    ser = nan(size(tecKay)); ser(tf) = fk(loc(tf));
    hF = plot(epochs, ser, '-', 'Color', F.color, 'LineWidth', 1.4); hh(end+1) = hF; leg{end+1} = F.name;
    e = ser - tecKay; out.kayseriRMSE(end+1, :) = {string(F.name), sqrt(mean(e.^2, 'omitnan'))};
    [pk, ipk] = max(abs(e)); fprintf('%s: RMSE %.2f TECU over window, peak error %.1f TECU at %s UT\n', ...
        F.name, out.kayseriRMSE.RMSE_TECU(end), pk, datestr(epochs(ipk), 'dd mmm HH:MM'));
end
if ~isempty(opts.station)
    O = opts.station; O = O(O.ok & O.el >= 30, :); tt = O.Properties.RowTimes; tt.TimeZone = 'UTC';
    st = nan(size(tecKay));
    for n = 1:numel(epochs)
        s = abs(minutes(tt - epochs(n))) <= 30; if nnz(s) >= 3, st(n) = mean(O.VTEC(s)); end
    end
    gimSt = interpKay(TEC, lat, lon, opts.stationLat, opts.stationLon);
    hGs = plot(epochs, gimSt, '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.9);
    hSt = plot(epochs, st, 'o', 'Color', [0.5 0.5 0.5], 'MarkerSize', 4, 'MarkerFaceColor', [0.6 0.6 0.6]);
    hh(end+1) = hGs; leg{end+1} = sprintf('CODE GIM, %s', opts.stationName);
    hh(end+1) = hSt; leg{end+1} = sprintf('%s dual-frequency VTEC', opts.stationName);
    out.stationGIM_RMS = sqrt(mean((st - gimSt).^2, 'omitnan'));
    fprintf('%s station vs CODE GIM at %s: %.2f TECU RMS over window\n', opts.stationName, opts.stationName, out.stationGIM_RMS);
end
shadeStorm(stormWin);
ylabel('VTEC over Kayseri (TECU)');
yyaxis right; hD = plot(epochs, drivers.Dst, '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.0); ylabel('Dst (nT)');
hh(end+1) = hD; leg{end+1} = 'Dst';
grid on; set(gca,'FontSize',10); legend(hh, leg, 'Location', 'northoutside', 'NumColumns', 3, 'FontSize', 8);
title(sprintf('Vertical TEC over Kayseri (%.2f\\circN, %.2f\\circE), 6-14 May 2024', kayLat, kayLon));
xlabel('Universal time');
exportgraphics(fC, [opts.outPrefix 'fig3.jpg'], 'Resolution', 300);

%% ---------------- key numbers (Sect. 4.4) ----------------
[out.kayMax, im] = max(tecKay); out.kayMaxTime = epochs(im); out.kayRefAtMax = tecRef(im);
out.enhancementPct = 100 * (out.kayMax / out.kayRefAtMax - 1);
d12 = day(epochs) == 12 & hour(epochs) >= 9 & hour(epochs) <= 15;
out.depletion12 = mean(tecRef(d12) - tecKay(d12)); out.depletion12Pct = 100 * mean(1 - tecKay(d12) ./ tecRef(d12));
[out.dstMin, id] = min(drivers.Dst); out.dstMinTime = epochs(id); out.kpMax = max(drivers.Kp);
out.F107range = [min(drivers.F107) max(drivers.F107)];
fprintf(['Dst min %d nT at %s UT | Kp max %.1f | F10.7 %.0f-%.0f sfu\n' ...
    'Kayseri max %.1f TECU at %s UT (ref %.1f, +%.0f %%) | 12 May depletion %.1f TECU (%.0f %%)\n'], ...
    round(out.dstMin), datestr(out.dstMinTime,'dd mmm HH:MM'), out.kpMax, out.F107range, ...
    out.kayMax, datestr(out.kayMaxTime,'dd mmm HH:MM'), out.kayRefAtMax, out.enhancementPct, out.depletion12, out.depletion12Pct);
disp(out.kayseriRMSE);
% daytime (09-15 UT) deviation from the quiet-day reference, per day (Sect. 4.4)
dl = unique(day(epochs)); dev = zeros(numel(dl), 1); pct = dev; pk = dev; pkUT = dev;
for i = 1:numel(dl)
    m = day(epochs) == dl(i) & hour(epochs) >= 9 & hour(epochs) <= 15;
    dev(i) = mean(tecKay(m) - tecRef(m)); pct(i) = 100 * mean(tecKay(m) ./ tecRef(m) - 1);
    [~, j] = max(tecKay(day(epochs) == dl(i))); e = epochs(day(epochs) == dl(i)); pkUT(i) = hour(e(j)); pk(i) = max(tecKay(day(epochs) == dl(i)));
end
out.dailyDaytime = table(dl, dev, pct, pk, pkUT, 'VariableNames', {'Day','DaytimeDev_TECU','DaytimeDev_pct','DailyMax_TECU','MaxUT'});
fprintf('Daytime (09-15 UT) deviation from quiet reference:\n'); disp(out.dailyDaytime);

%% ---------------- helpers ----------------
    function drawMap(M, cl, cmap, ttl, cbLabel)
        nexttile; imagesc(lon, lat, M, cl); axis xy tight; colormap(gca, cmap);
        cb = colorbar; cb.Label.String = cbLabel;
        hold on; plot(kayLon, kayLat, 'kp', 'MarkerFaceColor','w', 'MarkerSize',10);
        text(kayLon+6, kayLat, 'Kayseri', 'FontSize', 8, 'FontWeight','bold');
        xlabel('Longitude (\circ)'); ylabel('Latitude (\circ)'); title(ttl, 'FontSize', 10); set(gca,'FontSize',9,'Layer','top');
    end
end

function shadeStorm(w)
yl = ylim;
patch([w(1) w(2) w(2) w(1)], [yl(1) yl(1) yl(2) yl(2)], [0.85 0.33 0.10], 'FaceAlpha', 0.08, 'EdgeColor', 'none');
end

function v = interpKay(T3, lat, lon, sLat, sLon)
% bilinear interpolation of an H x W x K stack to a point (circular lon)
[H, W, K] = size(T3);
dlon = lon(2) - lon(1);
fj = (sLon - lon(1)) / dlon; j1 = floor(fj); b = fj - j1; j1 = mod(j1, W) + 1; j2 = mod(j1, W) + 1;
fi = (sLat - lat(1)) / (lat(2) - lat(1)); i1 = floor(fi); a = fi - i1; i1 = min(max(i1 + 1, 1), H - 1); i2 = i1 + 1;
v = squeeze((1-a)*(1-b)*T3(i1,j1,:) + (1-a)*b*T3(i1,j2,:) + a*(1-b)*T3(i2,j1,:) + a*b*T3(i2,j2,:));
v = v(:);
end
