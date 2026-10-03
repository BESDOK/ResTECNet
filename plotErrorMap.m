function fig = plotErrorMap(statsRef, ds, lat, lon, opts)
%PLOTERRORMAP Fig. 6 of the revised manuscript: spatial distribution of the
%   24-h prediction RMSE of Res-TECNet (hindcast, against the CODE final
%   GIM) over the test period, with the dipole magnetic equator, the +/-15
%   deg magnetic-latitude contours and the GNSS validation stations of
%   Sect. 2.3 overlaid.
%   fig = PLOTERRORMAP(statsRef, ds, lat, lon, 'stations', stationList('all'), ...)
%     'stations' : table with code, lat, lon, regime ([] = only Kayseri)
%     'outFile'  : PNG/JPG path ('' = no export; default 'fig0.jpg')
%   Centered-dipole approximation, IGRF-13 pole epoch 2020 (80.65 N, 72.68 W).

arguments
    statsRef struct
    ds struct
    lat (:,1) double
    lon (:,1) double
    opts.stations = stationList('all')
    opts.outFile (1,:) char = 'fig0.jpg'
    opts.coast (1,1) logical = true
    opts.title (1,:) char = ''
end

E = double(statsRef.Yhat) - double(ds.Y) * ds.Tmax;
rmseMap = sqrt(mean(E.^2, 4, 'omitnan'));
poleLat = 80.65; poleLon = -72.68;
[LON, LAT] = meshgrid(lon, lat);
magLat = asind(sind(LAT) .* sind(poleLat) + cosd(LAT) .* cosd(poleLat) .* cosd(LON - poleLon));

fig = figure('Color', 'w', 'Position', [100 100 940 460]);
imagesc(lon, lat, rmseMap); axis xy tight; colormap(turbo);
clim([0, ceil(prctile(rmseMap(:), 99.5))]);               % RMSE is non-negative
cb = colorbar; cb.Label.String = 'RMSE (TECU)'; cb.Label.FontSize = 11;
hold on
[~, hEq] = contour(LON, LAT, magLat, [0 0], 'LineColor', 'w', 'LineWidth', 1.6);
[~, hPM] = contour(LON, LAT, magLat, [-15 15], 'LineColor', 'w', 'LineWidth', 1.1, 'LineStyle', '--');
if opts.coast
    try, C = load('coastlines'); plot(C.coastlon, C.coastlat, '-', 'Color', [1 1 1]*0.35, 'LineWidth', 0.4); catch, end
end
hS = [];
if ~isempty(opts.stations)
    S = opts.stations;
    mk = struct('EIA', 'o', 'MID', 's', 'HIGH', '^', 'OCEAN', 'd', 'ANATOLIA', 'p');
    regs = unique(S.regime, 'stable'); hS = gobjects(numel(regs), 1);
    for r = 1:numel(regs)
        i = S.regime == regs(r);
        hS(r) = plot(S.lon(i), S.lat(i), mk.(char(regs(r))), 'MarkerFaceColor', 'w', ...
            'MarkerEdgeColor', 'k', 'MarkerSize', 8, 'LineWidth', 0.8, 'LineStyle', 'none');
    end
    % labels: alternate left/right for clustered stations (Anatolia) to reduce overlap
    for i = 1:height(S)
        dx = 4; ha = 'left';
        if S.regime(i) == "ANATOLIA" && mod(i, 2) == 0, dx = -4; ha = 'right'; end
        text(S.lon(i) + dx, S.lat(i) + (S.regime(i) == "ANATOLIA") * (mod(i, 3) - 1) * 3, S.code(i), ...
            'Color', 'w', 'FontSize', 7, 'FontWeight', 'bold', 'HorizontalAlignment', ha);
    end
    legNames = [{'Magnetic equator', '\pm15\circ magnetic latitude'}, cellstr(regs' + " stations")];
    legend([hEq, hPM, hS'], legNames, 'TextColor', 'k', 'Location', 'southwest', 'FontSize', 7, 'NumColumns', 2);
else
    hK = plot(35.49, 38.72, 'kp', 'MarkerFaceColor', 'w', 'MarkerSize', 12);
    legend([hEq, hPM, hK], {'Magnetic equator', '\pm15\circ magnetic latitude', 'Kayseri'}, 'Location', 'southwest');
end
xlabel('Geographic longitude (\circ)'); ylabel('Geographic latitude (\circ)');
xticks(-180:60:180); yticks(-80:20:80); set(gca, 'FontSize', 11, 'Layer', 'top');
if ~isempty(opts.title), title(opts.title); end
if ~isempty(opts.outFile), exportgraphics(fig, opts.outFile, 'Resolution', 300); fprintf('Figure written to %s\n', opts.outFile); end
end
