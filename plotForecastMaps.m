function fig = plotForecastMaps(statsRes, ds, lat, lon, epochsShow, opts)
%PLOTFORECASTMAPS Fig. 2 of the revised manuscript: reference CODE map at the
%   target epoch, Res-TECNet 24-h forecast and forecast minus reference, for
%   a quiet epoch and for the May 2024 storm main phase. Colour bars are
%   labelled with quantity and unit (reply to the DOCX review, A2.41).
%   fig = PLOTFORECASTMAPS(statsRes, ds, lat, lon, epochsShow, 'outFile', 'fig4.jpg')
%     epochsShow : datetime vector of target epochs (default 8 May 2024 18 UT
%                  and 10 May 2024 18 UT)

arguments
    statsRes struct
    ds struct
    lat (:,1) double
    lon (:,1) double
    epochsShow (:,1) datetime = [datetime(2024,5,8,18,0,0,'TimeZone','UTC'); datetime(2024,5,10,18,0,0,'TimeZone','UTC')]
    opts.outFile (1,:) char = 'fig4.jpg'
    opts.labels cell = {'quiet, 8 May 2024 18 UT', 'storm, 10 May 2024 18 UT'}
end
epochsShow.TimeZone = 'UTC';
[LON, LAT] = meshgrid(lon, lat);
magLat = asind(sind(LAT) .* sind(80.65) + cosd(LAT) .* cosd(80.65) .* cosd(LON + 72.68));
nR = numel(epochsShow);
fig = figure('Color', 'w', 'Position', [60 60 1250 330 * nR]);
tl = tiledlayout(nR, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for r = 1:nR
    [tf, n] = ismember(epochsShow(r), ds.epochsY);
    if ~tf, error('plotForecastMaps:epoch', 'Epoch %s not in ds.epochsY', datestr(epochsShow(r))); end
    ref = double(ds.Y(:, :, 1, n)) * ds.Tmax;
    fc  = double(statsRes.Yhat(:, :, 1, n));
    cmax = max([ref(:); fc(:)]);
    dmax = max(abs(fc(:) - ref(:)));
    panel(ref, [0 cmax], turbo, sprintf('(%s) CODE GIM (%s)', char('a' + 3*(r-1)), opts.labels{r}), 'Vertical TEC (TECU)');
    panel(fc,  [0 cmax], turbo, sprintf('(%s) Res-TECNet 24-h forecast', char('b' + 3*(r-1))), 'Vertical TEC (TECU)');
    panel(fc - ref, dmax * [-1 1], parula, sprintf('(%s) Forecast - reference (RMSE %.2f TECU)', ...
        char('c' + 3*(r-1)), sqrt(mean((fc(:) - ref(:)).^2))), 'TEC difference (TECU)');
end
xlabel(tl, 'Geographic longitude (\circ)'); ylabel(tl, 'Geographic latitude (\circ)');
if ~isempty(opts.outFile), exportgraphics(fig, opts.outFile, 'Resolution', 300); end

    function panel(M, cl, cmap, ttl, cbLabel)
        nexttile; imagesc(lon, lat, M, cl); axis xy tight; colormap(gca, cmap);
        cb = colorbar; cb.Label.String = cbLabel; cb.Label.FontSize = 9;
        hold on; contour(LON, LAT, magLat, [0 0], 'w-', 'LineWidth', 1);
        title(ttl, 'FontSize', 9); set(gca, 'FontSize', 8, 'Layer', 'top');
    end
end
