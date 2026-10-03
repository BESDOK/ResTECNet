%RUNFULLPAPERSIMS End-to-end reproduction of the 6-14 May 2024 window (Sect. 4.4).
%   Retrieves the CODE final maps and the drivers from the public archives,
%   assembles the dataset with the revised pipeline (72 unique meridians,
%   trailing 81-day F10.7 mean, gap screening) and produces Figs. 3-5. If a
%   trained network is available (resTECNet_hindcast.mat / _operational.mat
%   from runExperiments.m) the 24-h forecasts are overlaid on the Kayseri
%   series; otherwise persistence alone is shown. The dual-frequency VTEC of
%   the open IGS station nearest to Kayseri (MERS, ~260 km) is added, with
%   the GIM interpolated to MERS, when station_MERS.mat exists.

buildTECDatabase('2024-05-05', '2024-05-14', 'TECmay2024.mat', 'ionex_cache', 'product', 'CODG');
load('TECmay2024.mat', 'TEC', 'epochs', 'lat', 'lon', 'badMask');
drivers = buildDriverTable(epochs, 'index_cache', 'mode', 'hindcast', 'meanType', 'trailing');
save('TECmay2024.mat', 'drivers', '-append');

fc = struct('name', {}, 'epochs', {}, 'Yhat', {}, 'color', {});
if isfile('resTECNet_hindcast.mat')
    S = load('resTECNet_hindcast.mat', 'net', 'cfg');
    ds = buildTECDataset(TEC, epochs, drivers, 'L', S.cfg.L, 'Delta', S.cfg.Delta, 'Tmax', S.cfg.Tmax, ...
        'zmu', S.cfg.zmu, 'zsig', S.cfg.zsig, 'badMask', badMask);
    st = evaluateTEC(S.net, ds, lat);
    fc(end+1) = struct('name', 'Res-TECNet 24-h forecast (hindcast)', 'epochs', ds.epochsY, 'Yhat', st.Yhat, 'color', [0.85 0.1 0.1]);
    if isfile('resTECNet_operational.mat')
        O = load('resTECNet_operational.mat', 'netOp');
        drvOp = buildDriverTable(epochs, 'index_cache', 'mode', 'operational', 'meanType', 'trailing');
        dsOp = buildTECDataset(TEC, epochs, drvOp, 'L', S.cfg.L, 'Delta', S.cfg.Delta, 'Tmax', S.cfg.Tmax, ...
            'zmu', S.cfg.zmu, 'zsig', S.cfg.zsig, 'mode', 'operational', 'badMask', badMask);
        stOp = evaluateTEC(O.netOp, dsOp, lat);
        fc(end+1) = struct('name', 'Res-TECNet 24-h forecast (operational)', 'epochs', dsOp.epochsY, 'Yhat', stOp.Yhat, 'color', [1 0.55 0]);
    end
end
station = []; if isfile('station_MERS.mat'), L = load('station_MERS.mat'); station = L.obsAll; end
win = epochs >= datetime(2024,5,6,'TimeZone','UTC');
out = visualizeMay2024(TEC(:,:,win), epochs(win), lat, lon, drivers(win,:), 'forecasts', fc, 'station', station, 'stationName', 'MERS', 'stationLat', 36.56, 'stationLon', 34.26);
