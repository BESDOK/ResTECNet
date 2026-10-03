function ok = compareRinexReaders(file)
%COMPARERINEXREADERS Gate: fast GPS reader vs rinexread on a real RINEX-3 file.
%   ok = COMPARERINEXREADERS(file) reads the file with readRinexObsGPS and
%   with rinexread, matches the records by (time, satellite) and compares
%   the observables used by the processing chain. ok = true only if the
%   record sets are identical and all values agree to 1e-6 (equal NaN
%   pattern). Run once on a real daily file before the station run.
tic; Gf = readRinexObsGPS(file); tf = toc;
tic; D = rinexread(file); tr = toc; Gr = D.GPS;
fprintf('fast reader: %d records (%.1f s) | rinexread: %d GPS records (%.1f s) | speed-up x%.1f\n', height(Gf), tf, height(Gr), tr, tr / tf);
kf = [datenum(Gf.Time), double(Gf.SatelliteID)]; kr = [datenum(Gr.Time), double(Gr.SatelliteID)];
[kf, of] = sortrows(kf); [kr, or] = sortrows(kr);
ok = isequal(size(kf), size(kr)) && max(abs(kf(:) - kr(:))) < 1e-9;
if ~ok
    fprintf('FAIL: record sets differ (%d vs %d records)\n', size(kf, 1), size(kr, 1)); return
end
vars = intersect(Gf.Properties.VariableNames, Gr.Properties.VariableNames);
vars = vars(~ismember(vars, {'SatelliteID'}));
for j = 1:numel(vars)
    a = Gf.(vars{j})(of); b = double(Gr.(vars{j})(or));
    sameNaN = isequal(isnan(a), isnan(b));
    d = max(abs(a(~isnan(a) & ~isnan(b)) - b(~isnan(a) & ~isnan(b))));
    if isempty(d), d = 0; end
    pass = sameNaN && d < 1e-6;
    fprintf('  %-4s max |diff| %.2e, same NaN pattern %d -> %s\n', vars{j}, d, sameNaN, tern(pass, 'OK', 'FAIL'));
    ok = ok && pass;
end
fprintf('compareRinexReaders: %s\n', tern(ok, 'PASS', 'FAIL'));
end
function s = tern(c, a, b), if c, s = a; else, s = b; end, end
