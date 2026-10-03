function [T, regional] = checkStationAvailability(opts)
%CHECKSTATIONAVAILABILITY Which stations have 30-s RINEX-3 files on the open archives?
%   [T, regional] = CHECKSTATIONAVAILABILITY(days, netrc)
%   Reads the directory listings of CDDIS (yyd/ and yyo/) and BKG (IGS and
%   EUREF obs) for a few sample days and reports, for every station of
%   stationList('all'), on how many of the sample days a 30-s daily file
%   exists and on which source. 'regional' lists all stations found in the
%   listings with a country code of the Anatolia / eastern Mediterranean /
%   Black Sea / Middle East region, to choose replacements.
%   Default days: 24 Mar 2023, 5 Nov 2023, 10 May 2024, 11 Oct 2024 (storm days).
arguments
    opts.days (:,1) datetime = [datetime(2023,3,24); datetime(2023,11,5); datetime(2024,5,10); datetime(2024,10,11)]
    opts.netrc (1,:) char = fullfile(pwd, 'earthdata_netrc')
    opts.extra (:,1) string = strings(0, 1)      % extra candidate station codes to check
end
days = opts.days; netrc = opts.netrc;
S = stationList('all');
pat = '([A-Z0-9]{4})00([A-Z]{3})_R_\d{11}_01D_30S_MO\.(?:crx|rnx)\.gz';
found = containers.Map('KeyType', 'char', 'ValueType', 'any');   % 'STAT|dayidx' -> source
allNames = strings(0, 2);
for k = 1:numel(days)
    d = days(k); d.TimeZone = 'UTC'; yr = year(d); ddd = day(d, 'dayofyear'); yy = mod(yr, 100);
    src = {sprintf('https://cddis.nasa.gov/archive/gnss/data/daily/%d/%03d/%02dd/', yr, ddd, yy), 'CDDIS-d', true; ...
           sprintf('https://cddis.nasa.gov/archive/gnss/data/daily/%d/%03d/%02do/', yr, ddd, yy), 'CDDIS-o', true; ...
           sprintf('https://igs.bkg.bund.de/root_ftp/EUREF/obs/%d/%03d/', yr, ddd), 'BKG-EUREF', false; ...
           sprintf('https://igs.bkg.bund.de/root_ftp/IGS/obs/%d/%03d/', yr, ddd), 'BKG-IGS', false};
    for s = 1:size(src, 1)
        txt = listing(src{s, 1}, netrc, src{s, 3});
        tok = regexp(txt, pat, 'tokens');
        for i = 1:numel(tok)
            key = sprintf('%s|%d', tok{i}{1}, k);
            if ~isKey(found, key), found(key) = src{s, 2}; else, found(key) = [found(key) '+' src{s, 2}]; end
            allNames(end+1, :) = [string(tok{i}{1}), string(tok{i}{2})]; %#ok<AGROW>
        end
        fprintf('%s  %-10s: %d station files listed\n', datestr(d, 'yyyy-mm-dd'), src{s, 2}, numel(tok));
    end
end
allNames = unique(allNames, 'rows');
if ~isempty(opts.extra)                                   % append candidates; ISO code taken from the listings
    ex = table(upper(opts.extra), repmat("CANDIDATE", numel(opts.extra), 1), nan(numel(opts.extra), 1), nan(numel(opts.extra), 1), ...
        repmat("IGS", numel(opts.extra), 1), strings(numel(opts.extra), 1), 'VariableNames', S.Properties.VariableNames);
    for i = 1:height(ex), r = find(allNames(:, 1) == ex.code(i), 1); if ~isempty(r), ex.iso(i) = allNames(r, 2); end, end
    S = [S; ex];
end
nD = numel(days); avail = zeros(height(S), 1); where = strings(height(S), 1);
for i = 1:height(S)
    for k = 1:nD
        key = sprintf('%s|%d', S.code(i), k);
        if isKey(found, key), avail(i) = avail(i) + 1; where(i) = string(found(key)); end
    end
end
T = table(S.code, S.regime, S.iso, avail, repmat(nD, height(S), 1), where, ...
    'VariableNames', {'code', 'regime', 'iso', 'daysAvailable', 'daysChecked', 'lastSource'});
disp(T);
reg = ismember(allNames(:, 2), ["TUR","CYP","ISR","GRC","ARM","GEO","BGR","ROU","AZE","JOR","EGY","SAU","IRN","UKR","RUS","SRB","ALB","ITA","MLT","LBN","IRQ"]);
regional = array2table(allNames(reg, :), 'VariableNames', {'code', 'iso'});
fprintf('\nStations of the region present in the listings (%d):\n', height(regional)); disp(regional);
end

function txt = listing(url, netrc, useNetrc)
f = [tempname '.html'];
if useNetrc
    cmd = sprintf('curl -s -L --netrc-file "%s" -c "%s" -b "%s" -o "%s" "%s"', netrc, [f '.ck'], [f '.ck'], f, url);
else
    cmd = sprintf('curl -s -L -o "%s" "%s"', f, url);
end
system(cmd);
txt = ''; if isfile(f), txt = fileread(f); delete(f); end
end
