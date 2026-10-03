function buildTECDatabase(startDate, endDate, outFile, workDir, opts)
%BUILDTECDATABASE Download GIM products and assemble an hourly TEC database.
%   BUILDTECDATABASE(startDate, endDate, outFile, workDir, 'product', P)
%     startDate, endDate : datetime or date strings, e.g. '2024-05-01'
%     outFile            : output MAT-file (default 'TECdatabase.mat')
%     workDir            : download/cache folder (default 'ionex_cache')
%   Options
%     'product'  : 'CODG' (CODE final, default)
%                  'C1PG' (CODE one-day predicted GIM, baseline ii of the
%                          revised manuscript, Sect. 3.6)
%                  'UQRG' (UPC rapid 15-min GIM, independent reference,
%                          Sect. 2.4; subsampled to the hourly epochs)
%     'dropDuplicateLon' : drop the +180 deg column (default true, W = 72)
%     'fillGaps'         : call fillTECGaps (default true)
%     'cddisNetrc'       : path of a netrc file with the Earthdata credentials
%                          for CDDIS ("machine urs.earthdata.nasa.gov login U
%                          password P"); default: <pwd>/earthdata_netrc if it
%                          exists, else curl -n (~/.netrc or %HOME%\_netrc)
%
%   Resumable: every completed year is saved as <outFile>_yYYYY.mat and
%   reused on restart (delete a year file to rebuild it). TEC is stored as
%   single.
%   Output variables: TEC (H x W x K), epochs (K x 1 hourly UTC), lat,
%   lon, badMask (K x 1, true where a gap > 6 h could not be filled),
%   product.
%
%   Archive naming handled:
%     CODG  <=2022 CODGddd0.yyI.Z      >=2023 COD0OPSFIN_yyyyddd0000_01D_01H_GIM.INX.gz
%     C1PG  <=2022 C1PGddd0.yyI.Z      >=2023 COD0OPSPRD_yyyyddd0000_01D_01H_GIM.INX.gz
%           (both at http://ftp.aiub.unibe.ch/CODE/yyyy/)
%     UQRG  uqrgddd0.yyi.Z / UPC0OPSRAP_yyyyddd0000_01D_15M_GIM.INX.gz at
%           https://cddis.nasa.gov/archive/gnss/products/ionex/yyyy/ddd/
%           (Earthdata login; curl with .netrc), with fallback
%           http://chapman.upc.es/tomion/rapid/yyyy/ddd_yymmdd.15min/
%   The pre-19-Oct-2014 2-h CODG record is resampled to 1 h in a sun-fixed
%   frame (Sect. 2.1); its influence is quantified by interpolationInfluence.m.

arguments
    startDate
    endDate
    outFile (1,:) char = 'TECdatabase.mat'
    workDir (1,:) char = 'ionex_cache'
    opts.product (1,:) char {mustBeMember(opts.product, {'CODG','C1PG','UQRG'})} = 'CODG'
    opts.dropDuplicateLon (1,1) logical = true
    opts.fillGaps (1,1) logical = true
    opts.cddisNetrc (1,:) char = defaultNetrc()
end

d1 = datetime(startDate); d2 = datetime(endDate);
d1.TimeZone = 'UTC'; d2.TimeZone = 'UTC';
if ~isfolder(workDir), mkdir(workDir); end
if ~isempty(opts.cddisNetrc) && ~isfile(opts.cddisNetrc)
    error('buildTECDatabase:netrc', 'Earthdata credentials file not found: %s', opts.cddisNetrc);
end

% Resumable: each calendar year is assembled and saved as
% <outFile>_yYYYY.mat as soon as it is complete; a restart skips those.
[outDir, outStem] = fileparts(outFile); if isempty(outDir), outDir = '.'; end
years = year(d1):year(d2);
lat = []; lon = []; consecMiss = 0;
yearFiles = strings(numel(years), 1);
for iy = 1:numel(years)
    yr = years(iy);
    yFile = fullfile(outDir, sprintf('%s_y%04d.mat', outStem, yr)); yearFiles(iy) = yFile;
    if isfile(yFile)
        M = matfile(yFile); lat = M.lat; lon = M.lon;
        fprintf('%s: year %d found (%s)\n', opts.product, yr, yFile); continue
    end
    yTEC = []; yEp = datetime.empty(0,1); yEp.TimeZone = 'UTC';
    dA = max(d1, datetime(yr,1,1,'TimeZone','UTC')); dB = min(d2, datetime(yr,12,31,'TimeZone','UTC'));
    for d = dA:dB
        localFile = fetchIONEX(opts.product, d, workDir, opts.cddisNetrc);
        if isempty(localFile)
            warning('buildTECDatabase:missing', 'No %s file for %s -- skipped.', opts.product, datestr(d));
            consecMiss = consecMiss + 1;
            if consecMiss >= 30
                error('buildTECDatabase:stalled', ['30 consecutive days without data: check the Earthdata ' ...
                    'netrc file (''cddisNetrc''), the network, or the product name -- aborting.']);
            end
            continue
        end
        consecMiss = 0;
        [Td, ed, lat, lon, okRead] = readDay(localFile, opts.dropDuplicateLon);
        if ~okRead                                   % corrupt / truncated download: refetch once
            delete(localFile); zc = dir([localFile '.*']); for z = 1:numel(zc), delete(fullfile(zc(z).folder, zc(z).name)); end
            localFile = fetchIONEX(opts.product, d, workDir, opts.cddisNetrc);
            if ~isempty(localFile), [Td, ed, lat, lon, okRead] = readDay(localFile, opts.dropDuplicateLon); end
            if ~okRead
                warning('buildTECDatabase:corrupt', '%s %s unreadable after re-download -- skipped.', opts.product, datestr(d));
                consecMiss = consecMiss + 1; continue
            end
        end
        if strcmp(opts.product, 'UQRG')
            keep = minute(ed) == 0; Td = Td(:, :, keep); ed = ed(keep);
        end
        if numel(ed) > 1 && hours(ed(2) - ed(1)) == 2     % 2-h record -> 1 h (uses the 24:00 map when present)
            [Td, ed] = resampleSunFixed(Td, ed, lon);
        end
        if d < d2 && ~isempty(ed) && ed(end) == dateshift(d, 'start', 'day') + days(1)
            Td(:, :, end) = []; ed(end) = [];            % 24:00 == next day's 00:00
        end
        yTEC = cat(3, yTEC, single(Td)); yEp = [yEp; ed]; %#ok<AGROW>
        fprintf('%s %s : %d hourly maps (year total %d)\n', opts.product, datestr(d, 'yyyy-mm-dd'), numel(ed), numel(yEp));
    end
    S.TEC = yTEC; S.epochs = yEp; S.lat = lat; S.lon = lon; %#ok<STRNU>
    save(yFile, '-struct', 'S', '-v7.3'); clear S yTEC
    fprintf('%s: year %d checkpoint written (%s)\n', opts.product, yr, yFile);
end
if isempty(lat), error('buildTECDatabase:empty', 'No maps retrieved.'); end

% ---- memory-lean assembly on a preallocated hourly grid ----
e0 = dateshift(d1, 'start', 'day'); e1 = dateshift(d2, 'start', 'day') + hours(23);
epochs = (e0:hours(1):e1)'; Kh = numel(epochs);
H = numel(lat); W = numel(lon);
TEC = nan(H, W, Kh, 'single');
have = false(Kh, 1); nRej = 0; nTot = 0;
for iy = 1:numel(years)
    Y = load(yearFiles(iy), 'TEC', 'epochs');
    if isempty(Y.epochs), continue; end
    [tf, loc] = ismember(dateshift(Y.epochs, 'start', 'hour'), epochs);
    for k = 1:numel(Y.epochs)
        nTot = nTot + 1;
        if ~tf(k) || have(loc(k)), continue; end
        m = Y.TEC(:, :, k);
        if opts.fillGaps && (any(isnan(m), 'all') || any(m < 0, 'all') || any(m > 250, 'all')), nRej = nRej + 1; continue; end
        TEC(:, :, loc(k)) = m; have(loc(k)) = true;
    end
    clear Y
end
fprintf('%s: %d maps placed on %d hourly epochs, %d rejected by screening (%.3f %%)\n', ...
    opts.product, nnz(have), Kh, nRej, 100 * nRej / max(nTot, 1));
if opts.fillGaps
    [TEC, badMask] = fillGapsInPlace(TEC, have, lon, 6);
    fprintf('%s: %d epochs filled, %d left as gaps\n', opts.product, nnz(~have) - nnz(badMask), nnz(badMask));
else
    badMask = ~have;                                     % NaN maps flagged, no interpolation
end
product = opts.product; %#ok<NASGU>
save(outFile, 'TEC', 'epochs', 'lat', 'lon', 'badMask', 'product', '-v7.3');
fprintf('Saved %d maps (%s) to %s\n', numel(epochs), opts.product, outFile);
end

% =========================================================================
function localFile = fetchIONEX(product, d, workDir, netrc)
% Candidate (url, compressed name) pairs, tried in order. NASA CDDIS
% (Earthdata login via --netrc-file) holds the complete IGS ionex archive
% since 1998 and is tried first; AIUB (https, then http) is the fallback for
% the CODE products, chapman.upc.es for UQRG.
yr = year(d); ddd = day(d, 'dayofyear'); yy = mod(yr, 100);
cddis = @(f) {sprintf('https://cddis.nasa.gov/archive/gnss/products/ionex/%d/%03d/%s', yr, ddd, f), f};
aiub  = @(f) {sprintf('https://ftp.aiub.unibe.ch/CODE/%d/%s', yr, f), f; ...
              sprintf('http://ftp.aiub.unibe.ch/CODE/%d/%s',  yr, f), f};
switch product
    case 'CODG'
        longN = sprintf('COD0OPSFIN_%04d%03d0000_01D_01H_GIM.INX.gz', yr, ddd);
        shortN = sprintf('CODG%03d0.%02dI.Z', ddd, yy);
        cands = [cddis(longN); cddis(lowerStem(shortN)); aiub(longN); aiub(shortN)];
    case 'C1PG'
        longN = sprintf('COD0OPSPRD_%04d%03d0000_01D_01H_GIM.INX.gz', yr, ddd);
        shortN = sprintf('C1PG%03d0.%02dI.Z', ddd, yy);
        cands = [cddis(longN); cddis(lowerStem(shortN)); aiub(longN); aiub(shortN)];
    case 'UQRG'
        longN = sprintf('UPC0OPSRAP_%04d%03d0000_01D_15M_GIM.INX.gz', yr, ddd);
        shortN = sprintf('uqrg%03d0.%02di.Z', ddd, yy);
        cands = [cddis(longN); cddis(shortN); ...
                 {sprintf('http://chapman.upc.es/tomion/rapid/%d/%03d_%s.15min/%s', yr, ddd, datestr(d, 'yymmdd'), shortN), shortN}];
end
localFile = '';
for c = 1:size(cands, 1)
    zName = cands{c, 2};
    plainName = regexprep(zName, '\.(gz|Z)$', '');
    plainPath = fullfile(workDir, plainName);
    if isfile(plainPath), localFile = plainPath; return; end
    zPath = fullfile(workDir, zName);
    if isfile(zPath) && ~isCompressedFile(zPath), delete(zPath); end   % stale HTML / truncated
    if ~isfile(zPath)
        ok = downloadFile(cands{c, 1}, zPath, netrc);
        if ~ok, continue; end
    end
    if ~decompressFile(zPath, workDir)
        if isfile(zPath), delete(zPath); end
        continue
    end
    if isfile(plainPath), localFile = plainPath; return; end
end
end

function n = lowerStem(f)
% lower-case the file name but keep the extension (CDDIS: codg0010.17i.Z)
[~, stem, ext] = fileparts(f); n = [lower(stem) ext];
end

function ok = isCompressedFile(f)
% gzip magic 1F 8B, Unix compress magic 1F 9D
ok = false; fid = fopen(f, 'r'); if fid < 0, return; end
b = fread(fid, 2, 'uint8')'; fclose(fid);
ok = numel(b) == 2 && b(1) == 31 && (b(2) == 139 || b(2) == 157);
end

% -------------------------------------------------------------------------
function ok = downloadFile(url, dest, netrc)
% websave for open archives; curl with --netrc-file + cookie jar for CDDIS.
% A host that fails with a connection timeout is skipped for the rest of
% the session (persistent list), so an unreachable mirror costs one wait.
persistent deadHosts
if isempty(deadHosts), deadHosts = {}; end
host = regexp(url, '^\w+://([^/]+)', 'tokens', 'once'); host = host{1};
ok = false;
if any(strcmp(deadHosts, host)), return; end
if contains(url, 'cddis.nasa.gov')
    cookie = fullfile(fileparts(dest), 'cddis_cookies.txt');
    if isempty(netrc), nflag = '-n'; else, nflag = sprintf('--netrc-file "%s"', netrc); end
    cmd = sprintf('curl -s -S -L %s -c "%s" -b "%s" -w "%%{http_code}" -o "%s" "%s"', nflag, cookie, cookie, dest, url);
    [st, msg] = system(cmd);
    code = regexp(msg, '(\d{3})\s*$', 'tokens', 'once');
    ok = (st == 0) && ~isempty(code) && code{1}(1) == '2' && fileBytes(dest) > 1000 && isCompressedFile(dest);
    if ~ok
        if isfile(dest), delete(dest); end
        if st == 26 || st == 2 || contains(msg, 'netrc')
            error('buildTECDatabase:netrc', ['curl cannot read the Earthdata credentials file. Pass ' ...
                '''cddisNetrc'', <path to earthdata_netrc> or place earthdata_netrc in the working folder.']);
        end
        if ~isempty(code) && strcmp(code{1}, '404')
            % file simply not on CDDIS under this name: silent, try next candidate
        else
            fprintf('  CDDIS download failed: %s (curl status %d, http %s) %s\n', url, st, ...
                strjoin(code, ''), strtrim(regexprep(msg, '\d{3}\s*$', '')));
        end
    end
    return
end
try
    websave(dest, url, weboptions('Timeout', 60));
    ok = fileBytes(dest) > 1000 && isCompressedFile(dest);
    if ~ok && isfile(dest), delete(dest); end
catch ME
    if isfile(dest), delete(dest); end
    if contains(ME.message, 'timed out') || contains(ME.message, 'connect')
        deadHosts{end+1} = host;
        fprintf('  %s unreachable (%s) -- skipped for the rest of this session\n', host, strtok(ME.message, '.'));
    elseif ~contains(ME.message, '404')
        fprintf('  download failed: %s (%s)\n', url, ME.message);
    end
end
end

% -------------------------------------------------------------------------
function ok = decompressFile(zPath, workDir)
% .gz via gunzip; .Z (Unix compress) via uncompress / gzip / 7-Zip.
ok = false;
try
    if endsWith(zPath, '.gz')
        gunzip(zPath, workDir); ok = true; return
    end
    tools = {'uncompress -k -f "%FILE%"', 'gzip -d -k -f "%FILE%"', ...
             '7z x -y -o"%DIR%" "%FILE%" > nul 2>&1', ...
             '"C:\Program Files\7-Zip\7z.exe" x -y -o"%DIR%" "%FILE%" > nul 2>&1'};
    plainPath = regexprep(zPath, '\.Z$', '');
    for t = 1:numel(tools)
        cmd = strrep(strrep(tools{t}, '%FILE%', zPath), '%DIR%', workDir);
        [st, ~] = system(cmd);
        if st == 0 && isfile(plainPath), ok = true; return; end
    end
    warning('buildTECDatabase:Z', ['%s downloaded but not decompressed: install 7-Zip ' ...
        '(https://www.7-zip.org, default folder) or gzip/uncompress on the PATH.'], zPath);
catch ME
    warning('buildTECDatabase:decomp', '%s', ME.message);
end
end

% -------------------------------------------------------------------------
function [Th, eh] = resampleSunFixed(T2, e2, lon)
% 2-h -> 1-h interpolation in a sun-fixed frame (IONEX recommendation).
[H, W, K] = size(T2);
eh = (e2(1):hours(1):e2(end))';
Th = zeros(H, W, numel(eh), 'like', T2);
dlon = lon(2) - lon(1);
for n = 1:numel(eh)
    k = find(e2 <= eh(n), 1, 'last');
    if e2(k) == eh(n) || k == K
        Th(:, :, n) = T2(:, :, min(k, K));
        continue
    end
    a = hours(eh(n) - e2(k)) / hours(e2(k+1) - e2(k));
    A = circshift(T2(:, :, k),   -round(15 * hours(eh(n) - e2(k)) / dlon), 2);
    B = circshift(T2(:, :, k+1),  round(15 * hours(e2(k+1) - eh(n)) / dlon), 2);
    Th(:, :, n) = (1 - a) * A + a * B;
end
end

% -------------------------------------------------------------------------
function n = fileBytes(f)
n = 0; if isfile(f), d = dir(f); n = d(1).bytes; end
end

function f = defaultNetrc()
f = fullfile(pwd, 'earthdata_netrc'); if ~isfile(f), f = ''; end
end

function [Td, ed, lat, lon, ok] = readDay(localFile, dropDup)
% read one IONEX file; ok = false when it is unreadable or has too few maps
Td = []; ed = []; lat = []; lon = []; ok = false;
try
    [Td, ed, lat, lon] = readIONEXTEC(localFile, 'dropDuplicateLon', dropDup);
    ok = numel(ed) >= 12;                                 % 13 (2-h), 25 (1-h) or 97 (15-min) expected
    if ~ok, warning('buildTECDatabase:short', '%s: only %d maps.', localFile, numel(ed)); end
catch ME
    warning('buildTECDatabase:read', '%s: %s', localFile, ME.message);
end
end
