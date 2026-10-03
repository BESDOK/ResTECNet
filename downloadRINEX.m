function [obsFile, navFile] = downloadRINEX(station, d, workDir, opts)
%DOWNLOADRINEX Retrieve daily 30-s RINEX observation and broadcast navigation files.
%   [obsFile, navFile] = DOWNLOADRINEX(station, date, workDir, ...)
%     station : 4-char IGS/EPN code (e.g. 'ANKR'). All stations of
%               stationList are on the open CDDIS/BKG archives; files
%               already present in workDir (<code><ddd>0.<yy>o or RINEX-3
%               long names) are returned without downloading.
%     date    : datetime (UTC day)
%   Options
%     'cddisNetrc' netrc file with Earthdata credentials (default
%                  <pwd>/earthdata_netrc if present, else curl -n)
%     'crx2rnx'    command for Hatanaka decompression (default 'CRX2RNX')
%     'country'    ISO code (stationList.iso): RINEX-3 long names such as
%                  ANKR00TUR_R_yyyyddd0000_01D_30S_MO.crx.gz are tried first
%   Sources (tried in order):
%     obs : https://cddis.nasa.gov/archive/gnss/data/daily/yyyy/ddd/yyd/<ssss><ddd>0.<yy>d.gz  (or .Z)
%           https://cddis.nasa.gov/archive/gnss/data/daily/yyyy/ddd/yyd/<SSSS>00<CCC>_R_yyyyddd0000_01D_30S_MO.crx.gz
%           https://igs.bkg.bund.de/root_ftp/IGS/obs/yyyy/ddd/<ssss><ddd>0.<yy>d.gz
%     nav : https://cddis.nasa.gov/archive/gnss/data/daily/yyyy/brdc/BRDC00IGS_R_yyyyddd0000_01D_MN.rnx.gz
%           https://igs.bkg.bund.de/root_ftp/IGS/BRDC/yyyy/ddd/BRDC00IGS_R_yyyyddd0000_01D_MN.rnx.gz
%   Returns '' when a file could not be obtained.

arguments
    station (1,:) char
    d (1,1) datetime
    workDir (1,:) char = 'rinex_cache'
    opts.cddisNetrc (1,:) char = defaultNetrc()
    opts.crx2rnx (1,:) char = 'CRX2RNX'
    opts.sources (1,:) char {mustBeMember(opts.sources, {'both','cddis','bkg'})} = 'both'
    opts.preferBKG (1,1) logical = false   % try the BKG mirror first (spreads the load over two servers)
    opts.verbose (1,1) logical = false   % print every candidate URL with its HTTP status
    opts.country (1,:) char = ''      % ISO code for RINEX-3 long names (stationList.iso); tried first when given
end
if ~isfolder(workDir), mkdir(workDir); end
d.TimeZone = 'UTC';
yr = year(d); ddd = day(d, 'dayofyear'); yy = mod(yr, 100);
ss = lower(station); SS = upper(station);

% ---------------- observation file ----------------
obsPlain = {sprintf('%s%03d0.%02do', ss, ddd, yy), sprintf('%s%03d0.%02dO', SS, ddd, yy), ...
            sprintf('%s00*_R_%04d%03d0000_01D_30S_MO.rnx', SS, yr, ddd)};
obsFile = firstExisting(workDir, obsPlain);
if isempty(obsFile)
    crx = sprintf('%s%03d0.%02dd', ss, ddd, yy);
    cands = {};
    if ~isempty(opts.country)                                  % RINEX-3 long name (most IGS stations since ~2020)
        long = sprintf('%s00%s_R_%04d%03d0000_01D_30S_MO.crx', SS, upper(opts.country), yr, ddd);
        longRnx = strrep(long, '.crx', '.rnx');                % plain RINEX-3 (CDDIS yyo/ directory)
        cands = [cands; {sprintf('https://cddis.nasa.gov/archive/gnss/data/daily/%04d/%03d/%02dd/%s.gz', yr, ddd, yy, long), long; ...
                         sprintf('https://cddis.nasa.gov/archive/gnss/data/daily/%04d/%03d/%02do/%s.gz', yr, ddd, yy, longRnx), longRnx; ...
                         sprintf('https://igs.bkg.bund.de/root_ftp/EUREF/obs/%04d/%03d/%s.gz', yr, ddd, long), long; ...
                         sprintf('https://igs.bkg.bund.de/root_ftp/IGS/obs/%04d/%03d/%s.gz', yr, ddd, long), long}];
    end
    cands = [cands; { ...
        sprintf('https://cddis.nasa.gov/archive/gnss/data/daily/%04d/%03d/%02dd/%s.gz', yr, ddd, yy, crx), crx;
        sprintf('https://cddis.nasa.gov/archive/gnss/data/daily/%04d/%03d/%02dd/%s.Z',  yr, ddd, yy, crx), crx;
        sprintf('https://igs.bkg.bund.de/root_ftp/IGS/obs/%04d/%03d/%s.gz', yr, ddd, crx), crx;
        sprintf('https://igs.bkg.bund.de/root_ftp/IGS/obs/%04d/%03d/%s.Z',  yr, ddd, crx), crx}];
    cands = orderSources(cands, opts.sources, opts.preferBKG);
    for c = 1:size(cands, 1)
        [~, nm, ext] = fileparts(cands{c, 1});
        zPath = fullfile(workDir, [nm ext]);
        if ~isfile(zPath) && ~fetch(cands{c, 1}, zPath, opts.cddisNetrc, opts.verbose), continue; end
        if ~decompress(zPath, workDir), if isfile(zPath), delete(zPath); end, continue; end
        crxPath = fullfile(workDir, cands{c, 2});
        if endsWith(cands{c, 2}, '.rnx')                       % already plain RINEX
            obsFile = firstExisting(workDir, obsPlain);
            if ~isempty(obsFile), break; end
        elseif isfile(crxPath)
            st = system(sprintf('%s -f "%s"%s', opts.crx2rnx, crxPath, quiet()));   % Hatanaka -> RINEX
            if st ~= 0 && st ~= 2
                warning('downloadRINEX:crx', 'CRX2RNX failed or not available (https://terras.gsi.go.jp/ja/crx2rnx.html)');
                break
            end
            obsFile = firstExisting(workDir, obsPlain);
            if ~isempty(obsFile), break; end
        end
    end
end

% ---------------- navigation file ----------------
% RINEX-3 merged broadcast ephemeris of the IGS (rinexread does not read RINEX-2 nav files)
navLong = sprintf('BRDC00IGS_R_%04d%03d0000_01D_MN.rnx', yr, ddd);
navFile = firstExisting(workDir, {navLong});
if isempty(navFile)
    cands = {sprintf('https://cddis.nasa.gov/archive/gnss/data/daily/%04d/brdc/%s.gz', yr, navLong), ...
             sprintf('https://igs.bkg.bund.de/root_ftp/IGS/BRDC/%04d/%03d/%s.gz', yr, ddd, navLong)};
    cands = orderSources(reshape(cands, [], 1), opts.sources, opts.preferBKG)';
    for c = 1:numel(cands)
        zPath = fullfile(workDir, [navLong '.gz']);
        if ~isfile(zPath) && ~fetch(cands{c}, zPath, opts.cddisNetrc, opts.verbose), continue; end
        if decompress(zPath, workDir) && isfile(fullfile(workDir, navLong))
            navFile = fullfile(workDir, navLong); break
        end
    end
end
end

% -------------------------------------------------------------------------
function f = firstExisting(workDir, names)
f = '';
for i = 1:numel(names)
    L = dir(fullfile(workDir, names{i}));
    if ~isempty(L), f = fullfile(workDir, L(1).name); return; end
end
end

function ok = fetch(url, dest, netrc, verbose)
% Download with curl (same chain for CDDIS and BKG); success needs HTTP 2xx and a
% gzip / Unix-compress magic number (an HTML error page is rejected). Three
% attempts; a host that does not answer at all is skipped for 5 minutes (not
% for the whole session, so that a temporary overload does not disable a mirror).
persistent coolUntil
if isempty(coolUntil), coolUntil = containers.Map('KeyType', 'char', 'ValueType', 'double'); end
host = regexp(url, '^\w+://([^/]+)', 'tokens', 'once'); host = host{1};
ok = false; code = '';
if isKey(coolUntil, host) && now < coolUntil(host), return; end
if contains(url, 'cddis.nasa.gov')
    if isempty(netrc), nflag = '-n'; else, nflag = sprintf('--netrc-file "%s"', netrc); end
else
    nflag = '';
end
cookie = [tempname '.ck'];                                        % own cookie jar per call: safe with parallel workers
for attempt = 1:3
    [st, msg] = system(sprintf('curl -s -L %s --connect-timeout 30 --max-time 300 -c "%s" -b "%s" -w "%%{http_code}" -o "%s" "%s"', nflag, cookie, cookie, dest, url));
    code = ''; c = regexp(msg, '(\d{3})\s*$', 'tokens', 'once'); if ~isempty(c), code = c{1}; end
    ok = st == 0 && ~isempty(code) && code(1) == '2' && fileBytes(dest) > 1000 && isCompressedFile(dest);
    if ok || strcmp(code, '404'), break; end
    if isfile(dest), delete(dest); end
    pause(2 * attempt);
end
if isfile(cookie), delete(cookie); end
if ~ok && (isempty(code) || strcmp(code, '000'))
    coolUntil(host) = now + 5 / 1440;
    fprintf('  %s not answering -- skipped for 5 min\n', host);
end
if verbose, fprintf('  [%s] http %s  %s\n', tern(ok, 'OK  ', 'FAIL'), code, url); end
if ~ok && isfile(dest), delete(dest); end
end

function ok = isCompressedFile(f)
ok = false; fid = fopen(f, 'r'); if fid < 0, return; end
b = fread(fid, 2, 'uint8')'; fclose(fid);
ok = numel(b) == 2 && b(1) == 31 && (b(2) == 139 || b(2) == 157);
end

function s = tern(c, a, b), if c, s = a; else, s = b; end, end

function ok = decompress(zPath, workDir)
ok = false;
try
    if endsWith(zPath, '.gz'), gunzip(zPath, workDir); ok = true;
    else
        plain = regexprep(zPath, '\.Z$', '');
        tools = {'uncompress -k -f "%FILE%"', 'gzip -d -k -f "%FILE%"', '7z x -y -o"%DIR%" "%FILE%"', ...
                 '"C:\\Program Files\\7-Zip\\7z.exe" x -y -o"%DIR%" "%FILE%"'};
        for t = 1:numel(tools)
            st = system([strrep(strrep(tools{t}, '%FILE%', zPath), '%DIR%', workDir) quiet()]);
            if st == 0 && isfile(plain), ok = true; return; end
        end
    end
catch
end
end

function q = quiet()
% redirect tool chatter (e.g. "'uncompress' is not recognized")
if ispc, q = ' >nul 2>&1'; else, q = ' >/dev/null 2>&1'; end
end

function n = fileBytes(f)
n = 0; if isfile(f), d = dir(f); n = d(1).bytes; end
end

function c = orderSources(c, sources, preferBKG)
% filter / reorder candidate URLs (first column, or the cell itself for nav) by server
if iscell(c) && size(c, 2) >= 2, urls = c(:, 1); else, urls = c; end
isB = contains(urls, 'bkg.bund.de');
switch sources
    case 'cddis', keep = ~isB;
    case 'bkg',   keep = isB;
    otherwise,    keep = true(size(isB));
end
c = c(keep, :); isB = isB(keep);
if preferBKG, c = [c(isB, :); c(~isB, :)]; end
end
