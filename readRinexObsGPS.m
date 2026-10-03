function G = readRinexObsGPS(file)
%READRINEXOBSGPS Fast GPS-only reader for RINEX 3.x observation files.
%   G = READRINEXOBSGPS(file) returns a timetable with the same layout as
%   rinexread(file).GPS for the variables used by the VTEC / positioning
%   chain: row times (GPS time, no time zone), SatelliteID (double) and the
%   observables among C1C C1W C1P C2W C2P C2L C2X C2S L1C L1W L1P L2W L2P
%   L2L L2X L2S (code in m, phase in cycles; blank -> NaN).
%   Only the lines of GPS satellites are parsed, which is several times
%   faster than rinexread on a multi-GNSS daily file. Event records (epoch
%   flag > 1) are skipped. Errors for RINEX-2 files (the caller falls back
%   to rinexread, see readObsGPS). Validated against rinexread on real
%   files by compareRinexReaders and on a synthetic file by
%   verifyRevisionCoverage (B15).
txt = fileread(file);
k = strfind(txt, 'END OF HEADER');
if isempty(k), error('readRinexObsGPS:header', 'No END OF HEADER in %s', file); end
k = k(1);
nl = find(txt(k:min(k + 200, numel(txt))) == newline, 1);
hl = splitlines(string(txt(1:k + nl - 2)));
body = txt(k + nl:end); clear txt
% ---- header ----
first = strtrim(hl(1));
ver = str2double(extractBefore(first, ' '));
if ~contains(hl(1), 'OBSERVATION DATA') || isnan(ver) || ver < 3
    error('readRinexObsGPS:version', 'Only RINEX 3.x observation files are supported (%s).', file);
end
types = strings(1, 0); cur = ""; nExp = 0;
for i = 1:numel(hl)
    L = char(hl(i));
    if numel(L) < 61 || ~contains(L, 'SYS / # / OBS TYPES'), continue; end
    if L(1) ~= ' '
        cur = string(L(1));
        if cur == "G", nExp = str2double(L(4:6)); types = strings(1, 0); end
    end
    if cur == "G"
        tk = split(strtrim(string(L(8:60))))';
        types = [types, tk(tk ~= "")]; %#ok<AGROW>
    end
end
if isempty(types), error('readRinexObsGPS:types', 'No GPS observation types in %s', file); end
types = types(1:min(nExp, numel(types)));
% ---- body: epoch and GPS lines ----
lines = splitlines(string(body)); clear body
isEp = startsWith(lines, ">");
epIdx = find(isEp); nEp = numel(epIdx);
epT = NaT(nEp, 1); flag = zeros(nEp, 1); nsat = zeros(nEp, 1);
for e = 1:nEp
    v = sscanf(char(extractAfter(lines(epIdx(e)), 1)), '%f');
    if numel(v) >= 8
        epT(e) = datetime(v(1), v(2), v(3), v(4), v(5), 0) + seconds(v(6));
        flag(e) = v(7); nsat(e) = v(8);
    end
end
isG = startsWith(lines, "G") & ~isEp;
for q = find(flag > 1)'                                    % event records: skip the following lines
    isG(epIdx(q) + 1 : min(epIdx(q) + nsat(q), numel(lines))) = false;
end
epId = cumsum(isEp);
gl = lines(isG); epOfLine = epId(isG); clear lines
need = 3 + 16 * numel(types);
chars = char(gl);
if size(chars, 2) < need, chars(:, end + 1:need) = ' '; end
sat = str2double(cellstr(chars(:, 2:3)));
want = ["C1C","C1W","C1P","C2W","C2P","C2L","C2X","C2S","L1C","L1W","L1P","L2W","L2P","L2L","L2X","L2S"];
names = want(ismember(want, types));
vals = zeros(numel(sat), numel(names));
for j = 1:numel(names)
    it = find(types == names(j), 1); c0 = 3 + 16 * (it - 1);
    vals(:, j) = str2double(cellstr(chars(:, c0 + 1:c0 + 14)));
end
Time = epT(epOfLine);
G = array2timetable(vals, 'RowTimes', Time, 'VariableNames', cellstr(names));
G.SatelliteID = sat;
G = movevars(G, 'SatelliteID', 'Before', 1);
% ---- plausibility ----
if numel(unique(Time)) < 100, warning('readRinexObsGPS:epochs', '%s: only %d epochs.', file, numel(unique(Time))); end
for j = 1:numel(names)
    if startsWith(names(j), "C")
        v = vals(~isnan(vals(:, j)), j);
        if ~isempty(v) && mean(v < 1e7 | v > 6e7) > 0.01
            warning('readRinexObsGPS:range', '%s: %s has implausible pseudoranges (check the column layout).', file, names(j));
        end
    end
end
end
