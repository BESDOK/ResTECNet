function [TEC, epochs, lat, lon, RMS] = readIONEXTEC(fileName, opts)
%READIONEXTEC Read vertical TEC (and optionally RMS) maps from an IONEX file.
%   [TEC, epochs, lat, lon] = READIONEXTEC(fileName)
%   [TEC, epochs, lat, lon, RMS] = READIONEXTEC(fileName, 'readRMS', true)
%
%     TEC    : H x W x K array of vertical TEC values [TECU]
%     epochs : K x 1 datetime vector (UTC) of the map epochs
%     lat    : H x 1 latitudes  [deg]  (e.g. +87.5 : -2.5 : -87.5)
%     lon    : W x 1 longitudes [deg]
%     RMS    : H x W x K RMS maps [TECU] (NaN if the file has none)
%
%   Options
%     'dropDuplicateLon' (default true)
%         The IONEX grid lists both -180 and +180 deg, which describe the
%         same meridian. When true, the redundant last column is removed
%         so that W = 72 unique meridians (Sect. 3.1 of the revised
%         manuscript, Reviewer #2 minor comment 1). Use
%         DUPLICATELONCOLUMN to re-create the 73rd column when writing
%         IONEX output.
%     'readRMS' (default false)   also parse the RMS map blocks.
%
%   Works for CODG (final), C1PG/C2PG (predicted) and UQRG (UPC rapid,
%   15-min) products. 9999 -> NaN. MATLAB R2024a.

arguments
    fileName (1,:) char
    opts.dropDuplicateLon (1,1) logical = true
    opts.readRMS (1,1) logical = false
end

fid = fopen(fileName, 'r');
if fid < 0
    error('readIONEXTEC:openFail', 'Cannot open file: %s', fileName);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

expFactor = -1;
lat1 = []; lat2 = []; dlat = [];
lon1 = []; lon2 = []; dlon = [];

% ---------- header ----------
while true
    line = fgetl(fid);
    if ~ischar(line), error('readIONEXTEC:header', 'Unexpected end of header.'); end
    if numel(line) < 61, continue; end
    label = strtrim(line(61:end));
    switch label
        case 'EXPONENT'
            expFactor = sscanf(line(1:60), '%f', 1);
        case 'LAT1 / LAT2 / DLAT'
            v = sscanf(line(1:60), '%f', 3); lat1 = v(1); lat2 = v(2); dlat = v(3);
        case 'LON1 / LON2 / DLON'
            v = sscanf(line(1:60), '%f', 3); lon1 = v(1); lon2 = v(2); dlon = v(3);
        case 'END OF HEADER'
            break
    end
end
if isempty(lat1) || isempty(lon1)
    error('readIONEXTEC:grid', 'Grid definition missing in header.');
end

lat = (lat1:dlat:lat2)';
lonFull = (lon1:dlon:lon2)';
Wfull = numel(lonFull); H = numel(lat);

TECc = {}; RMSc = {}; epc = {};

% ---------- body ----------
while ~feof(fid)
    line = fgetl(fid);
    if ~ischar(line), break; end
    if numel(line) < 61, continue; end
    label = strtrim(line(61:end));
    isTEC = strcmp(label, 'START OF TEC MAP');
    isRMS = strcmp(label, 'START OF RMS MAP');
    if ~(isTEC || (isRMS && opts.readRMS)), continue; end

    line = fgetl(fid);                           % EPOCH OF CURRENT MAP
    t = sscanf(line(1:60), '%f', 6);
    thisEpoch = datetime(t(1), t(2), t(3), t(4), t(5), t(6), 'TimeZone', 'UTC');
    map = nan(H, Wfull);
    for i = 1:H
        fgetl(fid);                              % LAT/LON1/LON2/DLON/H
        row = zeros(1, 0);
        truncated = false;
        while numel(row) < Wfull
            dline = fgetl(fid);
            if ~ischar(dline), truncated = true; break; end
            row = [row, sscanf(dline, '%f')']; %#ok<AGROW>
        end
        if truncated
            warning('readIONEXTEC:truncated', '%s: file ends inside map %d (row %d) -- map discarded.', fileName, numel(TECc) + 1, i);
            map = []; break
        end
        row(row == 9999) = NaN;
        map(i, :) = row(1:Wfull);
    end
    if isempty(map), break; end
    map = map * 10^expFactor;
    if isTEC
        TECc{end+1} = map; epc{end+1} = thisEpoch; %#ok<AGROW>
    else
        RMSc{end+1} = map; %#ok<AGROW>
    end
end

K = numel(TECc);
TEC = zeros(H, Wfull, K);
for k = 1:K, TEC(:, :, k) = TECc{k}; end
epochs = datetime.empty(0, 1); epochs.TimeZone = 'UTC';
for k = 1:K, epochs(k, 1) = epc{k}; end

if opts.readRMS && numel(RMSc) == K
    RMS = zeros(H, Wfull, K);
    for k = 1:K, RMS(:, :, k) = RMSc{k}; end
else
    RMS = nan(H, Wfull, K);
end

lon = lonFull;
if opts.dropDuplicateLon && abs((lonFull(end) - lonFull(1)) - 360) < 1e-6
    TEC(:, end, :) = [];
    RMS(:, end, :) = [];
    lon(end) = [];
end
end
