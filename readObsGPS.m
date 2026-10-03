function G = readObsGPS(obsFile, reader)
%READOBSGPS GPS observations of a RINEX file: fast reader with rinexread fallback.
%   G = READOBSGPS(obsFile) uses rinexread (measured faster than the
%   alternative parser). READOBSGPS(obsFile, 'fast') uses readRinexObsGPS
%   (validated identical to rinexread; falls back to rinexread on failure).
if nargin < 2, reader = 'rinexread'; end
if strcmp(reader, 'fast')
    try
        G = readRinexObsGPS(obsFile); return
    catch ME
        warning('readObsGPS:fast', 'Fast RINEX reader failed (%s) -- using rinexread.', ME.message);
    end
end
data = rinexread(obsFile);
if ~isfield(data, 'GPS'), error('readObsGPS:nogps', 'No GPS observations in %s.', obsFile); end
G = data.GPS;
end
