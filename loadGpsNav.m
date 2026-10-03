function Eph = loadGpsNav(navFile)
%LOADGPSNAV GPS broadcast records of a RINEX-3 navigation file, cached per file.
%   The merged multi-GNSS file is parsed once; the small GPS table is stored
%   next to it as <navFile>.gps.mat and reused by every station of that day.
cacheFile = [navFile '.gps.mat'];
if isfile(cacheFile)
    S = load(cacheFile, 'Eph'); Eph = S.Eph; return
end
nav = rinexread(navFile); Eph = nav.GPS;
tmp = [cacheFile '.' char(java.util.UUID.randomUUID) '.tmp'];
try
    save(tmp, 'Eph'); movefile(tmp, cacheFile);
catch
    if isfile(tmp), delete(tmp); end
end
end
