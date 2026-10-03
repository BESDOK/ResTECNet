function [alpha, beta] = readNavIono(file)
%READNAVIONO GPS Klobuchar coefficients (alpha, beta; 1x4 each) from a broadcast navigation header.
%   Handles RINEX-2 ('ION ALPHA' / 'ION BETA') and RINEX-3 ('GPSA' / 'GPSB'
%   under 'IONOSPHERIC CORR') headers; Fortran 'D' exponents are accepted.
fid = fopen(file, 'r'); if fid < 0, error('readNavIono:open', 'Cannot open %s', file); end
cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
alpha = []; beta = [];
num = @(s) sscanf(strrep(strrep(s, 'D', 'E'), 'd', 'e'), '%f')';
for i = 1:600
    l = fgetl(fid);
    if ~ischar(l) || contains(l, 'END OF HEADER'), break; end
    if contains(l, 'ION ALPHA'), alpha = num(l(1:min(60, numel(l)))); end
    if contains(l, 'ION BETA'),  beta  = num(l(1:min(60, numel(l)))); end
    if startsWith(l, 'GPSA'), alpha = num(l(6:min(55, numel(l)))); end
    if startsWith(l, 'GPSB'), beta  = num(l(6:min(55, numel(l)))); end
end
if numel(alpha) < 4 || numel(beta) < 4
    error('readNavIono:missing', 'No Klobuchar coefficients in %s.', file);
end
alpha = alpha(1:4); beta = beta(1:4);
end
