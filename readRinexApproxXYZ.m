function xyz = readRinexApproxXYZ(file)
%READRINEXAPPROXXYZ Approximate receiver ECEF position [m] from a RINEX observation header.
%   xyz = READRINEXAPPROXXYZ(file) parses the 'APPROX POSITION XYZ' record
%   (rinexinfo does not return it). Errors if the record is missing or zero.
fid = fopen(file, 'r'); if fid < 0, error('readRinexApproxXYZ:open', 'Cannot open %s', file); end
cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
xyz = [];
for i = 1:600
    l = fgetl(fid);
    if ~ischar(l) || contains(l, 'END OF HEADER'), break; end
    if contains(l, 'APPROX POSITION XYZ')
        xyz = sscanf(l(1:min(60, numel(l))), '%f')'; break
    end
end
if numel(xyz) ~= 3 || all(xyz == 0)
    error('readRinexApproxXYZ:missing', 'No usable APPROX POSITION XYZ in %s (pass ''rxXYZ'' explicitly).', file);
end
end
