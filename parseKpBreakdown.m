function kp8 = parseKpBreakdown(fName, targetDay)
%PARSEKPBREAKDOWN Read the 8 three-hourly Kp values of one day from a NOAA
%   SWPC 3-day forecast text file ("NOAA Kp index breakdown" table).
%   kp8 = PARSEKPBREAKDOWN(fName, targetDay) returns an 8x1 vector (00-03,
%   ..., 21-00 UT) for the column whose header matches targetDay, or []
%   if the file/column is not found. Blank lines after the header and the
%   storm-level tags such as "(G1)" are tolerated.
kp8 = [];
if ~isfile(fName), return; end
txt = readlines(fName);
hdr = find(contains(txt, 'NOAA Kp index breakdown'), 1);
if isempty(hdr), return; end
keep = strtrim(txt(hdr+1:end)) ~= "";              % drop blank lines
body = txt(hdr+1:end); body = body(keep);
if numel(body) < 9, return; end
tok = regexp(char(body(1)), '([A-Z][a-z]{2})\s+(\d{2})', 'tokens');
col = 0;
for c = 1:numel(tok)
    if strcmpi(tok{c}{1}, datestr(targetDay, 'mmm')) && str2double(tok{c}{2}) == day(targetDay)
        col = c; break
    end
end
if col == 0, return; end
vals = zeros(8, 1);
for r = 1:8
    ln = regexprep(char(body(1 + r)), '\([A-Z]\d\)', '');
    nums = sscanf(regexprep(ln, '^\s*\d{2}-\d{2}UT', ''), '%f');
    if numel(nums) < col, kp8 = []; return; end
    vals(r) = nums(col);
end
kp8 = min(max(vals, 0), 9);
end
