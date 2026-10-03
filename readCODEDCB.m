function dcb = readCODEDCB(yr, mo, workDir, kind)
%READCODEDCB Download and parse the CODE monthly satellite DCB file.
%   dcb = READCODEDCB(yr, mo, workDir, kind)
%     kind : 'P1P2' (default) or 'P1C1'
%   Returns a table with variables prn ("G01" ...), dcb_ns, rms_ns.
%   Sources: https://ftp.aiub.unibe.ch/CODE/yyyy/P1P2yymm.DCB.Z, with the
%   BKG IGS mirror https://igs.bkg.bund.de/root_ftp/IGS/products/dcb/
%   (p1p2yymm.dcb.Z, 2006-2014 and 2017-present) as fallback. Long-name
%   BIA files are handled by the same parser (" DSB" lines).
%   Satellite DCBs are used in processStationVTEC (Sect. 2.3, step 2).

arguments
    yr (1,1) double
    mo (1,1) double
    workDir (1,:) char = 'dcb_cache'
    kind (1,:) char {mustBeMember(kind, {'P1P2','P1C1'})} = 'P1P2'
end
if ~isfolder(workDir), mkdir(workDir); end
short = sprintf('%s%02d%02d.DCB', kind, mod(yr, 100), mo);
plain = fullfile(workDir, short);
if ~isfile(plain)
    % sources: BKG IGS mirror first (https://igs.bkg.bund.de/root_ftp/IGS/products/dcb/,
    % lower-case names), then AIUB (https, http). A host that times out is skipped
    % for the rest of the session.
    persistent deadHosts
    if isempty(deadHosts), deadHosts = {}; end
    urls = {sprintf('https://igs.bkg.bund.de/root_ftp/IGS/products/dcb/%s.Z', lower(short)), ...
            sprintf('https://ftp.aiub.unibe.ch/CODE/%d/%s.Z', yr, short), ...
            sprintf('http://ftp.aiub.unibe.ch/CODE/%d/%s.Z',  yr, short)};
    zPath = [plain '.Z'];
    if ispc, q = ' >nul 2>&1'; else, q = ' >/dev/null 2>&1'; end
    tools = {['uncompress -k -f "%FILE%"' q], ['gzip -d -k -f "%FILE%"' q], ['7z x -y -o"%DIR%" "%FILE%"' q], ...
             ['"C:\\Program Files\\7-Zip\\7z.exe" x -y -o"%DIR%" "%FILE%"' q]};
    for u = 1:numel(urls)
        host = regexp(urls{u}, '^\w+://([^/]+)', 'tokens', 'once'); host = host{1};
        if any(strcmp(deadHosts, host)), continue; end
        try
            websave(zPath, urls{u}, weboptions('Timeout', 30));
        catch ME
            if contains(ME.message, 'timed out') || contains(ME.message, 'connect'), deadHosts{end+1} = host; end %#ok<AGROW>
            continue
        end
        for t = 1:numel(tools)
            st = system(strrep(strrep(tools{t}, '%FILE%', zPath), '%DIR%', workDir));
            if st == 0 && isfile(plain), break; end
        end
        if ~isfile(plain)                       % 7-Zip writes the name without .Z; also try lower-case name
            lc = fullfile(workDir, lower(short)); if isfile(lc), movefile(lc, plain); end
        end
        if isfile(plain), break; end
    end
end
if ~isfile(plain)
    error('readCODEDCB:missing', 'DCB file %s not available.', short);
end
txt = readlines(plain);
prn = strings(0, 1); val = []; rms = [];
for i = 1:numel(txt)
    s = char(txt(i));
    tok = regexp(s, '^\s*([GRECJ]\d{2})\s+(?:\S+\s+)?([-+]?\d+\.\d+)\s+(\d+\.\d+)', 'tokens', 'once');
    if isempty(tok)
        tok2 = regexp(s, '^\s*DSB\s+\S*\s*([GRECJ]\d{2})\s.*?([-+]?\d+\.\d+)\s+(\d+\.\d+)\s*$', 'tokens', 'once');
        tok = tok2;
    end
    if ~isempty(tok)
        prn(end+1, 1) = string(tok{1}); %#ok<AGROW>
        val(end+1, 1) = str2double(tok{2}); %#ok<AGROW>
        rms(end+1, 1) = str2double(tok{3}); %#ok<AGROW>
    end
end
dcb = table(prn, val, rms, 'VariableNames', {'prn', 'dcb_ns', 'rms_ns'});
fprintf('CODE %s DCB %04d-%02d: %d satellites\n', kind, yr, mo, height(dcb));
end
