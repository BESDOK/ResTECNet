function summarizeStruct(S, outFile)
%SUMMARIZESTRUCT Write every small numeric/text/table leaf of a (nested) struct to a text file.
%   Large arrays (> 200 kB, e.g. maps or networks) are only listed with class and size.
%   Usage:  S = load('results\interpolation_influence.mat', 'Rint'); summarizeStruct(S.Rint, 'Rint_summary.txt'); type Rint_summary.txt
if nargin < 2, outFile = 'struct_summary.txt'; end
fid = fopen(outFile, 'w'); c = onCleanup(@() fclose(fid)); %#ok<NASGU>
walk(S, 'Rint');
    function walk(v, name)
        if isstruct(v)
            f = fieldnames(v);
            for k = 1:numel(v)
                nm = name; if numel(v) > 1, nm = sprintf('%s(%d)', name, k); end
                for i = 1:numel(f), walk(v(k).(f{i}), [nm '.' f{i}]); end
            end
            return
        end
        if isa(v, 'gpuArray'), v = gather(v); end
        w = whos('v');
        if w.bytes > 2e5
            fprintf(fid, '%s : %s %s (%.1f MB, omitted)\n', name, class(v), mat2str(size(v)), w.bytes / 1e6);
        elseif isnumeric(v) || islogical(v)
            fprintf(fid, '%s = %s\n', name, mat2str(double(v), 6));
        elseif ischar(v) || isstring(v)
            fprintf(fid, '%s = %s\n', name, strjoin(string(v), ', '));
        elseif istable(v)
            fprintf(fid, '%s : table\n%s\n', name, evalc('disp(v)'));
        else
            fprintf(fid, '%s : %s %s\n', name, class(v), mat2str(size(v)));
        end
    end
end
