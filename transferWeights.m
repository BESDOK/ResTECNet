function netTo = transferWeights(netTo, netFrom)
%TRANSFERWEIGHTS Warm start: copy learnables with matching layer/parameter
%   names and sizes from netFrom into netTo (dlnetwork). Layers whose
%   parameter sizes differ (e.g. the 2-filter head of Res-TECNet-U) keep
%   their initialization. Used to fine-tune the operational and the
%   uncertainty networks from the trained hindcast network (compute budget).
Lt = netTo.Learnables; Lf = netFrom.Learnables;
n = 0;
for i = 1:height(Lt)
    j = find(strcmp(Lf.Layer, Lt.Layer(i)) & strcmp(Lf.Parameter, Lt.Parameter(i)), 1);
    if ~isempty(j) && isequal(size(Lf.Value{j}), size(Lt.Value{i}))
        Lt.Value{i} = Lf.Value{j}; n = n + 1;
    elseif ~isempty(j) && ndims(Lf.Value{j}) == 4 && isequal(size(Lf.Value{j}, [1 2 4]), size(Lt.Value{i}, [1 2 4])) ...
            && size(Lf.Value{j}, 3) < size(Lt.Value{i}, 3)
        V = Lt.Value{i}; cin = size(Lf.Value{j}, 3);                  % extra input channels keep their init
        V(:, :, 1:cin, :) = Lf.Value{j}; Lt.Value{i} = V; n = n + 1;
    end
end
netTo.Learnables = Lt;
% batch-norm statistics
St = netTo.State; Sf = netFrom.State;
for i = 1:height(St)
    j = find(strcmp(Sf.Layer, St.Layer(i)) & strcmp(Sf.Parameter, St.Parameter(i)), 1);
    if ~isempty(j) && isequal(size(Sf.Value{j}), size(St.Value{i})), St.Value{i} = Sf.Value{j}; end
end
netTo.State = St;
fprintf('transferWeights: %d of %d learnable tensors copied\n', n, height(Lt));
end
