function stats = persistenceBaseline(ds, lat)
%PERSISTENCEBASELINE 24-h persistence forecast (Sect. 3.4, baseline i).
%   stats = PERSISTENCEBASELINE(ds, lat) evaluates Yhat = T_t.
dsP = ds;
dsP.Yhat = ds.Tlast;         % scaled T_t
stats = evaluateTEC([], dsP, lat);
end
