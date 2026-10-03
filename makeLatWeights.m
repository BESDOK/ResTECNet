function wlat = makeLatWeights(lat)
%MAKELATWEIGHTS Cosine-latitude weights, normalized to unit mean.
%   wlat = MAKELATWEIGHTS(lat) with lat in degrees (H x 1).
wlat = cosd(lat(:));
wlat = wlat / mean(wlat);
end
