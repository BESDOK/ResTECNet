function out = geoTools(op, varargin)
%GEOTOOLS Small geodetic utilities used by the station processing chain.
%   lla = GEOTOOLS('ecef2lla', xyz)            -> [lat lon h] (deg, deg, m), WGS-84
%   [az, el] = GEOTOOLS('azel', rxyz, sxyz)    -> deg (N x 1 each) via cell output
%   [ipLat, ipLon] = GEOTOOLS('ipp', lat, lon, az, el, hShell_km)
%   xyz' = GEOTOOLS('sagnac', xyz, tau)    Earth-rotation (Sagnac) correction of satellite positions
%   M = GEOTOOLS('mslm', el, hShell_km, alpha) -> modified single-layer mapping
%        function, M = 1 / sqrt(1 - (R/(R+h) * sin(alpha * z))^2), z = 90 - el
%        (CODE convention: alpha = 0.9782; shell height per Sect. 2.3)
%   Outputs are returned as a cell array {..} for multi-output ops.

R = 6371.0;
switch op
    case 'ecef2lla'
        xyz = varargin{1};
        a = 6378137; f = 1/298.257223563; e2 = f*(2-f);
        x = xyz(:,1); y = xyz(:,2); z = xyz(:,3);
        lon = atan2(y, x); p = hypot(x, y);
        lat = atan2(z, p*(1-e2));
        for it = 1:5
            N = a ./ sqrt(1 - e2*sin(lat).^2);
            h = p ./ cos(lat) - N;
            lat = atan2(z, p .* (1 - e2 * N ./ (N + h)));
        end
        N = a ./ sqrt(1 - e2*sin(lat).^2); h = p ./ cos(lat) - N;
        out = [rad2deg(lat), rad2deg(lon), h];
    case 'azel'
        rxyz = varargin{1}; sxyz = varargin{2};
        lla = geoTools('ecef2lla', rxyz);
        lat = deg2rad(lla(1)); lon = deg2rad(lla(2));
        d = sxyz - rxyz;
        Rm = [-sin(lon) cos(lon) 0; -sin(lat)*cos(lon) -sin(lat)*sin(lon) cos(lat); ...
              cos(lat)*cos(lon) cos(lat)*sin(lon) sin(lat)];
        enu = d * Rm';
        az = mod(rad2deg(atan2(enu(:,1), enu(:,2))), 360);
        el = rad2deg(atan2(enu(:,3), hypot(enu(:,1), enu(:,2))));
        out = {az, el};
    case 'ipp'
        [lat, lon, az, el, h] = varargin{1:5};
        z = deg2rad(90 - el); zp = asin(R/(R+h) * sin(z));
        psi = z - zp;                                    % Earth central angle
        latr = deg2rad(lat); azr = deg2rad(az);
        ipLat = asin(sin(latr)*cos(psi) + cos(latr)*sin(psi).*cos(azr));
        ipLon = deg2rad(lon) + asin(sin(psi).*sin(azr)./cos(ipLat));
        out = {rad2deg(ipLat), wrapTo180(rad2deg(ipLon))};
    case 'sagnac'
        % express the emission-epoch satellite position(s) in the ECEF frame of
        % the reception epoch (rotated by omega*tau): x' = x cos + y sin, y' = -x sin + y cos
        xyz = varargin{1}; tau = varargin{2}; a = 7.2921151467e-5 * tau(:);
        out = [xyz(:,1) .* cos(a) + xyz(:,2) .* sin(a), -xyz(:,1) .* sin(a) + xyz(:,2) .* cos(a), xyz(:,3)];
    case 'mslm'
        [el, h, alpha] = varargin{1:3};
        z = deg2rad(90 - el);
        out = 1 ./ sqrt(1 - (R/(R+h) * sin(alpha * z)).^2);
    otherwise
        error('geoTools:op', 'Unknown operation %s', op);
end
end

function x = wrapTo180(x)
x = mod(x + 180, 360) - 180;
end
