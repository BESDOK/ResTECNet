function [xyz, dts] = gpsSatPos(eph, tow)
%GPSSATPOS GPS satellite ECEF position and clock from broadcast ephemeris (vectorized in tow).
%   [xyz, dts] = GPSSATPOS(eph, tow)
%     eph : ONE row of the GPS navigation table returned by rinexread
%           (Time = clock reference epoch Toc, Toe, sqrtA, Eccentricity, M0,
%            Delta_n, omega, OMEGA0, OMEGA_DOT, i0, IDOT, Cuc, Cus, Crc, Crs,
%            Cic, Cis, SVClockBias, SVClockDrift, SVClockDriftRate, TGD)
%     tow : N x 1 GPS time of week of signal transmission [s]
%   Returns ECEF position [m] (N x 3) in the frame of the transmission epoch
%   and the clock offset dts [s] (N x 1) incl. relativistic correction and
%   TGD (L1 users). IS-GPS-200 algorithm. Earth rotation during the signal
%   travel time must be applied by the caller (see rotateEarth in
%   processStationVTEC).
mu = 3.986005e14; OmE = 7.2921151467e-5; F = -4.442807633e-10;
tow = tow(:);
ecc = eph.Eccentricity; a = eph.sqrtA^2; n0 = sqrt(mu / a^3);
tk = tow - eph.Toe; tk = tk - 604800 * round(tk / 604800);
n = n0 + eph.Delta_n;
M = eph.M0 + n * tk;
E = M;
for it = 1:12, E = M + ecc * sin(E); end
nu = atan2(sqrt(1 - ecc^2) * sin(E), cos(E) - ecc);
phi = nu + eph.omega;
du = eph.Cus * sin(2*phi) + eph.Cuc * cos(2*phi);
dr = eph.Crs * sin(2*phi) + eph.Crc * cos(2*phi);
di = eph.Cis * sin(2*phi) + eph.Cic * cos(2*phi);
u = phi + du; r = a * (1 - ecc * cos(E)) + dr;
inc = eph.i0 + di + eph.IDOT * tk;
xp = r .* cos(u); yp = r .* sin(u);
Om = eph.OMEGA0 + (eph.OMEGA_DOT - OmE) * tk - OmE * eph.Toe;
xyz = [xp .* cos(Om) - yp .* cos(inc) .* sin(Om), ...
       xp .* sin(Om) + yp .* cos(inc) .* cos(Om), ...
       yp .* sin(inc)];
toc = seconds(eph.Time - dateshift(eph.Time, 'start', 'week'));   % clock reference epoch (TOW)
tc = tow - toc; tc = tc - 604800 * round(tc / 604800);
dts = eph.SVClockBias + eph.SVClockDrift * tc + eph.SVClockDriftRate * tc.^2 ...
    + F * ecc * eph.sqrtA * sin(E) - eph.TGD;
end
