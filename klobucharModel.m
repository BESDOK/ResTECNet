function dI = klobucharModel(alpha, beta, rxLat, rxLon, el, az, t)
%KLOBUCHARMODEL GPS broadcast ionospheric delay at L1 [m] (Klobuchar, 1987).
%   dI = KLOBUCHARMODEL(alpha, beta, rxLat, rxLon, el, az, t)
%     alpha, beta : 1x4 coefficients from the navigation header
%                   (rinexinfo(navFile).IonosphericCorr / 'ION ALPHA/BETA')
%     rxLat, rxLon: receiver latitude/longitude [deg]
%     el, az      : N x 1 elevation/azimuth [deg]
%     t           : N x 1 datetime (GPS/UTC)
c = 299792458;
E = el / 180; A = deg2rad(az);                      % semicircles / rad
phiU = rxLat / 180; lamU = rxLon / 180;
psi = 0.0137 ./ (E + 0.11) - 0.022;
phiI = phiU + psi .* cos(A); phiI = min(max(phiI, -0.416), 0.416);
lamI = lamU + psi .* sin(A) ./ cos(phiI * pi);
phiM = phiI + 0.064 * cos((lamI - 1.617) * pi);
tow = seconds(t - dateshift(t, 'start', 'day'));
tl = mod(4.32e4 * lamI + tow, 86400);
F = 1 + 16 * (0.53 - E).^3;
PER = beta(1) + beta(2) * phiM + beta(3) * phiM.^2 + beta(4) * phiM.^3; PER = max(PER, 72000);
AMP = alpha(1) + alpha(2) * phiM + alpha(3) * phiM.^2 + alpha(4) * phiM.^3; AMP = max(AMP, 0);
x = 2 * pi * (tl - 50400) ./ PER;
dT = F .* 5e-9;
i = abs(x) < 1.57;
dT(i) = F(i) .* (5e-9 + AMP(i) .* (1 - x(i).^2 / 2 + x(i).^4 / 24));
dI = c * dT;
end
