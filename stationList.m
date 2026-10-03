function S = stationList(set)
%STATIONLIST Station sets used for the independent validation (Sect. 2.3).
%   S = STATIONLIST('igs')      16 IGS stations in four regimes
%   S = STATIONLIST('anatolia') six open-access IGS/EPN stations of Anatolia and
%       the eastern Mediterranean that are actually archived at CDDIS/BKG
%       (verified with checkStationAvailability): ISTA, MERS (Turkiye), NICO
%       (Cyprus), BSHM (Israel), AUT1 (Greece), ARUC (Armenia). TUSAGA-Aktif data
%       need a licence and ANKR/TUBI are not archived, so they are not used.
%   S = STATIONLIST('all')
%   Table with variables: code, regime, lat, lon, source ('IGS'), iso
%   (ISO-3166 country code used in RINEX-3 long file names, e.g. ANKR00TUR).
%   Coordinates are approximate (deg, for maps and labels only); the exact ECEF position used in the
%   processing is taken from the RINEX header (APPROX POSITION XYZ) or from
%   the station log and overrides these values.

arguments
    set (1,:) char {mustBeMember(set, {'igs','anatolia','all'})} = 'all'
end
igs = {
 'BOGT' 'EIA' 4.64 -74.08 'IGS' 'COL'
 'IISC' 'EIA' 13.02 77.57 'IGS' 'IND'
 'NKLG' 'EIA' 0.35 9.67 'IGS' 'GAB'
 'GUAM' 'EIA' 13.59 144.87 'IGS' 'GUM'
 'ONSA' 'MID' 57.40 11.93 'IGS' 'SWE'
 'ALGO' 'MID' 45.96 -78.07 'IGS' 'CAN'
 'MIZU' 'MID' 39.14 141.13 'IGS' 'JPN'
 'STR1' 'MID' -35.32 149.01 'IGS' 'AUS'
 'KIR0' 'HIGH' 67.88 21.06 'IGS' 'SWE'
 'THU2' 'HIGH' 76.54 -68.83 'IGS' 'GRL'
 'MAC1' 'HIGH' -54.50 158.94 'IGS' 'AUS'
 'DAV1' 'HIGH' -68.58 77.97 'IGS' 'ATA'
 'KOKB' 'OCEAN' 22.13 -159.66 'IGS' 'USA'
 'ASCG' 'OCEAN' -7.92 -14.33 'IGS' 'SHN'
 'DGAR' 'OCEAN' -7.27 72.37 'IGS' 'GBR'
 'THTI' 'OCEAN' -17.58 -149.61 'IGS' 'PYF'};
ana = {
 'AUT1' 'ANATOLIA' 40.63 22.96 'IGS' 'GRC'
 'ISTA' 'ANATOLIA' 41.10 29.02 'IGS' 'TUR'
 'MERS' 'ANATOLIA' 36.56 34.26 'IGS' 'TUR'
 'NICO' 'ANATOLIA' 35.14 33.40 'IGS' 'CYP'
 'BSHM' 'ANATOLIA' 32.78 35.02 'IGS' 'ISR'
 'ARUC' 'ANATOLIA' 40.30 44.30 'IGS' 'ARM'};
switch set
    case 'igs',      c = igs;
    case 'anatolia', c = ana;
    case 'all',      c = [igs; ana];
end
S = cell2table(c, 'VariableNames', {'code', 'regime', 'lat', 'lon', 'source', 'iso'});
S.code = string(S.code); S.regime = string(S.regime); S.source = string(S.source); S.iso = string(S.iso);
end
