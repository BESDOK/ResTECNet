%RUNALL Whole unattended workflow in THIS MATLAB session (no external launcher needed).
%   Usage: make the code folder the current folder, then type   runAll
%   1. keeps the PC awake (AC power), logs everything to log_runAll.txt
%   2. station processing (runStationProcessing), repeated up to 4 times until every
%      station-month exists (finished months are skipped, so a repeat is cheap)
%   3. frees the workers and memory, then runs runExperiments (Tables 4-7, 1, figures)
%   An error in any stage is logged and the workflow continues with the next stage.
%   Leave MATLAB open (do not sign out); you may lock the screen. Progress: checkProgress.
diary('log_runAll.txt'); diary on
ra_t0 = tic;
try, system('powercfg /change standby-timeout-ac 0'); system('powercfg /change hibernate-timeout-ac 0'); catch, end
if isfile('STATIONS_DONE.txt') && ~isfile('stations_missing_months.csv'), delete('STATIONS_DONE.txt'); end   % stale marker
ra_pass = 0;
while ~isfile('STATIONS_DONE.txt') && ra_pass < 4
    ra_pass = ra_pass + 1;
    fprintf('\n##### %s  station processing, pass %d #####\n', datestr(now), ra_pass);
    try
        runStationProcessing
    catch ra_ME
        fprintf(2, '!!! pass %d failed: %s\n%s\n', ra_pass, ra_ME.message, getReport(ra_ME, 'basic'));
    end
    ra_nMiss = NaN; if isfile('stations_missing_months.csv'), ra_nMiss = height(readtable('stations_missing_months.csv')); end
    if exist('ra_prevMiss', 'var') && ra_nMiss >= ra_prevMiss
        fprintf('No further station-months recovered in pass %d (%d missing: data not available on the servers); stopping the repeats.\n', ra_pass, ra_nMiss); break
    end
    ra_prevMiss = ra_nMiss;
end
if ~isfile('STATIONS_DONE.txt')
    fid = fopen('STATIONS_DONE.txt', 'w'); fprintf(fid, 'incomplete after 4 passes: see stations_missing_months.csv\n'); fclose(fid);
    fprintf('Some station-months are still missing after 4 passes; continuing with the evaluation.\n');
end
fprintf('\n##### %s  stations done after %.1f h; freeing memory #####\n', datestr(now), toc(ra_t0) / 3600);
delete(gcp('nocreate'));                                       % release the 16 workers
clearvars -except ra_t0
fprintf('\n##### %s  evaluation (runExperiments) #####\n', datestr(now));
try
    runExperiments
catch ra_ME
    fprintf(2, '!!! runExperiments failed: %s\n%s\n', ra_ME.message, getReport(ra_ME, 'basic'));
end
fprintf('\n##### %s  FINISHED after %.1f h #####\n', datestr(now), toc(ra_t0) / 3600);
diary off
