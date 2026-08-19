function nDone = sb_sample_vpop(paramSpec, nSamples, baselineDay, outCsv, seed)
%SB_SAMPLE_VPOP Generate a virtual population by SAMPLING disease-driver parameters
%   and simulating each candidate to its untreated disease baseline, writing the
%   baseline DAS28-CRP per candidate. This is the Stage-3 (Vpop generation) knob:
%   the caller (an agent) chooses which parameters to vary and over what bounds,
%   and the resulting baseline DAS28-CRP distribution is scored against a clinical
%   target in Python.
%
%   paramSpec: ';'-separated entries "name,lo,hi,scale" where scale is 'lin' or
%   'log' (log = log-uniform between lo and hi). Each candidate patient draws every
%   listed parameter independently. No drug is given - we read the disease baseline.
%   seed makes the sampling reproducible.

    nDone = 0;
    origWarn = warning('off', 'all');
    cleanupWarn = onCleanup(@() warning(origWarn)); %#ok<NASGU>

    m  = evalin('base', 'sbmodel');
    cs = getconfigset(m);
    try, cs.RuntimeOptions.StatesToLog = 'all'; catch, end
    try, cs.StopTime = baselineDay; catch, end
    try, cs.SolverOptions.MaxStep = []; catch, end

    if nargin < 5 || isempty(seed), seed = 1; end
    rng(double(seed), 'twister');

    % parse paramSpec
    entries = strsplit(string(paramSpec), ';');
    names = {}; los = []; his = []; logs = [];
    for k = 1:numel(entries)
        e = strtrim(char(entries(k)));
        if isempty(e), continue; end
        parts = strsplit(e, ',');
        if numel(parts) < 3
            fprintf('WARNING: bad param spec "%s"\n', e); continue;
        end
        nm = strtrim(char(parts(1)));
        obj = sbioselect(m, 'Name', nm);
        if isempty(obj)
            fprintf('WARNING: parameter "%s" not in model - skipped\n', nm); continue;
        end
        names{end+1} = nm; %#ok<AGROW>
        los(end+1)  = str2double(parts(2)); %#ok<AGROW>
        his(end+1)  = str2double(parts(3)); %#ok<AGROW>
        isLog = numel(parts) >= 4 && strcmpi(strtrim(char(parts(4))), 'log');
        logs(end+1) = isLog; %#ok<AGROW>
    end
    nP = numel(names);
    if nP == 0, error('sb_sample_vpop:noparams', 'no valid parameters to sample'); end
    fprintf('sampling %d candidates over %d parameters to day %g (seed %d)\n', ...
            nSamples, nP, baselineDay, seed);

    fid = fopen(outCsv, 'w', 'n', 'UTF-8');
    hdr = 'sample';
    for j = 1:nP, hdr = [hdr ',' names{j}]; end %#ok<AGROW>
    fprintf(fid, '%s,DAS28_base\n', hdr);

    for i = 1:nSamples
        content = {}; vals = zeros(1, nP);
        for j = 1:nP
            if logs(j)
                vals(j) = 10 ^ (log10(los(j)) + rand() * (log10(his(j)) - log10(los(j))));
            else
                vals(j) = los(j) + rand() * (his(j) - los(j));
            end
            content{end+1} = {'parameter', names{j}, 'Value', vals(j)}; %#ok<AGROW>
        end
        v = sbiovariant(sprintf('cand%d', i)); v.Content = content;
        try
            sd = sbiosimulate(m, cs, v);           % NO dose: untreated baseline
            das = selectbyname(sd, 'DAS28_CRP').Data;
            [~, ib] = min(abs(sd.Time - baselineDay));
            fprintf(fid, '%d', i);
            for j = 1:nP, fprintf(fid, ',%g', vals(j)); end
            fprintf(fid, ',%g\n', das(ib));
            nDone = nDone + 1;
        catch ME
            fprintf('candidate %d FAILED: %s\n', i, ME.message);
        end
    end
    fclose(fid);
    fprintf('done: %d candidates written\n', nDone);
end
