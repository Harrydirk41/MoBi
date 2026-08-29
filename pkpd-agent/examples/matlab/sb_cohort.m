function nDone = sb_cohort(paramSpec, armsSpec, baselineDay, readoutDay, ...
                           nSamples, seed, stateSpec, outCsv, nExtra)
%SB_COHORT Sample a virtual COHORT and record each candidate's untreated baseline
%   severity AND its primary response flag under several therapy ARMS. This is the
%   richer cohort the multi-anchor Vpop selection needs: with per-candidate response
%   under MTX / ADA / TCZ etc., Python can then optimize selection weights so the
%   population matches SEVERAL clinical anchors at once (baseline distribution + each
%   arm's response rate) - the gQSPsim-style calibration the paper's GA implements.
%
%   Each candidate is simulated ONCE per arm here (in MATLAB, no Python bridge in the
%   loop); the weight optimization afterwards is cheap post-processing on this table.
%
%   Inputs (strings unless noted):
%     paramSpec   ';'-separated "name,lo,hi,scale" (scale 'lin'/'log') - drivers to vary.
%     armsSpec    arms separated by ';;'; each arm "label:dose1,dose2,..." (plain dose
%                 names applied from their shipped start). '' = no arms (baseline only).
%     baselineDay day to read the untreated baseline severity.
%     readoutDay  day to read each arm's primary first-line response flag.
%     nSamples    number of candidates.  seed  reproducible sampling.
%     stateSpec   ';'-joined 11 role-ordered state names (as sb_run_vpop); '' = RA
%                 defaults. The PRIMARY response flag is role 1; severity is role 10.
%     outCsv      output: sample, <params...>, sev_base, <arm labels...>.

    nDone = 0;
    origWarn = warning('off', 'all');
    cleanupWarn = onCleanup(@() warning(origWarn)); %#ok<NASGU>

    m  = evalin('base', 'sbmodel');
    cs = getconfigset(m);
    try, cs.RuntimeOptions.StatesToLog = 'all'; catch, end
    if nargin < 6 || isempty(seed), seed = 1; end
    if nargin < 9 || isempty(nExtra), nExtra = 2; end   % extra response roles per arm
    rng(double(seed), 'twister');

    % -- role-ordered state names (primary flag + severity) ---------------- %
    st = {'ACR20','ACR50','ACR70','Remission','MTX_NonResp', ...
          'MTX_NonResp_TCZ_ACR20','MTX_NonResp_TCZ_ACR50', ...
          'MTX_NonResp_TCZ_ACR70','MTX_NonResp_TCZ_Rem','DAS28_CRP','DAS28_BL'};
    if nargin >= 7 && ~isempty(stateSpec)
        parts = cellstr(strsplit(string(stateSpec), ';'));
        if numel(parts) == numel(st), st = parts; end
    end
    primaryFlag = char(st{1});
    trajState   = char(st{10});

    % -- parse paramSpec --------------------------------------------------- %
    entries = strsplit(string(paramSpec), ';');
    names = {}; los = []; his = []; logs = [];
    for k = 1:numel(entries)
        e = strtrim(char(entries(k)));
        if isempty(e), continue; end
        parts = strsplit(e, ',');
        nm = strtrim(char(parts(1)));
        if isempty(sbioselect(m, 'Name', nm))
            fprintf('WARNING: parameter "%s" not in model - skipped\n', nm); continue;
        end
        names{end+1} = nm; los(end+1) = str2double(parts(2)); %#ok<AGROW>
        his(end+1) = str2double(parts(3)); %#ok<AGROW>
        logs(end+1) = numel(parts) >= 4 && strcmpi(strtrim(char(parts(4))), 'log'); %#ok<AGROW>
    end
    nP = numel(names);
    if nP == 0, error('sb_cohort:noparams', 'no valid parameters to sample'); end

    % -- parse arms "label:doseA,doseB" ;; ... ----------------------------- %
    armLabels = {}; armDoses = {};
    if ~isempty(strtrim(char(string(armsSpec))))
        aEntries = strsplit(string(armsSpec), ';;');
        for k = 1:numel(aEntries)
            e = strtrim(char(aEntries(k)));
            if isempty(e), continue; end
            c = strfind(e, ':');
            if isempty(c), fprintf('WARNING: bad arm "%s"\n', e); continue; end
            armLabels{end+1} = strtrim(e(1:c(1)-1)); %#ok<AGROW>
            armDoses{end+1}  = strtrim(e(c(1)+1:end)); %#ok<AGROW>
        end
    end
    nArms = numel(armLabels);
    endDay = max([baselineDay, readoutDay, 1]) + 10;
    fprintf('cohort: %d candidates x %d params, %d arms, readout day %g\n', ...
            nSamples, nP, nArms, readoutDay);

    % extra first-line response roles (2,3 = ACR50/70) recorded per arm alongside the
    % primary (role 1), so Python can match the full response DISTRIBUTION, not just the
    % primary rate. Column '<label>' stays the primary (backward compatible); the extras
    % are '<label>__<state>'.
    extraRoles = {};
    for r = 2:min(1 + round(nExtra), numel(st)), extraRoles{end+1} = char(st{r}); end %#ok<AGROW>

    fid = fopen(outCsv, 'w', 'n', 'UTF-8');
    hdr = 'sample';
    for j = 1:nP, hdr = [hdr ',' names{j}]; end %#ok<AGROW>
    hdr = [hdr ',sev_base'];
    for a = 1:nArms
        hdr = [hdr ',' armLabels{a}]; %#ok<AGROW>
        for r = 1:numel(extraRoles)
            hdr = [hdr ',' armLabels{a} '__' extraRoles{r}]; %#ok<AGROW>
        end
    end
    fprintf(fid, '%s\n', hdr);

    for i = 1:nSamples
        vals = zeros(1, nP); content = {};
        for j = 1:nP
            if logs(j)
                vals(j) = 10 ^ (log10(los(j)) + rand() * (log10(his(j)) - log10(los(j))));
            else
                vals(j) = los(j) + rand() * (his(j) - los(j));
            end
            content{end+1} = {'parameter', names{j}, 'Value', vals(j)}; %#ok<AGROW>
        end
        v = sbiovariant(sprintf('cand%d', i)); v.Content = content;

        try, cs.StopTime = baselineDay; catch, end
        sevBase = NaN;
        try
            sd = sbiosimulate(m, cs, v);                 % untreated baseline
            das = selectbyname(sd, trajState).Data;
            [~, ib] = min(abs(sd.Time - baselineDay));
            sevBase = das(ib);
        catch ME
            fprintf('candidate %d baseline FAILED: %s\n', i, ME.message);
        end

        armResp  = nan(1, nArms);
        armExtra = nan(nArms, numel(extraRoles));
        for a = 1:nArms
            dn = strsplit(string(armDoses{a}), ',');
            d = [];
            for q = 1:numel(dn)
                nm = strtrim(char(dn(q)));
                if isempty(nm), continue; end
                dk = getdose(m, nm);
                if isempty(dk), fprintf('WARNING: dose "%s" not found\n', nm); continue; end
                d = [d, dk]; %#ok<AGROW>
            end
            try, cs.StopTime = endDay; catch, end
            try
                sd = sbiosimulate(m, cs, v, d);
                [~, ir] = min(abs(sd.Time - readoutDay));
                ir = max(1, ir);
                y = selectbyname(sd, primaryFlag).Data;
                if ~isempty(y), armResp(a) = y(min(ir, numel(y))); end
                for r = 1:numel(extraRoles)
                    try
                        yr = selectbyname(sd, extraRoles{r}).Data;
                        if ~isempty(yr), armExtra(a, r) = yr(min(ir, numel(yr))); end
                    catch
                    end
                end
            catch ME
                fprintf('candidate %d arm %s FAILED: %s\n', i, armLabels{a}, ME.message);
            end
        end

        fprintf(fid, '%d', i);
        for j = 1:nP, fprintf(fid, ',%g', vals(j)); end
        fprintf(fid, ',%g', sevBase);
        for a = 1:nArms
            fprintf(fid, ',%g', armResp(a));
            for r = 1:numel(extraRoles), fprintf(fid, ',%g', armExtra(a, r)); end
        end
        fprintf(fid, '\n');
        nDone = nDone + 1;
        if mod(i, max(1, round(nSamples / 20))) == 0
            fprintf('  cohort progress: %d/%d candidates\n', i, nSamples);
        end
    end
    fclose(fid);
    fprintf('cohort done: %d candidates written\n', nDone);
end
