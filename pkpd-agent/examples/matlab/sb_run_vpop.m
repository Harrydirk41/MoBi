function nDone = sb_run_vpop(vpopXlsx, doseNames, stopTime, baselineDay, readoutDay, outCsv, nLimit)
%SB_RUN_VPOP Run the loaded model across a virtual population (an .xlsx whose row
%   1 is parameter names and each following row is one patient's parameter set),
%   optionally under one or more named doses, and write per-patient the MODEL's
%   OWN clinical-response flags. The Vantage RA model encodes the entire trial as
%   events: it captures DAS28_BL at day 199, sets ACR20/ACR50/ACR70/Remission at
%   day 284 (week 12, first-line readout), and - for patients flagged MTX_NonResp -
%   sets MTX_NonResp_TCZ_ACR20/50/70/Rem at day 600 (second-line readout). So the
%   flagship validation (TCZ response in MTX-inadequate-responders) is a built-in
%   model output: we read the flags, we do not recompute ACR.
%
%   doseNames may be several dose names joined by ';' (a SimBiology dose array is
%   applied), e.g. 'MTX_15mg_Q1W_SC_t200;TCZ8mgkg_Q4W_IV_t200'. stopTime>0 forces
%   the simulation end time (must exceed the readout the flags fire at: >=285 for
%   first-line, >=601 for the second-line TCZ readout). MATLAB reads the Excel and
%   writes the CSV (no arrays cross to Python). nLimit>0 runs an evenly-spaced
%   representative subsample of nLimit patients (the Vpop rows are severity-ordered).

    nDone = 0;
    % The Vantage model's Hill/MM terms (X^n) go briefly complex when a state
    % overshoots slightly negative; SimBiology discards the imaginary part. This is
    % a benign, known property of the model, but it emits a warning on nearly every
    % solver step and buries the results. Suppress the flood for the batch run (our
    % own diagnostics use fprintf, not warning, so nothing useful is hidden); the
    % warning state is restored automatically when the function returns.
    origWarn = warning('off', 'all');
    cleanupWarn = onCleanup(@() warning(origWarn)); %#ok<NASGU>

    m  = evalin('base', 'sbmodel');
    cs = getconfigset(m);
    try
        cs.RuntimeOptions.StatesToLog = 'all';
    catch
    end
    if stopTime > 0
        try
            cs.StopTime = stopTime;
        catch ME
            fprintf('WARNING: could not set StopTime=%g: %s\n', stopTime, ME.message);
        end
    end
    % Force solver output onto a fixed daily grid through the event-readout points.
    % The first-line ACR/MTX_NonResp events have a NARROW trigger window
    % (time>=284 & time<285): the compound expression is 0 on both sides of the
    % window, so an adaptive solver that steps over [284,285) never sees the rising
    % edge and the flags silently read 0 - unless a dose (e.g. TCZ Q4W lands on day
    % 284) happens to force a step inside it. A daily grid puts an output point AT
    % 284 (and just inside each window) so the events fire regardless of the dosing
    % schedule. Duplicates from dose-event timestamps are handled by nearest-time
    % lookup downstream.
    try
        st  = cs.StopTime;
        pts = [0:1:st, 199, 284, 284.5, 600, 600.5];
        cs.SolverOptions.OutputTimes = unique(pts(pts <= st));
    catch ME
        fprintf('WARNING: could not set OutputTimes: %s\n', ME.message);
    end

    raw   = readcell(vpopXlsx, 'Sheet', 1);
    names = raw(1, :);
    keep  = cellfun(@(x) (ischar(x) || isstring(x)) && strlength(string(x)) > 0, names);
    names = names(keep);
    nCols = numel(names);
    nTot  = size(raw, 1) - 1;
    % Representative subsample: the Vpop rows are ORDERED by disease severity, so
    % the first N patients are the sickest slice. When limiting, take evenly-spaced
    % rows across the whole population (deterministic) instead of the first N.
    if nLimit > 0 && nLimit < nTot
        sel = unique(round(linspace(1, nTot, nLimit)));
    else
        sel = 1:nTot;
    end
    nP = numel(sel);
    fprintf('vpop: %d of %d patients (evenly spaced) x %d parameters; stopTime %g, baseline day %g, readout day %g\n', ...
            nP, nTot, nCols, stopTime, baselineDay, readoutDay);

    types   = cell(1, nCols);
    missing = {};
    for j = 1:nCols
        obj = sbioselect(m, 'Name', names{j});
        if isempty(obj)
            types{j} = ''; missing{end+1} = names{j}; %#ok<AGROW>
        elseif isa(obj(1), 'SimBiology.Species')
            types{j} = 'species';
        else
            types{j} = 'parameter';
        end
    end
    if ~isempty(missing)
        fprintf('WARNING: %d vpop names not in model, e.g.: %s\n', ...
                numel(missing), strjoin(missing(1:min(6, end)), ', '));
    end

    % assemble the dose array (';'-separated names). A "name@START" suffix clones
    % the dose with an overridden StartTime (days) - this builds a SEQUENTIAL
    % switch the model does not ship: e.g. give MTX from day 200 but start TCZ
    % after the day-284 first-line readout, so MTX_NonResp is a clean pure-MTX
    % classification and TCZ then rescues those non-responders through day 600.
    d = [];
    if ~isempty(doseNames)
        parts = strsplit(string(doseNames), ';');
        for k = 1:numel(parts)
            spec = strtrim(char(parts(k)));
            if isempty(spec), continue; end
            startOverride = NaN;
            at = strfind(spec, '@');
            if ~isempty(at)
                nm = strtrim(spec(1:at(1)-1));
                startOverride = str2double(spec(at(1)+1:end));
            else
                nm = spec;
            end
            dk = getdose(m, nm);
            if isempty(dk)
                fprintf('WARNING: dose "%s" not found - skipped\n', nm);
                continue;
            end
            if ~isnan(startOverride)
                dk = copyobj(dk);                 % detached copy we can retime
                try
                    dk.StartTime = startOverride;
                    fprintf('dose "%s" retimed to StartTime=%g\n', nm, startOverride);
                catch ME
                    fprintf('WARNING: could not retime "%s": %s\n', nm, ME.message);
                end
            end
            d = [d, dk]; %#ok<AGROW>
        end
    end

    % clinical endpoints computed by the model's own events (read, never recompute)
    flags = {'ACR20','ACR50','ACR70','Remission','MTX_NonResp', ...
             'MTX_NonResp_TCZ_ACR20','MTX_NonResp_TCZ_ACR50', ...
             'MTX_NonResp_TCZ_ACR70','MTX_NonResp_TCZ_Rem'};

    fid = fopen(outCsv, 'w', 'n', 'UTF-8');
    fprintf(fid, ['patient,DAS28_BL,DAS28_base,DAS28_read,DAS28_end,' ...
                  'ACR20,ACR50,ACR70,Rem,MTX_NonResp,' ...
                  'TCZ_ACR20,TCZ_ACR50,TCZ_ACR70,TCZ_Rem\n']);

    for ii = 1:nP
        i = sel(ii);                          % real Vpop row index (1..nTot)
        rowvals = raw(i + 1, keep);
        content = {};
        for j = 1:nCols
            if isempty(types{j}) || ~isnumeric(rowvals{j})
                continue;
            end
            if strcmp(types{j}, 'species')
                content{end+1} = {'species', names{j}, 'InitialAmount', rowvals{j}}; %#ok<AGROW>
            else
                content{end+1} = {'parameter', names{j}, 'Value', rowvals{j}}; %#ok<AGROW>
            end
        end
        v = sbiovariant(sprintf('vp%d', i));
        v.Content = content;

        try
            sd = sbiosimulate(m, cs, v, d);
        catch ME
            fprintf('patient %d sim FAILED: %s\n', i, ME.message);
            continue;
        end
        try
            das = selectbyname(sd, 'DAS28_CRP').Data;
            [~, ib] = min(abs(sd.Time - baselineDay));   % nearest time (dose events
            [~, ir] = min(abs(sd.Time - readoutDay));    % create duplicate stamps)
            fv = zeros(1, numel(flags));
            for f = 1:numel(flags)
                fv(f) = local_lastval(sd, flags{f});
            end
            bl = local_lastval(sd, 'DAS28_BL');
            fprintf(fid, '%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n', ...
                    i, bl, das(ib), das(ir), das(end), ...
                    fv(1), fv(2), fv(3), fv(4), fv(5), fv(6), fv(7), fv(8), fv(9));
            nDone = nDone + 1;
        catch ME
            fprintf('patient %d readout FAILED: %s\n', i, ME.message);
        end
    end
    fclose(fid);
    fprintf('done: %d patients written\n', nDone);
end

function val = local_lastval(sd, name)
%LOCAL_LASTVAL Final logged value of a state by name; NaN if the model does not
%   log it (e.g. a constant parameter or a name absent from this build).
    val = NaN;
    try
        y = selectbyname(sd, name).Data;
        if ~isempty(y)
            val = y(end);
        end
    catch
    end
end
