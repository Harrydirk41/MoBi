function nDone = sb_run_vpop(vpopXlsx, doseName, baselineDay, readoutDay, outCsv, nLimit)
%SB_RUN_VPOP Run the loaded model across a virtual population (an .xlsx whose row
%   1 is parameter names and each following row is one patient's parameter set),
%   optionally under a named dose, and write per-patient DAS28-CRP at the
%   treatment-start BASELINE day and at the READOUT day. The clinical response
%   (ACR20/50/70 = % DAS28 improvement from that patient's own baseline) is
%   computed in Python - this matches the model's own DAS28_BL definition and
%   needs only one arm. MATLAB reads the Excel and writes the CSV (no arrays cross
%   to Python). nLimit>0 runs only the first nLimit patients.

    nDone = 0;
    m  = evalin('base', 'sbmodel');
    cs = getconfigset(m);
    try
        cs.RuntimeOptions.StatesToLog = 'all';
    catch
    end

    raw   = readcell(vpopXlsx, 'Sheet', 1);
    names = raw(1, :);
    keep  = cellfun(@(x) (ischar(x) || isstring(x)) && strlength(string(x)) > 0, names);
    names = names(keep);
    nCols = numel(names);
    nP    = size(raw, 1) - 1;
    if nLimit > 0
        nP = min(nP, nLimit);
    end
    fprintf('vpop: %d patients x %d parameters; baseline day %g, readout day %g\n', ...
            nP, nCols, baselineDay, readoutDay);

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

    d = [];
    if ~isempty(doseName)
        d = getdose(m, doseName);
    end

    fid = fopen(outCsv, 'w', 'n', 'UTF-8');
    fprintf(fid, 'patient,DAS28_base,DAS28_read\n');

    for i = 1:nP
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
            y = selectbyname(sd, 'DAS28_CRP').Data;
            [~, ib] = min(abs(sd.Time - baselineDay));   % nearest time (dose events
            [~, ir] = min(abs(sd.Time - readoutDay));    % create duplicate stamps)
            fprintf(fid, '%d,%g,%g\n', i, y(ib), y(ir));
            nDone = i;
        catch ME
            fprintf('patient %d readout FAILED: %s\n', i, ME.message);
        end
    end
    fclose(fid);
    fprintf('done: %d patients written\n', nDone);
end
