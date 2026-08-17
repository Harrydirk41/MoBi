function nDone = sb_run_vpop(vpopXlsx, doseName, readoutTime, outCsv, nLimit)
%SB_RUN_VPOP Run the loaded model across a virtual population (an .xlsx whose row
%   1 is parameter names and each following row is one patient's parameter set),
%   optionally under a named dose, and write the clinical readouts per patient to
%   outCsv. MATLAB reads the Excel and writes the CSV, so no arrays cross to
%   Python. nLimit>0 runs only the first nLimit patients (use a few to test fast).

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
    fprintf('vpop: %d patients x %d parameters\n', nP, nCols);

    % resolve each parameter/species once (validate names, cache type)
    types   = cell(1, nCols);
    missing = {};
    for j = 1:nCols
        obj = sbioselect(m, 'Name', names{j});
        if isempty(obj)
            types{j} = '';
            missing{end+1} = names{j}; %#ok<AGROW>
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
    readouts = {'DAS28_CRP', 'ACR20', 'ACR50', 'ACR70', 'Remission', 'Response'};

    fid = fopen(outCsv, 'w', 'n', 'UTF-8');
    fprintf(fid, 'patient,%s\n', strjoin(readouts, ','));

    for i = 1:nP
        rowvals = raw(i + 1, keep);
        content = {};
        for j = 1:nCols
            if isempty(types{j}) || ~isnumeric(rowvals{j})
                continue;                       % skip unmatched / non-numeric
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

        vals = nan(1, numel(readouts));
        for k = 1:numel(readouts)
            try
                dd = selectbyname(sd, readouts{k});
                y  = dd.Data;
                if readoutTime > 0
                    % nearest time point - robust to the DUPLICATE timestamps
                    % SimBiology inserts at each dose event (interp1 errors on
                    % non-unique sample points).
                    [~, idx] = min(abs(sd.Time - readoutTime));
                    vals(k) = y(idx);
                else
                    vals(k) = y(end);
                end
            catch
                vals(k) = NaN;
            end
        end
        fprintf(fid, '%d,%s\n', i, strjoin(cellstr(compose('%g', vals)), ','));
        nDone = i;
    end
    fclose(fid);
    fprintf('done: %d patients written to CSV\n', nDone);
end
