function nCols = sb_simulate_csv(doseName, variantName, stopTime, outCsv)
%SB_SIMULATE_CSV Simulate the loaded model (optionally with a named dose and/or
%   named variant) and write Time + all state trajectories to outCsv (a header row
%   of exact names, then numeric rows). Returns the number of columns (a scalar,
%   which marshals fine). Empty doseName/variantName mean "none"; stopTime <= 0
%   means "use the model's configured StopTime".

    % The Vantage model's Hill/MM terms (X^n) go briefly complex when a state overshoots
    % slightly negative; SimBiology discards the imaginary part and warns on nearly every
    % solver step, burying real output. Suppress the flood for this call (our diagnostics
    % use fprintf, not warning); the state is restored automatically on return.
    origWarn = warning('off', 'all');
    cleanupWarn = onCleanup(@() warning(origWarn)); %#ok<NASGU>

    m  = evalin('base', 'sbmodel');
    cs = getconfigset(m);
    if stopTime > 0
        cs.StopTime = stopTime;
    end
    % the project may have been saved logging only a few states; log ALL so every
    % species (cytokines, ACR/DAS28 clinical endpoints) is returned.
    try
        cs.RuntimeOptions.StatesToLog = 'all';
    catch
    end

    v = [];
    if ~isempty(variantName)
        v = getvariant(m, variantName);
    end
    d = [];
    if ~isempty(doseName)
        d = getdose(m, doseName);
    end

    sd = sbiosimulate(m, cs, v, d);

    names = reshape(cellstr(sd.DataNames), 1, []);
    header = [{'Time'}, names];
    M = [sd.Time, sd.Data];

    fid = fopen(outCsv, 'w', 'n', 'UTF-8');
    if fid < 0
        error('sb_simulate_csv:open', 'cannot open %s for writing', outCsv);
    end
    fprintf(fid, '%s\n', strjoin(header, ','));
    fclose(fid);
    writematrix(M, outCsv, 'WriteMode', 'append');

    nCols = numel(header);
end
