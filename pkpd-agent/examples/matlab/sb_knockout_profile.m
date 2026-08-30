function nCols = sb_knockout_profile(paramCsv, readoutDay, stopTime, outCsv)
%SB_KNOCKOUT_PROFILE Freeze the rule-parameters named in paramCsv (one per line; empty/missing
%   file = freeze nothing = baseline), simulate to the disease steady state with NO drug, and
%   write every species' value at readoutDay to outCsv (header of names, one value row). Same
%   freezing mechanism as sb_knockout_readout (deactivate the assigning rule, hold constant),
%   restored before returning - this is the ablated-model observable profile used to form the
%   structure-discovery SYMPTOM (which species moved vs the intact model).

    m  = evalin('base', 'sbmodel');
    cs = getconfigset(m);
    try, cs.RuntimeOptions.StatesToLog = 'all'; catch, end
    if stopTime > 0, cs.StopTime = stopTime; end
    origWarn = warning('off', 'all');
    cleanupWarn = onCleanup(@() warning(origWarn)); %#ok<NASGU>

    names = {};
    if ~isempty(paramCsv) && exist(char(paramCsv), 'file') == 2
        txt = fileread(char(paramCsv));
        parts = regexp(strtrim(txt), '\r\n|\r|\n', 'split');
        names = parts(~cellfun(@isempty, strtrim(parts)));
    end

    savedRules = {}; savedParams = {};
    rules = m.Rules;
    for i = 1:numel(names)
        nm = strtrim(names{i});
        obj = sbioselect(m, 'Name', nm);
        if isempty(obj) || ~isa(obj(1), 'SimBiology.Parameter'), continue; end
        p = obj(1);
        for r = 1:numel(rules)
            ru = rules(r);
            lhs = regexp(ru.Rule, '^\s*([\w\.]+)\s*=', 'tokens', 'once');
            if ~isempty(lhs)
                tgt = lhs{1}; dot = strfind(tgt, '.');
                if ~isempty(dot), tgt = tgt(dot(end)+1:end); end
                if strcmp(tgt, nm) && ru.Active
                    ru.Active = false; savedRules{end+1} = ru; %#ok<AGROW>
                end
            end
        end
        savedParams{end+1} = {p, p.ConstantValue}; %#ok<AGROW>
        p.ConstantValue = true;
    end
    cleanupModel = onCleanup(@() i_restore(savedRules, savedParams));

    sd = sbiosimulate(m, cs);                          % no dose = disease steady state
    names2 = reshape(cellstr(sd.DataNames), 1, []);
    t = sd.Time; idx = find(t >= readoutDay, 1);
    if isempty(idx), idx = numel(t); end
    row = sd.Data(idx, :);

    fid = fopen(outCsv, 'w', 'n', 'UTF-8');
    fprintf(fid, '%s\n', strjoin(names2, ','));
    fclose(fid);
    writematrix(row, outCsv, 'WriteMode', 'append');
    nCols = numel(names2);
end

function i_restore(savedRules, savedParams)
    for r = 1:numel(savedRules), try, savedRules{r}.Active = true; catch, end, end
    for p = 1:numel(savedParams)
        try, savedParams{p}{1}.ConstantValue = savedParams{p}{2}; catch, end
    end
end
