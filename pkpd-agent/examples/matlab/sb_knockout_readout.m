function val = sb_knockout_readout(paramCsv, doseNames, stopTime, readoutDay, readoutState)
%SB_KNOCKOUT_READOUT Ablate regulatory edges and read the clinical readout.
%   Freezes the rule-parameters named in paramCsv (one name per line; an empty/missing
%   file = knock out nothing = the baseline run) by DEACTIVATING the rule that assigns
%   each and holding the parameter constant at its current value. This severs the
%   parameter's source-dependence - the regulatory edge - while leaving its nominal
%   value in place. Then simulates the loaded model under doseNames (one or more dose
%   names joined by ';', a SimBiology dose array) and returns readoutState (e.g.
%   'DAS28_CRP') at readoutDay. The model is RESTORED before returning (rules
%   reactivated, ConstantValue reset), so many knockouts can run against one loaded
%   project without reloading.
%
%   The Vantage model's Hill/MM terms go briefly complex on slight negative overshoot;
%   SimBiology discards the imaginary part and warns on nearly every step, so warnings
%   are silenced for the run (restored on return).

    m  = evalin('base', 'sbmodel');
    cs = getconfigset(m);
    try
        cs.RuntimeOptions.StatesToLog = 'all';
    catch
    end
    if stopTime > 0
        cs.StopTime = stopTime;
    end
    origWarn = warning('off', 'all');
    cleanupWarn = onCleanup(@() warning(origWarn)); %#ok<NASGU>

    % --- read the parameter names to freeze (empty file / missing -> baseline run) ---
    names = {};
    if ~isempty(paramCsv) && exist(char(paramCsv), 'file') == 2
        txt = fileread(char(paramCsv));
        parts = regexp(strtrim(txt), '\r\n|\r|\n', 'split');
        names = parts(~cellfun(@isempty, strtrim(parts)));
    end

    % --- freeze each named rule-parameter; remember what we changed, to restore ---
    savedRules  = {};   % {ruleObj}
    savedParams = {};   % {paramObj, oldConstantValue}
    rules = m.Rules;
    for i = 1:numel(names)
        nm  = strtrim(names{i});
        obj = sbioselect(m, 'Name', nm);
        if isempty(obj) || ~isa(obj(1), 'SimBiology.Parameter')
            continue
        end
        p = obj(1);
        % deactivate every rule that assigns this parameter (Target = ... form)
        for r = 1:numel(rules)
            ru = rules(r);
            lhs = regexp(ru.Rule, '^\s*([\w\.]+)\s*=', 'tokens', 'once');
            if ~isempty(lhs)
                tgt = lhs{1};
                dot = strfind(tgt, '.');
                if ~isempty(dot), tgt = tgt(dot(end)+1:end); end
                if strcmp(tgt, nm) && ru.Active
                    ru.Active = false;
                    savedRules{end+1} = ru; %#ok<AGROW>
                end
            end
        end
        savedParams{end+1} = {p, p.ConstantValue}; %#ok<AGROW>
        p.ConstantValue = true;                     % hold at its current value
    end
    cleanupModel = onCleanup(@() i_restore(savedRules, savedParams));

    % --- doses (';'-joined dose names -> dose array) ---
    d = [];
    dn = strtrim(regexp(char(doseNames), ';', 'split'));
    dn = dn(~cellfun(@isempty, dn));
    for k = 1:numel(dn)
        dk = getdose(m, dn{k});
        if ~isempty(dk), d = [d; dk]; end %#ok<AGROW>
    end

    sd = sbiosimulate(m, cs, [], d);

    % --- read readoutState at readoutDay (nearest logged time at/after it) ---
    dnames = reshape(cellstr(sd.DataNames), 1, []);
    col = find(strcmp(dnames, char(readoutState)), 1);
    if isempty(col)
        error('sb_knockout_readout:state', 'readout state %s not found', char(readoutState));
    end
    t = sd.Time;
    idx = find(t >= readoutDay, 1);
    if isempty(idx), idx = numel(t); end
    val = sd.Data(idx, col);
end

function i_restore(savedRules, savedParams)
    for r = 1:numel(savedRules)
        try, savedRules{r}.Active = true; catch, end
    end
    for p = 1:numel(savedParams)
        try, savedParams{p}{1}.ConstantValue = savedParams{p}{2}; catch, end
    end
end
