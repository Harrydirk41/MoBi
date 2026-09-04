function report = sb_drug_mechanism(drugPattern)
%SB_DRUG_MECHANISM Find HOW a drug acts in the loaded model (base 'sbmodel'): the reactions and
%   rules whose rate law / expression references the drug's species or parameters. Lets us see
%   whether a drug's pharmacodynamic effect rides on IMMUNE reactions (which the transplant removes,
%   killing the effect) or on a species-level binding a transplanted network still connects to.
%
%   drugPattern: a regexp matched against species AND parameter names, e.g. 'MTX', 'TCZ', 'ADA'.
%   Returns report with: drugQuantities (matching species/params), reactionsAffected (name +
%   whether the rate law also touches an immune species), rulesAffected.

    m = evalin('base', 'sbmodel');
    spNames  = string({m.Species.Name});
    prNames  = string({m.Parameters.Name});
    isDrug = @(nm) ~isempty(regexp(char(nm), char(drugPattern), 'once'));
    drugSp = spNames(arrayfun(isDrug, spNames));
    drugPr = prNames(arrayfun(isDrug, prNames));
    quants = [drugSp, drugPr];

    reacts = {};
    for i = 1:numel(m.Reactions)
        r = m.Reactions(i);
        rate = char(r.ReactionRate);
        if any(arrayfun(@(q) i_wordIn(rate, char(q)), quants))
            immune = any(arrayfun(@(s) i_wordIn(rate, char(s)), spNames));
            reacts{end+1} = sprintf('%s%s : %s', r.Name, ...
                ternary(immune, ' [touches a species]', ''), char(r.Reaction)); %#ok<AGROW>
        end
    end

    rules = {};
    for i = 1:numel(m.Rules)
        ex = char(m.Rules(i).Rule);
        if any(arrayfun(@(q) i_wordIn(ex, char(q)), quants))
            rules{end+1} = ex; %#ok<AGROW>
        end
    end

    report = struct();
    report.drugQuantities = cellstr(quants(:).');
    report.reactionsAffected = reacts;
    report.rulesAffected = rules;

    fprintf('== drug mechanism for /%s/ ==\n', char(drugPattern));
    fprintf('  drug species/params (%d): %s\n', numel(quants), strjoin(cellstr(quants), ', '));
    fprintf('  reactions whose rate law references the drug (%d):\n', numel(reacts));
    for i = 1:numel(reacts), fprintf('     %s\n', reacts{i}); end
    fprintf('  rules referencing the drug (%d):\n', numel(rules));
    for i = 1:numel(rules), fprintf('     %s\n', rules{i}); end
end

function out = ternary(c, a, b)
    if c, out = a; else, out = b; end
end

function tf = i_wordIn(str, word)
    tf = ~isempty(regexp(str, ['(?<![A-Za-z0-9_.])' regexptranslate('escape', word) ...
                               '(?![A-Za-z0-9_.])'], 'once'));
end
