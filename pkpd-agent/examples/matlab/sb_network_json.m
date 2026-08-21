function sb_network_json(outJson)
%SB_NETWORK_JSON Dump the loaded model's full STRUCTURE (the network answer-key)
%   to JSON: every species, reaction (with its reactant/product species and the
%   ReactionRate expression), and rule (RepeatedAssignment/Rate/Algebraic, with the
%   rule expression). This is the ground-truth topology a Stage-1 reconstruction
%   benchmark scores against - the regulatory edges live both in the reaction
%   stoichiometry AND in the rule expressions (e.g. Pro_IL6Sec_byMacro_effect).
%
%   Arrays cannot cross this engine build to Python, so the structure goes via file.
%   Load the project first (sb_load), then:  sb_network_json('network.json')

    m = evalin('base', 'sbmodel');
    s = struct();
    s.name = m.Name;

    % -- species ------------------------------------------------------- %
    sp = m.Species;
    spc = cell(1, numel(sp));
    for i = 1:numel(sp)
        r = struct();
        r.name = sp(i).Name;
        try, r.compartment = sp(i).Parent.Name; catch, r.compartment = ''; end
        try, r.initial = sp(i).InitialAmount;   catch, r.initial = NaN;   end
        spc{i} = r;
    end
    s.species = spc;

    % -- reactions (stoichiometry + rate law) -------------------------- %
    rx = m.Reactions;
    rxc = cell(1, numel(rx));
    for i = 1:numel(rx)
        r = struct();
        r.name      = rx(i).Name;
        r.reaction  = rx(i).Reaction;         % 'A + B -> C' string
        r.rate      = rx(i).ReactionRate;     % kinetic law expression
        r.reactants = local_names(rx(i).Reactants);
        r.products  = local_names(rx(i).Products);
        rxc{i} = r;
    end
    s.reactions = rxc;

    % -- rules (RepeatedAssignment etc. - the regulatory expressions) --- %
    ru = m.Rules;
    ruc = cell(1, numel(ru));
    for i = 1:numel(ru)
        r = struct();
        try, r.type = ru(i).RuleType; catch, r.type = ''; end
        try, r.rule = ru(i).Rule;     catch, r.rule = ''; end
        ruc{i} = r;
    end
    s.rules = ruc;

    % -- parameters (name/value/constant) ------------------------------ %
    pp = m.Parameters;
    ppc = cell(1, numel(pp));
    for i = 1:numel(pp)
        r = struct();
        r.name = pp(i).Name;
        try, r.value = pp(i).Value;            catch, r.value = NaN;   end
        try, r.units = pp(i).ValueUnits;       catch, r.units = '';    end
        try, r.constant = logical(pp(i).ConstantValue); catch, r.constant = true; end
        ppc{i} = r;
    end
    s.parameters = ppc;

    s.counts = struct('species', numel(sp), 'reactions', numel(rx), ...
                      'rules', numel(ru), 'parameters', numel(pp));

    txt = jsonencode(s);
    fid = fopen(outJson, 'w', 'n', 'UTF-8');
    if fid < 0, error('sb_network_json:open', 'cannot open %s', outJson); end
    fprintf(fid, '%s', txt);
    fclose(fid);
    fprintf('wrote %s: %d species, %d reactions, %d rules, %d parameters\n', ...
            outJson, numel(sp), numel(rx), numel(ru), numel(pp));
end

function names = local_names(objs)
    names = {};
    for i = 1:numel(objs)
        try, names{end+1} = objs(i).Name; catch, end %#ok<AGROW>
    end
end
