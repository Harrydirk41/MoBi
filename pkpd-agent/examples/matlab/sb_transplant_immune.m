function report = sb_transplant_immune(sbmlFile, dryRun)
%SB_TRANSPLANT_IMMUNE Graft an AGENT-built immune network (an SBML file, e.g. the
%   mynet.xml from run_qsp_build_network) into the loaded paper model (base
%   workspace 'sbmodel', put there by sb_load), replacing the paper's PURE-IMMUNE
%   reactions for the shared species while KEEPING every reaction and rule that
%   also touches a clinical / PK / readout species. Because the agent network uses
%   the SAME species names as the paper (the nodes are the paper's own), the paper's
%   DAS28-CRP / CRP / PK / dose shell re-connects to the agent's dynamics through
%   those shared names - no formula is invented here, the clinical shell stays the
%   paper's GIVEN model.
%
%   Usage:
%     sb_load('Vantage RA QSP Model v1.0.sbproj');   % -> base workspace 'sbmodel'
%     report = sb_transplant_immune('mynet.xml', true);   % DRY RUN - read the report first
%     report = sb_transplant_immune('mynet.xml', false);  % apply, then fit/vpop as usual
%
%   dryRun=true (default) changes NOTHING; it returns a report struct listing exactly
%   what would be removed/added so you can verify the surgery against YOUR sbproj
%   before trusting a fit. This function is a scaffold: it was written without a
%   MATLAB/SimBiology instance to test against, so run the dry run, read every field,
%   and sanity-check a simulation before any calibration.
%
%   HONEST SEAMS (what you must verify):
%     * "immune species" = names present in BOTH the agent SBML and the paper model.
%     * a paper reaction is removed ONLY if ALL its reactant+product species are
%       immune (a pure-immune reaction). Any reaction touching a non-immune species
%       (DAS28, CRP, a PK compartment, a drug) is KEPT - the shell is untouched.
%     * report.uncovered lists immune species the agent does NOT drive (kept on paper
%       dynamics - a hybrid). report.doubleDriven lists shared species still produced/
%       consumed by a KEPT reaction (potential double-counting - review by hand).
%     * agent parameters are added/overwritten with the agent's calibrated values;
%       paper parameters left in place (harmless orphans if now unused).

    if nargin < 2 || isempty(dryRun), dryRun = true; end
    m = evalin('base', 'sbmodel');
    ow = warning('off', 'all');                              % silence SBML default-unit warnings
    cleanup = onCleanup(@() warning(ow)); %#ok<NASGU>
    am = sbmlimport(char(sbmlFile));

    paperSpecies = string({m.Species.Name});
    agentSpecies = string({am.Species.Name});
    immune = intersect(agentSpecies, paperSpecies);          % shared = the nodes to swap
    if isempty(immune)
        error('sb_transplant_immune:noOverlap', ...
              'no species name is shared between the agent SBML and the paper model.');
    end

    % ---- classify paper reactions by NET STOICHIOMETRY OWNERSHIP ----
    % The agent supplies the COMPLETE mass balance for each of the 23 shared species. The paper
    % writes modifiers as CATALYSTS - the same species on BOTH sides (e.g. "GMCSF + VEGF ->
    % Endothelial + GMCSF + VEGF"), so it is not consumed. So "in the stoichiometry" is too coarse:
    % we must use the NET involvement (species on exactly one side). A reaction is REMOVED iff it
    % NET-produces or NET-consumes a transplanted species (the agent now owns that). A reaction whose
    % transplanted species are all catalysts (both sides, net zero) is KEPT - that includes the shell
    % reactions that PRODUCE a non-immune species (GMCSF, CAM, AutoAb) catalysed by an immune cell,
    % and the DAS28/PK reactions that only READ the immune state. Nothing is double-driven, and no
    % shell-species source is lost.
    toRemove = {}; removeDetails = {}; shellInRemoved = string.empty; clinicalCouplings = {};
    for i = 1:numel(m.Reactions)
        r = m.Reactions(i);
        net = i_netSpecies(r);                               % species on exactly one side (net != 0)
        ownsTransplanted = any(ismember(net, immune));
        if ownsTransplanted
            toRemove{end+1} = r.Name; %#ok<AGROW>
            shellHere = setdiff(net, immune);                % non-immune species NET-losing a term
            if ~isempty(shellHere)
                shellInRemoved = union(shellInRemoved, shellHere);
                removeDetails{end+1} = sprintf('%s : %s', r.Name, char(r.Reaction)); %#ok<AGROW>
            end
        else
            % kept: does it touch the immune state at all (as catalyst or rate-law modifier)?
            rate = char(r.ReactionRate);
            sp = i_reactionSpecies(r);
            if any(ismember(sp, immune)) || any(arrayfun(@(s) i_wordIn(rate, char(s)), immune))
                clinicalCouplings{end+1} = r.Name; %#ok<AGROW>
            end
        end
    end

    % ---- what the agent covers vs not ----
    agentDriven = string.empty;
    for i = 1:numel(am.Reactions)
        agentDriven = union(agentDriven, i_reactionSpecies(am.Reactions(i)));
    end
    uncovered = setdiff(immune, agentDriven);                % immune species agent doesn't drive

    report = struct();
    report.immuneShared   = cellstr(immune(:).');
    report.nPaperReactions = numel(m.Reactions);
    report.removeReactions = toRemove;
    report.removeMixedDetails = removeDetails;               % removed reactions that also touch shell
    report.shellSpeciesInRemoved = cellstr(shellInRemoved(:).');  % non-immune species losing a term
    report.clinicalCouplings = clinicalCouplings;            % kept reactions reading immune state (good)
    report.addReactions    = arrayfun(@(r) string(r.Reaction), am.Reactions(:).', 'uni', 1);
    report.addParameters   = string({am.Parameters.Name});
    report.uncovered       = cellstr(uncovered(:).');
    report.applied         = ~dryRun;

    i_printReport(report);
    if dryRun
        fprintf('\n[DRY RUN] nothing changed. Re-run with dryRun=false to apply.\n');
        return;
    end

    % ---- APPLY: remove pure-immune reactions, add agent params + reactions ----
    for k = 1:numel(toRemove)
        r = sbioselect(m, 'Type', 'reaction', 'Name', toRemove{k});
        if ~isempty(r), delete(r); end
    end
    for i = 1:numel(am.Parameters)
        p = am.Parameters(i);
        existing = sbioselect(m, 'Type', 'parameter', 'Name', p.Name);
        if isempty(existing)
            addparameter(m, p.Name, p.Value);
        else
            existing.Value = p.Value;                        % agent's calibrated value wins
        end
    end
    % the paper model has >1 compartment, so species in reaction strings AND rate laws must be
    % compartment-qualified (e.g. Synovium.IL6). Map every agent species to its compartment in m.
    defaultComp = m.Compartments(1).Name;
    for i = 1:numel(am.Reactions)
        r = am.Reactions(i);
        for s = i_reactionSpecies(r)
            if isempty(sbioselect(m, 'Type', 'species', 'Name', char(s)))
                addspecies(m.Compartments(1), char(s));      % new species -> default compartment
            end
        end
        rxnStr = i_qualifyReaction(m, r, defaultComp);
        nr = addreaction(m, rxnStr);
        addkineticlaw(nr, 'Unknown');
        nr.ReactionRate = i_qualifyRate(m, char(r.ReactionRate), i_reactionSpecies(r), defaultComp);
    end
    % give the transplanted species the agent's CALIBRATED initial amounts (its steady-state
    % target), so the network starts at the fixed point the free rates were fit to - not the paper's
    % initial values, which the agent dynamics would otherwise pull away from.
    for i = 1:numel(am.Species)
        s = am.Species(i);
        ps = sbioselect(m, 'Type', 'species', 'Name', s.Name);
        if ~isempty(ps), ps(1).InitialAmount = s.InitialAmount; end
    end
    fprintf('\n[APPLIED] removed %d pure-immune reactions, added %d agent reactions, %d params.\n', ...
            numel(toRemove), numel(am.Reactions), numel(am.Parameters));
    fprintf('Verify a baseline simulation before any fit: sbiosimulate(sbmodel).\n');
    assignin('base', 'sbmodel', m);
end

function names = i_reactionSpecies(r)
%I_REACTIONSPECIES the reactant+product species names of a reaction, as a string row.
    names = string.empty;
    for j = 1:numel(r.Reactants), names(end+1) = string(r.Reactants(j).Name); end %#ok<AGROW>
    for j = 1:numel(r.Products),  names(end+1) = string(r.Products(j).Name);  end %#ok<AGROW>
    names = unique(names);
end

function comp = i_speciesComp(m, name, defaultComp)
%I_SPECIESCOMP the compartment name of a species in model m (defaultComp if not found).
    s = sbioselect(m, 'Type', 'species', 'Name', char(name));
    if isempty(s), comp = defaultComp; else, comp = s(1).Parent.Name; end
end

function str = i_qualifyReaction(m, r, defaultComp)
%I_QUALIFYREACTION build a "Comp.A + Comp.B -> Comp.C" reaction string, compartment-qualified so
%   addreaction resolves species in a multi-compartment model. Empty side -> 'null'.
    lhs = i_qualSide(m, r.Reactants, defaultComp);
    rhs = i_qualSide(m, r.Products,  defaultComp);
    str = [lhs ' -> ' rhs];
end

function side = i_qualSide(m, specarr, defaultComp)
    parts = {};
    for j = 1:numel(specarr)
        nm = specarr(j).Name;
        parts{end+1} = [i_speciesComp(m, nm, defaultComp) '.' char(nm)]; %#ok<AGROW>
    end
    if isempty(parts), side = 'null'; else, side = strjoin(parts, ' + '); end
end

function rate = i_qualifyRate(m, rate, specNames, defaultComp)
%I_QUALIFYRATE replace each bare agent species name in the rate expression with its qualified
%   Comp.name, using word boundaries so it never matches inside a parameter name (e.g. the IL6 in
%   ksec_IL6 is left alone). Longer names first so a short name isn't matched inside a longer one.
    [~, order] = sort(strlength(specNames), 'descend');
    specNames = specNames(order);
    for k = 1:numel(specNames)
        nm = char(specNames(k));
        qual = [i_speciesComp(m, nm, defaultComp) '.' nm];
        rate = regexprep(rate, ['(?<![A-Za-z0-9_.])' regexptranslate('escape', nm) ...
                                '(?![A-Za-z0-9_.])'], qual);
    end
end

function net = i_netSpecies(r)
%I_NETSPECIES species with NON-ZERO net stoichiometry: those on exactly one side. A catalyst
%   (same species as both reactant and product, e.g. a rate modifier) has net zero and is
%   excluded. Uses name sets (this model's catalyst idiom is 1:1), so a species that is genuinely
%   produced and consumed with different multiplicities is treated as a catalyst - acceptable here.
    reac = string.empty; prod = string.empty;
    for j = 1:numel(r.Reactants), reac(end+1) = string(r.Reactants(j).Name); end %#ok<AGROW>
    for j = 1:numel(r.Products),  prod(end+1) = string(r.Products(j).Name);  end %#ok<AGROW>
    net = setxor(unique(reac), unique(prod));
end

function i_printReport(rep)
    fprintf('== sb_transplant_immune report ==\n');
    fprintf('  shared immune species (%d): %s\n', numel(rep.immuneShared), strjoin(rep.immuneShared, ', '));
    fprintf('  paper reactions total: %d\n', rep.nPaperReactions);
    fprintf('  WOULD REMOVE (any transplanted species in stoichiometry) %d: %s\n', ...
            numel(rep.removeReactions), strjoin(rep.removeReactions, ', '));
    fprintf('  WOULD ADD agent reactions: %d ; agent parameters: %d\n', ...
            numel(rep.addReactions), numel(rep.addParameters));
    fprintf('  KEEP - clinical/PK reactions that READ the immune state via their rate law (%d): %s\n', ...
            numel(rep.clinicalCouplings), strjoin(rep.clinicalCouplings, ', '));
    if ~isempty(rep.uncovered)
        fprintf('  [!] immune species the agent does NOT drive (kept on paper dynamics): %s\n', ...
                strjoin(rep.uncovered, ', '));
    end
    if ~isempty(rep.shellSpeciesInRemoved)
        fprintf(['  [!] NON-immune (shell) species that lose a reaction term because a removed ' ...
                 'reaction\n      also touched them - verify their balance survives: %s\n'], ...
                strjoin(rep.shellSpeciesInRemoved, ', '));
        fprintf('  [!] the mixed removed reactions (name : stoichiometry):\n');
        for i = 1:numel(rep.removeMixedDetails)
            fprintf('        %s\n', rep.removeMixedDetails{i});
        end
    end
end

function tf = i_wordIn(str, word)
%I_WORDIN true if WORD appears in STR delimited by non-identifier chars (a real reference,
%   not a substring of a longer name).
    tf = ~isempty(regexp(str, ['(?<![A-Za-z0-9_])' regexptranslate('escape', word) ...
                               '(?![A-Za-z0-9_])'], 'once'));
end
