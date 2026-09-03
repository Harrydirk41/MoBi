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
    am = sbmlimport(char(sbmlFile));

    paperSpecies = string({m.Species.Name});
    agentSpecies = string({am.Species.Name});
    immune = intersect(agentSpecies, paperSpecies);          % shared = the nodes to swap
    if isempty(immune)
        error('sb_transplant_immune:noOverlap', ...
              'no species name is shared between the agent SBML and the paper model.');
    end

    % ---- classify paper reactions: pure-immune (remove) vs shell-touching (keep) ----
    toRemove = {}; keptTouchingImmune = {};
    for i = 1:numel(m.Reactions)
        r = m.Reactions(i);
        sp = i_reactionSpecies(r);
        if isempty(sp), continue; end
        allImmune = all(ismember(sp, immune));
        if allImmune
            toRemove{end+1} = r.Name; %#ok<AGROW>
        elseif any(ismember(sp, immune))
            keptTouchingImmune{end+1} = r.Name; %#ok<AGROW>
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
    report.addReactions    = arrayfun(@(r) string(r.Reaction), am.Reactions(:).', 'uni', 1);
    report.addParameters   = string({am.Parameters.Name});
    report.uncovered       = cellstr(uncovered(:).');
    report.doubleDriven    = keptTouchingImmune;             % review these for double-counting
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
    for i = 1:numel(am.Reactions)
        r = am.Reactions(i);
        % ensure every species the agent reaction needs exists in the paper model
        for s = i_reactionSpecies(r)
            if isempty(sbioselect(m, 'Type', 'species', 'Name', char(s)))
                addspecies(m.Compartments(1), char(s));
            end
        end
        nr = addreaction(m, char(r.Reaction));
        addkineticlaw(nr, 'Unknown');
        nr.ReactionRate = char(r.ReactionRate);
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

function i_printReport(rep)
    fprintf('== sb_transplant_immune report ==\n');
    fprintf('  shared immune species (%d): %s\n', numel(rep.immuneShared), strjoin(rep.immuneShared, ', '));
    fprintf('  paper reactions total: %d\n', rep.nPaperReactions);
    fprintf('  WOULD REMOVE (pure-immune) %d: %s\n', numel(rep.removeReactions), ...
            strjoin(rep.removeReactions, ', '));
    fprintf('  WOULD ADD agent reactions: %d ; agent parameters: %d\n', ...
            numel(rep.addReactions), numel(rep.addParameters));
    if ~isempty(rep.uncovered)
        fprintf('  [!] immune species the agent does NOT drive (kept on paper dynamics): %s\n', ...
                strjoin(rep.uncovered, ', '));
    end
    if ~isempty(rep.doubleDriven)
        fprintf(['  [!] kept reactions that still touch a transplanted species ' ...
                 '(review for double-counting): %s\n'], strjoin(rep.doubleDriven, ', '));
    end
end
