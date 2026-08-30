function val = sb_perturb_response(species, highValue, readoutState, readoutDay, stopTime, ...
                                   clampSpecies, clampValue)
%SB_PERTURB_RESPONSE Isolating single-species perturbation: clamp one species (a cytokine)
%   at highValue, simulate, and return readoutState (a cell) at readoutDay. Everything else
%   stays at the model's current parameter values, EXCEPT any species named in clampSpecies
%   (a ';'-joined list), which are held at clampValue - use this to DECOUPLE the experiment
%   (clamp the other cytokines to ~0, reproducing an in-vitro dish of the cell + only this one
%   cytokine, so the response isolates this cytokine's regulator with no network feedback).
%   Clamping is by BoundaryCondition; the model is restored before returning.
%
%   clampSpecies/clampValue are optional (nargin < 6 -> no extra clamping).

    if nargin < 6, clampSpecies = ''; end
    if nargin < 7, clampValue = 0; end

    m  = evalin('base', 'sbmodel');
    cs = getconfigset(m);
    try, cs.RuntimeOptions.StatesToLog = 'all'; catch, end
    if stopTime > 0, cs.StopTime = stopTime; end
    origWarn = warning('off', 'all');
    cleanupWarn = onCleanup(@() warning(origWarn)); %#ok<NASGU>

    % list of (species, value) to clamp: the elevated target + any decoupling species
    targets = [{char(species), double(highValue)}];
    others = strtrim(regexp(char(clampSpecies), ';', 'split'));
    others = others(~cellfun(@isempty, others));
    for k = 1:numel(others)
        targets(end+1, :) = {others{k}, double(clampValue)}; %#ok<AGROW>
    end

    saved = {};                                    % {speciesObj, oldInit, oldBC}
    restore = onCleanup(@() i_restore(saved));
    for k = 1:size(targets, 1)
        s = sbioselect(m, 'Type', 'species', 'Name', targets{k, 1});
        if isempty(s)
            if k == 1
                error('sb_perturb_response:nospecies', 'species "%s" not found', targets{k, 1});
            end
            continue                               % a decoupling species not in model: skip
        end
        s = s(1);
        saved{end+1} = {s, s.InitialAmount, s.BoundaryCondition}; %#ok<AGROW>
        s.InitialAmount   = targets{k, 2};
        s.BoundaryCondition = true;                % hold it fixed against reaction flux
    end

    sd = sbiosimulate(m, cs);
    dnames = reshape(cellstr(sd.DataNames), 1, []);
    col = find(strcmp(dnames, char(readoutState)), 1);
    if isempty(col)
        error('sb_perturb_response:noreadout', 'readout "%s" not found', char(readoutState));
    end
    t = sd.Time;
    idx = find(t >= readoutDay, 1);
    if isempty(idx), idx = numel(t); end
    val = sd.Data(idx, col);
end

function i_restore(saved)
    for k = 1:numel(saved)
        try, saved{k}{1}.InitialAmount = saved{k}{2}; catch, end
        try, saved{k}{1}.BoundaryCondition = saved{k}{3}; catch, end
    end
end
