function val = sb_perturb_response(species, highValue, readoutState, readoutDay, stopTime)
%SB_PERTURB_RESPONSE Isolating single-species perturbation: clamp one species (a cytokine)
%   at highValue, simulate, and return readoutState (a cell) at readoutDay. Everything else
%   stays at the model's current parameter values. Because only one cytokine is elevated, the
%   readout's response isolates that cytokine's regulator - the experiment that pins one
%   coupled parameter. The species is clamped by BoundaryCondition (held at its value, not
%   consumed by reactions); the model is restored before returning.

    m  = evalin('base', 'sbmodel');
    cs = getconfigset(m);
    try, cs.RuntimeOptions.StatesToLog = 'all'; catch, end
    if stopTime > 0, cs.StopTime = stopTime; end
    origWarn = warning('off', 'all');
    cleanupWarn = onCleanup(@() warning(origWarn)); %#ok<NASGU>

    s = sbioselect(m, 'Type', 'species', 'Name', char(species));
    if isempty(s)
        error('sb_perturb_response:nospecies', 'species "%s" not found', char(species));
    end
    s = s(1);
    oldVal = s.InitialAmount;
    % hold the cytokine elevated: BoundaryCondition keeps it fixed against reaction flux
    oldBC = s.BoundaryCondition;
    restore = onCleanup(@() i_restore(s, oldVal, oldBC));
    s.InitialAmount   = double(highValue);
    s.BoundaryCondition = true;

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

function i_restore(s, oldVal, oldBC)
    try, s.InitialAmount = oldVal; catch, end
    try, s.BoundaryCondition = oldBC; catch, end
end
