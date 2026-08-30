function val = sb_import_simulate(sbmlFile, readoutState, stopTime)
%SB_IMPORT_SIMULATE Import an SBML file (an assembled model), simulate it, and return
%   readoutState at the end of the run. Used to run a from-scratch assembled subsystem and
%   compare it to the real model. Independent of the loaded project (imports its own model).

    m = sbmlimport(char(sbmlFile));
    cs = getconfigset(m);
    if stopTime > 0, cs.StopTime = stopTime; end
    try, cs.RuntimeOptions.StatesToLog = 'all'; catch, end
    origWarn = warning('off', 'all');
    cleanupWarn = onCleanup(@() warning(origWarn)); %#ok<NASGU>

    sd = sbiosimulate(m, cs);
    names = reshape(cellstr(sd.DataNames), 1, []);
    col = find(strcmp(names, char(readoutState)), 1);
    if isempty(col)
        error('sb_import_simulate:noreadout', 'readout "%s" not found (have: %s)', ...
              char(readoutState), strjoin(names, ', '));
    end
    val = sd.Data(end, col);
end
