function oldVal = sb_set_param(name, value)
%SB_SET_PARAM Set a model parameter's Value, returning its previous value.
%   A minimal, general primitive for perturb/restore workflows (e.g. the calibration
%   recover demo: perturb a rate, fit it back). Errors if the parameter is not found.

    m = evalin('base', 'sbmodel');
    p = sbioselect(m, 'Type', 'parameter', 'Name', char(name));
    if isempty(p)
        error('sb_set_param:notfound', 'parameter "%s" not in model', char(name));
    end
    oldVal = p(1).Value;
    p(1).Value = double(value);
end
