function sb_add_drug(targetParam, efficacy, startDay)
%SB_ADD_DRUG Add a designed anti-cytokine biologic to the loaded model. At startDay
%   the drug suppresses the chosen disease-driver parameter to (1-efficacy) of its
%   baseline value, modeling an agent that inhibits that cytokine pathway (efficacy
%   0 = no drug, 1 = full blockade). This is the Stage-2 structural edit: it adds an
%   event to 'sbmodel' in the base workspace. Reload the project (sb_load) to reset
%   before designing a different drug.
%
%   targetParam : a disease-driver parameter name, e.g. 'F_IL6', 'F_TNFa', 'F_IL17'.
%   efficacy    : fractional suppression in [0,1] (the drug's effect size).
%   startDay    : day the drug takes effect (e.g. 200, treatment start).

    m = evalin('base', 'sbmodel');
    p = sbioselect(m, 'Name', targetParam);
    if isempty(p)
        error('sb_add_drug:notfound', 'parameter "%s" not in model', targetParam);
    end
    if ~isa(p(1), 'SimBiology.Parameter')
        error('sb_add_drug:notparam', '"%s" is not a parameter', targetParam);
    end
    base = p(1).Value;
    supp = base * (1 - efficacy);
    % the parameter must be non-constant for an event to reassign it
    set(p(1), 'ConstantValue', false);
    ev = addevent(m, sprintf('time >= %g', startDay), ...
                  sprintf('%s = %.12g', targetParam, supp));
    ev.Active = true;
    fprintf('designed drug: suppress %s  %.6g -> %.6g  (efficacy %g) from day %g\n', ...
            targetParam, base, supp, efficacy, startDay);
    assignin('base', 'sbmodel', m);
end
