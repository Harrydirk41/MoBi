function [t, c] = sb_smoke()
%SB_SMOKE Build a trivial 1-compartment IV-bolus PK model in SimBiology,
%   simulate it, and return time (t) and drug concentration (c). Used to prove
%   the Python -> MATLAB Engine -> SimBiology toolchain end to end, entirely in
%   code (no GUI). Returns two column vectors.

    m = sbiomodel('smoke');
    central = addcompartment(m, 'central', 1.0);            % volume 1 L
    addspecies(central, 'drug', 10.0);                      % initial amount 10 (mg)
    addparameter(m, 'ke', 0.5);                             % first-order rate (1/hr)

    r  = addreaction(m, 'drug -> null');                    % first-order elimination
    kl = addkineticlaw(r, 'MassAction');
    kl.ParameterVariableNames = {'ke'};                     % rate = ke * drug

    cs = getconfigset(m);
    cs.StopTime = 24;
    cs.SolverType = 'ode15s';

    sd = sbiosimulate(m);
    t  = sd.Time;
    dd = selectbyname(sd, 'drug');
    c  = dd.Data;
end
