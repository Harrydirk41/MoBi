function nSel = sb_select_ga(cohortCsv, anchorSpec, popTarget, outCsv)
%SB_SELECT_GA Select a virtual population from a cohort with a NATIVE genetic algorithm.
%   This is the paper's Vpop method: given an already-simulated virtual cohort, choose a
%   SUBSET (include/exclude each candidate) whose aggregate statistics match several
%   clinical anchors at once - the baseline-severity distribution AND each therapy arm's
%   response rate. Unlike prevalence weighting it returns an ACTUAL population of real
%   virtual patients (a binary selection), which is what a Vpop is.
%
%   The optimisation is MATLAB-native: ga (Global Optimization Toolbox) searches over a
%   binary inclusion vector; the fitness is pure arithmetic on the pre-simulated cohort
%   (no re-simulation), so a generation is cheap. Nothing here is model-specific - the
%   anchors name columns of the cohort CSV, whatever the disease.
%
%   Inputs (strings unless noted; the Python engine builds them):
%     cohortCsv   path to a cohort CSV (from sb_cohort): a header row then one row per
%                 candidate, with a severity column and one column per arm (0/1 flag).
%     anchorSpec  ';'-separated anchors, each one of:
%                    "moment:COL:MEAN:SD"  - selected COL mean->MEAN, sd->SD (SD '' = skip)
%                    "rate:COL:TARGET"     - % of selected with COL>=0.5 -> TARGET
%     popTarget   desired population size (numeric; 0 = unconstrained, else a soft pull).
%     outCsv      where to write the SELECTED rows (same columns as the cohort).
%
%   Writes outCsv and returns the number selected. Requires the Global Optimization Tbx.

    nSel = 0;
    if nargin < 3 || isempty(popTarget), popTarget = 0; end
    popTarget = double(popTarget);

    T = readtable(cohortCsv);
    n = height(T);
    if n == 0, error('sb_select_ga:empty', 'cohort CSV has no rows'); end

    % -- parse anchors into a struct array bound to column vectors ------------ %
    A = struct('type', {}, 'col', {}, 'a', {}, 'b', {});
    entries = strsplit(string(anchorSpec), ';');
    for k = 1:numel(entries)
        e = strtrim(char(entries(k)));
        if isempty(e), continue; end
        parts = strsplit(e, ':');
        kind = strtrim(char(parts(1)));
        colnm = strtrim(char(parts(2)));
        if ~ismember(colnm, T.Properties.VariableNames)
            fprintf('WARNING: anchor column "%s" not in cohort - skipped\n', colnm); continue;
        end
        v = T.(colnm);
        if strcmpi(kind, 'moment')
            A(end+1) = struct('type','moment','col',v, ...
                'a', str2double(parts(3)), ...
                'b', (numel(parts)>=4)*str2double(parts(4)) + (numel(parts)<4)*NaN); %#ok<AGROW>
            if numel(parts) < 4 || isempty(strtrim(char(parts(4)))), A(end).b = NaN; end
        elseif strcmpi(kind, 'rate')
            A(end+1) = struct('type','rate','col', double(v>=0.5), ...
                'a', str2double(parts(3)), 'b', NaN); %#ok<AGROW>
        else
            fprintf('WARNING: unknown anchor kind "%s" - skipped\n', kind);
        end
    end
    if isempty(A), error('sb_select_ga:noanchors', 'no valid anchors'); end

    minSel = max(5, round(0.02*n));   % refuse near-empty selections

    function err = fitness(x)
        sel = x > 0.5;
        ns = sum(sel);
        if ns < minSel
            err = 1e6 * (minSel - ns + 1); return;
        end
        err = 0;
        for j = 1:numel(A)
            c = A(j).col(sel);
            if strcmp(A(j).type, 'rate')
                r = 100 * sum(c) / ns;                 % c is already 0/1
                scale = max(A(j).a, 5);
                err = err + ((r - A(j).a)/scale)^2;
            else
                m = mean(c);
                sdref = A(j).a; if sdref == 0, sdref = 1; end
                err = err + ((m - A(j).a)/sdref)^2;
                if ~isnan(A(j).b)
                    s = std(c);
                    err = err + ((s - A(j).b)/max(A(j).b,1e-6))^2;
                end
            end
        end
        if popTarget > 0
            err = err + 0.1 * ((ns - popTarget)/popTarget)^2;
        end
    end

    % -- native genetic algorithm over the binary inclusion vector ----------- %
    if isempty(which('ga'))
        error('sb_select_ga:noGA', ['ga is unavailable (Global Optimization Toolbox ' ...
            'not installed). Use the weighting selector (select_multi_anchor) instead.']);
    end
    opts = optimoptions('ga', 'PopulationType', 'bitstring', ...
                        'PopulationSize', min(200, max(50, 10*round(sqrt(n)))), ...
                        'MaxGenerations', 300, 'MaxStallGenerations', 40, ...
                        'FunctionTolerance', 1e-6, 'Display', 'off');
    [x, fval] = ga(@fitness, n, [], [], [], [], [], [], [], opts);

    sel = x > 0.5;
    nSel = sum(sel);
    writetable(T(sel,:), outCsv);
    fprintf('sb_select_ga: selected %d / %d candidates (fitness %.4g)\n', nSel, n, fval);
end
