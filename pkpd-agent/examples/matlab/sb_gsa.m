function nOk = sb_gsa(paramSpec, observable, readoutDay, nSamples, outCsv)
%SB_GSA Global sensitivity RANKING of parameters via SimBiology's native Sobol method.
%   This is the NUMERICAL half of "select which parameters to vary": rank candidate
%   parameters by how much they drive the clinical readout, computed by MATLAB's own
%   variance-based global sensitivity analysis (sbiosobol) - not a stored list. An agent
%   then applies biological reasoning on top of this ranking to choose the varied set.
%   Nothing here is model-specific: the parameters and the observable are passed in.
%
%   Inputs (strings unless noted; the Python engine builds them):
%     paramSpec   ';'-separated "name,lo,hi" - the candidate parameters to screen and the
%                 physiological bounds to sample them over.
%     observable  the model quantity to compute sensitivity of (e.g. "DAS28_CRP").
%     readoutDay  numeric; the time at which the observable is evaluated (StopTime).
%     nSamples    numeric; Sobol sample size (default 1000).
%     outCsv      where to write "parameter,first_order,total_order", ranked by total order.
%
%   Writes outCsv and returns 1 on success. Requires SimBiology (sbiosobol).

    nOk = 0;
    origWarn = warning('off', 'all');
    cleanupWarn = onCleanup(@() warning(origWarn)); %#ok<NASGU>

    m = evalin('base', 'sbmodel');
    if nargin < 4 || isempty(nSamples), nSamples = 1000; end
    nSamples = double(nSamples);

    % -- candidate parameters + physiological bounds ------------------------- %
    entries = strsplit(string(paramSpec), ';');
    names = {}; lo = []; hi = [];
    for k = 1:numel(entries)
        e = strtrim(char(entries(k)));
        if isempty(e), continue; end
        p = strsplit(e, ',');
        nm = strtrim(char(p(1)));
        if isempty(sbioselect(m, 'Name', nm))
            fprintf('WARNING: parameter "%s" not in model - skipped\n', nm); continue;
        end
        names{end+1} = nm; %#ok<AGROW>
        lo(end+1) = str2double(p(2)); %#ok<AGROW>
        hi(end+1) = str2double(p(3)); %#ok<AGROW>
    end
    if isempty(names), error('sb_gsa:noparams', 'no valid parameters to screen'); end
    bounds = [lo(:) hi(:)];

    % evaluate the observable at the readout day
    cs = getconfigset(m);
    set(cs, 'StopTime', double(readoutDay));

    fprintf('sb_gsa: Sobol over %d parameters, %d samples, observable %s ...\n', ...
            numel(names), nSamples, char(observable));
    results = sbiosobol(m, names, {char(observable)}, ...
                        'Bounds', bounds, 'NumberSamples', nSamples);

    % -- reduce the time-resolved indices to the value at the final time ----- %
    SI = results.SobolIndices;
    fo = zeros(numel(names), 1);
    to = zeros(numel(names), 1);
    for i = 1:height(SI)
        pnm = char(string(SI.Parameter(i)));
        idx = find(strcmp(names, pnm), 1);
        if isempty(idx), continue; end
        fval = SI.FirstOrder(i, :); if iscell(fval), fval = fval{1}; end
        tval = SI.TotalOrder(i, :); if iscell(tval), tval = tval{1}; end
        fo(idx) = fval(end);
        to(idx) = tval(end);
    end

    % rank by total-order sensitivity, descending
    [~, order] = sort(to, 'descend');
    fid = fopen(outCsv, 'w', 'n', 'UTF-8');
    fprintf(fid, 'parameter,first_order,total_order\n');
    for i = order(:)'
        fprintf(fid, '%s,%g,%g\n', names{i}, fo(i), to(i));
    end
    fclose(fid);
    fprintf('sb_gsa: wrote %d parameter sensitivities to %s\n', numel(names), outCsv);
    nOk = 1;
end
