function nOk = sb_fit(paramSpec, dataCsv, responseMap, method, doseNames, outCsv)
%SB_FIT Multi-parameter calibration using SimBiology's NATIVE estimator (sbiofit).
%   The numerical optimization runs INSIDE MATLAB (Optimization / Global Optimization
%   Toolbox via sbiofit), not over the Python bridge - so a joint fit of many
%   parameters to a data curve is one native call, fast and battle-tested. The agent's
%   job shrinks to SETTING UP the problem (which parameters, bounds, which data,
%   which method) and INTERPRETING the result; the heavy lifting is MATLAB's.
%
%   Inputs (all strings except where noted; the Python engine builds them):
%     paramSpec   ';'-separated "name,lo,hi,scale" (scale 'lin'/'log') - the parameters
%                 to estimate and their bounds. log => the estimate is on a log scale.
%     dataCsv     path to a CSV of observed data with columns: Time, then one column per
%                 observed model component (the response). One group (single patient/arm).
%     responseMap ';'-separated "ModelComponent = DataColumn" (e.g. "DAS28_CRP = das28").
%     method      sbiofit method: 'lsqnonlin' (default), 'fmincon', 'particleswarm',
%                 'ga', 'nlinfit', ... - the native optimizer to use.
%     doseNames   optional ';'-joined dose names applied during the fit (or '').
%     outCsv      where to write the result table: parameter,estimate,lo,hi.
%
%   Writes outCsv and returns 1 on success, 0 on failure. Requires the Optimization
%   Toolbox (and Global Optimization Toolbox for 'ga'/'particleswarm').

    nOk = 0;
    origWarn = warning('off', 'all');
    cleanupWarn = onCleanup(@() warning(origWarn)); %#ok<NASGU>

    m = evalin('base', 'sbmodel');
    if nargin < 4 || isempty(method), method = 'lsqnonlin'; end
    if nargin < 5, doseNames = ''; end

    % -- parameters to estimate -> array of estimatedInfo ------------------- %
    entries = strsplit(string(paramSpec), ';');
    est = [];
    for k = 1:numel(entries)
        e = strtrim(char(entries(k)));
        if isempty(e), continue; end
        parts = strsplit(e, ',');
        nm = strtrim(char(parts(1)));
        if isempty(sbioselect(m, 'Name', nm))
            fprintf('WARNING: parameter "%s" not in model - skipped\n', nm); continue;
        end
        lo = str2double(parts(2)); hi = str2double(parts(3));
        isLog = numel(parts) >= 4 && strcmpi(strtrim(char(parts(4))), 'log');
        ei = estimatedInfo(nm, 'Bounds', [lo hi]);
        if isLog, ei.Transform = 'log'; end
        est = [est ei]; %#ok<AGROW>
    end
    if isempty(est), error('sb_fit:noparams', 'no valid parameters to estimate'); end

    % -- observed data -> groupedData -------------------------------------- %
    T = readtable(dataCsv);
    gd = groupedData(T);

    % -- response map "ModelComponent = DataColumn" ------------------------ %
    rmEntries = strsplit(string(responseMap), ';');
    rmap = {};
    for k = 1:numel(rmEntries)
        e = strtrim(char(rmEntries(k)));
        if isempty(e), continue; end
        kv = strsplit(e, '=');
        rmap{end+1} = sprintf('%s = %s', strtrim(char(kv(1))), strtrim(char(kv(2)))); %#ok<AGROW>
    end

    % -- doses applied during the fit (optional) --------------------------- %
    doseObjs = [];
    if ~isempty(strtrim(char(string(doseNames))))
        dn = strsplit(string(doseNames), ';');
        for k = 1:numel(dn)
            nm = strtrim(char(dn(k)));
            if isempty(nm), continue; end
            d = sbioselect(m, 'Name', nm);
            if ~isempty(d), doseObjs = [doseObjs d(1)]; end %#ok<AGROW>
        end
    end

    fprintf('sb_fit: estimating %d parameter(s) with %s ...\n', numel(est), method);
    if isempty(doseObjs)
        res = sbiofit(m, gd, rmap, est, [], method);
    else
        res = sbiofit(m, gd, rmap, est, doseObjs, method);
    end

    % -- write the estimates ----------------------------------------------- %
    pe = res.ParameterEstimates;      % table: Name, Estimate, StandardError, ...
    fid = fopen(outCsv, 'w', 'n', 'UTF-8');
    fprintf(fid, 'parameter,estimate\n');
    for i = 1:height(pe)
        fprintf(fid, '%s,%g\n', pe.Name{i}, pe.Estimate(i));
    end
    fclose(fid);
    fprintf('sb_fit: wrote %d estimate(s) to %s\n', height(pe), outCsv);
    nOk = 1;
end
