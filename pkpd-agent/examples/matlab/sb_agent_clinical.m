function sb_agent_clinical(sbprojPath, agentXml, outSbproj, checkReadout, checkDose, ...
                           secFactorSpec, influxFactor)
%SB_AGENT_CLINICAL Build a clinical sbproj from the AGENT-built immune network, then
%   hand off to the existing train/test/simulate helpers.
%
%   It does the model surgery ONCE and saves the result, so the whole downstream
%   (sb_fit for training, sb_run_vpop for the held-out qualification and the therapy
%   switch) runs on the AGENT's immune mechanism wearing the paper's GIVEN clinical
%   shell (DAS28-CRP / PK / doses), which stays untouched.
%
%   Steps:
%     1. sb_load(sbprojPath)                 - the paper's given clinical model
%     2. sb_transplant_immune(agentXml,false)- swap in the agent's immune reactions
%     3. baseline sanity simulation          - the readout must be finite & sensible
%     4. sbsaveproject(outSbproj, sbmodel)   - the agent-based clinical sbproj
%
%   Then run YOUR usual train / test / simulate on outSbproj, e.g. (from Python via the
%   existing engine, or directly):
%     TRAIN     : sb_fit('kg_...,lo,hi,log; ...', 'trial_das28.csv', ...
%                        'DAS28_CRP = das28', 'lsqnonlin', '<doses>', 'fit_out.csv')
%     TEST      : sb_run_vpop('Vpop1.xlsx', '<first-line doses>', 285, 199, 284, ...
%                             'vpop_first.csv', 300, '', '')
%     SIMULATE  : sb_run_vpop('Vpop1.xlsx', '<MTX;TCZ switch doses>', 601, 199, 600, ...
%                             'vpop_second.csv', 300, '', '')
%     COMPARE   : sb_paper_compare(...)      - agent-model response vs the trial
%
%   HONEST NOTE: this pairing (agent immune mechanism + given clinical shell) is the
%   real from-scratch->clinical test. Because the agent network still over-includes
%   edges (lower precision), expect the clinical fit to be WORSE than the paper's own
%   model - that gap is the point of the experiment, not a bug. Written without a
%   MATLAB instance to test on: run step 3's baseline sim and read the transplant
%   report before trusting any number.

    if nargin < 4 || isempty(checkReadout), checkReadout = 'DAS28_CRP'; end
    if nargin < 5, checkDose = ''; end
    if nargin < 6, secFactorSpec = ''; end                   % optional drug re-wire (e.g. MTX)
    if nargin < 7, influxFactor = ''; end
    if nargin < 3 || isempty(outSbproj), outSbproj = 'agent_clinical.sbproj'; end

    fprintf('== 1. load the paper clinical model ==\n');
    sb_load(sbprojPath);

    fprintf('\n== 2. transplant the agent immune network (dry run first) ==\n');
    sb_transplant_immune(agentXml, true, secFactorSpec, influxFactor);      % report only
    fprintf('\n--- applying the transplant ---\n');
    if ~isempty(secFactorSpec) || ~isempty(influxFactor)
        fprintf('  re-wiring a drug PD onto the agent secretion/influx reactions\n');
    end
    sb_transplant_immune(agentXml, false, secFactorSpec, influxFactor);

    fprintf('\n== 3. baseline sanity simulation ==\n');
    m = evalin('base', 'sbmodel');
    cs = getconfigset(m);
    try, cs.RuntimeOptions.StatesToLog = 'all'; catch, end
    % Diagnostic across the clinical protocol's horizons. The response READOUTS (ACR_Perc, the ACR/
    % Remission/NonResp flags, DAS28_BL) are set by EVENTS at day 199/284/600, so they are
    % legitimately non-finite/zero before their event fires - they are excluded from the
    % divergence check, which watches only the CONTINUOUS states (immune species + DAS28_CRP).
    readoutPat = '^ACR|Remission|Response|NonResp|_BL$|^delta_';
    for T = [50, 210, 285]
        cs.StopTime = T;
        try
            if ~isempty(checkDose)
                sd = sbiosimulate(m, cs, [], i_getdoses(m, checkDose));
            else
                sd = sbiosimulate(m, cs);
            end
        catch e
            fprintf('  t=%-4g simulation FAILED: %s\n', T, e.message); break;
        end
        names = reshape(cellstr(sd.DataNames), 1, []);
        last = sd.Data(end, :);
        isReadout = ~cellfun('isempty', regexp(names, readoutPat, 'once'));
        bad = names(~isfinite(last) & ~isReadout);           % non-finite CONTINUOUS states
        col = find(strcmp(names, checkReadout), 1);
        v = NaN; if ~isempty(col), v = last(col); end
        fprintf('  t=%-4g  %s=%g (finite=%d)   diverged continuous states: %d %s\n', T, ...
                checkReadout, v, isfinite(v), numel(bad), i_head(bad, 6));
        if ~isempty(bad), break; end
    end

    fprintf('\n== 4. save the agent-based clinical sbproj ==\n');
    sbmodel = m; %#ok<NASGU>   % sbiosaveproject reads the caller workspace by variable name
    sbiosaveproject(char(outSbproj), 'sbmodel');
    fprintf('  saved -> %s\n', outSbproj);
    fprintf(['\nNext: run sb_fit (train) / sb_run_vpop (test + simulate switch) / ' ...
             'sb_paper_compare on this project, exactly as for the paper model.\n']);
end

function s = i_head(c, n)
%I_HEAD first n entries of a cellstr joined by commas (for a compact NaN-species preview).
    if isempty(c), s = ''; return; end
    k = min(n, numel(c));
    s = ['[' strjoin(c(1:k), ', ')];
    if numel(c) > k, s = [s ', ...']; end
    s = [s ']'];
end

function d = i_getdoses(m, doseSpec)
%I_GETDOSES resolve a ';'-joined list of dose names to a SimBiology dose array.
    parts = strsplit(char(doseSpec), ';');
    d = [];
    for i = 1:numel(parts)
        nm = strtrim(parts{i});
        if isempty(nm), continue; end
        dd = sbioselect(m, 'Type', 'dose', 'Name', nm);
        if isempty(dd), error('dose "%s" not found in the model.', nm); end
        d = [d, dd]; %#ok<AGROW>
    end
end
