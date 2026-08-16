function info = sb_project_info(projPath)
%SB_PROJECT_INFO Load a SimBiology project (.sbproj) headlessly and report its
%   structure (species / parameters / reactions / variants / doses), then run a
%   baseline simulation. Proves an existing QSP model (e.g. the Vantage RA model)
%   can be loaded and driven from code, no GUI. Returns a struct summary.

    % the model's reactions may call custom functions shipped alongside the
    % project (e.g. the Vantage RA model's MM.m). Put the project folder on the
    % path so those resolve during load + simulation.
    projDir = fileparts(projPath);
    if ~isempty(projDir)
        addpath(projDir);
    end

    data = sbioloadproject(projPath);      % loads project variables into a struct

    % find the SimBiology model among the loaded variables (its variable name
    % varies by how the project was saved)
    m = [];
    f = fieldnames(data);
    for i = 1:numel(f)
        v = data.(f{i});
        if isa(v, 'SimBiology.Model')
            m = v;
            break;
        end
    end
    if isempty(m)
        error('sb_project_info:noModel', ...
              'No SimBiology.Model found in %s (fields: %s)', ...
              projPath, strjoin(f, ', '));
    end

    variants = getvariant(m);
    doses    = getdose(m);

    info = struct();
    info.name        = m.Name;
    info.nSpecies    = numel(m.Species);
    info.nParameters = numel(m.Parameters);
    info.nReactions  = numel(m.Reactions);
    info.nVariants   = numel(variants);
    info.nDoses      = numel(doses);
    info.variantNames = local_names(variants);
    info.doseNames    = local_names(doses);

    fprintf('model: %s\n', info.name);
    fprintf('  species: %d  parameters: %d  reactions: %d\n', ...
            info.nSpecies, info.nParameters, info.nReactions);
    fprintf('  variants: %d  doses: %d\n', info.nVariants, info.nDoses);
    local_print_names('  doses', info.doseNames, 20);
    local_print_names('  variants (first 8)', info.variantNames, 8);

    % baseline simulation smoke (no dose/variant) - just confirm it integrates
    try
        sd = sbiosimulate(m);
        info.baselineOk = true;
        info.baselineStopTime = sd.Time(end);
        fprintf('  baseline simulation: OK (t_end = %g, %d states logged)\n', ...
                sd.Time(end), size(sd.Data, 2));
    catch ME
        info.baselineOk = false;
        info.baselineError = ME.message;
        fprintf('  baseline simulation FAILED: %s\n', ME.message);
    end
end

function names = local_names(objs)
    names = {};
    for i = 1:numel(objs)
        try
            names{end+1} = objs(i).Name; %#ok<AGROW>
        catch
        end
    end
end

function local_print_names(label, names, maxn)
    if isempty(names)
        fprintf('%s: (none)\n', label);
        return;
    end
    n = min(numel(names), maxn);
    fprintf('%s: %s\n', label, strjoin(names(1:n), ' | '));
end
