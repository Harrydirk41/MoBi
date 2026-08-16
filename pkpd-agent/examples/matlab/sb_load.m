function sb_load(projPath)
%SB_LOAD Load a SimBiology project headlessly and stash the model in the base
%   workspace as 'sbmodel' for later operations. Also puts the project folder on
%   the path so model-referenced functions (e.g. the RA model's MM.m) resolve.

    projDir = fileparts(projPath);
    if ~isempty(projDir)
        addpath(projDir);
    end
    data = sbioloadproject(projPath);
    m = [];
    f = fieldnames(data);
    for i = 1:numel(f)
        if isa(data.(f{i}), 'SimBiology.Model')
            m = data.(f{i});
            break;
        end
    end
    if isempty(m)
        error('sb_load:noModel', 'No SimBiology.Model in %s (fields: %s)', ...
              projPath, strjoin(f, ', '));
    end
    assignin('base', 'sbmodel', m);
end
