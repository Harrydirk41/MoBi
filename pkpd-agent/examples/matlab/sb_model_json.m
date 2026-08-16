function sb_model_json(outJson)
%SB_MODEL_JSON Write the loaded model's structure (names of species, parameters,
%   reactions, doses, variants) as JSON to outJson. Arrays cannot be marshaled
%   back to Python by this engine build, so structured data goes via a file.

    m = evalin('base', 'sbmodel');
    s = struct();
    s.name        = m.Name;
    s.nSpecies    = numel(m.Species);
    s.nParameters = numel(m.Parameters);
    s.nReactions  = numel(m.Reactions);
    s.species     = local_names(m.Species);
    s.parameters  = local_names(m.Parameters);
    s.reactions   = local_names(m.Reactions);
    s.doses       = local_names(getdose(m));
    s.variants    = local_names(getvariant(m));

    txt = jsonencode(s);
    fid = fopen(outJson, 'w', 'n', 'UTF-8');
    if fid < 0
        error('sb_model_json:open', 'cannot open %s for writing', outJson);
    end
    fprintf(fid, '%s', txt);
    fclose(fid);
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
