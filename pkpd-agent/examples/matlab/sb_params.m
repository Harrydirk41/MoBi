function nP = sb_params(outCsv)
%SB_PARAMS Enumerate EVERY parameter in the loaded model (name, value, constant) to CSV.
%   The full candidate set for driver selection, taken straight from the model itself -
%   so the agent starts from all parameters, not a pre-curated list. Model-agnostic.
%
%   Writes outCsv with columns name,value,constant and returns the parameter count.

    m = evalin('base', 'sbmodel');
    p = sbioselect(m, 'Type', 'parameter');
    fid = fopen(outCsv, 'w', 'n', 'UTF-8');
    fprintf(fid, 'name,value,constant\n');
    for i = 1:numel(p)
        nm = p(i).Name;
        v = p(i).Value;
        try
            c = double(p(i).ConstantValue);
        catch
            c = -1;                       % older releases: property absent
        end
        fprintf(fid, '"%s",%g,%d\n', nm, v, c);
    end
    fclose(fid);
    nP = numel(p);
    fprintf('sb_params: wrote %d parameters to %s\n', nP, outCsv);
end
