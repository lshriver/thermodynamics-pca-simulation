% ---------- thin wrapper for invisible figures (CI-friendly) -----------
function fh = newFigure(name)
    isCI = ismember(getenv('CI'), {'1','true','TRUE'});
    fh   = figure('Name',name,'Visible', ternary(isCI,'off','on'));
end

function out = ternary(cond,a,b)
    out = b;
    if cond
        out = a;
    end
end
