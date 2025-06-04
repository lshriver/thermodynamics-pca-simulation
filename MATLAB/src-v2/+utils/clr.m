% ---------- centred-log-ratio transform --------------------------------
function Xclr = clr(X)
    Xsafe = X + eps;                % avoid log(0)
    gm    = geomean(Xsafe,2);
    Xclr  = log(Xsafe ./ gm);
end
