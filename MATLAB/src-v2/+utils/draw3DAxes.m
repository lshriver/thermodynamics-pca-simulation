% ---------- draw3DAxes --------------------------------------------------
function h = draw3DAxes(ax,c)
    if nargin<1||isempty(ax), ax=gca; end
    if nargin<2||isempty(c),  c=0.4*[1 1 1]; end
    xl = xlim(ax); yl = ylim(ax); zl = zlim(ax); hold(ax,'on');
    hx = line(ax,xl,[0 0],[0 0],'Color',c);
    hy = line(ax,[0 0],yl,[0 0],'Color',c);
    hz = line(ax,[0 0],[0 0],zl,'Color',c);
    h  = [hx hy hz];
end
