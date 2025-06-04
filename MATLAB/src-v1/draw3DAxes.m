function h = draw3DAxes(ax, color)
%DRAW3DAXES  Draw X-, Y-, and Z-axes through the origin on a 3-D plot.
%
%   h = DRAW3DAXES()                draws grey axes in the current axes.
%   h = DRAW3DAXES(ax)              draws in the specified axes object.
%   h = DRAW3DAXES(ax, color)       uses RGB color [r g b].
%
%   OUTPUT:  h – 3-element vector of line handles [hx, hy, hz].

    % Defaults -----------------------------------------------------------
    if nargin < 1 || isempty(ax),    ax = gca;                end
    if nargin < 2 || isempty(color), color = [0.3 0.3 0.3];   end

    % Ensure we don’t overwrite existing hold state
    holdState = ishold(ax);
    hold(ax, 'on');

    % Current axis limits
    xl = xlim(ax);
    yl = ylim(ax);
    zl = zlim(ax);

    % Draw lines
    hx = line(ax, [xl(1) xl(2)], [0 0],    [0 0],    'Color', color);
    hy = line(ax, [0 0],         [yl(1) yl(2)], [0 0],    'Color', color);
    hz = line(ax, [0 0],         [0 0],    [zl(1) zl(2)], 'Color', color);

    % Restore previous hold state
    if ~holdState,  hold(ax,'off');  end

    if nargout,  h = [hx hy hz];  end
end
