function hTxt = annotateSpecies(score, feature, nLabel, ax, varargin)
%ANNOTATESPECIES  Put text labels on selected 3-D points.
%
%   h = ANNOTATESPECIES(score, feature) labels *all* points
%   using their row index (S1, S2, …).  Labels are drawn in
%   the current axes (gca).
%
%   h = ANNOTATESPECIES(score, feature, nLabel) labels only
%   the nLabel highest-feature points.
%
%   h = ANNOTATESPECIES(score, feature, nLabel, ax) draws in
%   the axes handle ax.
%
%   Extra name–value pairs after ax are passed directly to TEXT,
%   e.g. 'FontSize',12,'Color','r'.
%
%   INPUTS
%     score     n-by-3 matrix of coordinates (PC1,PC2,PC3)
%     variable  n-by-1 vector of scalar scores (larger = more important)
%     nLabel    scalar, how many points to label  (default = n)
%     ax        axes handle                       (default = gca)
%
%   OUTPUT
%     hTxt    vector of handles to the created text objects
%
%   Example
%     annotateSpecies(score_P(:,1:3), feature, 5, gca, ...
%                     'FontSize',10,'Color','k','FontWeight','bold');

    %---------------- defaults ------------------------------------------
    if nargin < 3 || isempty(nLabel), nLabel = size(score,1); end
    if nargin < 4 || isempty(ax),     ax     = gca;           end

    %---------------- find top-entropy indices --------------------------
    nLabel = min(nLabel, numel(feature));           % safety
    [~, idx] = maxk(feature(:), nLabel);            % row indices to label

    %---------------- create labels -------------------------------------
    holdState = ishold(ax);   % remember so we can restore at the end
    hold(ax,'on');

    hTxt = gobjects(nLabel,1);
    for k = 1:nLabel
        i   = idx(k);                       % species index
        txt = sprintf('S%d', i);            % label text

        hTxt(k) = text(ax, ...
            score(i,1), score(i,2), score(i,3), ...
            txt, ...
            'HorizontalAlignment','left', ...
            'VerticalAlignment','bottom', ...
            varargin{:});                   % user-supplied opts
    end

    if ~holdState, hold(ax,'off'); end
end


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