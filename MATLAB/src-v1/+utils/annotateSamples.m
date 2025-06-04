function hTxt = annotateSamples(score, feature, nLabel, ax, varargin)
%ANNOTATESAMPLES  Put text labels on selected 3-D points.
%
%   h = ANNOTATESAMPLES(score, feature) labels *all* points
%   using their row index (N1, N2, …).  Labels are drawn in
%   the current axes (gca).
%
%   h = ANNOTATESAMPLES(score, feature, nLabel) labels only
%   the nLabel highest-feature points.
%
%   h = ANNOTATESAMPLES(score, feature, nLabel, ax) draws in
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
%     annotateSamples(score_P(:,1:3), feature, 5, gca, ...
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
        i   = idx(k);                       % sample index
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
