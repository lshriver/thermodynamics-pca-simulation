function [hArrow,hTxt] = drawLoadings3D(coeff,score,pcIdx,labels,varargin)
%DRAWLOADINGS3D  Plot loading vectors of selected principal components.
%
%   utils.drawLoadings3D(coeff,score)
%   utils.drawLoadings3D(coeff,score,pcIdx,labels)
%   utils.drawLoadings3D(___,Name,Value)
%
%   INPUTS
%     coeff    p-by-m  PCA coefficients  (columns = PCs)
%     score    n-by-m  PCA scores        (used only for scaling)
%     pcIdx    [1×3]   which PCs to show             (default = 1:3)
%     labels   {p×1}   cellstr / string array labels (default = 'Var#')
%
%   Name–Value pairs are forwarded to both QUIVER3 and TEXT.
%
%   OUTPUTS
%     hArrow   p×1 vector of quiver handles
%     hTxt     p×1 vector of text   handles
%
%   Author:  Your Name – yyyy-mm
% -------------------------------------------------------------------------

    if nargin < 3 || isempty(pcIdx), pcIdx = 1:3; end
    if nargin < 4 || isempty(labels)
        labels = compose("Var%d", 1:size(coeff,1));
    end

    ax = newplot;           % current axes (creates one if none)
    hold(ax,'on');

    coeff     = coeff(:,pcIdx);                 % keep wanted PCs
    scale     = 0.4 * max(abs(score(:)));       % nice heuristic scale
    coeff     = coeff * scale;                  % rescale vectors
    nVar      = size(coeff,1);

    extraArgs = varargin;                       % *** key line! ***
    hArrow    = gobjects(nVar,1);
    hTxt      = gobjects(nVar,1);

    for k = 1:nVar
        hArrow(k) = quiver3(ax, 0,0,0, ...
                            coeff(k,1), coeff(k,2), coeff(k,3), ...
                            'LineWidth', 2, 'Color', [0 0 0], ...
                            extraArgs{:});

        hTxt(k)   = text(ax, coeff(k,1), coeff(k,2), coeff(k,3), ...
                         labels{k}, ...
                         'Interpreter','latex', ...
                         'HorizontalAlignment','center', ...
                         'FontSize', 11, ...
                         extraArgs{:});
    end
end
