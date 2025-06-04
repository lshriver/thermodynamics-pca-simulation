function [hArrow,hText] = drawLoadings3D(coeff,score,pcIdx,varargin)
    % INPUTS
        % drawLoadings3D(coeff,score) % default 1:3 PCs
        % drawLoadings3D(coeff,score,[1 2 4]) % pick PCs
        % drawLoadings3D(coeff,score,[1 2 3],ax) % give axes
        % drawLoadings3D(coeff,score,[1 2 3],labels) % give labels
        % drawLoadings3D(coeff,score,[1 2 3],ax,labels) % both
        % drawLoadings3D(coeff,score,[1 2 3],ax,'Color','r')% Name-Value pairs
        %
    % OUTPUTS
        % hArrow – N×1 quiver handles
        % hText – N×1 text handles

    % ------------------------------------------------ defaults -----------
    if nargin < 3 || isempty(pcIdx), pcIdx = 1:3; end
    validateattributes(pcIdx,{'numeric'},{'numel',3,'integer','>=',1});

    % -------------------- parse the variable inputs ---------------------
    ax = []; % default: create later if needed
    labels = []; % default labels will be generated
    nvPairs = {}; % any remaining Name-Value pairs

    if ~isempty(varargin)

        % If first extra arg is an axes handle …
        if isgraphics(varargin{1},'axes')
            ax = varargin{1};
            varargin = varargin(2:end);
        end


        % If first remaining extra arg is a cell-array of strings …
        if ~isempty(varargin) && iscell(varargin{1})
            labels   = varargin{1};
            varargin = varargin(2:end);
        end

        % anything left is treated as Name-Value for quiver/text
        nvPairs = varargin;
    end

    if isempty(ax), ax = gca; end

    % Keep only the requested principal components
    coeff = coeff(:,pcIdx);

    % ------------------------------------------------ graphics options ---
    pQuiv = {'LineWidth',2,'MaxHeadSize',0.5, 'Color', 'k'};
    pText = {'Interpreter','latex','HorizontalAlignment','center', ...
                'VerticalAlignment','middle','FontSize',12};

    pQuiv = [pQuiv nvPairs];
    pText = [pText nvPairs];

    % ------------------------------------------ scale + prepare ----
    scaleFactor = 0.4 * max(abs(score(:))); % 40 % of extremer score
    nVar = size(coeff,1);
    origin = [0 0 0];
    hArrow = gobjects(nVar,1);
    hText = gobjects(nVar,1);

    % default labels if none supplied
    if isempty(labels)
        labels = arrayfun(@(k)sprintf('$\text{feature}_{%d}$',k),1:nVar,'uni',false);
    end
    assert(numel(labels)>=nVar,'Not enough entries in ''labels''.');

    % ----------------------------------------- plotting -------
    keepHold = ishold(ax); hold(ax,'on');

    for j = 1:nVar
        vec = scaleFactor * coeff(j,:);


        hArrow(j) = quiver3(ax,origin(1),origin(2),origin(3), ...
                        vec(1),vec(2),vec(3),pQuiv{:});

        hText(j)  = text(ax,1*vec(1),1*vec(2),1*vec(3), ...
                     labels{j},pText{:});
    end

    if ~keepHold, hold(ax,'off'); end
end
