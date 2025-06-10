function main_demo
%MAIN_DEMO  Toy thermodynamics + PCA example.
%
%   This is a self-contained demonstration.  All physical constants are
%   hard-coded.  Helper functions are defined at the end of the file – you
%   can of course move them to their own *.m files if you wish.
%
%   OUTPUT  – the script prints two tables and opens six figures.

    %#ok<*NASGU> % silence lint for variables used only for plotting
    close all;  clear; clc;

    %% ------------------------------------------------------------------
    %  1.  User-set parameters
    % -------------------------------------------------------------------
    nSpecies   = 10;   % number of chemical species
    nLevels    = 5;    % energy levels per species
    %TEMP_K     = 300;  % temperature (K)
    THRESH_P   = 1e-6; % “inaccessible” probability threshold
    PC_TO_PLOT = 1:3;  % principal components to show in 3-D plots
    N_LABEL    = 5;    % how many species to annotate in scatter plots

    % --- Temperature selection ---
    temp_options = [100, 200, 300, 400, 500]; % Available choices (K)
    [temp_idx, ok_temp] = listdlg('PromptString', 'Select Temperature (K):', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(temp_options)));

    if ~ok_temp
        error('Temperature selection cancelled.');
    end
    TEMP_K = temp_options(temp_idx);

    % --- Color mapping variable selection ---
    color_options = {'AvgE', 'Entropy', 'Z', 'F', 'PctInaccess'};
    [color_idx, ok_color] = listdlg('PromptString', 'Select colormap variable:', ...
        'SelectionMode', 'single', ...
        'ListString', color_options);

    if ~ok_color
        error('Color variable selection cancelled.');
    end
    color_choice = color_options{color_idx};

    %% ------------------------------------------------------------------
    %  2.  Generate random energy landscape
    % -------------------------------------------------------------------
    E = rand(nSpecies, nLevels);                       % 0–1 eV uniform
    levelNames = "E" + (1:nLevels);                    % string array

    T_E = array2table(E, VariableNames = levelNames);
    T_E.Species = (1:nSpecies).';                      % add index column
    T_E = movevars(T_E, "Species", "Before", 1);
    disp("=== Random energy levels (eV) ===");  disp(T_E);

    %% ------------------------------------------------------------------
    %  3.  Boltzmann statistics
    % -------------------------------------------------------------------
    kB_eVK = 8.617E-5;                     % eV/K
    beta   = 1 / (kB_eVK * TEMP_K);        % 1/eV

    Z          = sum(exp(-beta*E), 2);     % partition function (vectorised)
    P          = exp(-beta*E) ./ Z;        % normalised Boltzmann probs
    meanE      = sum(P .* E, 2);           % ⟨E⟩_i for each species
    S          = -kB_eVK * sum(P .* log(max(P, eps)), 2);   % entropy (eV/K)
    F          = -kB_eVK * TEMP_K * log(Z);                 % Helmholtz (eV)
    pctNA      = 100 * sum(P < THRESH_P, 2) / nLevels;      % inaccessible %

    T_thermo = table((1:nSpecies).', meanE, S, Z, F, pctNA, ...
                     'VariableNames', ...
                     {'Species','AvgE','Entropy','Z','F','PctInaccess'});
    disp("=== Derived thermodynamic features ===");  disp(T_thermo);

    %% ------------------------------------------------------------------
    %  4.  Visualise probability matrix
    % -------------------------------------------------------------------
    newFigure("Boltzmann probability heat-map");
    imagesc(P);   axis square;  colormap(flipud(summer(256))); colorbar;
    set(gca,'FontSize',12,'TickDir','out');
    xlabel('$j$  (energy-level index)','Interpreter','latex');
    ylabel('$i$  (species index)',     'Interpreter','latex');
    title('Boltzmann probability matrix  $P_{ij}$','Interpreter','latex');

    %% ------------------------------------------------------------------
    %  5.  PCA on probability matrix  (compositional data)
    % -------------------------------------------------------------------
    P_clr = clr(P);                         % centred log-ratio transform
    [coeffP,scoreP,~,~,explP] = pca(zscore(P_clr)); %#ok<ASGLU>

    listLoadings(coeffP(:,1:3), levelNames, "Probabilities – loadings");

    newFigure("PCA – probabilities (3-D)");
    scatter3(scoreP(:,1),scoreP(:,2),scoreP(:,3),60,S,'filled'); hold on;
    xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
    title('PCA of Boltzmann distributions');  
    grid on;  colorbar;  colormap cool;
    annotateSpecies(scoreP(:,PC_TO_PLOT),S,N_LABEL);
    drawLoadings3D(coeffP,scoreP,PC_TO_PLOT,levelNames);
    draw3DAxes(gca);  view(135,30);  hold off;

    newFigure("PCA – probabilities (variance)");
    bar(explP(1:5));  xlabel('PC'); ylabel('% var. explained'); grid on;

    %% ------------------------------------------------------------------
    %  6.  PCA on thermodynamic feature matrix
    % -------------------------------------------------------------------
    X = zscore([meanE,S,Z,F,pctNA]);                   % standardised
    featNames = {'⟨E⟩','S','Z','F','%NA'};
    [coeffX,scoreX,~,~,explX] = pca(X);

    listLoadings(coeffX(:,1:3), featNames, "Thermo features – loadings");

    newFigure("PCA – thermo (3-D)");
    scatter3(scoreX(:,1),scoreX(:,2),scoreX(:,3),60,S,'filled'); hold on;
    xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
    title('PCA of thermodynamic features'); 
    grid on; colorbar;  colormap cool;
    annotateSpecies(scoreX(:,PC_TO_PLOT),S,N_LABEL);
    drawLoadings3D(coeffX,scoreX,PC_TO_PLOT,featNames);
    draw3DAxes(gca);  view(135,30);  hold off;

    newFigure("PCA – thermo (variance)");
    bar(explX(1:5)); xlabel('PC'); ylabel('% var. explained'); grid on;


    %% ------------------------------------------------------------------
    %  7.  Correlation heatmaps
    % -------------------------------------------------------------------
    plotCorrelationHeatmap(scoreX, X, featNames, "X features vs PCA scores");
    plotCorrelationHeatmap(scoreP, P_clr, levelNames, "P(clr) vs PCA scores");

end

% ---------- centred-log-ratio transform --------------------------------
function Xclr = clr(X)
    Xsafe = X + eps;                      % avoid log(0)
    gm    = geomean(Xsafe,2);
    Xclr  = log(Xsafe ./ gm);
end

% ---------- thin wrapper for invisible figures (CI-friendly) -----------
function fh = newFigure(name)
    isCI = ismember(getenv('CI'), {'1','true','TRUE'});
    fh   = figure('Name',name,'Visible', ternary(isCI,'off','on'));
end
function out = ternary(cond,a,b),  out = b; if cond, out = a; end, end

% ---------- pretty loading-table to Command Window ---------------------
function listLoadings(coeff, rowNames, hdr)
    fprintf('\n--- %s ---\n', hdr);
    disp(array2table(coeff, RowNames=rowNames, ...
                     VariableNames={'PC1','PC2','PC3'}));
end

% ---------- drawLoadings3D (same API, simplified) ----------------------
function [hArrow,hTxt] = drawLoadings3D(coeff,score,pc,labels,varargin)
%DRAWLOADINGS3D  Quiver/label principal-component loading vectors in 3-D.
    if nargin<3||isempty(pc), pc=1:3; end
    if nargin<4||isempty(labels)
        labels = "Var"+(1:size(coeff,1));
    end
    ax = newplot;  hold on;
    scale = 0.4*max(abs(score(:)));
    nVar   = size(coeff,1);               % number of loading vectors
    hArrow = quiver3( zeros(nVar,1), ...  % X
                  zeros(nVar,1),     ...  % Y  ← fixed
                  zeros(nVar,1),     ...  % Z
                  coeff(:,1), coeff(:,2), coeff(:,3), ...
                  'LineWidth',2, varargin{:}, 'Color', 'k');
    hTxt   = text(coeff(:,1),coeff(:,2),coeff(:,3),labels,...
                  'Interpreter','latex','FontSize',11,'Color',[0 0 0],varargin{:});
end

% ---------- annotateSpecies --------------------------------------------
function hTxt = annotateSpecies(V,metric,n,ax)
% annotateSpecies(V,metric,n,ax) – label the n highest-metric points.
    if nargin<4||isempty(ax), ax = gca; end
    if nargin<3||isempty(n),  n = size(V,1); end
    [~,idx] = maxk(metric(:),n);
    hold(ax,'on');
    hTxt = text(ax,V(idx,1),V(idx,2),V(idx,3), compose('S%d',idx), ...
                'FontWeight','bold','FontSize',9,'Color','k');
end

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

% ---------- Correlation Heatmap --------------------------------------------------
function plotCorrelationHeatmap(scores, data, varNames, figTitle)
    nPCs = size(scores, 2);
    nVars = size(data, 2);
    corrMatrix = zeros(nVars, nPCs);
    for i = 1:nVars
        for j = 1:nPCs
            corrMatrix(i, j) = corr(scores(:, j), data(:, i));
        end
    end
    newFigure(figTitle);
    imagesc(corrMatrix);
    colormap(parula); colorbar; axis square;
    xticks(1:nPCs); yticks(1:nVars);
    xticklabels("PC" + (1:nPCs));
    yticklabels(varNames);
    set(gca,'FontSize',12,'TickDir','out');
    title(figTitle, 'Interpreter', 'none');
end
