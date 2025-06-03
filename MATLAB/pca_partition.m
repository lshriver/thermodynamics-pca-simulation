clc, clear, close all

% UNITS:
%   Energy levels (E)               -> eV
%   Temperature (T)                 -> K
%   Boltzmann constant (kB)         -> 8.617e-5 eV/K
%   Partition function (Z)          -> unitless
%   Probabilities (Z)               -> unitless
%   Average energy (avgE)           -> eV
%   Entropy (S)                     -> eV/K

%% Create Data
num_species = 10;   % number of chemical species (n)
num_levels = 5;     % number of energy levels (p)

% Generate random energy levels between 0 and 1 eV
E = rand(num_species, num_levels);

% Create variable names for energy levels
energy_level_titles = strcat('E', string(1:num_levels));

% Convert the matrix to a table with specified variable names
E_table = array2table(E, VariableNames=energy_level_titles);

% Add species index as a new variable
E_table.SpeciesIndex = (1:num_species)';
E_table = movevars(E_table, "SpeciesIndex", "Before", 1);
disp(E_table);

%% Calculate Partition Function at T=300K
% Constants
kB = 8.617e-5;      % Boltzmann constant (eV/K)
T = 300;            % Temperature (K)
beta = 1 / (kB*T);  % Thermodynamic beta/coldness (1/eV)

% Boltzmann distribution
Z = zeros(num_species, 1);              % Partition function
P = zeros(num_species, num_levels);     % Boltzmann probabilities

% Compute Z and probabilities 
for i = 1:num_species
    Z(i) = sum(exp(-E(i,:)*beta));
    P(i, :) = exp(-E(i,:)*beta) / Z(i);
end

% Plot heat map for Boltzmann probabilities
figure();
imagesc(P);                     % Rows: species, Cols: energy levels
colormap(flipud(summer(100)));  % perceptually uniform
colorbar;

label_opts_axes = {'Fontsize', 14, 'Interpreter', 'latex'};
xlabel('Energy Level Index $j$', label_opts_axes{:});
ylabel('Species Index $i$', label_opts_axes{:});
title('Boltzmann Probability Matrix ($P_{i,j}$)', 'FontSize', 16, 'Interpreter', 'latex')

xticks(1:num_levels);
yticks(1:num_species);
axis square;

%% Calculate the Average Energy, Entropy, Helmholts Free Energy, and Inaccessible Energies
avgE = sum(P .* E,2);               % dot product along each row

entropy = zeros(num_species, 1);     % eV/K
for i = 1:num_species
    p_i = P(i,:);
    mask = p_i > 0;
    entropy(i) = -kB * sum(p_i(mask) .* log(p_i(mask)));
end

F = -kB * T * log(Z);   % Helmholtz free energy

% Find percent inaccessible energies
threshold = 1e-6;   % threshold for inaccessibility
percent_inaccessible = (sum(P < threshold, 2) ./ (num_levels))*100;

% Create a table for feature thermodynamic variables
species_indices = (1:num_species)';
feature_thermo_table = table(species_indices, avgE, entropy, Z, F, percent_inaccessible, ...
    'VariableNames', {'Species', 'AvgEnergy', 'Entropy', 'PartitionFunction', ...
    'HelmholtzEnergy', 'PercentInaccessible'});
disp(feature_thermo_table);

%% PCA Analysis - Boltzmann Probabilities
% Step 0 :  pseudo-count to avoid log(0)
eps      = 1e-12;
P_safe   = P + eps;               % same size, strictly positive

% Step 1 :  centred log-ratio transform
gMean    = geomean(P_safe, 2);    % geometric mean of each row
P_clr    = log( P_safe ./ gMean );% rows now sum to 0

% (optional) z-score by column – keeps the spirit of your original code
muP      = mean(P_clr, 1);
sigP     = std (P_clr, 0, 1);
P_z      = (P_clr - muP) ./ sigP;

% Step 2 :  PCA
[coeff_P, score_P, latent_P, ~, explained_P] = pca(P_z);

% Step 2: Run PCA
[coeff_P, score_P, latent_P, ~, explained_P] = pca(P_z);

% Step 3: Visualize
% --- Create Table ---
disp(array2table(coeff_P(:,1:3), 'RowNames', energy_level_titles, 'VariableNames', {'PC1', 'PC2', 'PC3'}))

% --- Create a 3D scatter plot whose axes correspopnd with first three loadings ---
figure();
scatter3(score_P(:,1), score_P(:,2), score_P(:,3), 60, entropy, 'filled');

label_opts_axes = {'Fontsize', 14, 'Interpreter', 'latex'};
xlabel('PC 1', label_opts_axes{:});
ylabel('PC 2', label_opts_axes{:});
zlabel('PC 3', label_opts_axes{:});
title('PCA of Boltzmann Distributions (Colored by Entropy)', 'FontSize', 16, 'Interpreter', 'latex');
colorbar; colormap cool;
grid on; hold on;

label_opts_data = {'FontSize', 10, 'Color', 'k', 'FontWeight', 'bold'}
annotateSpecies(score_P(:,1:3), entropy, 5, gca, label_opts_data{:});

energy_variable_names = arrayfun(@(k)sprintf('$E_{%d}$',k),1:num_levels,'uni',false);
drawLoadings3D(coeff_P, score_P, [1 2 3], energy_variable_names);

draw3DAxes(gca) % draw x, y, and z axes
view(135, 30);  % set a nice 3D viewing angle
hold off

% --- Create a bar plot of the explained variances for the first 5
% principal components ---
figure();
bar(explained_P(1:5));
xlabel('Principal Components', label_opts_axes{:});
ylabel('Variance Explained', label_opts_axes{:});
title('Explained Variance by Principal Components - Boltzmann Probabilities', ...
    'FontSize', 16, 'Interpreter', 'latex');
grid on;


%% PCA Analysis - Thermodynamic Feature Matrix
% Step 1: Build the Feature Matrix
% --- Raw feature matrix  % Rows: species, Cols: features ---
X = [avgE, entropy, Z, F, percent_inaccessible];

% Step 2: Normalize the Data (subtract mean and divide by standard
% deviation -- i.e., standardization -- same as last time)
mu = mean(X, 1);
X_centered = X - mu;
X_std = std(X, 0, 1);   % sample standard deviation for each column (energy level)

% --- Avoid dividing by zero and avoid inflating data ---
tol = 1e-12;
X_std_safe = X_std;
X_std_safe(X_std < tol) = 1;    % any positive constant works (zero variance -> no contribution to PCA)
 
% --- Standardize Thermodynamic Features ---
X_z = X_centered ./ X_std_safe;

% Step 3: Run PCA
[coeff_X, score_X, latent_X, ~, explained_X, mu_pca] = pca(X_z);

% Step 4: Visualize
% Print table of principal component loadings
feature_names = {'AvgEnergy', 'Entropy', 'PartitionFunction', ...
    'HelmholtzEnergy', 'PercentInaccessible'}

disp(array2table(coeff_X(:,1:3), 'RowNames', feature_names, 'VariableNames', {'PC1', 'PC2', 'PC3'}))

% --- Create a 3D scatter plot whose axes correspopnd with first three PC vectors ---
figure(); 
scatter3(score_X(:,1), score_X(:,2), score_X(:,3), 60, entropy, 'filled');

label_opts_axes = {'Fontsize', 14, 'Interpreter', 'latex'};
xlabel('PC 1', label_opts_axes{:});
ylabel('PC 2', label_opts_axes{:});
zlabel('PC 3', label_opts_axes{:});
title('PCA of Thermodynamic Feature Matrix', 'FontSize', 16, 'Interpreter', 'latex')
colorbar; colormap cool;
grid on; hold on;

label_opts_data = {'FontSize', 10, 'Color', 'k', 'FontWeight', 'bold'}
annotateSpecies(score_X(:,1:3), entropy, 5, gca, label_opts_data{:});

feature_variable_names = {'$\langle E \rangle$', '$S$', '$Z$', '$F$', '$\%\mathrm{NA}$'};
drawLoadings3D(coeff_X, score_X, [1 2 3], feature_variable_names);
draw3DAxes(gca) % draw x, y, and z axes

view(135, 30);  % set a nice 3D viewing angle
hold off

% --- Create a bar plot of the explained variances for the first 5
% principal components ---
figure();
bar(explained_X(1:5));
xlabel('Principal Components', label_opts_axes{:});
ylabel('Variance Explained', label_opts_axes{:});
title('Explained Variance by Principal Components - Thermodynamic Feature Matrix', ...
    'FontSize', 16, 'Interpreter', 'latex');
grid on;

% =====================================================================

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
