function thermo_pca_demo
%MAIN_DEMO Toy thermodynamics & PCA example.
%
%   This is a self-contained demonstration. All physical constants are
%   hard-coded. Helper functions are defined in their own *.m files.
%
%   OUTPUT - the script prints two tables and opens six figures.
%
%   UNITS:
%       Energy levels (E)               -> eV
%       Temperature (T)                 -> K
%       Boltzmann constant (kB)         -> 8.617e-5 eV/K
%       Partition function (Z)          -> unitless
%       Probabilities (Z)               -> unitless
%       Average energy (avgE)           -> eV
%       Entropy (S)                     -> eV/K

clc; clear; close all;

%% ------------------------------------------------------------------
%  1.  User-set parameters
% -------------------------------------------------------------------
nSpecies    = 10;   % number of chemical species
nLevels     = 5;    % energy levels per species
TEMP        = 300;  % temperature (K)
THRESH_P    = 1e-6; % "inaccessible" probability threshold
PC_TO_PLOT  = 1:3;  % principal components to show in 3-D plots
N_LABEL     = 5;    % how many species to annotate in scatter plots

%% ------------------------------------------------------------------
%  2.   Generate random energy landscape
% -------------------------------------------------------------------
E = rand(nSpecies, nLevels);                % 0-1 eV uniform
levelNames = "E" + (1:nLevels);             % string array

table_E = array2table(E, VariableNames = levelNames);
table_E.Species = (1:nSpecies).';                       % add index column
table_E = movevars(table_E, "Species", "Before", 1);
disp(table_E);

%% ------------------------------------------------------------------
%  3.   Boltzmann statistics
% -------------------------------------------------------------------
kB = 8.617e-5;              % Boltzmann constant (eV/K)
beta = 1 / (kB * TEMP);     % thermodynamic beta (1/eV)

Z       = sum(exp(-beta*E), 2);     % partition function for each species
P       = exp(-beta*E) ./ Z;         % normalized Boltzmann probabilities 
avgE    = sum(P .* E, 2);           % average energy for each species
S       = -kB * sum(P .* log(max(P, eps)), 2);      % entropy (eV/K)
F       = -kB * TEMP  * log(Z);                     % Helmholtz (eV)
pctNA   = 100 * sum(P < THRESH_P, 2) / nLevels;     % inaccessible %

table_thermo = table((1:nSpecies).', avgE, S, Z, F, pctNA, ...
    'VariableNames', ...
    {'Species', 'AvgE', 'Entropy', 'Z', 'F', 'PcntInaccess'});
disp(table_thermo)

%% ------------------------------------------------------------------
%  4.   Visualize probability matrix
% -------------------------------------------------------------------
utils.newFigure("Boltzmann probability heat-map");
imagesc(P); axis square; colormap(flipud(summer(256))); colorbar;
set(gca,'FontSize',12,'TickDir','out');
xlabel('$j$ (energy-level index)', 'Interpreter','latex');
ylabel('$i$ (species index)',      'Interpreter','latex');
title('Boltzmann Probability matrix $P_{ij}$','Interpreter','latex')

%% ------------------------------------------------------------------
%  5.   PCA on probability matrix (compositional data)
% -------------------------------------------------------------------
P_clr = utils.clr(P);           % centered log-ratio transfrom
[coeff_P,score_P,latent_P,explained_P] = pca(zscore(P_clr));

utils.listLoadings(coeff_P(:,1:3), levelNames, "Probabilities - loadings");

utils.newFigure("PCA - Boltzmann probabilities (3-D)");
scatter3(score_P(:,1),score_P(:,2),score_P(:,3),60,S,'filled');

hold on;
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title('PCA of Boltzmann distributions'); 
grid on; colorbar; colormap cool;
utils.annotateSpecies(score_P(:,PC_TO_PLOT),S,N_LABEL);
utils.drawLoadings3D(coeff_P,score_P,PC_TO_PLOT,levelNames);
utils.draw3DAxes(gca); view(135,30); 
hold off;

utils.newFigure("PCA - Boltzmann Probabilities (variance)");
bar(explained_P(1:5));
xlabel('PC');
ylabel('% variance explained');
grid on;

%% ------------------------------------------------------------------
%  6.   PCA on thermoddynamic feature matrix
% -------------------------------------------------------------------
X = zscore([avgE,S,Z,F,pctNA]);                   % standardised
    featNames = {'⟨E⟩','S','Z','F','%NA'};
    [coeffX,scoreX,~,~,explX] = pca(X);

    utils.listLoadings(coeffX(:,1:3), featNames, "Thermo features – loadings");

    utils.newFigure("PCA – thermo (3-D)");
    scatter3(scoreX(:,1),scoreX(:,2),scoreX(:,3),60,S,'filled'); hold on;
    xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
    title('PCA of thermodynamic features'); 
    grid on; colorbar; colormap cool
    utils.annotateSpecies(scoreX(:,PC_TO_PLOT),S,N_LABEL);
    utils.drawLoadings3D(coeffX,scoreX,PC_TO_PLOT,featNames);
    utils.draw3DAxes(gca);  view(135,30);  hold off;

    utils.newFigure("PCA – thermo (variance)");
    bar(explX(1:5)); xlabel('PC'); ylabel('% var. explained'); grid on;
end
