function Main_beta_diversity_analysis_metatrasncriptomics
% Run beta-diversity PCoA and correlation biplot
%
%   This script loads metatranscriptomics composition (Species) data, computes Bray–Curtis distances,
%   runs PCoA-based beta diversity analyses for different groupings, and
%   overlays correlations between taxa relative abundances and the first two
%   PCoA axes.

clc;clear;close all

% ---------------------------
% Load data and basic variables
% ---------------------------
load('../data/data_metatranscriptomics_composition.mat');
tax = x_tbl.tax; 
Rel = x_tbl.rel;

% ---------------------------
% Define color palettes
% ---------------------------
FaceColor = [51  187 238
             170 220 120
             17  119 51
             255 200 210
             204 51  17] / 255;   % 5-color palette (RGB normalized)

% Create a second palette (blend specific colors)
FaceColor2 = [ FaceColor(1,:)
               (FaceColor(2,:)+FaceColor(3,:))/2
               (FaceColor(4,:)+FaceColor(5,:))/2 ];

FaceColor(FaceColor>1) = 1;   % Clamp any values >1

% ---------------------------
% Distance matrix (samples x samples)
% ---------------------------
Dvec = f_dis(Rel' , 'bc');   % Bray–Curtis distance

% ---------------------------
% Prepare grouping variables
% ---------------------------
type1 = catLabel(meta,4);  % grouping by serum
type2 = catLabel(meta,3);  % grouping by serum distinguishing initial and shift

groups_y = type1.y;
header = type1.legend;

gy2 = type2.y;
gh2 = type2.legend;

% ---------------------------
% Beta diversity analysis (PCoA) settings
% ---------------------------

dis = Dvec;
pcoa = f_pcoa(dis,1);
mapped = pcoa.scores(:,1:2)';

% ---------------------------
% Correlate taxa with PCoA axes
% ---------------------------
% For components 1..2 compute Spearman correlation between taxon relative
% abundance (across samples) and mapped PCoA axis scores.

for i=1:2
    for d=1:length(tax)
        [rho(d,i),p(d,i)]=corr(mapped(i,:)',Rel(d,:)','type','spearman');
    end
end

% ---------------------------
% Multiple testing correction (simple Bonferroni-like)
% ---------------------------
adjusted_pvals =  f_Bonferroni(p); 
pp = (adjusted_pvals(:,1)) < 0.05;   % taxa significant with axis 1 after correction

% Prepare color array for biplot (default zeros => black)
colorArray = zeros(size(adjusted_pvals,1), 3);

% ---------------------------
% Final plotting: biplot overlay
% ---------------------------
betaAnalysis_bi(dis, gy2, gh2, 2, 0, FaceColor, groups_y, FaceColor2);
hold on
% Overlay taxa that are significant on the biplot using f_biplot
f_biplot(rho(pp,:), 1, 0, tax(pp), colorArray(pp,:));
legend(gh2);
axis([-1.6 1.6 -1.2 1.2]);
axis equal;
set(gca, 'FontSize', 20);
end

