function Main_differential_analysis_metatranscriptomics
% Run differential (group) plotting for metatranscriptomics taxa
%
%   This entry-point function loads precomputed metatranscriptomics data (species), 
% prepares color palettes and group labels, and calls a plotting helper to visualize
%   relative abundances (percent) grouped by a metadata category.
%
%   Notes:
%     - This function uses relative abundances (0-1) from x_tbl.rel and
%       converts them to percentages for plotting (multiply by 100).
%

clc;clear;close all;

% ---------------------------
% Load data and basic variables
% ---------------------------
load('../data/data_metatranscriptomics_composition.mat');
tax = x_tbl.tax;
Rel = x_tbl.rel;

% ---------------------------
% Grouping labels
% ---------------------------
des = catLabel(meta,4); % grouping by serum

% ---------------------------
% Color palettes and plotting
% ---------------------------
% Define color matrix (normalized RGB).
FaceColor = [51	68	17	238	204	43	221
    187	170	119	119	51	144	85
    238	153	51	51	17	102	34]'/255;

% Build a secondary palette by blending selected colors
FaceColor2 = [FaceColor(1,:)
    (FaceColor(2,:)+FaceColor(3,:))/2
    (FaceColor(4,:)+FaceColor(5,:))/2
    ];
FaceColor(FaceColor>1)  =1;

% Plot grouped taxa abundances as percentages
fun_plotGroup(Rel*100,tax,des.y,des.legend,FaceColor2);
xtickangle(90)
ytickangle(90)
pbaspect([4 1.5 1]);

% Determine the significance of difference between groups
figure
logRel = log10(Rel+10^-6); % add a small number (10^-6) to avoid log10(0)
for i=1:length(tax)
    [anova_p(i,1),tbl,stats] = anova1(logRel(i,:)',des.y,"off");
    c = multcompare(stats,'CriticalValueType','bonferroni');
    pairwise_p(i,:) = c(:,6)';
end
end
