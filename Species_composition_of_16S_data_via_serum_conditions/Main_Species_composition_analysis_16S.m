function Main_Species_composition_analysis_16S
%   This entry-point function loads precomputed 16S data, prepares color
%   palettes and group labels, and calls a plotting helper to visualize
%   relative abundances (percent) grouped by a metadata category.
%
%   Notes:
%     - This function uses relative abundances (0-1) from x_tbl.rel and
%       converts them to percentages for plotting (multiply by 100).

clc;clear;close all;

% ---------------------------
% Load data
% ---------------------------
load('../data/data_16S.mat');
tax = x_tbl.tax;
Rel = x_tbl.rel;

% ---------------------------
% Grouping labels
% ---------------------------
des = catLabel(meta,7); % grouping by serum

% ---------------------------
% Color palettes and plotting
% ---------------------------
% Define color matrix (normalized RGB).
FaceColor = [51	68	17	238	204	43	221
    187	170	119	119	51	144	85
    238	153	51	51	17	102	34]'/255;
FaceColor(FaceColor>1)  =1;

% Build a secondary palette by blending selected colors
FaceColor2 = [FaceColor(1,:)
    (FaceColor(2,:)+FaceColor(3,:))/2
    (FaceColor(4,:)+FaceColor(5,:))/2
    ];

% Plot grouped taxa abundances as percentages
fun_plotGroup(Rel*100,tax,des.y,des.legend,FaceColor2);
xtickangle(90)
ytickangle(90)
pbaspect([4 1.5 1]);

end
