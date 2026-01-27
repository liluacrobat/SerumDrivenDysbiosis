function Main_beta_diversity_metabolite_ESI_pos
% This function performs PCA on the normalized LC-MS data of ESI Positive
%
%
clc;clear;close all
rng default

% Read normalized positive-mode metabolomics table (tab-delimited)
tbl = readtable('../data/pos_normalized.txt','delimiter','\t');
% Convert feature columns (from column 6 onward) to a matrix:
mx = table2array(tbl(:,6:end))';

% Extract categorical labels for plotting
cat = catLabel(tbl,4);% grouping by serum distinguishing initial or shift
g = catLabel(tbl,5);  % grouping by serum
group =cat.y;

% Define primary color palette (rows = colors)
FaceColor = [51	187	238
    170 220 120
    17	119	51
    255 200 210
    204	51	17]/255;
% Build a secondary palette by blending selected colors
FaceColor2 = [FaceColor(1,:)
    (FaceColor(2,:)+FaceColor(3,:))/2
    (FaceColor(4,:)+FaceColor(5,:))/2
    ];
FaceColor(FaceColor>1)  =1;

% Plot PCA for the primary grouping and the secondary grouping (with ellipses)
plotPCA(mx,group,FaceColor,2,0);
hold on
plotPCA2(mx,g.y,FaceColor2,2,1);

% Legend and aesthetics
header = {'0% serum'
    '5% serum (initial)'
    '5% serum (shift)'
    '50% serum (initial)'
    '50% serum (shift)'};
legend(header,'location','eastoutside');
set(gca,'fontsize',18);
pbaspect([1 1 1]);
end
function [mapped_data,power,fig,FaceColor]=plotPCA2(X,Y,FaceColor,d,flag)
% X: m*n data, m rows n samples
% Y: group label, 1,2, 3, 4, 5, 6 ...
% FaceColor: color used for each group
% flag: plot ellipses

if nargin<2
    Y = ones(size(X,2),1);
end
if nargin<4
    d = 2;
end
if nargin<5
    flag = 1;
end
U=sort(unique(Y));
str=cell(1,length(U));

if nargin<3
    FaceColor = defaultColor(length(unique(Y)));

end

[mapped_data,~,power]=compute_mapping(X','PCA',d);
mapped_data=mapped_data';


hold on
if d==3
    for i=1:length(U)
        plot3(mapped_data(1,Y==U(i)),mapped_data(2,Y==U(i)),mapped_data(3,Y==U(i)),'o','MarkerFaceColor',FaceColor(i,:),'MarkerSize',12,'MarkerEdgeColor','k');
    end
    xlabel(['PC1 (' num2str(round(power(1)*10000)/100) '%)']);
    ylabel(['PC2 (' num2str(round(power(2)*10000)/100) '%)']);
    zlabel(['PC3 (' num2str(round(power(3)*10000)/100) '%)']);
else
    if flag == 1
        for i=1:length(U)
            if sum(Y==U(i))>2
                plotEllipse(mapped_data(1:2,Y==U(i)),FaceColor(i,:))
            end
        end
    end
end
set(gca,'FontSize',14);
pbaspect([1 1 1]);
box on

end
function [mapped_data,power,fig,FaceColor]=plotPCA(X,Y,FaceColor,d,flag)
% X: m*n data, m rows n samples
% Y: group label, 1,2, 3, 4, 5, 6 ...
% FaceColor: color used for each group
% flag: plot ellipses
if nargin<2
    Y = ones(size(X,2),1);
end
if nargin<4
    d = 2;
end
if nargin<5
    flag = 1;
end
U=sort(unique(Y));
str=cell(1,length(U));

if nargin<3
    FaceColor = defaultColor(length(unique(Y)));

end

[mapped_data,~,power]=compute_mapping(X','PCA',d);
mapped_data=mapped_data';


hold on
if d==3
    for i=1:length(U)
        plot3(mapped_data(1,Y==U(i)),mapped_data(2,Y==U(i)),mapped_data(3,Y==U(i)),'o','MarkerFaceColor',FaceColor(i,:),'MarkerSize',8,'MarkerEdgeColor','k');
    end
    xlabel(['PC1 (' num2str(round(power(1)*10000)/100) '%)']);
    ylabel(['PC2 (' num2str(round(power(2)*10000)/100) '%)']);
    zlabel(['PC3 (' num2str(round(power(3)*10000)/100) '%)']);
else
    for i=1:length(U)
        plot(mapped_data(1,Y==U(i)),mapped_data(2,Y==U(i)),'o','MarkerFaceColor',FaceColor(i,:),'MarkerSize',10,'MarkerEdgeColor','k');
    end
    xlabel(['PC1 (' num2str(round(power(1)*10000)/100) '%)']);
    ylabel(['PC2 (' num2str(round(power(2)*10000)/100) '%)']);
    if flag == 1
        for i=1:length(U)
            if sum(Y==U(i))>2
                plotEllipse(mapped_data(1:2,Y==U(i)),FaceColor(i,:))
            end
        end
    end
end
set(gca,'FontSize',14);
pbaspect([1 1 1]);
box on
end
function plotEllipse(X,FaceColor)
data = X';
% Calculate the eigenvectors and eigenvalues
covariance = cov(data);
[eigenvec, eigenval ] = eig(covariance);

% Get the index of the largest eigenvector
[largest_eigenvec_ind_c, r] = find(eigenval == max(max(eigenval)));
largest_eigenvec = eigenvec(:, largest_eigenvec_ind_c);

% Get the largest eigenvalue
largest_eigenval = max(max(eigenval));

% Get the smallest eigenvector and eigenvalue
if(largest_eigenvec_ind_c == 1)
    smallest_eigenval = max(eigenval(:,2));
    smallest_eigenvec = eigenvec(:,2);
else
    smallest_eigenval = max(eigenval(:,1));
    smallest_eigenvec = eigenvec(1,:);
end

% Calculate the angle between the x-axis and the largest eigenvector
angle = atan2(largest_eigenvec(2), largest_eigenvec(1));

% This angle is between -pi and pi.
% Let's shift it such that the angle is between 0 and 2pi
if(angle < 0)
    angle = angle + 2*pi;
end

% Get the coordinates of the data mean
avg = mean(data);

% Get the 95% confidence interval error ellipse
chisquare_val = 2.4477;
theta_grid = linspace(0,2*pi);
phi = angle;
X0=avg(1);
Y0=avg(2);
a=chisquare_val*sqrt(largest_eigenval);
b=chisquare_val*sqrt(smallest_eigenval);

% the ellipse in x and y coordinates
ellipse_x_r  = a*cos( theta_grid );
ellipse_y_r  = b*sin( theta_grid );

%Define a rotation matrix
R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];

%let's rotate the ellipse to some angle phi
r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;

% Draw the error ellipse
plot(r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0,'-','color',FaceColor,'linewidth',1);
hold on;
fill(r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0,FaceColor,'EdgeColor','none','FaceColor',FaceColor,'FaceAlpha',0.15);
%
% % Plot the original data
% plot(data(:,1), data(:,2), '.');
% mindata = min(min(data));
% maxdata = max(max(data));
% Xlim([mindata-3, maxdata+3]);
% Ylim([mindata-3, maxdata+3]);
% hold on;

% % Plot the eigenvectors
% quiver(X0, Y0, largest_eigenvec(1)*sqrt(largest_eigenval), largest_eigenvec(2)*sqrt(largest_eigenval), '-m', 'LineWidth',2);
% quiver(X0, Y0, smallest_eigenvec(1)*sqrt(smallest_eigenval), smallest_eigenvec(2)*sqrt(smallest_eigenval), '-g', 'LineWidth',2);

end
