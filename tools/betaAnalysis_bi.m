function betaAnalysis_bi(measure,groups_y,header,d,flag_pca,facecolor,groups_y2,facecolor2)
% facecolor = defaultColor(length(header));
d = min(d,3);
if flag_pca==0
    dis = measure;
    pcoa = f_pcoa(dis,1);
    mapped_data = pcoa.scores(:,1:d)'*2;
    power = pcoa.expl(:,1);
    plotFIGURE_bi(mapped_data,groups_y,facecolor);
    hold on
    for i=1:3
        plotEllipse_bi(mapped_data(1:2,groups_y2==i),facecolor2(i,:))
    end
    set(gca,'FontSize',14);
    xlabel(['PCo1 (' num2str(round(power(1)*10)/10) '%)']);
    ylabel(['PCo2 (' num2str(round(power(2)*10)/10) '%)']);
    if d==3
        ylabel(['PCo3 (' num2str(round(power(2)*10)/10) '%)']);
    end
else
    plotPCA(measure,groups_y,facecolor,d);
end
legend(header);
set(gca,'FontSize',14);
pbaspect([1 1 1]);
box on
end