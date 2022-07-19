function distributionFigure(Ndrug,NN,EN,binIncrements,maxbin)
    bins=-1-binIncrements/2:binIncrements:maxbin;
    [count_NN,xcenters]=hist(NN,bins); 
    [count_EN,~]=hist(EN,bins); 
    size(count_NN)
    size(count_EN)
    xtick=-1:0.5:maxbin;
    figure; hold on; set(gca,'fontsize',24,'FontName','Arial');
    
    b=bar(xcenters,[count_NN',count_EN'],'BarWidth', 1); hold on;
    b(1).EdgeColor = 'k'; b(1).FaceColor = 'w'; %Net Colors
    b(2).EdgeColor = 'k'; b(2).FaceColor = 'k'; %Emergent Colors

    ymax=max([count_NN,count_EN]+1); ymax=ymax+ymax/5;
    xlim([-1-0.2,maxbin+0.2]); ylim([0, ymax]); set(gca,'ygrid','on');
    set(gca,'XTick',xtick,'XTickLabel',[num2cell(xtick(1:end-1)),'>3'])
    ylabel('Drug combinations');
    xlabel([num2str(Ndrug),'-drug interaction metric']); box on;

end